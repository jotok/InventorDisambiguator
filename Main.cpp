#include <cassert>
#include <cstdio>
#include <ctime>
#include <set>
#include <string>
#include <vector>

#include <Eigen/SparseCore>

#ifdef _WIN32
#define PATHSEP '\\'
#else
#define PATHSEP '/'
#endif

#define WORDBUFSIZE 1000
#define PATHBUFSIZE 1000

char pathbuffer[PATHBUFSIZE];

typedef Eigen::Triplet<bool> Tripletb;
typedef std::vector<Eigen::Triplet<bool>> TripletbList;

/**
 * A bare-bones sparse boolean matrix type.
 */
struct Sp {
    template <typename T>
    Sp(Eigen::SparseMatrix<T>& mat);

    int m;             // number of rows
    int n;             // number of columns
    const int* colptr; // column i is in colptr[i]:(colptr[i+1]-1)
    const int* rowval; // row values of nonzeros
};

/**
 * Construct a Sp struct which is a view to the data in an Eigen::SparseMatrix.
 */
template <typename T>
Sp::Sp(Eigen::SparseMatrix<T>& mat) 
    : m{mat.rows()}
    , n{mat.cols()}
    , colptr{mat.outerIndexPtr()}
    , rowval{mat.innerIndexPtr()}
{}

/**
 * A view to a single column in a sparse boolean matrix.
 */
struct ColumnView {
    ColumnView(int, const int*);
    int nnz;
    const int* rowval;
};

ColumnView::ColumnView(int nnz, const int* rowval)
    : nnz{nnz}
    , rowval{rowval}
{}

/**
 * Get a view of the data in column `col` in the sparse boolean matrix `mat`.
 */
const ColumnView getColumnView(const Sp& mat, int col) {
    int start = mat.colptr[col];
    int nnz = mat.colptr[col + 1] - start;
    const int* rowval = mat.rowval + start;

    return ColumnView(nnz, rowval);
}

/**
 * Search the array `ary` for an element greater than or equal to `value` between 
 * the indices `start` and `end` (inclusive). If such an element is found, write
 * the index of the element to `result`.
 *
 * \param ary A sorted array of int values to be searched.
 * \param start The starting index into the array where we will search.
 * \param end The last index in the array where we will search (inclusive).
 * \param value The value we are searching for.
 * \param result A pointer to an inter where we will store the index of the found
 * element.
 *
 * \return 1 if an element is found, 0 otherwise.
 */
int binarySearch(const int* ary, int start, int end, int value, int* result) {
    assert(start <= end);
    
    if (end - start <= 1) {
        if (ary[start] >= value) {
            *result = start;
            return 1;
        }

        if (ary[end] >= value) {
            *result = end;
            return 1;
        }

        return 0;
    }
    
    int mid = (start + end) / 2;

    if (ary[mid] >= value)
        return binarySearch(ary, start, mid, value, result);
    else
        return binarySearch(ary, mid, end, value, result);
}

/**
 * Finds the nonzero indices of `step & fn & ln & as & ct & cl` and adds them to `index`.
 */
void makeLumpIndex1(const std::vector<bool>& step,
                    const ColumnView& fn,
                    const ColumnView& ln,
                    const ColumnView& as,
                    const ColumnView& ct,
                    const ColumnView& cl,
                    std::vector<int>& index)
{
    if (fn.nnz == 0 || ln.nnz == 0 || as.nnz == 0 || ct.nnz == 0 || cl.nnz == 0)
        return;

    int ifn = 0;
    int iln = 0;
    int ias = 0;
    int ict = 0;
    int icl = 0;

    int r = fn.rowval[ifn];
    if (ln.rowval[iln] > r) r = ln.rowval[iln];
    if (as.rowval[ias] > r) r = as.rowval[ias];
    if (ct.rowval[ict] > r) r = ct.rowval[ict];
    if (cl.rowval[icl] > r) r = cl.rowval[icl];

    int rfn, rln, ras, rct, rcl, temp;

    while (true) {
        if (!binarySearch(ln.rowval, iln, ln.nnz - 1, r, &temp))
            break;
        iln = temp;
        rln = ln.rowval[iln];

        if (!binarySearch(fn.rowval, ifn, fn.nnz - 1, r, &temp))
            break;
        ifn = temp;
        rfn = fn.rowval[ifn];

        if (!binarySearch(as.rowval, ias, as.nnz - 1, r, &temp))
            break;
        ias = temp;
        ras = as.rowval[ias];

        if (!binarySearch(ct.rowval, ict, ct.nnz - 1, r, &temp))
            break;
        ict = temp;
        rct = ct.rowval[ict];

        if (!binarySearch(cl.rowval, icl, cl.nnz - 1, r, &temp))
            break;
        icl = temp;
        rcl = cl.rowval[icl];

        if (rfn == rln && rln == ras && ras == rct && rct == rcl) {
            if (step[rfn])
                index.push_back(rfn);

            r++;
        }

        if (rfn > r) r = rfn;
        if (rln > r) r = rln;
        if (ras > r) r = ras;
        if (rct > r) r = rct;
        if (rcl > r) r = rcl;
    }
}

/**
 * Finds the nonzero indices of `step & nm & (as | ct | cl)` and adds them to `index`.
 */
void makeLumpIndex2(const std::vector<bool>& step,
                    const ColumnView& nm,
                    const ColumnView& as,
                    const ColumnView& ct,
                    const ColumnView& cl,
                    std::vector<int>& index)
{
    if (nm.nnz == 0)
        return;

    int inm = 0;
    int ias = 0;
    int ict = 0;
    int icl = 0;
    int foundas = 1;
    int foundct = 1; 
    int foundcl = 1;
    int temp;

    int rnm;

    while (inm < nm.nnz) {
        rnm = nm.rowval[inm];

        if (!step[rnm]) {
            inm++;
            continue;
        }

        if (foundas && (foundas = binarySearch(as.rowval, ias, as.nnz - 1, rnm, &temp)))
            ias = temp;

        if (foundct && (foundct = binarySearch(ct.rowval, ict, ct.nnz - 1, rnm, &temp)))
            ict = temp;

        if (foundcl && (foundcl = binarySearch(cl.rowval, icl, cl.nnz - 1, rnm, &temp)))
            icl = temp;

        if (!(foundas || foundct || foundcl))
            break;

        if (as.rowval[ias] == rnm || ct.rowval[ict] == rnm || cl.rowval[icl] == rnm) {
            index.push_back(rnm);
        }

        inm++;
    }
}

/**
 * Read in the attribute matrix file and create a list of triplets that will be used to
 * set the nonzero values in an Eigen::SparseMatrix object.
 */
std::pair<int, int> loadAttributeMatrixTriplets(const char* directory, TripletbList& list) {

    *pathbuffer = '\0';
    snprintf(pathbuffer, PATHBUFSIZE, "%s%c%s", directory, PATHSEP, "_attribute_matrix");
    printf("Loading attribute matrix from %s\n", pathbuffer);

    FILE* in = fopen(pathbuffer, "r");

    int i, j;
    std::pair<int, int> size{0, 0};

    while (fscanf(in, "%d,%d,%*d", &i, &j) == 2) {

        if (i > size.first)
            size.first = i;

        if (j > size.second)
            size.second = j;

        list.push_back(Tripletb(i - 1, j - 1, 1));
    }

    fclose(in);
    return size;
}

/**
 * Load the attribute dictionary.
 */
void loadAttributeDictionary(const char* directory, std::vector<std::string>& words) {

    *pathbuffer = '\0';
    snprintf(pathbuffer, PATHBUFSIZE, "%s%c%s", directory, PATHSEP, "_attribute_dictionary");
    printf("Loading attribute dictionary from %s\n", pathbuffer);

    FILE* in = fopen(pathbuffer, "r");
    char word[WORDBUFSIZE];

    while (fscanf(in, "%*d %s %*d", word) == 1) {
        words.push_back(std::string(word));
    }

    fclose(in);
}

/**
 * Load the attribute metric file. This function assumes that `patno` and `inventorID` 
 * are already malloc'd.
 */
void loadAttributeMetric(const char* directory, 
                         char** patno,
                         int* patnoCount,
                         char** inventorID,
                         int* inventorIDCount) 
{
    *pathbuffer = '\0';
    snprintf(pathbuffer, PATHBUFSIZE, "%s%c%s", directory, PATHSEP, "_attribute_metric");
    printf("Loading attribute metrics from %s\n", pathbuffer);

    FILE* in = fopen(pathbuffer, "r");
    char c, pn[WORDBUFSIZE], id[WORDBUFSIZE];
    *patnoCount = *inventorIDCount = 0;
    int i, j;

    c = fgetc(in);
    while (c != EOF) {
        // read the patent number
        i = 0;
        while (c != '\t') {
            pn[i++] = c;
            c = fgetc(in);
        } 
        pn[i++] = '\0';

        // skip 5 tabs
        j = 0;
        while (j < 5) {
            while ((c = fgetc(in)) != '\t')
                ;

            j++;
        }

        // read the ID
        j = 0;
        while ((c = fgetc(in)) != '\n') {
            id[j++] = c;
        }
        id[j++] = '\0';

        char* thispn = (char*)malloc(i * sizeof(char));
        strncpy(thispn, pn, i);
        patno[(*patnoCount)++] = thispn;

        char* thisid = (char*)malloc(j * sizeof(char));
        strncpy(thisid, id, j);
        inventorID[(*inventorIDCount)++] = thisid;

        c = fgetc(in);
    }

    fclose(in);
}

/**
 * Load keys from the disambiguator input file. This function assumes `key` is already malloc'd
 */
void loadDisambiguatorInput(const char* directory, char** key, int* keyCount) {
    *pathbuffer = '\0';
    snprintf(pathbuffer, PATHBUFSIZE, "%s%c%s", directory, PATHSEP, "_disambiguator_input.csv");
    printf("Loading disambiguator input from %s\n", pathbuffer);

    FILE* in = fopen(pathbuffer, "r");
    char c, ky[WORDBUFSIZE];
    *keyCount = 0;
    int i;

    c = fgetc(in);
    while (c != EOF) {
        i = 0;
        while (c != '\t') {
            ky[i++] = c;
            c = fgetc(in);
        }
        ky[i++] = '\0';

        char* thisky = (char*)malloc(i * sizeof(char));
        strncpy(thisky, ky, i);
        key[(*keyCount)++] = thisky;

        while ((c = fgetc(in)) != '\n')
            ;

        c = fgetc(in);
    }
}

/**
 * Given the list of words loaded from the attribute dictionary, record the indices
 * that begin with the given prefix.
 */
void findAttributeIndices(const std::vector<std::string>& words,
                          const std::string& prefix,
                          std::vector<int>& keep)
{
    keep.clear();
    int n = prefix.size();
    int i = 0;

    for (auto& w: words) {
        if (w.compare(0, n, prefix) == 0)
            keep.push_back(i);

        i += 1;
    }
}

/**
 * Given a sparse matrix, return a new sparse matrix by extracting a list of columns.
 *
 * \param mat A sparse matrix object.
 * \param keep A list of column indices that we want to keep from `mat`
 *
 * \return A `mat.rows() x keep.size()` sparse matrix consisting of the specified columns
 * from `mat`.  
 */
template <typename T>
Eigen::SparseMatrix<T> extractColumns(Eigen::SparseMatrix<T> mat, std::vector<int> keep) {
    Eigen::SparseMatrix<T> result(mat.rows(), keep.size());
    TripletbList list;

    const int* colptr = mat.outerIndexPtr();
    const int* rowval = mat.innerIndexPtr();

    for (int i = 0; i < keep.size(); i++) {
        int col = keep[i];
        int start = colptr[col];
        int nnz = colptr[col + 1] - start;
        const int* rowval_ = rowval + start;

        for (int j = 0; j < nnz; j++) {
            list.push_back(Tripletb(rowval_[j], i, 1));
        }
    }

    result.setFromTriplets(list.begin(), list.end());
    return result;
}

/**
 * Return the first row index containing a nonzero value in column `col` of `mat`.
 */
inline int find(const Sp& mat, int col) { return mat.rowval[mat.colptr[col]]; }

int main(int argc, char* argv[]) {
    if (argc == 1) {
        printf("Usage: disambig DIRECTORY\n");
        return 1;
    }

    char* directory = argv[1];

    TripletbList attributeTriplets;
    std::pair<int, int> size = loadAttributeMatrixTriplets(directory, attributeTriplets);

    Eigen::SparseMatrix<bool> X(size.first, size.second);
    X.setFromTriplets(attributeTriplets.begin(), attributeTriplets.end());

    printf("Attribute matrix dimensions (%d, %d)\n", X.rows(), X.cols());
    printf("Nonzero entries: %d\n", X.nonZeros());

    std::vector<std::string> words;
    loadAttributeDictionary(directory, words);

    printf("Words: %d\n", words.size());

    char** patno = (char**)malloc(X.rows() * sizeof(char*));
    char** inventorID = (char**)malloc(X.rows() * sizeof(char*));
    int patnoCount, inventorIDCount;
    loadAttributeMetric(directory, patno, &patnoCount, inventorID, &inventorIDCount);

    printf("Patent numbers: %d\n", patnoCount);
    printf("Inventor IDs: %d\n", inventorIDCount);

    char** key = (char**)malloc(X.rows() * sizeof(char*));
    int keyCount;
    loadDisambiguatorInput(directory, key, &keyCount);

    printf("Keys: %d\n", keyCount);

    Eigen::SparseMatrix<bool> Xt = X.transpose();
    X.makeCompressed();
    Xt.makeCompressed();

    printf("Starting clock...\n");
    clock_t startTime = clock();

    std::vector<int> keep;
    findAttributeIndices(words, "LN", keep);
    Eigen::SparseMatrix<bool> XLN = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XLNt = XLN.transpose();
    XLN.makeCompressed();
    XLNt.makeCompressed();
    Sp XXLN(XLN);
    Sp XXLNt(XLNt);

    findAttributeIndices(words, "FN", keep);
    Eigen::SparseMatrix<bool> XFN = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XFNt = XFN.transpose();
    XFN.makeCompressed();
    XFNt.makeCompressed();
    Sp XXFN(XFN);
    Sp XXFNt(XFNt);

    findAttributeIndices(words, "NM", keep);
    Eigen::SparseMatrix<bool> XNM = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XNMt = XNM.transpose();
    XNM.makeCompressed();
    XNMt.makeCompressed();
    Sp XXNM(XNM);
    Sp XXNMt(XNMt);

    findAttributeIndices(words, "AS", keep);
    Eigen::SparseMatrix<bool> XAS = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XASt = XAS.transpose();
    XAS.makeCompressed();
    XASt.makeCompressed();
    Sp XXAS(XAS);
    Sp XXASt(XASt);

    findAttributeIndices(words, "CT", keep);
    Eigen::SparseMatrix<bool> XCT = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XCTt = XCT.transpose();
    XCT.makeCompressed();
    XCTt.makeCompressed();
    Sp XXCT(XCT);
    Sp XXCTt(XCTt);

    findAttributeIndices(words, "CL", keep);
    Eigen::SparseMatrix<bool> XCL = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XCLt = XCL.transpose();
    XCL.makeCompressed();
    XCLt.makeCompressed();
    Sp XXCL(XCL);
    Sp XXCLt(XCLt);

    clock_t thisTime = clock();
    double thisSeconds = double(thisTime - startTime) / CLOCKS_PER_SEC;
    printf("elapsed time: %f seconds (building data structures)\n", thisSeconds);

    std::vector<bool> step(X.rows(), true);

    std::vector<int> lump_index_1;
    std::vector<int> lump_index_2;
    std::set<std::string> lump_patno_2;

    std::set<int> lump;
    std::set<int> lump_initial;

    int colFN, colLN, colNM, colAS, colCT, colCL;

    for (int index = 0; index < X.rows(); index++) {
        if (!step[index])
            continue;

        step[index] = 0;

        lump_index_1.clear();
        lump_index_2.clear();
        lump_patno_2.clear();
        lump.clear();
        lump_initial.clear();

        colFN = find(XXFNt, index);
        colLN = find(XXLNt, index);
        colNM = find(XXNMt, index);
        colAS = find(XXASt, index);
        colCT = find(XXCTt, index);
        colCL = find(XXCLt, index);

        ColumnView fn = getColumnView(XXFN, colFN);
        ColumnView ln = getColumnView(XXLN, colLN);
        ColumnView nm = getColumnView(XXNM, colNM);
        ColumnView as = getColumnView(XXAS, colAS);
        ColumnView ct = getColumnView(XXCT, colCT);
        ColumnView cl = getColumnView(XXCL, colCL);

        makeLumpIndex2(step, nm, as, ct, cl, lump_index_2);

        lump_patno_2.insert(patno[index]);

        for (auto ix: lump_index_2) {
            step[ix] = 0;
            lump.insert(ix);
            lump_initial.insert(ix);
            lump_patno_2.insert(patno[ix]);
        }

        makeLumpIndex1(step, fn, ln, as, ct, cl, lump_index_1);

        for (auto ix: lump_index_1) {
            auto it = lump_patno_2.find(patno[ix]);
            
            if (it == lump_patno_2.end()) {
                step[ix] = 0;
                lump.insert(ix);
                lump_initial.insert(ix);
            }
        }

        for (int indexy: lump_initial) {

            lump_index_1.clear();
            lump_index_2.clear();

            colFN = find(XXFNt, indexy);
            colLN = find(XXLNt, indexy);
            colNM = find(XXNMt, indexy);
            colAS = find(XXASt, indexy);
            colCT = find(XXCTt, indexy);
            colCL = find(XXCLt, indexy);

            ColumnView fn_ = getColumnView(XXFN, colFN);
            ColumnView ln_ = getColumnView(XXLN, colLN);
            ColumnView nm_ = getColumnView(XXNM, colNM);
            ColumnView as_ = getColumnView(XXAS, colAS);
            ColumnView ct_ = getColumnView(XXCT, colCT);
            ColumnView cl_ = getColumnView(XXCL, colCL);

            makeLumpIndex2(step, nm_, as_, ct_, cl_, lump_index_2);

            for (int ix: lump_index_2) {
                lump.insert(ix);
                lump_patno_2.insert(patno[ix]);
            }
            makeLumpIndex1(step, fn_, ln_, as_, ct_, cl_, lump_index_1);

            for (int ix: lump_index_1) {
                auto it = lump_patno_2.find(patno[ix]);

                if (it == lump_patno_2.end()) {
                    lump.insert(ix);
                }
            }
        }

        // We could be freeing memory of inventorIDs that are reassigned,
        // but we can get away without it

        for (auto it = lump.begin(); it != lump.end(); ++it) {
            step[*it] = 0;
            inventorID[*it] = inventorID[index];
        }
    }

    *pathbuffer = '\0';
    snprintf(pathbuffer, PATHBUFSIZE, "%s%c%s", directory, PATHSEP, "_disambiguator_output_cpp.tsv");

    FILE *out = fopen(pathbuffer, "w");

    for (int i = 0; i < keyCount; i++) {
        fprintf(out, "%s\t%s\n", key[i], inventorID[i]);
    }
    fclose(out);

    clock_t endTime = clock();
    double elapsedSeconds = double(endTime - startTime) / CLOCKS_PER_SEC;
    printf("elapsed time: %f seconds (total)\n", elapsedSeconds);
}
