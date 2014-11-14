#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <string>
#include <vector>

#include <Eigen/SparseCore>

#ifdef _WIN32
#define PATHSEP '\\'
#else
#define PATHSEP '/'
#endif

typedef Eigen::Triplet<bool> Tripletb;
typedef std::vector<Eigen::Triplet<bool>> TripletbList;

struct ColumnView {
    ColumnView(int, const int*);
    int nnz;
    const int* rowval;
};

ColumnView::ColumnView(int nnz, const int* rowval)
    : nnz{nnz}
    , rowval{rowval}
{}

template <typename T>
const ColumnView getColumnView(const Eigen::SparseMatrix<T>& mat, int col) {
    const int* colptr = mat.outerIndexPtr();
    int start = colptr[col];
    int nnz = colptr[col + 1] - start;
    const int* rowval = mat.innerIndexPtr() + start;

    return ColumnView(nnz, rowval);
}

void makeLumpIndex1(const std::vector<bool>& step,
                    const ColumnView& fn,
                    const ColumnView& ln,
                    const ColumnView& as,
                    const ColumnView& ct,
                    const ColumnView& cl,
                    std::vector<int>& index)
{
    if (fn.nnz == 0 || ln.nnz == 0 || as.nnz == 0 || ct.nnz == 0 || cl.nnz == 0)
    {
        return;
    }

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

    int rfn, rln, ras, rct, rcl;

    while (ifn < fn.nnz && iln < ln.nnz && ias < as.nnz && ict < ct.nnz && icl < cl.nnz) {

        rfn = fn.rowval[ifn];
        while (rfn < r && ifn < fn.nnz) {
            ifn++;
            rfn = fn.rowval[ifn];
        }

        rln = ln.rowval[iln];
        while (rln < r && iln < ln.nnz) {
            iln++;
            rln = ln.rowval[iln];
        }

        ras = as.rowval[ias];
        while (ras < r && ias < as.nnz) {
            ias++;
            ras = as.rowval[ias];
        }

        rct = ct.rowval[ict];
        while (rct < r && ict < ct.nnz) {
            ict++;
            rct = ct.rowval[ict];
        }

        rcl = cl.rowval[icl];
        while (rcl < r && icl < cl.nnz) {
            icl++;
            rcl = cl.rowval[icl];
        }

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

        if (rfn < r) ifn++;
        if (rln < r) iln++;
        if (ras < r) ias++;
        if (rct < r) ict++;
        if (rcl < r) icl++;
    }
}

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

    int rnm;

    while (inm < nm.nnz) {
        rnm = nm.rowval[inm];

        if (!step[rnm]) {
            inm++;
            continue;
        }

        while (ias < as.nnz && as.rowval[ias] < rnm)
            ias++;

        while (ict < ct.nnz && ct.rowval[ict] < rnm)
            ict++;

        while (icl < cl.nnz && cl.rowval[icl] < rnm)
            icl++;

        if (ias >= as.nnz && ict >= ct.nnz && icl >= cl.nnz)
            break;

        if (as.rowval[ias] == rnm || ct.rowval[ict] == rnm || cl.rowval[icl] == rnm) {
            index.push_back(rnm);
        }

        inm++;
    }
}

std::pair<int, int> loadAttributeMatrixTriplets(const std::string& directory, TripletbList& list) {

    std::string path = directory + PATHSEP + "_attribute_matrix";
    std::cout << "Loading attribute matrix from '" << path << "'" << std::endl;

    std::ifstream in;
    in.open(path);

    int i, j;
    bool b;
    char c;
    std::string line;

    std::pair<int, int> size{0, 0};

    while (std::getline(in, line)) {
        std::istringstream iss(line);
        iss >> i >> c;
        iss >> j >> c;
        iss >> b;

        if (i > size.first)
            size.first = i;

        if (j > size.second)
            size.second = j;

        list.push_back(Tripletb(i - 1, j - 1, b));
    }

    return size;
}

void loadAttributeDictionary(const std::string& directory, std::vector<std::string>& words) {
    std::string path = directory + PATHSEP + "_attribute_dictionary";
    std::cout << "Loading attribute dictionary from '" << path << "'" << std::endl;

    std::ifstream in;
    in.open(path);

    int _x, _y;

    std::string line;
    while (std::getline(in, line)) {
        std::string s;
        std::istringstream iss(line);
        iss >> _x >> s >> _y;
        words.push_back(s);
    }
}

void loadAttributeMetric(const std::string& directory, 
                         std::vector<std::string>& patno,
                         std::vector<std::string>& inventorID) 
{
    std::string path = directory + PATHSEP + "_attribute_metric";
    std::cout << "Loading attribute metrics from '" << path << "'" << std::endl;

    std::ifstream in;
    in.open(path);

    std::string line;
    std::string _s;

    while (std::getline(in, line)) {
        std::string pn;
        std::string id;
        std::istringstream iss(line);

        std::getline(iss, pn, '\t');
        std::getline(iss, _s, '\t');
        std::getline(iss, _s, '\t'); // inventor
        std::getline(iss, _s, '\t'); // assignee
        std::getline(iss, _s, '\t'); // class 
        std::getline(iss, _s, '\t'); // lastname
        std::getline(iss, id, '\t');

        patno.push_back(pn);
        inventorID.push_back(id);
    }
}

void loadDisambiguatorInput(const std::string& directory, std::vector<std::string>& key) {
    std::string path = directory + PATHSEP + "_disambiguator_input.csv";
    std::cout << "Loading disambiguator input from '" << path << "'" << std::endl;

    std::ifstream in;
    in.open(path);

    std::string line;

    while (std::getline(in, line)) {
        std::string k;
        std::istringstream iss(line);
        iss >> k;
        key.push_back(k);
    }
}

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

template <typename T>
Eigen::SparseMatrix<T> extractColumns(Eigen::SparseMatrix<T> mat, std::vector<int> keep) {
    Eigen::SparseMatrix<T> result(mat.rows(), keep.size());
    for (int i = 0; i < keep.size(); i++)
        result.col(i) = mat.col(keep[i]);

    return result;
}

template <typename T>
int find(Eigen::SparseMatrix<T> mat, int col) {
    typename Eigen::SparseMatrix<T>::InnerIterator it(mat, col);
    return it.index();
    // const int* outer = mat.outerIndexPtr();
    // const int* inner = mat.innerIndexPtr();
    // return inner[outer[col]];
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cout << "Usage: disambig DIRECTORY" << std::endl;
        return 1;
    }

    std::string directory(argv[1]);

    TripletbList attributeTriplets;
    std::pair<int, int> size = loadAttributeMatrixTriplets(directory, attributeTriplets);

    Eigen::SparseMatrix<bool> X(size.first, size.second);
    X.setFromTriplets(attributeTriplets.begin(), attributeTriplets.end());

    std::cout << "Attribute matrix dimensions: (" << X.rows() << ", " << X.cols() << ")" << std::endl;
    std::cout << "Nonzero entries: " << X.nonZeros() << std::endl;

    std::vector<std::string> words;
    loadAttributeDictionary(directory, words);

    std::cout << "Words: " << words.size() << std::endl;

    std::vector<std::string> patno;
    std::vector<std::string> inventorID;
    loadAttributeMetric(directory, patno, inventorID);

    std::cout << "Patent Numbers: " << patno.size() << std::endl;
    std::cout << "Inventor IDs: " << inventorID.size() << std::endl;

    std::vector<std::string> key;
    loadDisambiguatorInput(directory, key);

    std::cout << "Keys: " << key.size() << std::endl;

    Eigen::SparseMatrix<bool> Xt = X.transpose();
    X.makeCompressed();
    Xt.makeCompressed();

    clock_t startTime = clock();

    std::vector<int> keep;
    findAttributeIndices(words, "LN", keep);
    Eigen::SparseMatrix<bool> XLN = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XLNt = XLN.transpose();
    XLN.makeCompressed();
    XLNt.makeCompressed();

    findAttributeIndices(words, "FN", keep);
    Eigen::SparseMatrix<bool> XFN = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XFNt = XFN.transpose();
    XFN.makeCompressed();
    XFNt.makeCompressed();

    findAttributeIndices(words, "NM", keep);
    Eigen::SparseMatrix<bool> XNM = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XNMt = XNM.transpose();
    XNM.makeCompressed();
    XNMt.makeCompressed();

    findAttributeIndices(words, "AS", keep);
    Eigen::SparseMatrix<bool> XAS = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XASt = XAS.transpose();
    XAS.makeCompressed();
    XASt.makeCompressed();

    findAttributeIndices(words, "CT", keep);
    Eigen::SparseMatrix<bool> XCT = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XCTt = XCT.transpose();
    XCT.makeCompressed();
    XCTt.makeCompressed();

    findAttributeIndices(words, "CL", keep);
    Eigen::SparseMatrix<bool> XCL = extractColumns<bool>(X, keep);
    Eigen::SparseMatrix<bool> XCLt = XCL.transpose();
    XCL.makeCompressed();
    XCLt.makeCompressed();

    std::vector<bool> step(X.rows(), true);

    std::vector<int> lump_index_1;
    std::vector<int> lump_index_2;
    std::vector<std::string> lump_patno_1;
    std::vector<std::string> lump_patno_2;

    std::set<int> lump;
    std::set<int> lump_initial;

    int colFN, colLN, colNM, colAS, colCT, colCL;

    for (int index = 0; index < X.rows(); index++) {
        if (!step[index])
            continue;

        lump_index_1.clear();
        lump_index_2.clear();
        lump_patno_1.clear();
        lump_patno_2.clear();
        lump.clear();
        lump_initial.clear();

        colFN = find(XFNt, index);
        colLN = find(XLNt, index);
        colNM = find(XNMt, index);
        colAS = find(XASt, index);
        colCT = find(XCTt, index);
        colCL = find(XCLt, index);

        ColumnView fn = getColumnView(XFN, colFN);
        ColumnView ln = getColumnView(XLN, colLN);
        ColumnView nm = getColumnView(XNM, colNM);
        ColumnView as = getColumnView(XAS, colAS);
        ColumnView ct = getColumnView(XCT, colCT);
        ColumnView cl = getColumnView(XCL, colCL);

        makeLumpIndex1(step, fn, ln, as, ct, cl, lump_index_1);
        makeLumpIndex2(step, nm, as, ct, cl, lump_index_2);

        for (int ix: lump_index_1) {
            lump_patno_1.push_back(patno[ix]);
        }

        for (int ix: lump_index_2) {
            lump_patno_2.push_back(patno[ix]);
        }

        // TODO test to make sure this matches previous output

        for (auto ix: lump_index_1) {
            lump.insert(ix);
            lump_initial.insert(ix);
        }

        for (auto ix: lump_index_2) {
            lump.insert(ix);
            lump_initial.insert(ix);
        }

        for (int indexy: lump_initial) {
            lump_index_1.clear();
            lump_index_2.clear();
            lump_patno_1.clear();
            lump_patno_2.clear();

            colFN = find(XFNt, indexy);
            colLN = find(XLNt, indexy);
            colNM = find(XNMt, indexy);
            colAS = find(XASt, indexy);
            colCT = find(XCTt, indexy);
            colCL = find(XCLt, indexy);

            ColumnView fn = getColumnView(XFN, colFN);
            ColumnView ln = getColumnView(XLN, colLN);
            ColumnView nm = getColumnView(XNM, colNM);
            ColumnView as = getColumnView(XAS, colAS);
            ColumnView ct = getColumnView(XCT, colCT);
            ColumnView cl = getColumnView(XCL, colCL);

            makeLumpIndex1(step, fn, ln, as, ct, cl, lump_index_1);
            makeLumpIndex2(step, nm, as, ct, cl, lump_index_2);

            for (int ix: lump_index_1) {
                lump_patno_1.push_back(patno[ix]);
                lump.insert(ix);
            }

            for (int ix: lump_index_2) {
                lump_patno_2.push_back(patno[ix]);
                lump.insert(ix);
            }
        }

        if (lump.size() == 0) {
            step[index] = 0;
            continue;
        }

        for (auto it = lump.rbegin(); it != lump.rend(); ++it) {
            step[*it] = 0;
            inventorID[*it] = inventorID[index];
        }
    }

    std::string path = directory + PATHSEP + "_disambiguator_output_cpp.tsv";
    std::ofstream out;
    out.open(path);

    for (int i = 0; i < key.size(); i++) {
        out << key[i] << "\t" << inventorID[i] << std::endl;
    }


    clock_t endTime = clock();
    double elpasedSeconds = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "elapsed time: " << elpasedSeconds << " seconds" << std::endl;
}
