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

    Eigen::SparseVector<bool> C_buf_1(X.rows());
    Eigen::SparseVector<bool> C_buf_2(X.rows());

    std::vector<int> lump_index_1;
    std::vector<int> lump_index_2;
    std::vector<std::string> lump_patno_1;
    std::vector<std::string> lump_patno_2;

    std::set<int> lump;
    std::set<int> lump_initial;

    for (int index = 0; index < X.rows(); index++) {
        if (!step[index])
            continue;

        lump_index_1.clear();
        lump_index_2.clear();
        lump_patno_1.clear();
        lump_patno_2.clear();
        lump.clear();
        lump_initial.clear();

        // compute C_buf_1 = C_firstname & C_lastname & C_assignee & C_city & C_class
        // compute C_buf_2 = C_name & (C_assignee | C_city | C_class)

        int target = find(XCLt, index);
        C_buf_1 = XLN.col(target);
        C_buf_2 = XLN.col(target);

        target = find(XCTt, index);
        C_buf_1 = C_buf_1.cwiseProduct(XCT.col(target));
        C_buf_2 += XCT.col(target);

        target = find(XASt, index);
        C_buf_1 = C_buf_1.cwiseProduct(XAS.col(target));
        C_buf_2 += XAS.col(target);

        target = find(XLNt, index);
        C_buf_1 = C_buf_1.cwiseProduct(XLN.col(target));

        target = find(XFNt, index);
        C_buf_1 = C_buf_1.cwiseProduct(XFN.col(target));

        target = find(XNMt, index);
        C_buf_2 = C_buf_2.cwiseProduct(XNM.col(target));

        for (Eigen::SparseVector<bool>::InnerIterator it(C_buf_1); it; ++it) {
            if (step[it.index()]) {
                lump_index_1.push_back(it.index());
                lump_patno_1.push_back(patno[it.index()]);
            }
        }

        for (Eigen::SparseVector<bool>::InnerIterator it(C_buf_2); it; ++it) {
            if (step[it.index()]) {
                lump_index_2.push_back(it.index());
                lump_patno_2.push_back(patno[it.index()]);
            }
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

            // compute C_buf_1 = C_firstname & C_lastname & C_assignee & C_city & C_class
            // compute C_buf_2 = C_name & (C_assignee | C_city | C_class)

            int target = find(XCLt, indexy);
            C_buf_1 = XLN.col(target);
            C_buf_2 = XLN.col(target);

            target = find(XCTt, index);
            C_buf_1 = C_buf_1.cwiseProduct(XCT.col(target));
            C_buf_2 += XCT.col(target);

            target = find(XASt, index);
            C_buf_1 = C_buf_1.cwiseProduct(XAS.col(target));
            C_buf_2 += XAS.col(target);

            target = find(XLNt, index);
            C_buf_1 = C_buf_1.cwiseProduct(XLN.col(target));

            target = find(XFNt, index);
            C_buf_1 = C_buf_1.cwiseProduct(XFN.col(target));

            target = find(XNMt, index);
            C_buf_2 = C_buf_2.cwiseProduct(XNM.col(target));

            for (Eigen::SparseVector<bool>::InnerIterator it(C_buf_1); it; ++it) {
                if (step[it.index()]) {
                    lump_index_1.push_back(it.index());
                    lump_patno_1.push_back(patno[it.index()]);
                }
            }

            for (Eigen::SparseVector<bool>::InnerIterator it(C_buf_2); it; ++it) {
                if (step[it.index()]) {
                    lump_index_2.push_back(it.index());
                    lump_patno_2.push_back(patno[it.index()]);
                }
            }

            for (auto ix: lump_index_1) {
                lump.insert(ix);
            }

            for (auto ix: lump_index_2) {
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
