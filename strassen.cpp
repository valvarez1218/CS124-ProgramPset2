#include <vector>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <utility>
#include <string>
#include <fstream>

using namespace std;

const int gBASE_N = 5;

// structure used to represent a matrix
struct Matrix {
    // how many rows in matrix (height)
    int rowDim;
    // how many columns in matrix (width)
    int colDim;
    // boolean indicating whether matrix is padded or not
    bool padded;

    // 2D vector containing values
    vector<vector<int>> Mat;

    // constructor function
    Matrix(int r, int c) {
        rowDim = r;
        colDim = c;

        // initialize matrix with all 0s
        Mat.resize(r, vector<int> {});
        for (int i = 0; i < r; i++) {
            Mat[i].resize(c, 0);
        }
        padded = false;
    }

    vector<int>& operator[] (const int i) {return Mat[i];}

    // method to print matrix
    void printMat() {
        for (vector<int> row : Mat) {
            for (int entry : row) {
                cout << entry  << "  " ;
            }
            cout << endl;
        }
    }

    void printDiagonal() {
        for (int i = 0; i < rowDim; i++) {
            cout << Mat[i][i] << endl;
        }
        cout << endl;
    }

    // add zero row to end
    void padRow() {
        vector<int> zeros(colDim, 0);
        Mat.push_back(zeros);
        rowDim += 1;
    }

    // add zero column to end
    void padCol() {
        // for each row, add a zero to the end
        for (int i = 0; i < rowDim; i++) {
            Mat[i].push_back(0);
        }
        colDim += 1;
    }

    // add column and row of zeros
    void pad() {
        padRow();
        padCol();
        padded = true;
    }

    // removes 0 padding of matrix
    void depad() {
        if (!padded) {
            cout << "Warning: Attempted to depad a non-padded matrix" << endl;
            return;
        }
        Mat.pop_back();
        rowDim -= 1;
        for (int i = 0; i < rowDim; i++) {
            Mat[i].pop_back();
        }
        colDim -= 1;
        padded = false;
    }
};

// regular matrix multiplication algorithm, used after we hit base case
//      assumes input matrices have appropriate dimensions; (n x m) and (m x n)
Matrix NaiveMatMult(Matrix M1, Matrix M2) {
    Matrix toRet(M1.rowDim, M1.colDim);

    for (int i = 0; i < M1.rowDim; i++) {
        for (int j = 0; j < M1.colDim; j++) {
            toRet[i][j] = 0;
            for (int k = 0; k < M1.colDim; k++) {
                toRet[i][j] += M1[i][k]*M2[k][j];
            }
        }
    }

    return toRet;
}

// for submatrices of the form (A + B) calcualte their sum and return as a matrix
Matrix sumSectors(Matrix M, pair<int, int> start1, pair<int, int> start2, bool subtracting = false) {
    int rowDim = M.rowDim/2;
    int colDim = M.colDim/2;
    Matrix toRet(rowDim, colDim);

    for (int i=0; i < rowDim; i++) {
        for (int j = 0; j < colDim; j++) {
            // if we are subtracting matrices
            if (subtracting) {
                toRet[i][j] = M[start1.first + i][start1.second + j] - M[start2.first + i][start2.second + j];                
            }
            // otherwise add as usual
            else {
                toRet[i][j] = M[start1.first + i][start1.second + j] + M[start2.first + i][start2.second + j];
            }
        }
    }

    return toRet;
}

// for submatrices of the form A simply calcualte them and return
Matrix getSubMat(Matrix M, pair<int, int> start) {
    int rowDim = M.rowDim/2;
    int colDim = M.colDim/2;
    Matrix toRet(rowDim, colDim);

    for (int i = 0; i < rowDim; i++) {
        for (int j = 0; j < colDim; j++) {
            toRet[i][j] = M[start.first + i][start.second + j];
        }
    }

    return toRet;
}

// given the 7 products from Strassen's combine into final matrix
//      manipulate Prod in place
void combineMats(Matrix &Prod, vector<Matrix> subProducts) {
    // we should only be receiving the 7 products from Strassen's
    assert(subProducts.size() == 7);

    // combine into one large matrix
    for (int i = 0; i < Prod.rowDim; i++) {
        for (int j = 0; j < Prod.colDim; j++) {
            // perform correct sum depending on which entry we are in
            // if we're in upper left portion
            if (i < Prod.rowDim/2 && j < Prod.colDim/2) {
                // upper left is -P2 + P4 + P5 + P6
                Prod[i][j] = (-subProducts[1][i][j]) + subProducts[3][i][j] + 
                                    subProducts[4][i][j] + subProducts[5][i][j];
            }
            // if we're in upper right portion
            else if (i < Prod.rowDim/2 && j >= Prod.colDim/2) {
                // upper right is P1 + P2
                int subJ = j - Prod.colDim/2;
                Prod[i][j] = subProducts[0][i][subJ] + subProducts[1][i][subJ];
            }
            // if we're in lower left portion
            else if (i >= Prod.rowDim/2 && j < Prod.colDim/2) {
                // lower left is P3 + P4
                int subI = i - Prod.rowDim/2;
                Prod[i][j] = subProducts[2][subI][j] + subProducts[3][subI][j];
            }
            // if we're in lower right portion
            else if (i >= Prod.rowDim/2 && j >= Prod.colDim/2) {
                // lower right is P1 - P3 + P5 + P7
                int subI = i - Prod.rowDim/2;
                int subJ = j - Prod.colDim/2;
                Prod[i][j] = subProducts[0][subI][subJ] - subProducts[2][subI][subJ] + subProducts[4][subI][subJ] + subProducts[6][subI][subJ];
            }
        }
    }
}

// Strassen's Matrix Multiplication Algorithm
Matrix StrassMult(Matrix M1, Matrix M2) {
    // arguments MUST be square matrices
    assert(M1.colDim == M1.rowDim);
    assert(M2.colDim == M2.rowDim);
    assert(M1.rowDim == M2.rowDim);

    // if we have base case we use naive matrix multiplication
    if (M1.colDim <= gBASE_N) {
        return NaiveMatMult(M1, M2);
    }

    // matrix where we will store final value
    Matrix Product(M1.rowDim, M1.colDim);

    // if matrices are of odd dimension pad them with zeros
    if (M1.colDim % 2 == 1) {
        M1.pad();
        M2.pad();
        Product.pad();
    }

    // calculate the 14 unique matrices needed in the 7 multiplications
    // calculates A(F-H)
    Matrix A = getSubMat(M1, make_pair(0, 0));
    Matrix FH = sumSectors(M2, make_pair(0, M2.colDim/2), make_pair(M2.rowDim/2, M2.colDim/2), true);
    Matrix P1 = StrassMult(A, FH);

    // calculates (A+B)H
    Matrix AB = sumSectors(M1, make_pair(0, 0), make_pair(0, M1.colDim/2));
    Matrix H = getSubMat(M2, make_pair(M2.rowDim/2, M2.colDim/2));
    Matrix P2 = StrassMult(AB, H);

    // calculates (C+D)E
    Matrix CD = sumSectors(M1, make_pair(M1.rowDim/2, 0), make_pair(M1.rowDim/2, M1.rowDim/2));
    Matrix E = getSubMat(M2, make_pair(0, 0));
    Matrix P3 = StrassMult(CD, E);
    
    // calculates D(G-E)
    Matrix D = getSubMat(M1, make_pair(M1.rowDim/2, M1.colDim/2));
    Matrix GE = sumSectors(M2, make_pair(M2.rowDim/2, 0), make_pair(0, 0), true);
    Matrix P4 = StrassMult(D, GE);

    // calculates (A+D)(E+H)
    Matrix AD = sumSectors(M1, make_pair(0, 0), make_pair(M1.rowDim/2, M1.colDim/2));
    Matrix EH = sumSectors(M2, make_pair(0, 0), make_pair(M2.rowDim/2, M2.colDim/2));
    Matrix P5 = StrassMult(AD, EH);

    // calculates (B-D)(G+H)
    Matrix BD = sumSectors(M1, make_pair(0, M1.colDim/2), make_pair(M1.rowDim/2, M1.colDim/2), true);
    Matrix GH = sumSectors(M2, make_pair(M2.rowDim/2, 0), make_pair(M2.rowDim/2, M2.colDim/2));
    Matrix P6 = StrassMult(BD, GH);

    // calculates (C-A)(E+F)
    Matrix CA = sumSectors(M1, make_pair(M1.rowDim/2, 0), make_pair(0, 0), true);
    Matrix EF = sumSectors(M2, make_pair(0, 0), make_pair(0, M2.colDim/2));
    Matrix P7 = StrassMult(CA, EF);

    combineMats(Product, vector<Matrix>{P1, P2, P3, P4, P5, P6, P7});

    if (Product.padded) {
        Product.depad();
    }

    return Product;
}

Matrix generateRandMat(int dimension) {

}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::printf("Usage: ./strassen 0 dimension inputfile\n");
        return -1;
    }

    // get dimension from command line
    int dimension = strtol(argv[2], nullptr, 0);
    Matrix M1(dimension, dimension);
    Matrix M2(dimension, dimension);

    // if testing, generate two random matrices and calculate product with both Strassen's and naive algorithm
    if (strtol(argv[1], nullptr, 0) == 1) {

    }
    // if flag is '2' then calculate optimal n0
    else if (strtol(argv[1], nullptr, 0) == 2) {

    }
    // otherwise read inputs from text file
    else {
        bool enteringM2 = false;;
        string filename(argv[3]);
        ifstream infile(filename);
        int i = 0;
        int j = 0;
        int entry;
        while (infile >> entry) {
            // enter values into correct matrix
            if (enteringM2) {
                M2[i][j] = entry;
            }
            else {
                M1[i][j] = entry;
            }

            j++;
            if (j >= dimension) {
                j = 0; 
                i++;
            }
            if (i >= dimension) {
                enteringM2 = true;
                i = 0;
                j=0;
            }
        }
    }

    Matrix Prod = StrassMult(M1, M2);
    // Prod.printMat();
    Prod.printDiagonal();

    return 0;
}