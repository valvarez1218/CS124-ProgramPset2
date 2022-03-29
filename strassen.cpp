#include <vector>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <utility>
#include <string>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <bits/stdc++.h>

using namespace std;

int N_0 = 0;
const int gBASE_N_ODD = 37;
const int gBASE_N_EVEN = 15;
// structure used to represent a matrix
struct Matrix {
    // how many rows in matrix (height)
    int rowDim;
    // how many columns in matrix (width)
    int colDim;
    // boolean indicating whether matrix is padded or not
    bool padded;

    // 2D vector containing values
    vector<vector<long>> Mat;

    // constructor function
    Matrix(int r, int c) {
        rowDim = r;
        colDim = c;

        // initialize matrix with all 0s
        Mat.resize(r, vector<long> {});
        for (int i = 0; i < r; i++) {
            Mat[i].resize(c, 0);
        }
        padded = false;
    }

    vector<long>& operator[] (const int i) {return Mat[i];}

    // method to print matrix
    void printMat() {
        for (vector<long> row : Mat) {
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
        vector<long> zeros(colDim, 0);
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
    if (N_0 != 0) {
        if (M1.colDim <= N_0) {
            return NaiveMatMult(M1, M2);
        }
    } else if (M1.colDim <= gBASE_N_ODD) {
        if (M1.colDim % 2 == 1) {
            // dimension is odd
            return NaiveMatMult(M1, M2);
        } else if (M1.colDim <= gBASE_N_EVEN) {
            // dimension is even
            return NaiveMatMult(M1, M2);
        }
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

// generate random matrices
Matrix generateRandMat(int dimension) {
    Matrix M = Matrix(dimension, dimension);
    srand(time(NULL));

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            // random integers between 0 and 100
            M[i][j] = rand() % 100;
            // make 10% of values negative
            if (rand() % 10 == 1) {
                M[i][j] *= -1;
            }
        }
    }

    return M;
}

// generate random graph
Matrix generateRandGraph(int dimension, float p) {
    Matrix M = Matrix(dimension, dimension);
    srand(time(NULL));

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            if (i == j) {
                M[i][j] = 0;
            } else {
                int edge = rand() % 100 < (int)(p * 100) ? 1 : 0;
                M[i][j] = edge;
                M[j][i] = edge;
            }
        }
    }

    return M;
}

int CountTriangles(Matrix MTri) {
    // multiply matrix by itself 3 times
    Matrix Prod = StrassMult(StrassMult(MTri, MTri), MTri);

    int rowDim = Prod.rowDim;
    int sum = 0;
    for (int i = 0; i < rowDim; i++) {
        // only add values where i = j
        sum += Prod[i][i];
    }

    return sum / 6;
}

int main(int argc, char** argv) {
    // get dimension from command line
    int dimension = strtol(argv[2], nullptr, 0);
    Matrix M1(dimension, dimension);
    Matrix M2(dimension, dimension);

    // if testing, generate two random matrices and calculate product with both Strassen's and naive algorithm
    if (strtol(argv[1], nullptr, 0) == 1) {
        if (argc != 3) {
            std::printf("Usage: ./strassen 1 dimension\n");
            return -1;
        }
        
        M1 = generateRandMat(dimension);
        M2 = generateRandMat(dimension);

        Matrix ProdNaive = NaiveMatMult(M1, M2);
        Matrix ProdStrass = StrassMult(M1, M2);
        cout << "Naive method:" << endl;
        ProdNaive.printDiagonal();
        cout << "\n Strassen's method:" << endl;
        ProdStrass.printDiagonal();
        return 0;
    }
    // if flag is '2' then calculate optimal n0
    else if (strtol(argv[1], nullptr, 0) == 2) {
        if (argc != 3) {
            std::printf("Usage: ./strassen 2 dimension\n");
            return -1;
        }
        
        clock_t start, end;
        double naive_time = 0;
        double strassens_time = INT_MAX;

        M1 = generateRandMat(dimension);
        M2 = generateRandMat(dimension);

        while (naive_time < strassens_time) {
            N_0 += 1;

            start = clock();
            NaiveMatMult(M1, M2);
            end = clock();
            naive_time = ((double)(end - start)) / (CLOCKS_PER_SEC / 1000);

            start = clock();
            StrassMult(M1, M2);
            end = clock();
            strassens_time = ((double)(end - start)) / (CLOCKS_PER_SEC / 1000);
        }
        cout << N_0 << endl;
        cout << "Naive finished in time " << naive_time << endl;
        cout << "Strassen's finished in time " << strassens_time << endl;

        return 0;
    }
    // if flag is '3' then calculate triangles
    else if (strtol(argv[1], nullptr, 0) == 3) {
        if (argc != 5) {
            std::printf("Usage: ./strassen 3 vertices probability trials\n");
            return -1;
        }
        // Usage: ./strassen 3 nodes probability
        float prob = strtod(argv[3], nullptr);
        int vertices = strtol(argv[2], nullptr, 0);
        int trials = strtol(argv[4], nullptr, 0);

        // perform requested number of trials on graphs of specified vertices and edge probability
        int totalTriangles = 0;
        for (int i = 0; i < trials; i++) {
            Matrix MRand = generateRandGraph(vertices, prob);
            int triangles = CountTriangles(MRand);
            cout << "P: " << prob << " | Triangles: " << triangles << endl;
            totalTriangles += triangles;
        }

        cout << "P: " << prob << " | AVERAGE: " << (float)totalTriangles/trials << endl;
        // end program
        return 0;
    }
    // otherwise read inputs from text file
    else {
        if (argc != 4) {
            std::printf("Usage: ./strassen 0 dimension inputfile\n");
            return -1;
        }
        
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
                j = 0;
            }
        }
    }

    Matrix Prod = StrassMult(M1, M2);
    Prod.printDiagonal();

    return 0;
}