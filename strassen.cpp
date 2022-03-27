#include <vector>
#include <iostream>
#include <assert.h>
#include <cmath>

using namespace std;


// structure used to represent a matrix
struct Matrix {
    // how many rows in matrix (height)
    int rowDim;
    // how many columns in matrix (width)
    int colDim;

    // 2D vector containing values
    vector<vector<int>> Mat;

    // constructor function
    Matrix(int r, int c) {
         rowDim = r;
         colDim = c;
        Mat.resize(r);
        for (vector<int> row : Mat) {
            row.resize(c);
        }
    }

    vector<int> &operator[] (const int i) {return Mat[i];}

    // method to print matrix
    void printMat() {
        for (vector<int> row : Mat) {
            for (int entry : row) {
                cout << entry  << "  " << endl;
            }
            cout << endl;
        }
    }

    // add zero row to end
    void padRow() {
        vector<int> zeros(colDim, 0);
        Mat.push_back(zeros);
    }

    // add zero column to end
    void padCol() {
        // for each row, add a zero to the end
        for (int i = 0; i < rowDim; i++) {
            Mat[i].push_back(0);
        }
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


// split an (n x n) matrix into sub-matrices 
vector<Matrix> split(Matrix M) {
    // argument MUST be a square matrix
    assert(M.colDim == M.rowDim);

    // create sub-matrices of correct dimension even if n is odd
    Matrix A(ceil(M.colDim), ceil(M.colDim));
    Matrix B(ceil(M.colDim), floor(M.colDim));
    Matrix C(floor(M.colDim), ceil(M.colDim));
    Matrix D(floor(M.colDim), floor(M.colDim));

    // copy values into matrices
    for (int i = 0; i < M.rowDim; i++) {
        for (int j = 0; j < M.colDim; j++) {
            // if in upper left add to A
            if (i < A.rowDim-1 && j < A.colDim-1) {
                A[i][j] = M[i][j];
            }
            // if in upper right add to B
            if (i < A.rowDim-1 && j >= A.colDim-1) {
                B[i][j - A.colDim] = M[i][j];
            }
            // if in lower left add to C
            if (i >= A.rowDim-1 && j < A.colDim-1) {
                C[i - A.rowDim][j] = M[i][j];
            }
            // if in lower right add to D
            if (i >= A.rowDim-1 && j >= A.colDim-1) {
                D[i - A.rowDim][j - A.colDim] = M[i][j];
            }
        }
    }

    return vector<Matrix> {A, B, C, D};
}

// This helper function performs the matrix multiplications for Strassen's Algorithm
Matrix StrassMult(Matrix M1, Matrix M2) {
    return M1;
}

// Strassen's matrix multiplication algorithm
Matrix Strassen(Matrix M1, Matrix M2) {
    return M1;
}


int main(void) {
    return 0;
}