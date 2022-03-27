#include <vector>
#include <iostream>

using namespace std;

struct Matrix {
    int dimension;
    Matrix(int d) {
        dimension = d;
    }
};

// regular matrix multiplication algorithm, used after we hit base case
vector<vector<int>> NaiveMatMult(vector<vector<int>> M1, vector<vector<int>> M2) {
    return M1;
}

// This helper function performs the matrix multiplications for Strassen's Algorithm
vector<vector<int>> StrassMult(vector<vector<int>> M1, vector<vector<int>> M2) {
    return M1;
}

// Strassen's matrix multiplication algorithm
vector<vector<int>> Strassen(vector<vector<int>> M1, vector<vector<int>> M2) {
    return M1;
}


int main(void) {
    return 0;
}