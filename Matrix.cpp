#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <math.h>
#include "Matrix.h"
using namespace std;

Matrix Matrix::get_detrm_subm(int split) {
    // cout << "Getting sub-matrix" << endl;
    Matrix new_m(size - 1);
    for (int i = 0; i < this->size; i++) {
        if (i == split) {
            split = -1;
            continue;
        }
        int x = split == -1 ? i - 1 : i;  // if already split, shift left
        for (int j = 1; j < this->size; j++) {
            new_m(x, j - 1) = (*this)(i, j);
        }
    }
    // new_m.print();
    return new_m;
}

int Matrix::detrm_helper(Matrix& m) {
    if (m.size == 1) {
        return m(0, 0);
    }
    int op = 1;
    int detrm = 0;
    for (int i = 0; i < m.size; i++) {
        // get sub-matrix without column i
        Matrix sub = m.get_detrm_subm(i);
        // recursively find determinant of sub-matrix
        int rec = detrm_helper(sub);
        // get product
        int result = m(i, 0) * rec;
        // add to or subtract from current total
        if (op) {
            detrm += result;
            op--;
        } else {
            detrm -= result;
            op++;
        }
    }
    return detrm;
}

Matrix::Matrix() {
    // cout << "Default Constructor Called" << endl;
    size = 0;
}

Matrix::Matrix(int n) {
    // cout << "Creating matrix" << endl;
    size = n;
    matrix = new int*[size];
    arr = new int[size * size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new int[size];
        for (int j = 0; j < size; j++) {
            matrix[i][j] = 0;
            arr[i + j * size] = 0;
        }
    }
}

Matrix::Matrix(int* buf, int size) {
    this->size = size;
    int len = size * size;
    matrix = new int*[size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new int[size];
    }
    arr = new int[len];
    int i = 0;
    int j = 0;
    for (int ind = 0; ind < len; ind++) {
        // cout << i << " " << j << ind << " " << buff[ind] << endl;
        matrix[i][j] = buf[ind];
        arr[i + j * size] = buf[ind];
        if (i == size - 1) {
            j++;
            i = 0;
        } else {
            i++;
        }
    }
}

Matrix::Matrix(const Matrix& other) {
    // cout << "Copy Constructor Called" << endl;
    size = other.size;
    matrix = new int*[size];
    arr = new int[size * size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new int[size];
        for (int j = 0; j < size; j++) {
            matrix[i][j] = other.matrix[i][j];
            arr[i + j * size] = other.matrix[i][j];
        }
    }
}

Matrix::~Matrix() {
    // cout << "Deleting matrix" << endl;
    // cout << "SIze " << size << endl;
    // print();

    if (size == 0) return;
    for (int i = 0; i < size; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] arr;
}

void Matrix::fill_rand(int seed = -1) {
    // cout << "Filling with random values" << endl;
    srand(seed == -1 ? time(NULL) : seed);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int r = 3 - rand() % 4;
            matrix[i][j] = r;
            arr[i + j * size] = r;
        }
    }
}

Matrix Matrix::operator*(const Matrix& other) {
    Matrix new_m(size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int sum = 0;
            for (int k = 0; k < size; k++) {
                sum += (*this)(i, k) * other.matrix[k][j];
            }
            new_m(i, j) = sum;
        }
    }
    return new_m;
}

Matrix Matrix::operator+(const Matrix& other) {
    Matrix new_m(size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            new_m(i, j) = (*this)(i, j) + other.matrix[i][j];
        }
    }
    return new_m;
}

Matrix& Matrix::operator=(const Matrix& other) {
    // cout << "Copy Operator Called" << endl;
    if (size != 0) {
        for (int i = 0; i < size; i++) {
            delete[] matrix[i];
        }
        delete[] matrix;
        delete[] arr;
    }

    size = other.size;
    matrix = new int*[size];
    arr = new int[size * size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new int[size];
        for (int j = 0; j < size; j++) {
            matrix[i][j] = other.matrix[i][j];
            arr[i + j * size] = other.matrix[i][j];
        }
    }
    return *this;
}

Matrix Matrix::get_subm(int len, int x, int y) {
    Matrix new_m(len);
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            new_m(i, j) = (*this)(i + x, j + y);
            new_m.arr[i + j * len] = (*this)(i + x, j + y);
        }
    }
    return new_m;
}

void Matrix::add_subm(Matrix& sub, int len, int x, int y) {
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            (*this)(i + x, j + y) = sub(i, j);
            this->arr[(i + x) + (j + y) * len] = sub(i, j);
        }
    }
}

int* Matrix::get_1d() {
    return arr;
}

int Matrix::determinant() {
    int result = detrm_helper(*this);
    return result;
}

int& Matrix::operator()(int x, int y) {
    return matrix[x][y];
}

void Matrix::print() {
    // cout << "Printing matrix" << endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            cout << (*this)(i, j) << " ";
        }
        cout << endl;
    }
    cout << endl;
}