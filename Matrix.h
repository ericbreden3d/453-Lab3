class Matrix {
    int** matrix;
    int* arr;  // also maintain 1-D array for sending and receiving
    int size;
    Matrix get_detrm_subm(int split);   // get sub-matrix without split column for determinant
    int detrm_helper(Matrix& m);   // for determinant()
 public:
    Matrix();
    Matrix(int n);
    Matrix(int* buf, int size);   // build matrix from 1D array
    Matrix(const Matrix& other);
    ~Matrix();
    void fill_rand(int seed);   // fill matrix with random values
    Matrix get_subm(int len, int x, int y);   // return sub-matrix of length len at location x,y
    void add_subm(Matrix& sub, int len, int x, int y);   // add sub-matrix to matrix at location x,y 
    int* get_1d();  // return 1D array representation of matrix
    int determinant();  // use detrm_helper to recursvely find determinant of matrix
    Matrix operator*(const Matrix& other);   // basic O(n^3) matrix multiplication algorithm
    Matrix operator+(const Matrix& other);
    Matrix& operator=(const Matrix& other);
    int& operator()(int x, int y);
    void print();
};

