#include <iostream>
#include <string>
#include <math.h>
#include <mpi.h>
#include "Matrix.h"
using namespace std;

// get # of processors in each dim - dims_arr is out param
void get_dim_counts(int m, MPI_Comm comm, int* dims_arr) {
    int dim_num;
    int periods[m];
    int coords[m];
    MPI_Cart_get(comm, m, dims_arr, periods, coords);
}

int main(int argc, char** argv) {    
    int this_rank;
    int num_procs;
    int dims[2] = {0, 0};
    int dim_counts[2];
    int periods[2] = {true, true};
    int this_coord[2];
    int neighbors[4] = {};
    int n = stoi(argv[1]);
    int k = stoi(argv[2]);
    int sub_n;
    int serial_result;
    double start;
    MPI_Comm cart_comm;
    MPI_Request req;
    MPI_Status stat;
    Matrix X;
    Matrix multA;
    Matrix multB;
    Matrix A;
    Matrix B;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    sub_n = n / sqrt(num_procs);

    MPI_Dims_create(num_procs, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &cart_comm);

    // get new rank, cart coords, and amount of procs in each dim
    MPI_Comm_rank(cart_comm, &this_rank);
    MPI_Cart_coords(cart_comm, this_rank, 2, this_coord);
    get_dim_counts(2, cart_comm, dim_counts);

    if (this_rank == 0) {
        // create initial randomized matrix of size n
        X = Matrix(n);
        X.fill_rand(1);

        cout << "n = " << n << endl;
        cout << "k = " << k << endl;

        if (num_procs == 1) {
            start = MPI_Wtime();
            Matrix result = X;
            for (int i = 0; i < k; i++) {
                result = result * X;
            }
            double serial_runtime = MPI_Wtime() - start;
            serial_result = result.determinant();
            cout << "Serial result: " << serial_result << endl;
            cout << "Serial runtime: " << serial_runtime << endl;
            cout << endl;

            return 0;
        }
    
        start = MPI_Wtime();
    }

    for (int i = 0; i < k; i++) {
        // Cannon's Algorithm
        if (this_rank == 0) {
            // If first iteration then current matrix is initial matrix.
            // Otherwise use last iteration result currently stored in multA.
            if (i == 0) {
                multA = X;
            }
            multB = X;

            // cout << "Disassembling A" << endl;
            Matrix partsA[num_procs] = {};
            int ind = 0;
            for (int i = 0; i < n; i+=sub_n) {
                for (int j = 0; j < n; j+=sub_n) {
                    partsA[ind++] = multA.get_subm(sub_n, i, j);
                }
            }
            // cout << "Disassembling B" << endl;
            Matrix partsB[num_procs] = {};
            ind = 0;
            for (int i = 0; i < n; i+=sub_n) {
                for (int j = 0; j < n; j+=sub_n) {
                    partsB[ind++] = multB.get_subm(sub_n, i, j);
                }
            }

            // cout << "Distributing submatrices" << endl;
            ind = 1;
            for (int i = 0; i < dims[0]; i++) {
                for (int j = 0; j < dims[1]; j++) {
                    if (i == 0 && j == 0) continue;
                    int targ_rank;
                    int coord[2] = {i, j};
                    MPI_Cart_rank(cart_comm, coord, &targ_rank);
                    MPI_Send(partsA[ind].get_1d(), sub_n*sub_n, MPI_INT, targ_rank, 0, cart_comm);
                    MPI_Send(partsB[ind].get_1d(), sub_n*sub_n, MPI_INT, targ_rank, 0, cart_comm);
                    ind++;
                }
            }
            
            // root doesn't ned to send/recv to itself
            A = partsA[0];
            B = partsB[0];
        } else {
            // cout << "Other processes receiving" << endl;
            int buf[sub_n * sub_n];
            MPI_Recv(buf, sub_n * sub_n, MPI_INT, 0, 0, cart_comm, &stat);
            A = Matrix(buf, sub_n);
            MPI_Recv(buf, sub_n * sub_n, MPI_INT, 0, 0, cart_comm, &stat);
            B = Matrix(buf, sub_n);
        }

        // MPI_Barrier(cart_comm);
        // cout << "past barrier 0" << endl;

        // Initial Alignment
        int A_src;
        int B_src;
        int A_dest;
        int B_dest;
        int buf[sub_n * sub_n];
        MPI_Cart_shift(cart_comm, 1, -this_coord[0], &A_src, &A_dest);
        MPI_Cart_shift(cart_comm, 0, -this_coord[1], &B_src, &B_dest);
        if (this_coord[0] != 0) {
            MPI_Sendrecv(A.get_1d(), sub_n * sub_n, MPI_INT, A_dest, 0, buf, sub_n * sub_n, MPI_INT, A_src, 0, cart_comm, &stat);
            A = Matrix(buf, sub_n);
        }
        if (this_coord[1] != 0) {
            MPI_Sendrecv(B.get_1d(), sub_n * sub_n, MPI_INT, B_dest, 0, buf, sub_n * sub_n, MPI_INT, B_src, 0, cart_comm, &stat);
            B = Matrix(buf, sub_n);
        }

        MPI_Barrier(cart_comm);

        // Calculate and Shift loop
        Matrix sum(sub_n);
        sum = sum + (A * B);
        for (int i = 1; i < dims[0]; i++) {
            MPI_Cart_shift(cart_comm, 1, -1, &A_src, &A_dest);
            MPI_Cart_shift(cart_comm, 0, -1, &B_src, &B_dest);
            MPI_Sendrecv(A.get_1d(), sub_n * sub_n, MPI_INT, A_dest, 0, buf, sub_n * sub_n, MPI_INT, A_src, 0, cart_comm, &stat);
            A = Matrix(buf, sub_n);
            MPI_Sendrecv(B.get_1d(), sub_n * sub_n, MPI_INT, B_dest, 0, buf, sub_n * sub_n, MPI_INT, B_src, 0, cart_comm, &stat);
            B = Matrix(buf, sub_n);
            sum = sum + (A * B);
        }

        MPI_Barrier(cart_comm);

        // Collect submatrices at root and assemble matrix
        if (this_rank != 0) {
            // Send to root
            MPI_Send(sum.get_1d(), sub_n * sub_n, MPI_INT, 0, 0, cart_comm);
        } else {
            // Receive and store
            Matrix parts[num_procs] = {};
            parts[0] = sum;
            int ind = 0;
            for (int i = 1; i < num_procs; i++) {
                int buf[sub_n * sub_n];
                MPI_Recv(buf, sub_n * sub_n, MPI_INT, i, 0, cart_comm, &stat);
                int coord[2];
                MPI_Cart_coords(cart_comm, i, 2, coord);
                parts[coord[1] + coord[0] * dims[0]] = Matrix(buf, sub_n);
            }
            
            // Assemble
            ind = 0;
            for (int i = 0; i < n; i+=sub_n) {
                for (int j = 0; j < n; j+=sub_n) {
                    multA.add_subm(parts[ind++], sub_n, i, j);
                }
            }
            
            // Result
            if (i == k - 1) {
                double parallel_runtime = MPI_Wtime() - start;
                cout << "Parallel result: " << multA.determinant() << endl; 
                cout << "Parallel runtime: " << parallel_runtime << endl;
                cout << endl;
            }
        }

        MPI_Barrier(cart_comm);
    }

    MPI_Finalize();
}