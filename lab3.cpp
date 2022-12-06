#include <iostream>
#include <string>
#include <math.h>
#include <mpi.h>
#include "Matrix.h"
using namespace std;

// get # of processors in each dim - dims_arr is out param
void calc_row(int i, int n, float* base_row, float* cur_row) {
    cur_row[i] = cur_row[i] / base_row[i];
    for (int j = i + 1; j < n; j++) {
        cur_row[j] = cur_row[j] - base_row[j] * cur_row[i];
    }
}

// update L and A (which becomes U) with received row
void update_row(int i, int k, int n, float* row, Matrix& L, Matrix& A) {
    L(k, i) = row[i];
    row[i] = 0;
    for (int j = 0; j < n; j++) {
        A(k, j) = row[j];
    }
}

int main(int argc, char** argv) {    
    int this_rank;
    int num_procs;
    int n = stoi(argv[1]);
    int sub_n;
    double start;
    MPI_Status stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (this_rank == 0) {
        // create initial randomized matrix of size n
        cout << "n = " << n << endl;
        Matrix A(n);
        A.fill_rand(-1);
        A.print();
        float serial_result = A.determinant();
        start = MPI_Wtime();

        // algo
        Matrix L(n);
        int root_rows[n - 1] = {};
        for (int i = 0; i < n - 1; i++) {
            // get base row i for this iter
            float base_buf[n];
            A.get_row(i, base_buf);

            // distribute
            int root_ind = 0;
            for (int k = i + 1; k < n; k++) {
                if (A(i, k) == 0) {
                    continue;
                } 

                int dest = k % num_procs - 1;

                // if root's responsibility, store k
                if (dest == 0) {
                    root_rows[root_ind++] = k;
                    continue;
                }

                float cur_buf[n];
                A.get_row(k, cur_buf);

                // send row i if haven't sent already then current row k
                MPI_Request req1, req2;
                if (dest < num_procs) {
                    MPI_Isend(base_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &req1);
                }
                MPI_Isend(cur_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &req2);
                
                // ensure buffers copied before they go out of scope
                MPI_Wait(&req1, &stat);
                MPI_Wait(&req2, &stat);
            }

            // root calculations
            for (int j = 0; j < root_ind; j++) {
                int k = root_rows[j];
                float cur_buf[n];
                A.get_row(k, cur_buf);
                calc_row(i, n, base_buf, cur_buf);
                update_row(i, k, n, cur_buf, L, A);
            }

            // loop -> recv and update
            for (int k = i + 1; k < n; k++) {
                if (A(i, k) == 0) {
                    continue;
                }

                int dest = k % num_procs - 1;
                
                // receive modified row
                float cur_buf[n];
                MPI_Recv(cur_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &stat);

                // update L with multiplier stored at row[i].
                // Then set to 0 in row and add row to U (A becomes U)
                update_row(i, k, n, cur_buf, L, A);
                
                L.print();
                A.print();
            }
        }

        // add 1s to diagonal of L
        for (int i = 0; i < n; i++) {
            L(i, i) = 1;
        }

        L.print();
        A.print();
        float det_L = L.determinant();
        float det_U = A.determinant();
        cout << "Serial Result: " << serial_result << endl;
        cout << "Parallel Result: " << det_L * det_U << endl;
        cout << "Determinant of L: " << det_L << endl;
        cout << "Determinant of U: " << det_U << endl;
        cout << "Parallel runtime: " << MPI_Wtime() - start << endl;

    } else {
        // child logic
        for (int i = 0; i < n - 1; i++) {
            float base_buf[n];
            for (int k = i + 1; k < n; k++) {
                // dont't proceed unless this proc is needed
                int recv_proc = k % num_procs - 1;
                if (recv_proc != this_rank) {
                    continue;
                }
                            
                // receive base row for iteration i if haven't already
                if (recv_proc < num_procs) {
                    MPI_Recv(base_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stat);
                }

                // receicve cur row k
                if (recv_proc == this_rank) {
                    float cur_buf[n];
                    MPI_Recv(cur_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stat);

                    // calculate multiplier, subtract row, and send back
                    calc_row(i, n, base_buf, cur_buf);
                    
                    MPI_Request req;
                    MPI_Isend(cur_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &req);
                    MPI_Wait(&req, &stat);
                }
            }
        }
    }

    MPI_Finalize();
}