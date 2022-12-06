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
        // A.print();
        Matrix U = A;
        start = MPI_Wtime();

        // algo
        Matrix L(n);
        for (int i = 0; i < n - 1; i++) {
            // get base row i and broadcast
            float base_buf[n];
            A.get_row(i, base_buf);
            MPI_Bcast(base_buf, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

            // if (base_buf[i] == 0) {
            //     cout << "BAD ZERO" << endl;
            // }

            // distribute
            int base_sent[num_procs] = {};
            int root_rows[n - 1] = {};
            int root_ind = 0;
            int child_rows[n - 1] = {};
            int child_ind = 0;
            float child_data[n - 1][n] = {};
            MPI_Request reqs[n - 1] = {};
            for (int k = i + 1; k < n; k++) {
                int dest = (k - 1) % num_procs;

                // if root's responsibility, store k
                if (dest == 0) {
                    root_rows[root_ind++] = k;
                    continue;
                } else {
                    child_rows[child_ind++] = k;
                }

                float cur_buf[n];
                A.get_row(k, cur_buf);

                // send row i if haven't sent already then current row k
                // if (!base_sent[dest]) {
                //     MPI_Isend(base_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &req1);
                //     base_sent[dest] = 1;
                // }
                MPI_Request req;
                MPI_Isend(cur_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &stat);
                
                // ensure buffers copied before they go out of scope
                // MPI_Wait(&req1, &stat);
            }

            // async recv
            for (int j = 0; j < child_ind; j++) {
                int k = child_rows[child_ind];
                int src = (k - 1) % num_procs;

                // first elem of row sent to child already 0, so it will not respond
                if (A(k, i) != 0) {
                    MPI_Irecv(child_data[k], n, MPI_FLOAT, src, 0, MPI_COMM_WORLD, &reqs[k]);
                }
            }

            // update with child data
            for (int j = 0; j < child_ind; j++) {
                int k = child_rows[child_ind];
                if (A(k, i) != 0) {
                    MPI_Wait(&reqs[k], &stat);
                    update_row(i, k, n, child_data[k], L, A);
                }
            }

            // update (after calc) with root data
            for (int j = 0; j < root_ind; j++) {
                int k = root_rows[j];
                float cur_buf[n];
                A.get_row(k, cur_buf);
                if (cur_buf[i] == 0) {
                    continue;
                }
                calc_row(i, n, base_buf, cur_buf);
                update_row(i, k, n, cur_buf, L, A);
            }
        }

        // add 1s to diagonal of L
        // for (int i = 0; i < n; i++) {
        //     L(i, i) = 1;
        // }

        // L.print();
        // A.print();

        // float serial_result = U.determinant();
        // float parallel_result = A.determinant();
        // cout << "Serial Result: " << serial_result << endl;
        // cout << "Parallel Result: " << parallel_result << endl;
        cout << "Parallel runtime: " << MPI_Wtime() - start << endl;

    } else {
        // child logic
        for (int i = 0; i < n - 1; i++) {
            float base_buf[n];
            int base_recvd = 0;

            MPI_Bcast(base_buf, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

            for (int k = i + 1; k < n; k++) {
                // dont't proceed unless this proc is needed
                int recv_proc = (k - 1) % num_procs;
                if (recv_proc != this_rank) {
                    continue;
                }
                            
                // receive base row for iteration i if haven't already
                // if (!base_recvd) {
                //     // cout << "child waiting " << this_rank << " on base " << i << " for " << i << ", " << k << endl;
                //     MPI_Recv(base_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stat);
                //     // cout << "child " << this_rank << " received " << k << endl;
                //     base_recvd = 1;
                // }

                // receicve cur row k
                float cur_buf[n];
                // cout << "child " << this_rank << " waiting " << i << ", " << k << endl;
                MPI_Recv(cur_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stat);
                // cout << "child " << this_rank << " received " << i << ", " << k << endl;

                // received already zeroed row, ignore
                if (cur_buf[i] == 0) {
                    cout << "child " << this_rank << ": already zeroed at " << i << ", " << k << endl;
                    continue;
                }

                // calculate multiplier, subtract row, and send back
                calc_row(i, n, base_buf, cur_buf);
                
                MPI_Request req;
                MPI_Isend(cur_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &stat);
            }
        }
    }

    MPI_Finalize();
}