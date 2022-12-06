#include <iostream>
#include <string>
#include <math.h>
#include <mpi.h>
#include "Matrix.h"
using namespace std;

// calculate and store multiplier at index i, then subtract rows
void calc_row(int i, int n, float* base_row, float* cur_row) {
    cur_row[i] = cur_row[i] / base_row[i];
    for (int j = i + 1; j < n; j++) {
        cur_row[j] = cur_row[j] - base_row[j] * cur_row[i];
    }
}

// update L and A (which becomes U) with received row.
// Also set index i to 0 to get echelon form
void update_row(int i, int k, int n, float* row, Matrix& U) {
    // L(k, i) = row[i];
    row[i] = 0;
    for (int j = 0; j < n; j++) {
        U(k, j) = row[j];
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
        // create initial randomized matrix of size n * n
        cout << "n = " << n << endl;
        Matrix A(n);
        A.fill_rand(-1);
        Matrix U = A;  // original A is modified into U
        U.print();
        // Matrix L(n);   // matrix L not needed for finding determinant

        MPI_Barrier(MPI_COMM_WORLD);

        // start timer
        start = MPI_Wtime();

        // algo
        for (int i = 0; i < n - 1; i++) {
            // get base row i and broadcast
            float base_buf[n];
            U.get_row(i, base_buf);
            MPI_Bcast(base_buf, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

            // distribute responsibilities to children
            int root_rows[n - 1] = {};        // k vals that root is responsible for
            int root_ind = 0;                 // cur index of root_rows
            int child_rows[n - 1] = {};       // k vals that children are responsible for
            int child_ind = 0;                // cur index of child_rows
            float child_data[n - 1][n] = {};  // row data for each val of k to allow async recvs
            MPI_Request reqs[n - 1] = {};     // array of request objs to allow for async recvs
            for (int k = i + 1; k < n; k++) {
                // store k val in either root or child arr
                int dest = (k - 1) % num_procs;
                if (dest == 0) {
                    root_rows[root_ind++] = k;
                    continue;  // if root, then don't need to send
                } else {
                    child_rows[child_ind++] = k;
                }

                // get row k and send if child's responsibility
                float cur_buf[n];
                U.get_row(k, cur_buf);
                MPI_Request req;
                MPI_Isend(cur_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &stat);
            }

            // async recv using child_rows gathered in last block
            for (int j = 0; j < child_ind; j++) {
                int k = child_rows[j];
                int src = (k - 1) % num_procs;

                // receive response from child with subtracted row and multiplier at i.
                // if ith elem was already 0, child is set to not respond
                if (U(k, i) != 0) {
                    MPI_Irecv(child_data[k], n, MPI_FLOAT, src, 0, MPI_COMM_WORLD, &reqs[k]);
                }
            }

            // update U with child data
            for (int j = 0; j < child_ind; j++) {
                int k = child_rows[j];
                if (U(k, i) != 0) {
                    MPI_Wait(&reqs[k], &stat);   // ensure data received
                    update_row(i, k, n, child_data[k], U);
                }
            }

            // calculate root's rows and update U
            for (int j = 0; j < root_ind; j++) {
                int k = root_rows[j];
                float cur_buf[n];
                U.get_row(k, cur_buf);
                if (cur_buf[i] == 0) {
                    continue;
                }
                calc_row(i, n, base_buf, cur_buf);
                update_row(i, k, n, cur_buf, U);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        //
        // Results
        //
        // float serial_result = A.determinant();
        // float parallel_result = U.determinant();
        // cout << "Serial Result: " << serial_result << endl;
        // cout << "Parallel Result: " << parallel_result << endl;
        cout << "Parallel runtime: " << MPI_Wtime() - start << endl;

    } else {

        MPI_Barrier(MPI_COMM_WORLD);

        // child logic
        for (int i = 0; i < n - 1; i++) {

            // debug
            cout << "ITER: " << i << endl;

            // receive broadcast base row for this iter i
            float base_buf[n];
            MPI_Bcast(base_buf, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

            int my_rows[n - 1] = {};
            float my_data[n - 1][n] = {};
            int my_ind = 0;
            MPI_Request reqs[n - 1] = {};
            for (int k = i + 1; k < n; k++) {
                // dont't proceed unless this proc is needed
                int recv_proc = (k - 1) % num_procs;
                if (recv_proc == this_rank) {
                    my_rows[my_ind++] = k;
                } else {
                    continue;
                }

                // async receicve cur row k from root
                float cur_buf[n];
                MPI_Recv(cur_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stat);

                if (cur_buf[i] == 0) {
                    continue;
                }
                
                calc_row(i, n, base_buf, cur_buf);
                
                // send back to root
                MPI_Request req;
                MPI_Isend(my_data[k], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &stat);
            }

            // for (int j = 0; j < my_ind; j++) {
            //     int k = my_rows[j];

            //     // calc multiplier and store at ind i ->  Rk[i] = Rk[i] / Rb[i], 
            //     // subtract multiplied base from row k -> Rk - Rb*multiplier
            //     // if row already zeroed out, ignore
            //     MPI_Wait(&reqs[k], &stat);  // ensure received
            //     if (my_data[k][i] == 0) {
            //         continue;
            //     }
            //     calc_row(i, n, base_buf, my_data[k]);
                
            //     // send back to root
            //     MPI_Request req;
            //     MPI_Isend(my_data[k], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &req);
            //     MPI_Wait(&req, &stat);
            // }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
}