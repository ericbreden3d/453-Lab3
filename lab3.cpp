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
void update_row(int i, int k, int n, float* row, Matrix& A) {
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

    // root logic
    if (this_rank == 0) {
        // create initial randomized matrix of size n
        cout << "n = " << n << endl;
        Matrix A(n);
        A.fill_rand(-1);
        Matrix U = A;

        // start timer
        start = MPI_Wtime();

        // Outermost loop: n - 1 rounds
        for (int i = 0; i < n - 1; i++) {
            
            // get base row i for this iter and broadcast
            float base_buf[n];
            U.get_row(i, base_buf);
            MPI_Bcast(base_buf, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

            // variables
            int root_rows[n - 1] = {};   // rows that root is responsible for in this iter
            int root_ind = 0;            // index for root_rows
            int child_rows[n - 1];       // rows that children are responsible for in this iter
            int child_ind = 0;           // index for child_rows
            float child_data[n][n];      // stores child data indexed by n possible values of k
            MPI_Request reqs[n];         // request structures to allow for async receive/wait calls

            // distribute rows for round + 1 to n
            for (int k = i + 1; k < n; k++) {
                // map k to processor
                int dest = (k - 1) % num_procs;

                // determine who is responsible for kth row and store val
                if (dest == 0) {
                    root_rows[root_ind++] = k;
                    continue;
                } else {
                    child_rows[child_ind++] = k;
                }

                // send kth row to child
                float cur_buf[n];
                U.get_row(k, cur_buf);
                MPI_Send(cur_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
            }

            // async recv using child_rows gathered in last block
            for (int j = 0; j < child_ind; j++) {
                // calculate source based ons tored k val
                int k = child_rows[j];
                int src = (k - 1) % num_procs;

                // async receive response from child with subtracted row and multiplier at i.
                // if ith elem was already 0, child is set to not respond
                if (U(k, i) != 0) {
                    MPI_Irecv(child_data[k], n, MPI_FLOAT, src, 0, MPI_COMM_WORLD, &reqs[k]);
                }
            }

            // update U based on child response
            for (int j = 0; j < child_ind; j++) {
                int k = child_rows[j];
                MPI_Wait(&reqs[k], &stat);
                update_row(i, k, n, child_data[k], U);
            }

            // calc root rows and update U
            for (int j = 0; j < root_ind; j++) {
                int k = root_rows[j];
                float cur_buf[n];
                A.get_row(k, cur_buf);
                if (cur_buf[i] != 0) {
                    calc_row(i, n, base_buf, cur_buf);
                    update_row(i, k, n, cur_buf, U);
                }
            }
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
        // child logic
        for (int i = 0; i < n - 1; i++) {
            // receive broadcast base row i for this iter
            float base_buf[n];
            MPI_Bcast(base_buf, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

            // variables
            int my_rows[n - 1];   // this child is responsible for these rows
            int my_ind = 0;       // index for my_rows
            float my_data[n][n];  // holds data indexable by k
            MPI_Request reqs[n];  // request structures for async receive/wait

            // receive from root if k values maps to this processor
            for (int k = i + 1; k < n; k++) {
                // store k val and calc after receiving all rows for this round
                int recv_proc = (k - 1) % num_procs;
                if (recv_proc == this_rank) {
                    my_rows[my_ind++] = k;
                } else {
                    continue;
                }

                // async receive row
                MPI_Irecv(my_data[k], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &reqs[k]);
            }

            // 
            for (int j = 0; j < my_ind; j++) {
                // wait for non-blocking recv to complete
                int k = my_rows[j];
                MPI_Wait(&reqs[k], &stat);

                // received already zeroed row, ignore
                if (my_data[k][i] == 0) {
                    continue;
                }

                // calculate multiplier, subtract row, and send back to root
                calc_row(i, n, base_buf, my_data[k]);
                MPI_Send(my_data[k], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    MPI_Finalize();
}