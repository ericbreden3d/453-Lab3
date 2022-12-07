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
        // Matrix U = A;
        // A.print();
        start = MPI_Wtime();

        // algo
        Matrix L(n);
        for (int i = 0; i < n - 1; i++) {
            
            // get base row i for this iter
            float base_buf[n];
            A.get_row(i, base_buf);
            MPI_Bcast(base_buf, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

            // distribute
            int root_rows[n - 1] = {};
            int root_ind = 0;
            int child_rows[n - 1];
            int child_ind = 0;
            float child_data[n][n];
            for (int k = i + 1; k < n; k++) {

                int dest = (k - 1) % num_procs;

                // if root's responsibility, store k
                if (dest == 0) {
                    root_rows[root_ind++] = k;
                    continue;
                }

                float cur_buf[n];
                A.get_row(k, cur_buf);

                // send row i if haven't sent already then current row k
                MPI_Request req;
                MPI_Send(cur_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
                
                // ensure buffers copied before they go out of scope
                // MPI_Wait(&req, &stat);
            }

            // root calculations (maybe move above receives)
            for (int j = 0; j < root_ind; j++) {
                int k = root_rows[j];
                float cur_buf[n];
                A.get_row(k, cur_buf);
                if (cur_buf[i] != 0) {
                    calc_row(i, n, base_buf, cur_buf);
                    update_row(i, k, n, cur_buf, L, A);
                }
            }

            // loop -> recv and update
            for (int k = i + 1; k < n; k++) {
                int dest = (k - 1) % num_procs;
                if (A(k, i) == 0 || dest == 0) {
                    continue;
                }
                
                // receive modified row
                float cur_buf[n];
                // cout << "root waiting on " << dest << endl;
                MPI_Recv(cur_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &stat);
                // cout << "root received " << dest << endl;

                // update L with multiplier stored at row[i].
                // Then set to 0 in row and add row to U (A becomes U)
                update_row(i, k, n, cur_buf, L, A);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }

        // float serial_result = U.determinant();
        // float parallel_result = A.determinant();
        // cout << "Serial Result: " << serial_result << endl;
        // cout << "Parallel Result: " << parallel_result << endl;
        cout << "Parallel runtime: " << MPI_Wtime() - start << endl;

    } else {
        // child logic
        for (int i = 0; i < n - 1; i++) {
            float base_buf[n];
            MPI_Bcast(base_buf, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

            int my_rows[n - 1];
            int my_ind = 0;
            float my_data[n][n];
            MPI_Request reqs[n];
            for (int k = i + 1; k < n; k++) {
                // dont't proceed unless this proc is needed
                int recv_proc = (k - 1) % num_procs;
                if (recv_proc == this_rank) {
                    my_rows[my_ind++] = k;
                } else {
                    continue;
                }

                // receicve cur row k
                // float cur_buf[n];
                // cout << "child " << this_rank << " waiting " << i << ",s " << k << endl;
                // cout << "here" << endl;
                // reqs[k] = MPI_Request();
                MPI_Irecv(my_data[k], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &reqs[k]);
                MPI_Wait(&reqs[k], &stat);

                // cout << "child " << this_rank << " received " << i << ", " << k << endl;

                // received already zeroed row, ignore
                if (my_data[k][i] == 0) {
                    continue;
                }

                // calculate multiplier, subtract row, and send back
                calc_row(i, n, base_buf, my_data[k]);
                
                MPI_Request req;
                MPI_Send(my_data[k], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
                // MPI_Wait(&req, ?&stat);

            }

            // for (int j = 0; j < my_ind; j++) {
            //     int k = my_rows[j];
            //     cout << "waiting"<< endl;
            //     MPI_Wait(&reqs[k], &stat);
            //     cout << "waited" << endl;

            //     // received already zeroed row, ignore
            //     if (my_data[k][i] == 0) {
            //         continue;
            //     }

            //     // calculate multiplier, subtract row, and send back
            //     calc_row(i, n, base_buf, my_data[k]);
                
            //     MPI_Request req;
            //     MPI_Isend(my_data[k], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &req);
            //     MPI_Wait(&req, &stat);
            // }

            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
}