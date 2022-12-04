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
    int n = stoi(argv[1]);
    int sub_n;
    double start;
    MPI_Status stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (this_rank == 0) {
        // create initial randomized matrix of size n
        Matrix A(n);
        A.fill_rand(1);
        A.print();
        float serial_result = A.determinant();
        start = MPI_Wtime();

        // algo
        Matrix L(n);
        for (int i = 0; i < n - 1; i++) {
            // get base row i for this iter
            float base_buf[n];
            A.get_row(i, base_buf);

            // cout << "-> ";
            // for (int j = 0; j < n; j++) {
            //     cout << base_buf[j] << " ";
            // }
            // cout << endl;

            // send
            for (int k = i + 1; k < n; k++) {
                if (A(i, k) == 0) {
                    continue;
                } 

                int dest = (k - 1) % (num_procs - 1) + 1;
                cout << "k: " << k << ", np: " << num_procs << ", dest: " << dest << endl;

                float cur_buf[n];
                A.get_row(k, cur_buf);

                // for (int j = 0; j < n; j++) {
                //     cout << cur_buf[j] << " ";
                // }
                // cout << endl;

                // send row i if haven't sent already then current row k
                MPI_Request req1, req2;
                if (k < num_procs) {
                    MPI_Isend(base_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &req1);
                }
                MPI_Isend(cur_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &req2);
                
                // ensure buffers copied before they go out of scope
                MPI_Wait(&req1, &stat);
                MPI_Wait(&req2, &stat);
            }
            cout << endl;

            // loop -> recv and update
            for (int k = i + 1; k < n; k++) {
                if (A(i, k) == 0) {
                    continue;
                }

                int dest = (k - 1) % num_procs + 1;
                
                // receive modified row
                float cur_buf[n];
                MPI_Recv(cur_buf, n, MPI_FLOAT, dest, 0, MPI_COMM_WORLD, &stat);

                // update L with multiplier stored at row[i].
                // Then set to 0 in row and add row to U (A becomes U)
                L(k, i) = cur_buf[i];
                cur_buf[i] = 0;
                for (int j = 0; j < n; j++) {
                    A(k, j) = cur_buf[j];
                }
                
                // if (dest == 2) {
                //     for (int j = 0; j < n; j++) {
                //         cout << cur_buf[j] << " ";
                //     }
                //     cout << endl;
                // }

                // L.print();
                // A.print();

                
            }
        }

        // add 1s to diagonal of L
        for (int i = 0; i < n; i++) {
            L(i, i) = 1;
        }

        L.print();
        A.print();
        int det_L = L.determinant();
        int det_U = A.determinant();
        cout << "Serial Result: " << serial_result << endl;
        cout << "Parallel Result: " << det_L * det_U << endl;
        cout << "Determinant of L: " << det_L << endl;
        cout << "Determinant of U: " << det_U << endl;

    } else {
        // child logic
        for (int i = 0; i < n - 1; i++) {
            float base_buf[n];
            for (int k = i + 1; k < n; k++) {
                // dont't proceed unless this proc is needed
                int recv_proc = (k - 1) % num_procs + 1;
                if (recv_proc != this_rank) {
                    continue;
                }
                            
                // receive base row for iteration i if haven't already
                if (k < num_procs) {
                    MPI_Recv(base_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stat);
                    // if (this_rank == 2) {
                    //     cout << "->";
                    //     for (int j = 0; j < n; j++) {
                    //         cout << base_buf[j] << " ";
                    //     }
                    //     cout << endl;
                    // }
                }

                // receicve cur row k
                if (recv_proc == this_rank) {
                    float cur_buf[n];
                    MPI_Recv(cur_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stat);

                    // calc and send back
                    cur_buf[i] = cur_buf[i] / base_buf[i];
                    for (int j = i + 1; j < n; j++) {
                        cur_buf[j] = cur_buf[j] - base_buf[j] * cur_buf[i];
                    }
                    
                    MPI_Request req;
                    MPI_Isend(cur_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &req);
                    MPI_Wait(&req, &stat);

                    // debug
                    // if (this_rank == 2) {
                    //     for (int j = 0; j < n; j++) {
                    //         cout << cur_buf[j] << " ";
                    //     }
                    //     cout << endl;
                    // }
                }
            }
        }
    }

    MPI_Finalize();
}