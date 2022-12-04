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
    // int dims[2] = {0, 0};
    // int dim_counts[2];
    // int periods[2] = {true, true};
    // int this_coord[2];
    // int neighbors[4] = {};
    int n = stoi(argv[1]);
    // int k = stoi(argv[2]);
    int sub_n;
    int serial_result;
    double start;
    MPI_Request req;
    MPI_Status stat;
    Matrix X;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // MPI_Dims_create(num_procs, 2, dims);
    // MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &cart_comm);

    // get new rank, cart coords, and amount of procs in each dim
    // MPI_Comm_rank(cart_comm, &this_rank);
    // MPI_Cart_coords(cart_comm, this_rank, 2, this_coord);
    // get_dim_counts(2, cart_comm, dim_counts);

    if (this_rank == 0) {
        // create initial randomized matrix of size n
        X = Matrix(n);
        X.fill_rand(1);
        X.print();
        start = MPI_Wtime();

        // algo
        for (int i = 0; i < n - 1; i++) {
            // get base row i for this iter
            float base_buf[n];
            X.get_row(i, base_buf);

            // cout << "-> ";
            // for (int j = 0; j < n; j++) {
            //     cout << base_buf[j] << " ";
            // }
            // cout << endl;

            // send
            for (int k = i + 1; k < n; k++) {
                if (X(i, k) == 0) {
                    continue;
                }

                int dest = (k - 1) % num_procs + 1;

                float cur_buf[n];
                X.get_row(k, cur_buf);

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

                // loop -> recv and update
            }
            cout << endl;

            // receive and update LU
            // for k = i + 1
        }
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
                    cout << "->";
                    if (this_rank == 3)
                        for (int j = 0; j < n; j++) {
                            cout << base_buf[j] << " ";
                        }
                        cout << endl;
                    }
                }

                // receicve cur row k
                if (recv_proc == this_rank) {
                    float cur_buf[n];
                    MPI_Recv(cur_buf, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &stat);

                    // calc and send back

                    // debug
                    if (this_rank == 3)
                        for (int j = 0; j < n; j++) {
                            cout << cur_buf[j] << " ";
                        }
                        cout << endl;
                    }
                }
            }
        }
    }

    MPI_Finalize();
}