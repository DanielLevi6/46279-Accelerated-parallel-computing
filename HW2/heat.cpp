/*
 * A parallel implementation of a stencil computation to solve the
 * 2-dimensional heat diffusion problem by Jacobi iteration using TBB.
 *
 * Copyright (c) 2013 Australian National University. All rights reserved.
 */

 /* Updated by Dr. Yariv Aridor, 2022 */

#include <iostream>
#include <algorithm>
#include <mutex>
#include <chrono>
#include <limits>
#include <tbb/tbb.h>

using namespace std;

    int main(int argc, char* argv[]) {
        size_t n_x, n_y;
        int max_iter;
        double t_edge, converge;
        size_t i, j;
        size_t num_threads;

        // configuration parameters 
        n_x = n_y = 2000;
        t_edge = 100.0;
        max_iter = 1000;
        converge = 0.01;
        int mode;

        // command-line parameters 
        if (argc < 3) {
            cout << "Argument missing! Usage: ./heat -mode <n>" << endl;
            exit(1);
        }

        if (strcmp(argv[1], "-mode") == 0) {
            mode = atoi(argv[2]);
        }
        else
            exit(1);

                
        cout << endl << "Stencil Computation" << " (mode=" << mode << ")" << endl;
        for (i = 0; i < 30; i++) cout << "-";
        cout << endl;

        // allocate data array
        double* t_old = new double[n_x * n_y]();
        double* t_new = new double[n_x * n_y](); 

        // fix boundary values
        j = 0;     for (i = 0; i < n_x; i++) t_new[j * n_x + i] = t_old[j * n_x + i] = t_edge;
        j = n_y - 1; for (i = 0; i < n_x; i++) t_new[j * n_x + i] = t_old[j * n_x + i] = t_edge;
        i = 0;     for (j = 0; j < n_y; j++) t_new[j * n_x + i] = t_old[j * n_x + i] = t_edge;
        i = n_x - 1; for (j = 0; j < n_y; j++) t_new[j * n_x + i] = t_old[j * n_x + i] = t_edge;


        int iter = 0;
        double max_diff;
        std::mutex m;

        auto start = chrono::steady_clock::now();

        while (iter < max_iter) {
            iter++;
            max_diff = 0.f;

            // sequential execution 
            if (mode == 0) {
                for (j = 1; j < n_y - 1; j++) {
                    for (i = 1; i < n_x - 1; i++) {
                        t_new[j * n_x + i] = 0.25 * (t_old[j * n_x + i + 1] + t_old[j * n_x + i - 1] +
                            t_old[(j + 1) * n_x + i] + t_old[(j - 1) * n_x + i]);
                        double tdiff = fabs(t_old[j * n_x + i] - t_new[j * n_x + i]);
                        max_diff = max(max_diff, tdiff);
                    }
                }
            }

            if (mode == 1) {
                // add your code for straightforward parallel_for implementation    
                tbb::parallel_for(tbb::blocked_range<size_t>(1, n_y - 1), 
                                        [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t j = r.begin(); j < r.end(); j++) {
                        for (int i = 1; i < n_x - 1; i++) {
                            t_new[j * n_x + i] = 0.25 * (t_old[j * n_x + i + 1] + t_old[j * n_x + i - 1] +
                                t_old[(j + 1) * n_x + i] + t_old[(j - 1) * n_x + i]);
                            double tdiff = fabs(t_old[j * n_x + i] - t_new[j * n_x + i]);
                            lock_guard<mutex> max_lock(m); 
                            max_diff = max(max_diff, tdiff);
                        }
                    }
                }); 
            }


            if (mode == 2) { 
                // add your code for optimized parallel_for implementation
                tbb::parallel_for(tbb::blocked_range<size_t>(1, n_y - 1), 
                                        [&](const tbb::blocked_range<size_t>& r) {
                    double local_max_diff=0.f;
                    for (size_t j = r.begin(); j < r.end(); j++) {
                        for (int i = 1; i < n_x - 1; i++) {
                            t_new[j * n_x + i] = 0.25 * (t_old[j * n_x + i + 1] + t_old[j * n_x + i - 1] +
                                t_old[(j + 1) * n_x + i] + t_old[(j - 1) * n_x + i]);
                            double tdiff = fabs(t_old[j * n_x + i] - t_new[j * n_x + i]);
                            local_max_diff = max(local_max_diff, tdiff);
                        }
                    }
                    lock_guard<mutex> max_lock(m); 
                    max_diff = max(local_max_diff, max_diff);
                }); 
            }


            if (mode == 3) {
                // add your code here for parallel_reduce implementation
                max_diff = tbb::parallel_reduce(
                                tbb::blocked_range2d<size_t, size_t>(1, n_x - 1, 1, n_y - 1),
                                0.0,
                                [&](const tbb::blocked_range2d<size_t, size_t> &r, double max_diff_local) -> double {
                                    for (size_t j = r.cols().begin(); j < r.cols().end(); j++) {
                                        for (size_t i = r.rows().begin(); i < r.rows().end(); i++) {
                                            t_new[j * n_x + i] = 0.25 * (t_old[j * n_x + i + 1] + t_old[j * n_x + i - 1] +
                                                t_old[(j + 1) * n_x + i] + t_old[(j - 1) * n_x + i]);
                                            double tdiff = fabs(t_old[j * n_x + i] - t_new[j * n_x + i]);
                                            max_diff_local = max(max_diff_local, tdiff);
                                        }
                                    }
                                    return max_diff_local;
                                },
                                /* reduction */
                                [](double x, double y) -> double {
                                    return max(x, y);
                                }
                );
            }

            // swap array pointers
            double* temp = t_new; t_new = t_old; t_old = temp;
            if (max_diff < converge) break;
        }

        cout << "iteration: " << iter << " convergence:" << max_diff << endl;
        auto end = chrono::steady_clock::now();
        cout << "Elapsed time in milliseconds: "
            << chrono::duration_cast<chrono::milliseconds>(end - start).count()
            << " ms" << endl;

        /* Printout the result */
        const size_t kMaxPrint = 5;
        cout << " final sample (20) grid values" << endl;
        for (j = 1; j < kMaxPrint; j++) {
            for (i = 1; i < kMaxPrint; i++) {
                cout << t_new[j * n_x + i];
                if (i % 5 == 4) cout << endl;
            }
        }
        cout << endl;

        delete[] t_old;
        delete[] t_new;

        return 0;
}


