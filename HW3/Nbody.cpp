#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <tbb/tbb.h>
#include <random>

using namespace std;

#define VECTOR_LENGTH 256
#define BITS_PER_BYTE 8
#define FLOAT_SIZE sizeof(float)

const int N = 16384 * 2;
double performance[2] = { 1,1 };

/**************************************************/
/*
/*             serial implementation 
/*  
/* 
/**************************************************/

struct ParticleType {
    float x, y, z;
    float vx, vy, vz;
};

// We work on 256 bits vectors = 8 floats

struct AllParticlesParallel{
    float x[N];
    float y[N];
    float z[N];
    float vx[N];
    float vy[N];
    float vz[N];
};

ParticleType particles[N];
AllParticlesParallel all_Particles_Parallel;

void  init_particles_serial(ParticleType* particles, const int nParticles) {
    for (unsigned int i = 0; i < nParticles; i++) {
        particles[i].x = (float)(i % 15);
        particles[i].y = (float)((i * i) % 15);
        particles[i].z = (float)((i * i * i) % 15);
        particles[i].vx = 0;
        particles[i].vy = 0;
        particles[i].vz = 0;
    }
}


void move_particles_serial(const int nParticles, ParticleType* const particles, const float dt) {

    // Loop over particles that experience force
    for (int i = 0; i < nParticles; i++) {

        // Components of the gravity force on particle i
        float Fx = 0, Fy = 0, Fz = 0;

        // Loop over particles that exert force: vectorization expected here
        for (int j = 0; j < nParticles; j++) {

            // Avoid singularity and interaction with self
            const float softening = 1e-20f;

            // Newton's law of universal gravity
            const float dx = particles[j].x - particles[i].x;
            const float dy = particles[j].y - particles[i].y;
            const float dz = particles[j].z - particles[i].z;

            const float rr1 = 1.0f / sqrt(dx * dx + dy * dy + dz * dz + softening);
            const float drPowerN32 = rr1 * rr1 * rr1;
            //Calculate the net force
            Fx += dx * drPowerN32;
            Fy += dy * drPowerN32;
            Fz += dz * drPowerN32;
        }
        // if(i%8000 == 0) printf("Fx = %f\n", Fx);

        // Accelerate particles in response to the gravitational force
        particles[i].vx += dt * Fx;
        particles[i].vy += dt * Fy;
        particles[i].vz += dt * Fz;
    }

    // Move particles according to their velocities
    // O(N) work, so using a serial loop
    for (int i = 0; i < nParticles; i++) {
        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;
    }
}

/**************************************************/
/*
/*             parallel implementation 
/*  
/* 
/**************************************************/

int nthreads;

// TODO: add missing function 
void init_particles_parallel(AllParticlesParallel* particles_parallel, const int nParticles) {
    for (unsigned int i = 0; i < nParticles; i++) {
        particles_parallel->x[i] = (float)(i % 15);
        particles_parallel->y[i] = (float)((i * i) % 15);
        particles_parallel->z[i] = (float)((i * i * i) % 15);
        particles_parallel->vx[i] = 0;
        particles_parallel->vy[i] = 0;
        particles_parallel->vz[i] = 0;
    }
}


float accumulate_vector(__m256 vec) {
    float tmp[8];
    _mm256_store_ps(tmp, vec);
    return tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
}


// TODO: add missing function
void move_particles_parallel(const int nParticles, AllParticlesParallel* const particles_parallel, const float dt) {
    int number_of_elements_per_op = (VECTOR_LENGTH/BITS_PER_BYTE)/FLOAT_SIZE;
    int loops = ceil(nParticles / number_of_elements_per_op);
    
    __m256 softening = _mm256_set1_ps(1e-20f);
    // Loop over particles that experience force
    tbb::parallel_for(tbb::blocked_range2d<int, int>(0, loops, 0, number_of_elements_per_op),
                            [&](const tbb::blocked_range2d<int, int>& r) {
                                for (int i = r.rows().begin(); i < r.rows().end(); i++) {
                                    for (int l = r.cols().begin(); l < r.cols().end(); l++) {
                                    // Components of the gravity force on particle i
                                        __m256 Fx_vec = _mm256_set1_ps(0);
                                        __m256 Fy_vec = _mm256_set1_ps(0);
                                        __m256 Fz_vec = _mm256_set1_ps(0);
                                        __m256 x_vec = _mm256_set1_ps(particles_parallel->x[i * number_of_elements_per_op + l]);
                                        __m256 y_vec = _mm256_set1_ps(particles_parallel->y[i * number_of_elements_per_op + l]);
                                        __m256 z_vec = _mm256_set1_ps(particles_parallel->z[i * number_of_elements_per_op + l]);
                                    // Loop over particles that exert force: vectorization expected here
                                        for (int j = 0; j < loops; j++) {
                                            // Newton's law of universal gravity
                                            __m256 dx = _mm256_sub_ps(_mm256_load_ps(&(particles_parallel->x[j*number_of_elements_per_op])), x_vec);
                                            __m256 dy = _mm256_sub_ps(_mm256_load_ps(&(particles_parallel->y[j*number_of_elements_per_op])), y_vec);
                                            __m256 dz = _mm256_sub_ps(_mm256_load_ps(&(particles_parallel->z[j*number_of_elements_per_op])), z_vec);

                                            __m256 rr1 = _mm256_div_ps(_mm256_set1_ps(1.0f), _mm256_sqrt_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(dx, dx), _mm256_mul_ps(dy, dy)), _mm256_add_ps(_mm256_mul_ps(dz, dz), softening))));
                                            __m256 drPowerN32 = _mm256_mul_ps(_mm256_mul_ps(rr1, rr1), rr1);
                                            
                                            //Calculate the net force
                                            Fx_vec = _mm256_add_ps(Fx_vec, _mm256_mul_ps(dx, drPowerN32));
                                            Fy_vec = _mm256_add_ps(Fy_vec, _mm256_mul_ps(dy, drPowerN32));
                                            Fz_vec = _mm256_add_ps(Fz_vec, _mm256_mul_ps(dz, drPowerN32));
                                        }
                                        float Fx = accumulate_vector(Fx_vec);
                                        float Fy = accumulate_vector(Fy_vec);
                                        float Fz = accumulate_vector(Fz_vec);
                                        particles_parallel->vx[i * number_of_elements_per_op + l] += dt * Fx;
                                        particles_parallel->vy[i * number_of_elements_per_op + l] += dt * Fy;
                                        particles_parallel->vz[i * number_of_elements_per_op + l] += dt * Fz;
                                    }
                                }
                            }
    );

    // Second loop
    tbb::parallel_for(tbb::blocked_range<size_t>(0, 3 * loops),
                    [=](const tbb::blocked_range<size_t>& r) {
                            for(int i = r.begin(); i!=r.end(); i++) {
                                __m256 xyz = _mm256_load_ps(&(particles_parallel->x[i*number_of_elements_per_op]));
                                __m256 v = _mm256_load_ps(&(particles_parallel->vx[i*number_of_elements_per_op]));
                                xyz = _mm256_add_ps(xyz, _mm256_mul_ps(v, _mm256_set1_ps(dt)));
                                _mm256_store_ps(&(particles_parallel->x[i*number_of_elements_per_op]), xyz);
                            }
                    }
    );
}

/**************************************************/
/*
/*             run simulations                    
/*
/**************************************************/

bool validate_results(float f1, float f2, float f3) {
    const float epsilon = 0.003f;

    std::cout << "f1 = " << f1 << "(2.351)" << endl;
    std::cout << "f2 = " << f2 << "(2.435)" << endl;
    std::cout << "f3 = " << f3 << "(9.546)" << endl;

    if ((fabs(f1 - 2.351) > epsilon) || (fabs(f2 - 2.435) > epsilon) || (fabs(f3 - 9.546) > epsilon)) {
        return false;
    }
    return true;
}

bool run_simulation(int mode) {  // mode 0 - serial ; mode 1 - parallel
    
    const int nParticles = N; // Problem size and other parameters
    const int nSteps = 10;  // number of simulation steps
    const float dt = 0.01f;

    double rate = 0, dRate = 0; 
    const float HztoInts = float(nParticles) * float(nParticles - 1);
    const int skipSteps = 3; // Skip first iterations as they are warm-up o
    double total_time = 0;

    if ((mode != 0) && (mode != 1)) {
        cout << "ERROR!";
        return false;
    }
        
    if (mode == 0) {
        init_particles_serial(particles, nParticles);
    }
    else {
        init_particles_parallel(&all_Particles_Parallel, nParticles);
    }
    cout << "Step" << "\t" << "Time(ms)" << "\t" << "Interact/s" << endl;

    for (int step = 1; step <= nSteps; step++) {
        auto tStart = chrono::steady_clock::now(); // Start timing
        if (mode == 0) {
            move_particles_serial(nParticles, particles, dt);
        }
        else {
            //TODO:  fill in invocation of move_particles_parallel()
            move_particles_parallel(nParticles, &all_Particles_Parallel, dt);
        }
        auto tEnd = chrono::steady_clock::now(); // End timing

        const float HztoInts = float(nParticles) * float(nParticles - 1);
        long long time = chrono::duration_cast<chrono::milliseconds>(tEnd - tStart).count();

        // Collect statistics
        if (step > skipSteps) {
            rate += (HztoInts / time);
            total_time += (double)time;
        }
        cout.precision(4);
        cout << step << "\t" << time << "\t\t" << 1000 * HztoInts / time << "  " << (step <= skipSteps ? "*" : "") << endl;
    }

    // Summarize results
    rate /= (double)(nSteps - skipSteps);  // only average the last (nSteps - skipSteps) steps
    total_time /= (nSteps - skipSteps);

    {
        string str="";
        if (mode == 0)
            str = "Average serial performance";
        else
            str = "Average parallel performance";

        cout << endl << str << ": " << total_time << "ms  " << "rate:" << rate << endl;
        cout << "-----------------------------------------------------" << endl;
        cout << "* - warm-up, not included in average" << endl << endl;
    }

    performance[mode] = total_time;

    // validation
    if (mode == 0)
       return validate_results(particles[0].x, particles[4].y, particles[2].z);
    else
       return validate_results(all_Particles_Parallel.x[0], all_Particles_Parallel.y[4], all_Particles_Parallel.z[2]); 
}


int main() {
    bool v1, v2;
    v1 = v2 = false;
    
    // run serial 
    cout << "\nPropagating " << N << " particles using serial implementation on CPU" << endl << endl;
    v1 = run_simulation(0);

    if (v1 == false)
        cout << "serial validation failed" << endl;

    // TODO: change number of threads
    nthreads = 1;
	
    cout << "\n\nPropagating " << N << " particles using " << nthreads << " threads on CPU" << endl << endl;
    v2 = run_simulation(1);
    
    if (v2 == false)
        cout << "parallel validation failed" << endl;
    
    if (v1 && v2) 
        cout << endl << "Speedup=" << ((performance[1]!=1) ? performance[0]/performance[1]: 1) << endl;
}
