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


// TODO: add missing function
void move_particles_parallel(const int nParticles, AllParticlesParallel* const particles_parallel, const float dt) {
    int number_of_elements_per_op = (VECTOR_LENGTH/BITS_PER_BYTE)/FLOAT_SIZE;
    int loops = ceil(nParticles / number_of_elements_per_op);
    
    // Loop over particles that experience force
    tbb::parallel_for(tbb::blocked_range<int>(0, nParticles),
                            [=](const tbb::blocked_range<int>& r) {
                                for (int i = r.begin(); i < r.end(); i++) {
                                    //printf("i = %d\n", i);
                                    // Components of the gravity force on particle i
                                    float Fx = 0, Fy = 0, Fz = 0;

                                    // Loop over particles that exert force: vectorization expected here
                                    for (int j = 0; j < loops; j++) {
                                        //printf("j = %d\n", j);
                                        // Avoid singularity and interaction with self
                                        __m256 softening = _mm256_set1_ps(1e-20f);
                                        
                                        // Newton's law of universal gravity
                                        __m256 dx = _mm256_sub_ps(_mm256_load_ps(&(particles_parallel->x[j*number_of_elements_per_op])), _mm256_set1_ps(particles_parallel->x[i]));
                                        __m256 dy = _mm256_sub_ps(_mm256_load_ps(&(particles_parallel->y[j*number_of_elements_per_op])), _mm256_set1_ps(particles_parallel->y[i]));
                                        __m256 dz = _mm256_sub_ps(_mm256_load_ps(&(particles_parallel->z[j*number_of_elements_per_op])), _mm256_set1_ps(particles_parallel->z[i]));

                                        __m256 rr1 = _mm256_div_ps(_mm256_set1_ps(1.0f), _mm256_sqrt_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(dx, dx), _mm256_mul_ps(dy, dy)), _mm256_add_ps(_mm256_mul_ps(dz, dz), softening))));
                                        __m256 drPowerN32 = _mm256_mul_ps(_mm256_mul_ps(rr1, rr1), rr1);
                                        
                                        //Calculate the net force
                                        for(int k=0; k<number_of_elements_per_op; k++) {
                                            Fx += dx[k] * drPowerN32[k];
                                            Fy += dy[k] * drPowerN32[k];
                                            Fz += dz[k] * drPowerN32[k];
                                        }
                                    }
                                    // if(i%1000 == 0) printf("Fx = %f\n", Fx);

                                    particles_parallel->vx[i] += dt * Fx;
                                    particles_parallel->vy[i] += dt * Fy;
                                    particles_parallel->vz[i] += dt * Fz;
                                }
                            }
    );

    // Second loop
    // for(int i=0; i<loops; i++) {
    //     __m256 x = _mm256_load_ps(&(particles_parallel->x[i*number_of_elements_per_op]));
    //     __m256 y = _mm256_load_ps(&(particles_parallel->y[i*number_of_elements_per_op]));
    //     __m256 z = _mm256_load_ps(&(particles_parallel->z[i*number_of_elements_per_op]));
    //     __m256 vx = _mm256_load_ps(&(particles_parallel->vx[i*number_of_elements_per_op]));
    //     __m256 vy = _mm256_load_ps(&(particles_parallel->vy[i*number_of_elements_per_op]));
    //     __m256 vz = _mm256_load_ps(&(particles_parallel->vz[i*number_of_elements_per_op]));
    //     x = _mm256_add_ps(x, _mm256_mul_ps(vx, _mm256_set1_ps(dt)));
    //     y = _mm256_add_ps(y, _mm256_mul_ps(vy, _mm256_set1_ps(dt)));
    //     z = _mm256_add_ps(z, _mm256_mul_ps(vz, _mm256_set1_ps(dt)));
    //     _mm256_store_ps(&(particles_parallel->x[i*number_of_elements_per_op]), x);
    //     _mm256_store_ps(&(particles_parallel->y[i*number_of_elements_per_op]), y);
    //     _mm256_store_ps(&(particles_parallel->z[i*number_of_elements_per_op]), z);
    // }
    // for (int i = 0; i < nParticles; i++) {
    //     particles_parallel->x[i] += particles_parallel->vx[i] * dt;
    //     particles_parallel->y[i] += particles_parallel->vy[i] * dt;
    //     particles_parallel->z[i] += particles_parallel->vz[i] * dt;
    // }
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nParticles),
                        [&](const tbb::blocked_range<size_t>& r) {
                            for(size_t i = r.begin(); i!=r.end(); i++) {
                                particles_parallel->x[i] += particles_parallel->vx[i] * dt;
                                particles_parallel->y[i] += particles_parallel->vy[i] * dt;
                                particles_parallel->z[i] += particles_parallel->vz[i] * dt;
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
    const float epsilon = 0.001f;

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
