#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <tbb/tbb.h>
#include <random>

using namespace std;

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

ParticleType particles[N];


void  init_particles_serial(ParticleType* particles, const int nParticles) {
    for (unsigned long long i = 0; i < nParticles; i++) {
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
void init_particles_parallel(ParticleType* particles, const int nParticles) {
    // There are N particles
    unsigned int internal_loops = (nParticles / 15) + 1;
    for (unsigned long long i=0; i<15; i++) {
        float i_rest = (float)(i);
        float i_power_two_rest = (float)((i * i) % 15);
        float i_power_three_rest = (float)((i * i * i) % 15);
        for (unsigned int j = 0; j < internal_loops; j++) {
            unsigned int curr_index = j * 15 + i;
            if(curr_index < nParticles) {
                particles[curr_index].x = i_rest;
                particles[curr_index].y = i_power_two_rest;
                particles[curr_index].z = i_power_three_rest;
                particles[curr_index].vx = 0;
                particles[curr_index].vy = 0;
                particles[curr_index].vz = 0;
            }
        }
    }
}


// TODO: add missing function
void move_particles_parallel(const int nParticles, ParticleType* const particles, const float dt) {
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
/*             run simulations                    
/*
/**************************************************/

bool validate_results(float f1, float f2, float f3) {
    const float epsilon = 0.001f;

    std::cout << "f1 = " << f1 << "(-0.5315)" << endl;
    std::cout << "f2 = " << f2 << "(2.263)" << endl;
    std::cout << "f3 = " << f3 << "(7.335)" << endl;

    if ((fabs(f1 - (-0.5315)) > epsilon) || (fabs(f2 - 2.263) > epsilon) || (fabs(f3 - 7.335) > epsilon)) {
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
        //TODO: fill in invocation of init_particles_parallel()
        init_particles_parallel(particles, nParticles);
    }
    cout << "Step" << "\t" << "Time(ms)" << "\t" << "Interact/s" << endl;

    for (int step = 1; step <= nSteps; step++) {
        auto tStart = chrono::steady_clock::now(); // Start timing
        if (mode == 0) {
            move_particles_serial(nParticles, particles, dt);
        }
        else {
            //TODO:  fill in invocation of move_particles_parallel()
            move_particles_parallel(nParticles, particles, dt);
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
       return validate_results(particles[0].x, particles[4].y, particles[2].z); // TODO: replace with invocation of validate_results(). 
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
