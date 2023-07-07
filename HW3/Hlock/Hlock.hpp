#ifndef HLOCK_HPP
#define HLOCK_HPP

#include <thread>
#include <mutex>
#include <iostream>
#include <stack>

class HierarchicalMutex {
    std::mutex mtx;
    int lvl;
public:

    static thread_local int threadLevel;

    // lvl indicate the order/number of the mutex. This class checks the order to prevent a deadlock
    HierarchicalMutex(int lvl);
    
    void lock();

    void unlock();

    bool try_lock();
};

#endif // HLOCK_HPP