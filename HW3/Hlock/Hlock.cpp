#include "Hlock.hpp"

static thread_local std::stack<int> hierarcyLocks;
thread_local int HierarchicalMutex::threadLevel = 0;

HierarchicalMutex::HierarchicalMutex(int lvl) : lvl(lvl) {}

void HierarchicalMutex::lock() {
    if (this->lvl <= threadLevel) {
        throw std::logic_error("mutex hierarchy violated");
    }
    this->mtx.lock();
    threadLevel = this->lvl;
    hierarcyLocks.push(this->lvl);
}

void HierarchicalMutex::unlock() {
    if(hierarcyLocks.empty()) {
        throw std::logic_error("mutex hierarchy violated");
    }
    hierarcyLocks.pop();
    threadLevel = hierarcyLocks.empty()? 0 : hierarcyLocks.top();
    this->mtx.unlock();
}

bool HierarchicalMutex::try_lock() {
    if (this->lvl <= threadLevel) {
        throw std::logic_error("mutex hierarchy violated");
    }
    if (mtx.try_lock()) {
        hierarcyLocks.push(this->lvl);
        threadLevel = this->lvl;
        return true;
    }
    return false;
}
