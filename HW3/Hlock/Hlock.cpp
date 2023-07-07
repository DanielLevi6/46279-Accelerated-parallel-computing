#include "Hlock.hpp"

static thread_local std::stack<int> hierarvcyLocks;

HierarchicalMutex::HierarchicalMutex(int lvl) : lvl(lvl) {}

void HierarchicalMutex::lock() {
    if (this->lvl <= threadLevel) {
        throw std::logic_error("mutex hierarchy violated");
    }
    this->mtx.lock();
    threadLevel = this->lvl;
    hierarvcyLocks.push(this->lvl);
}

void HierarchicalMutex::unlock() {
    hierarvcyLocks.pop();
    threadLevel = hierarvcyLocks.top();
    this->mtx.unlock();
}

bool HierarchicalMutex::try_lock() {
    if (this->lvl <= threadLevel) {
        throw std::logic_error("mutex hierarchy violated");
    }
    if (mtx.try_lock()) {
        hierarvcyLocks.push(this->lvl);
        threadLevel = this->lvl;
        return true;
    }
    return false;
}
