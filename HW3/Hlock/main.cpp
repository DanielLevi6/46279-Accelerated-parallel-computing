#include "Hlock.hpp"

void lock_task1(HierarchicalMutex& m1, HierarchicalMutex& m2, HierarchicalMutex& m3, HierarchicalMutex& m4, HierarchicalMutex& m5) {
    std::cout << "task1 - locks m1" << std::endl;
    m1.lock();
    std::this_thread::sleep_for(std::chrono::seconds(2));
    std::cout << "task1 - locks m2" << std::endl;
    m3.lock();
    std::cout << "task1 - locks m3" << std::endl;
    m2.lock();
    std::cout << "task1 - unlocks m3" << std::endl;
    m3.unlock();
    std::cout << "task1 - unlocks m2" << std::endl;
    m1.unlock();
}

// Checks the 
void lock_task2(HierarchicalMutex& m1, HierarchicalMutex& m2, HierarchicalMutex& m3, HierarchicalMutex& m4, HierarchicalMutex& m5) {
    std::cout << "task2 - locks m4" << std::endl;
    m4.lock();
    std::cout << "task2 - locks m5" << std::endl;
    m5.lock();
    try {
        std::cout << "task2 - locks m1" << std::endl;
        m1.lock();
    }
    catch(std::logic_error& e) {
        std::cout << "task2 - caught exception: " << e.what() << std::endl;
    }
    std::cout << "task2 - unlocks m5" << std::endl;
    m5.unlock();
    std::cout << "task2 - unlocks m4" << std::endl;
    m4.unlock();
}

void lock_task3(HierarchicalMutex& m1, HierarchicalMutex& m2, HierarchicalMutex& m3, HierarchicalMutex& m4, HierarchicalMutex& m5) {
    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "task3 - trys locking m1" << std::endl;
    if(m1.try_lock()) {
        std::cout << "task3 - succeded locking m1" << std::endl;
        m1.unlock();
    }
    else {
        std::cout << "task3 - failed locking m1" << std::endl;
    }
    m2.unlock();
    m1.unlock();
}

int main()
{
    HierarchicalMutex m1(1);
    HierarchicalMutex m2(2);
    HierarchicalMutex m3(3);
    HierarchicalMutex m4(4);
    HierarchicalMutex m5(5);

    std::thread t1(lock_task1, std::ref(m1), std::ref(m2), std::ref(m3), std::ref(m4), std::ref(m5));
    std::thread t2(lock_task2, std::ref(m1), std::ref(m2), std::ref(m3), std::ref(m4), std::ref(m5));
    std::thread t3(lock_task3, std::ref(m1), std::ref(m2), std::ref(m3), std::ref(m4), std::ref(m5));

    return 0;
}