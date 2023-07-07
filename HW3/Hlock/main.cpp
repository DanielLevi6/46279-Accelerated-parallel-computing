/**
 * Compilation command: g++ -std=c++11 -O3 -o Hlock ./Hlock.cpp ./Hlock.hpp ./main.cpp -lpthread
*/
#include "Hlock.hpp"

void lock_task1(HierarchicalMutex& m1, HierarchicalMutex& m2, HierarchicalMutex& m3, HierarchicalMutex& m4, HierarchicalMutex& m5) {
    m1.lock();
    std::cout << "task1 - locked m1" << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(2));
    m1.unlock();
    std::cout << "task1 - unlocked m1" << std::endl;
    std::cout << "end of task1" << std::endl;
}

void lock_task2(HierarchicalMutex& m1, HierarchicalMutex& m2, HierarchicalMutex& m3, HierarchicalMutex& m4, HierarchicalMutex& m5) {
    std::cout << "task2 - tring to lock m1" << std::endl;
    if(m1.try_lock()) {
        std::cout << "task2 - succeded locking m1" << std::endl;
        m1.unlock();
        std::cout << "task2 - unlocked m1" << std::endl;
    }
    else {
        std::cout << "task2 - failed locking m1" << std::endl;
    }
    std::cout << "end of task2" << std::endl;
}

void lock_task3(HierarchicalMutex& m1, HierarchicalMutex& m2, HierarchicalMutex& m3, HierarchicalMutex& m4, HierarchicalMutex& m5) {
    m1.lock();
    std::cout << "task3 - locked m1" << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    m5.lock();
    std::cout << "task3 - locked m5" << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    m5.unlock();
    std::cout << "task3 - unlocked m5" << std::endl;
    m1.unlock();
    std::cout << "task3 - unlocked m1" << std::endl;
    std::cout << "end of task3" << std::endl;
}

void lock_task4(HierarchicalMutex& m1, HierarchicalMutex& m2, HierarchicalMutex& m3, HierarchicalMutex& m4, HierarchicalMutex& m5) {
    m1.lock();
    std::cout << "task4 - locked m1" << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    m2.lock();
    std::cout << "task4 - locked m2" << std::endl;
    m3.lock();
    std::cout << "task4 - locked m3" << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    m3.unlock();
    std::cout << "task4 - unlocked m3" << std::endl;
    m2.unlock();
    std::cout << "task4 - unlocked m2" << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    m1.unlock();
    std::cout << "task4 - unlocked m1" << std::endl;
    std::cout << "end of task4" << std::endl;
}

void lock_task5(HierarchicalMutex& m1, HierarchicalMutex& m2, HierarchicalMutex& m3, HierarchicalMutex& m4, HierarchicalMutex& m5) {
    m2.lock();
    std::cout << "task5 - locked m2" << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    m3.lock();
    std::cout << "task5 - locked m3" << std::endl;
    m4.lock();
    std::cout << "task5 - locked m4" << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    m5.lock();
    std::cout << "task5 - locked m5" << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    m5.unlock();
    std::cout << "task5 - unlocked m5" << std::endl;
    m4.unlock();
    std::cout << "task5 - unlocked m4" << std::endl;
    m3.unlock();
    std::cout << "task5 - unlocked m3" << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    m2.unlock();
    std::cout << "task5 - unlocked m2" << std::endl;
    std::cout << "end of task5" << std::endl;
}


int main()
{
    HierarchicalMutex m1(1);
    HierarchicalMutex m2(2);
    HierarchicalMutex m3(3);
    HierarchicalMutex m4(4);
    HierarchicalMutex m5(5);

    std::cout << "############ main - thread1 starts ##############" << std::endl << std::endl;
    std::cout << "Test 1-try_lock fail" << std::endl;
    std::thread t1(lock_task1, std::ref(m1), std::ref(m2), std::ref(m3), std::ref(m4), std::ref(m5));
    std::thread t2(lock_task2, std::ref(m1), std::ref(m2), std::ref(m3), std::ref(m4), std::ref(m5));
    t1.join();
    t2.join();
    
    std::cout << std::endl << "Test 2-try_lock success" << std::endl;
    std::thread t3(lock_task2, std::ref(m1), std::ref(m2), std::ref(m3), std::ref(m4), std::ref(m5));
    t3.join();

    std::cout << std::endl << "Test 3-hirerchy violation" << std::endl;
    m1.lock();
    std::cout << "Locked m1" << std::endl;
    m5.lock();
    std::cout << "Locked m5" << std::endl;
    try {
        std::cout << "tring to lock m3(expecting to exception)" << std::endl;
        m3.lock();
    }
    catch(std::logic_error& e) {
        std::cout << "Exception caught - " << e.what() << std::endl;
    }
    m1.unlock();
    std::cout << "Unlocked m1" << std::endl;
    m5.unlock();
    std::cout << "Unlocked m5" << std::endl;

    std::cout << std::endl << "Test 4-three threads in parallel" << std::endl;
    std::thread t4(lock_task3, std::ref(m1), std::ref(m2), std::ref(m3), std::ref(m4), std::ref(m5));
    std::thread t5(lock_task4, std::ref(m1), std::ref(m2), std::ref(m3), std::ref(m4), std::ref(m5));
    std::thread t6(lock_task5, std::ref(m1), std::ref(m2), std::ref(m3), std::ref(m4), std::ref(m5));
    t4.join();
    t5.join();
    t6.join();

    std::cout << "############ main - end ##############" << std::endl;
    return 0;
}