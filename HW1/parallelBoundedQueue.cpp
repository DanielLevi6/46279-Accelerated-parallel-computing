#include <mutex>
#include <vector>
#include <condition_variable>
#include <iostream>
#include <thread>


class parallelBoundedQueue {
    size_t capacity;
    size_t used_counter;
    std::mutex queue_mutex;
    std::vector<int> buffer;
    std::condition_variable is_queue_empty;
    std::condition_variable is_queue_full;

    int tail; // Used as the tail index of the allocated block. decreased by 1 in each push
    int head; // Used as the head index of the allocated block. increased by 1 in each pop

public:
    parallelBoundedQueue(size_t new_queue_capacity) : capacity(new_queue_capacity), used_counter(0), tail(0), head(0) {
        buffer.resize(this->capacity);
    }

    ~parallelBoundedQueue() {}

    // Pop the next element from the queue.
    // If the buffer is empty, the calling thread waits until being notified of
    // new elements in the queue
    int pop() {
        std::unique_lock<std::mutex> lock(queue_mutex);
        is_queue_empty.wait(lock, [this](){return used_counter > 0;});
        
        int pop_v = buffer[tail];
        tail = (tail + 1) % capacity;
        used_counter--;

        is_queue_full.notify_one();

        return pop_v;
    }

    // Push a new integer to the queue.
    // If the buffer is full, the calling thread should wait until being notified 
    // that the queue is not full anymore
    // v - the new integer to push into the queue.
    void push(int v) {
        std::unique_lock<std::mutex> lock(queue_mutex);
        is_queue_full.wait(lock, [this](){return used_counter < capacity;});

        buffer[head] = v;
        head = (head + 1) % capacity;
        used_counter++;

        is_queue_empty.notify_one();
    }
};

void producer_func(parallelBoundedQueue& queue, int thread_id) {
    for(int i=0; i<10; i++) {
        queue.push(i);
        std::cout << "Producer " << thread_id << ": The producer pushed the value - " << i << std::endl;
    }
}

void consumer_func(parallelBoundedQueue& queue, int thread_id) {
    for(int i=0; i<10; i++) {
        int pop_v = queue.pop();
        std::cout << "Consumer " << thread_id << ": The consumer popped the value - " << pop_v << std::endl;
    }
}

int main() {
    parallelBoundedQueue queue(5);

    // two producer threads
    std::thread producer1(producer_func, std::ref(queue), 1);
    std::thread producer2(producer_func, std::ref(queue), 2);

    // two consumer threads
    std::thread consumer1(consumer_func, std::ref(queue), 1);
    std::thread consumer2(consumer_func, std::ref(queue), 2);

    producer1.join();
    producer2.join();
    consumer1.join();
    consumer2.join();

    return 0;
}
















