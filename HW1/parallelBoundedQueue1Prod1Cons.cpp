#include <mutex>
#include <vector>
#include <condition_variable>
#include <iostream>
#include <thread>


class parallelBoundedQueue1Prod1Cons {
    size_t capacity;
    size_t used_counter;
    std::mutex queue_mutex;
    std::vector<int> buffer;

    int tail; // Used as the tail index of the allocated block. decreased by 1 in each push
    int head; // Used as the head index of the allocated block. increased by 1 in each pop

public:
    parallelBoundedQueue1Prod1Cons(size_t new_queue_capacity) : capacity(new_queue_capacity), used_counter(0), tail(0), head(0) {
        buffer.resize(this->capacity);
    }

    ~parallelBoundedQueue1Prod1Cons() {}

    // Pop the next element (integer value) from the queue.
    // if the buffer is empty, ruturn_false. Otherwise, true.
    bool pop(int &val) {
        std::unique_lock<std::mutex> lock(queue_mutex);
        if(used_counter == 0) {
            return false;
        }
        
        val = buffer[tail];
        tail = (tail + 1) % capacity;
        used_counter--;

        return true;
    }

    // Push a new integer to the queue.
    // if queue is full, return false. Otherwise, true.
    // v - the new integer to push into the queue.
    bool push(int v) {
        std::unique_lock<std::mutex> lock(queue_mutex);
        if(used_counter == capacity) {
            return false;
        }

        buffer[head] = v;
        head = (head + 1) % capacity;
        used_counter++;

        return true;
    }
};

void producer_func(parallelBoundedQueue1Prod1Cons& queue) {
    for(int i=0; i<10; i++) {
        bool res = queue.push(i);
        if(res) {
            std::cout << "Producer : The producer pushed the value - " << i << std::endl;
        } else {
            std::cout << "Producer couldn't push a value" << std::endl;
            i--;
        }
    }
}

void consumer_func(parallelBoundedQueue1Prod1Cons& queue) {
    for(int i=0; i<10; i++) {
        int pop_v;
        bool res = queue.pop(pop_v);
        if(res) {
            std::cout << "Consumer : The consumer popped the value - " << pop_v << std::endl;
        } else {
            std::cout << "Consumer couldn't pop a value" << std::endl;
            i--;
        }
    }
}

int main() {
    parallelBoundedQueue1Prod1Cons queue(5);

    // sungle producer thread
    std::thread producer(producer_func, std::ref(queue));

    // single consumer thread
    std::thread consumer(consumer_func, std::ref(queue));

    producer.join();
    consumer.join();

    return 0;
}
















