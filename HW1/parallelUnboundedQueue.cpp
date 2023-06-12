#include <mutex>
#include <condition_variable>
#include <iostream>
#include <thread>

class Node {
public:
    int value;
    Node* next;

    Node(int val) : value(val), next(nullptr) {}
    ~Node() {}
};

class parallelUnboundedQueue {
    std::mutex queue_mutex;

    Node* tail;
    Node* head;

    int used_counter;

public:
    parallelUnboundedQueue() : tail(nullptr), head(nullptr), used_counter(0) {}

    ~parallelUnboundedQueue() {
        while(tail) {
            Node* to_delete = tail;
            tail = tail->next;
            delete to_delete;
        }
    }

    // returns false if queue is empty 
    bool pop(int &value) {
        queue_mutex.lock();

        bool res;
        if(used_counter == 0) {
            res = false;
        } else {
            value = tail->value;
            Node* to_delete = tail;
            tail = tail->next;
            head = nullptr;
            delete to_delete;

            used_counter--;
            res = true;
        }
        queue_mutex.unlock();
        return res;
    }

    // v - the new integer to push into the queue.
    void push(int val) {
        queue_mutex.lock();

        Node* new_node = new Node(val);
        if(used_counter == 0) {
            head = new_node;
            tail = new_node;
            used_counter++;
        } else {
            head->next = new_node;
            head = head->next;
            used_counter++;
        }

        queue_mutex.unlock();
    }
};

void producer_func(parallelUnboundedQueue& queue) {
    for(int i=0; i<10; i++) {
        queue.push(i);
        std::cout << "Producer : The producer pushed the value - " << i << std::endl;
    }
}

void consumer_func(parallelUnboundedQueue& queue) {
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
    parallelUnboundedQueue queue;

    // sungle producer thread
    std::thread producer(producer_func, std::ref(queue));

    // single consumer thread
    std::thread consumer(consumer_func, std::ref(queue));

    producer.join();
    consumer.join();

    return 0;
}
















