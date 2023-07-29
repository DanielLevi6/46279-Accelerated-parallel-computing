#pragma once

template <typename T>
class SmartPointer
{
private:
    T* ptr; // Is used for getting access to the memory address
public:
    SmartPointer(T* ptrToStore) : ptr(ptrToStore) { }

    ~SmartPointer()
    {
        delete pData;
    }

    T& operator*()
    {
        return *ptr;
    }

    T* operator->()
    {
        return ptr;
    }
};



