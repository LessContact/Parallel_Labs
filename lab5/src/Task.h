#pragma once

#include <cstdint>
#include <thread>
#include <chrono>

class Task{
public:
    Task(uint32_t time) : sleepMilli(time) {}
    ~Task() = default;

    void executeTask(){
        std::this_thread::sleep_for(std::chrono::milliseconds(sleepMilli));
    }

    uint32_t sleepMilli;
};