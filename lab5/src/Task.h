#pragma once

#include <cstdint>
#include <thread>
#include <chrono>

class Task{
public:
    Task(uint32_t time) : sleepMili(time) {}
    ~Task() = default;

    void executeTask(){
        std::this_thread::sleep_for(std::chrono::milliseconds(sleepMili));
    }
private:
    uint32_t sleepMili;
};