#pragma once
#include <condition_variable>
#include <optional>
#include <queue>


template<class T>
class LockingQueue {
public:
    LockingQueue() = default;
    ~LockingQueue() = default;

    void Push(const T& value)
    {
        {
            std::unique_lock lock(Mutex);
            queue.emplace(value);
        }
        CondVar.notify_all();
    }

    std::optional<T> TryPop()
    {
        std::unique_lock lock(Mutex);

        if (queue.empty()) {
            return std::nullopt;
        }

        T ret = queue.front();
        queue.pop();
        return ret;
    }

    std::optional<T> WaitAndPop()
    {
        std::unique_lock lock(Mutex);
        while (queue.empty() && !Interrupted) {
            CondVar.wait(lock);
        }

        if (Interrupted) {
            Interrupted = false;
            return std::nullopt;
        }

        T ret = queue.front();
        queue.pop();
        return ret;
    }

    std::optional<T> WaitAndFront()
    {
        std::unique_lock lock(Mutex);
        while (queue.empty() && !Interrupted) {
            CondVar.wait(lock);
        }

        if (Interrupted) {
            Interrupted = false;
            return std::nullopt;
        }

        return queue.front();
    }

    void InterruptWaiting()
    {
        CondVar.notify_all();
        Interrupted = true;
    }

    std::size_t Size() const
    {
//        std::lock_guard lock(Mutex);
        return queue.size();
    }

    bool IsEmpty() const
    {
//        std::lock_guard lock(Mutex);
        return queue.empty();
    }

private:
    std::queue<T> queue;
    std::condition_variable CondVar;
    std::mutex Mutex;

    bool Interrupted = false;
};
