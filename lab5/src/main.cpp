#include <iostream>
#include <thread>
#include <mpi.h>
#include "LockingQueue.h"
#include "Task.h"

#define ITERATIONS_TASK_COUNT 30
#define BASE_TASK_COUNT 2000
#define WORKER_FINISHED (-1)
//#define BALANCE

static int32_t ProcNum = 0;
static int32_t ClusterSize = 0;
static double TotalDisbalance = 0;
static bool IsFinished = false;
std::thread WorkerThread = {};
std::thread ReceiverThread = {};
auto TaskQueue = LockingQueue<Task>();

void GenerateTasks(int32_t iter, int32_t clusterSize) {
    for (uint32_t i = 0; i < BASE_TASK_COUNT; ++i) {
        TaskQueue.Push((std::abs(ProcNum - (iter % clusterSize))));
    }
}

void executeTasks() {
    while (!TaskQueue.IsEmpty()) {
        std::optional<Task> task;
        task = TaskQueue.TryPop();
        task->executeTask();
    }
}

void Receiver() {
    MPI_Barrier(MPI_COMM_WORLD);

    int32_t requester = 0;
    std::vector<Task> tasks = {};
}

void LaunchReceiverThread() {
    ReceiverThread = std::thread(Receiver);
}

void Work(int iter) {
    GenerateTasks(iter, ClusterSize);
    executeTasks();
}

void Worker(int iter) {
    double iterationDuration = 0.0, shortestDuration = 0.0, longestDuration = 0.0;
    double iterationStartTime = MPI::Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    Work(iter);
    iterationDuration = MPI::Wtime() - iterationStartTime;

    MPI_Allreduce(&iterationDuration, &longestDuration, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&iterationDuration, &shortestDuration, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    TotalDisbalance += (longestDuration - shortestDuration) / longestDuration;

    if (ProcNum == 0) {
        std::cout << "Disbalance is " << ((longestDuration - shortestDuration) / longestDuration) * 100.f << "%"
                  << std::endl;
        std::cout << "Worker " << ProcNum << ": "
                  << "Current iteration: " << iter << std::endl;
    }
}

void WorkerBalanced(int iter){

}

void LaunchWorkerThread(int iter) {
    WorkerThread = std::thread(Worker, iter);
}

void LaunchBalancedWorkerThread(int iter) {
    WorkerThread = std::thread(WorkerBalanced, iter);
}

int main(int argc, char **argv) {
    int a = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &a);
    if (a != MPI_THREAD_MULTIPLE) {
        std::cout << "MPI_THREAD_MULTIPLE not supported!" << std::endl;
        MPI::Finalize();
        return 0;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_size(MPI_COMM_WORLD, &ClusterSize);

    if (ProcNum == 0) {
        std::cout << "ClusterSize: " << ClusterSize << std::endl;
        std::cout << "Iterations: " << ITERATIONS_TASK_COUNT << std::endl;
        std::cout << "Tasks: " << BASE_TASK_COUNT << std::endl;
    }

    LaunchReceiverThread();

    const double tasksBeginTime = MPI::Wtime();

    for (int curIter = 0; curIter < ITERATIONS_TASK_COUNT; ++curIter) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (ProcNum == 0) {
            std::cout << "Began " << curIter << " iteration" << std::endl;
        }

        double beginTime = MPI::Wtime();
#ifndef BALANCE
        LaunchWorkerThread(curIter);
#endif
#ifdef BALANCE
        LaunchBalancedWorkerThread(curIter);
#endif
        WorkerThread.join();
        double endTime = MPI::Wtime();
        double timeElapsed = endTime - beginTime;
        std::cout << ProcNum << ": at " << curIter << "\nthis iteration took it " << timeElapsed << std::endl;
        std::cout << "it slept for " << std::abs(ProcNum - (curIter % ClusterSize)) << std::endl;
    }
    ReceiverThread.join();

    const double tasksEndTime = MPI::Wtime() - tasksBeginTime;
    double tasksMaxDuration = 0.0;
    MPI_Allreduce(&tasksEndTime, &tasksMaxDuration, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (ProcNum == 0) {
//        std::cout << "Disbalance: " << s_DisbalanceSum / ITERATIONS_TASK_COUNT * 100.f << "%" << std::endl;
        std::cout << "total time: " << tasksMaxDuration << std::endl;
    }

    MPI::Finalize();
    return 0;
}

