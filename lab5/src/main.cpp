#include <iostream>
#include <thread>
#include <mpi.h>
#include "LockingQueue.h"
#include "Task.h"

#define ITERATIONS_TASK_COUNT 5
#define BASE_TASK_COUNT 10000
#define WORKER_FINISHED (-1)
#define NO_TASKS_TO_SHARE (-2)
#define MIN_THRESHOLD 2000
#define BALANCE

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
//        std::cout << "Worker " << ProcNum << ": "
//                  << "Current iteration: " << iter << std::endl;
    }
}

void Receiver() {
    MPI_Barrier(MPI_COMM_WORLD);

    int32_t Msg = 0;
    std::vector<uint32_t> tasks = {};
    std::cout << ProcNum <<"'s receiver is up" << std::endl;
    while (!IsFinished) {
        MPI_Status status;
        MPI_Request request;
        uint32_t asd = 0;
        int flag = 0;
        std::cout <<  "b4 irecv on " << ProcNum << std::endl;
        MPI_Irecv(&Msg, 1, MPI_INT32_T, MPI_ANY_SOURCE, 255, MPI_COMM_WORLD, &request);
        while (!IsFinished && !flag) {
            MPI_Test(&request, &flag, &status);
//            if(!(asd % 10000000)) {
//                std::cout << "asdad " << ProcNum <<" " << status.MPI_ERROR << " "<< status.MPI_SOURCE << " " << status.MPI_TAG << " " << flag << std::endl;
//            }
//            asd++;
        }
        if (IsFinished) {
            std::cout << ProcNum << "'s receiver returned" << std::endl;
            return;
        }
        std::cout << "request from " << Msg << " received on " << ProcNum << std::endl;

//        if (Msg == WORKER_FINISHED) {
////            IsFinished = true;
//            continue;
//        }

        int32_t taskNum = 0;
//        if (TaskQueue.Size() > MIN_THRESHOLD) {
        if (TaskQueue.Size() > BASE_TASK_COUNT/ClusterSize) {
//            uint32_t tasksToShare = TaskQueue.Size()/2;
            uint32_t tasksToShare = BASE_TASK_COUNT / ClusterSize;
            tasks.resize(tasksToShare);
            for (uint32_t i = 0; i < tasksToShare; ++i) {
                std::optional<Task> task = TaskQueue.TryPop();
                if (task == std::nullopt) break;
                tasks[i] = task.value().sleepMilli;
                taskNum++;
            }
        }
        if (taskNum == 0) {
            taskNum = NO_TASKS_TO_SHARE;
        }
        MPI_Send(&taskNum, 1, MPI_INT32_T, Msg, ProcNum, MPI_COMM_WORLD);
        std::cout << taskNum << " sent from " << ProcNum << " to " << Msg << std::endl;
        if (taskNum != NO_TASKS_TO_SHARE) {
            MPI_Send(tasks.data(), taskNum, MPI_INT32_T, Msg, ProcNum, MPI_COMM_WORLD);
            std::cout << "tasks sent from " << ProcNum << " to " << Msg << std::endl;
        }
        Msg = 0;
        tasks.clear();
    }
}


void LaunchReceiverThread() {
    ReceiverThread = std::thread(Receiver);
}

void WorkerBalanced(int iter) {
    double iterationDuration = 0.0, shortestDuration = 0.0, longestDuration = 0.0;
    int32_t taskAmountRecv = 0;
    double iterationStartTime = MPI::Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    Work(iter);
    std::cout << "tasks done on " << ProcNum << std::endl;
    {
        for (int32_t curNum = 0; curNum < ClusterSize; ++curNum) {
            MPI_Request request;
            MPI_Status  status;
            if (curNum == ProcNum) continue;
            std::cout << ProcNum << " requested from " << curNum << std::endl;
            MPI_Send(&ProcNum, 1, MPI_INT32_T, curNum, 255, MPI_COMM_WORLD);
//            MPI_Isend(&ProcNum, 1, MPI_INT32_T, curNum, 255, MPI_COMM_WORLD, &request);
//            int flag = 0;
//            while(!flag){
//                MPI_Test(&request, &flag, &status);
//            }

            // Get task count
//            MPI_Status status;
//            std::cout << ProcNum << " receiving from " << curNum << std::endl;
            std::cout << "b4 response on " << ProcNum << std::endl;
            MPI_Recv(&taskAmountRecv, 1, MPI_INT32_T, curNum, curNum, MPI_COMM_WORLD, &status);
            std::cout << "after response on " << ProcNum << std::endl;
            if (taskAmountRecv != NO_TASKS_TO_SHARE) {
                auto receivedTaskCount = taskAmountRecv;
                std::vector<uint32_t> receiveTaskArray(receivedTaskCount);
                MPI_Recv(receiveTaskArray.data(), receivedTaskCount, MPI_UNSIGNED, curNum, curNum, MPI_COMM_WORLD,
                         &status);

                for (auto &val: receiveTaskArray)
                    TaskQueue.Push(val);

                executeTasks();
                std::cout << ProcNum << " finished requested tasks from " << curNum << std::endl;
            }
            else std::cout << ProcNum << " request yielded no tasks from " << curNum << std::endl;
        }
    }

    iterationDuration = MPI::Wtime() - iterationStartTime;

    MPI_Allreduce(&iterationDuration, &longestDuration, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&iterationDuration, &shortestDuration, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    TotalDisbalance += (longestDuration - shortestDuration) / longestDuration;

    if (ProcNum == 0) {
        std::cout << "\033[31mDisbalance is " << ((longestDuration - shortestDuration) / longestDuration) * 100.f << "%\033[0m"
                  << std::endl;
        std::cout << "Worker " << ProcNum << ": "
                  << "Current iteration: " << iter << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
//    if (ProcNum == 0) {
//        uint32_t executionFinishedMessage = WORKER_FINISHED;
//        for (int32_t curNum = 0; curNum < ClusterSize; ++curNum) {
//            if (curNum == ProcNum) continue;
//
//            MPI::COMM_WORLD.Send(&executionFinishedMessage, 1, MPI::UNSIGNED, procNum, s_ClusterSize + 1);
//            MPI_Send(&executionFinishedMessage, 1, MPI_UINT32_T, curNum, 255, MPI_COMM_WORLD);
//            std::cout << "SENT WORK_FIN" << std::endl;
//        }
//    }
    std::cout << "SET WORK_FIN" << std::endl;

    IsFinished = true;
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

    const double tasksBeginTime = MPI::Wtime();

//    LaunchReceiverThread();

    for (int curIter = 0; curIter <= ITERATIONS_TASK_COUNT; ++curIter) {
        MPI_Barrier(MPI_COMM_WORLD);
        IsFinished = false;
        if (ProcNum == 0) {
            std::cout << "\033[32mBegan " << curIter << " iteration\033[0m" << std::endl;
        }

        LaunchReceiverThread();
        MPI_Barrier(MPI_COMM_WORLD);
        double beginTime = MPI::Wtime();
#ifndef BALANCE
        LaunchWorkerThread(curIter);
#endif
#ifdef BALANCE
        LaunchBalancedWorkerThread(curIter);
#endif
        std::cout << " before join" << std::endl;
        WorkerThread.join();
        ReceiverThread.join();
        double endTime = MPI::Wtime();
        double timeElapsed = endTime - beginTime;
        std::cout << ProcNum << ": at " << curIter << "\tthis iteration took it " << timeElapsed << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        sleep(1);
//        std::cout << "it slept for " << std::abs(ProcNum - (curIter % ClusterSize)) << std::endl;
    }
//    ReceiverThread.join();

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

