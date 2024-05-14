#include <iostream>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <fstream>
#include <vector>
#define SCHEDULETYPE guided

#define ISPROCBIND /*proc_bind(close)*/

float g_normVal = 0;

    inline float NormOf(float *vec, size_t size) {
#pragma omp single
    g_normVal = 0;
#pragma omp for reduction(+:g_normVal) schedule(SCHEDULETYPE)
    for (size_t i = 0; i < size; i++) {
        g_normVal += vec[i] * vec[i];
    }
    return sqrtf(g_normVal);
}

void MatVecMult(float *mat, float *vec, float *vecEq, size_t N) {
#pragma omp for schedule(SCHEDULETYPE)
    for (std::size_t i = 0; i < N; ++i) // row
    {
        for (std::size_t j = 0; j < N; ++j) // col
        {
            vecEq[i] += mat[i * N + j] * vec[j];
        }
    }
}

void VecVecSub(float *vecA, float *vecB, size_t N) {
#pragma omp for schedule(SCHEDULETYPE)
    for (std::size_t i = 0; i < N; ++i) {
        vecA[i] -= vecB[i];
    }
}

void VecScalMul(float *vec, float scalar, size_t N) {
#pragma omp for schedule(SCHEDULETYPE)
    for (std::size_t i = 0; i < N; ++i) {
        vec[i] *= scalar;
    }
}

int main(int argc, char *argv[]) {
    const size_t N = 50 * 50;
    const float tau = -0.05;
    const float epsilon = 0.00001;
    const int maxThreads = omp_get_max_threads();

    size_t iterations = 0;
    std::chrono::time_point<std::chrono::high_resolution_clock> start;

    std::vector<float> matA(N * N);
    std::ifstream input("../src/testerPy/matA.bin", std::ios::binary);
    if (!input.is_open()) {
        throw std::runtime_error("Failed to load data file.");
    }
    input.read(reinterpret_cast<char *>(matA.data()), sizeof(float) * N * N);
    input.close();

    std::vector<float> vecB(N);
    std::ifstream inputVec("../src/testerPy/vecB.bin", std::ios::binary);
    if (!inputVec.is_open()) {
        throw std::runtime_error("Failed to load data file.");
    }
    inputVec.read(reinterpret_cast<char *>(vecB.data()), sizeof(float) * N);
    inputVec.close();

    std::vector<float> foundVecX(N);
    std::vector<float> intermediate(N);
    std::ofstream csv("../src/testerPy/PerfScheduleGUIDEDPROCWPlaces.csv", std::ios::out);
    float vecBNorm = NormOf(vecB.data(), N);

    for (int numThreads = 1; numThreads <= maxThreads; ++numThreads) {
    #pragma omp single
        {
            std::cout << "doing " << numThreads << " threads" << std::endl;

            omp_set_num_threads(numThreads);

            std::fill_n(foundVecX.data(), N, 0);

            start = std::chrono::high_resolution_clock::now();
            iterations = 0;
        }
        //calc
    #pragma omp parallel ISPROCBIND
        {
            while (true) {
            #pragma omp single
                std::fill_n(intermediate.data(), N, 0);

                MatVecMult(matA.data(), foundVecX.data(), intermediate.data(), N);

                VecVecSub(intermediate.data(), vecB.data(), N);

//            std::cout << (NormOf(intermediate.data(), N) / vecBNorm) << std::endl;

                if ((NormOf(intermediate.data(), N) / vecBNorm) < epsilon) {
                    break;
                }

                VecScalMul(intermediate.data(), tau, N);

                VecVecSub(foundVecX.data(), intermediate.data(), N);

            #pragma omp single
                iterations++;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        double seconds = duration.count();

        std::cout << "Execution time: " << seconds << " seconds" << std::endl;
        std::cout << "iterations: " << iterations << std::endl;

        csv << numThreads << ',' << seconds << std::endl;

        std::ofstream output("../src/testerPy/Outvec.bin", std::ios::binary);
        output.write(reinterpret_cast<const char *>(foundVecX.data()), sizeof(float) * foundVecX.size());
        output.close();
//#pragma omp barrier
    }
//    system("python ../src/testerPy/visualize_copy.py");
}
