#include <iostream>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <fstream>
#include <vector>

float NormOf(float* vec, size_t size){
    float norm = 0;
#pragma omp parallel for reduction(+:norm)
    for(size_t i = 0; i < size; i++){
        norm += vec[i] * vec[i];
    }
    return sqrtf(norm);
}

void MatVecMult(float* mat, float* vec, float* vecEq, size_t N){
#pragma omp parallel for
    for (std::size_t i = 0; i < N; ++i) // row
    {
        for (std::size_t j = 0; j < N; ++j) // col
        {
            vecEq[i] += mat[i * N + j] * vec[j];
        }
    }
}

void VecVecSub(float* vecA, float* vecB, size_t N){
#pragma omp parallel for
    for (std::size_t i = 0; i < N; ++i)
    {
        vecA[i] -= vecB[i];
    }
}

void VecScalMul(float* vec, float scalar, size_t N){
#pragma omp parallel for
    for (std::size_t i = 0; i < N; ++i)
    {
        vec[i] *= scalar;
    }
}

int main(int argc, char *argv[]) {
    const size_t N = 50 * 50;
    const float tau = -0.05;
    const float epsilon = 0.00001;
    const int maxThreads = omp_get_max_threads();

    std::vector<float> matA(N * N);
    std::ifstream input("../src/testerPy/matA.bin", std::ios::binary);
    if (!input.is_open())
    {
        throw std::runtime_error("Failed to load data file.");
    }
    input.read(reinterpret_cast<char*>(matA.data()), sizeof(float) * N * N);
    input.close();

    std::vector<float> vecB(N);
    std::ifstream inputVec("../src/testerPy/vecB.bin", std::ios::binary);
    if (!inputVec.is_open())
    {
        throw std::runtime_error("Failed to load data file.");
    }
    inputVec.read(reinterpret_cast<char*>(vecB.data()), sizeof(float) * N);
    inputVec.close();


    std::vector<float> foundVecX(N);
    std::vector<float> intermediate(N);
    std::ofstream csv("../src/testerPy/PerfEach.csv", std::ios::out);
    float vecBNorm = NormOf(vecB.data(), N);

    for (int numThreads = 1; numThreads <= maxThreads; ++numThreads)
    {
        std::cout << "doing "<< numThreads << " threads" << std::endl;

        omp_set_num_threads(numThreads);

//#pragma omp single
        std::fill_n(foundVecX.data(), N, 0);

        auto start = std::chrono::high_resolution_clock::now();
        size_t iterations = 0;

        //calc
        while(true){
//#pragma omp single
            std::fill_n(intermediate.data(), N, 0);

            MatVecMult(matA.data(), foundVecX.data(), intermediate.data(), N);

            VecVecSub(intermediate.data(), vecB.data(), N);

            float normofIntermediate;
//            std::cout << (NormOf(intermediate.data(), N) / vecBNorm) << std::endl;

            if((NormOf(intermediate.data(), N) / vecBNorm) < epsilon){
                break;
            }

            VecScalMul(intermediate.data(), tau, N);

            VecVecSub(foundVecX.data(), intermediate.data(), N);

#pragma omp single
            iterations++;
        }


        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        double seconds = duration.count();

        std::cout << "Execution time: " << seconds << " seconds" << std::endl;
        std::cout << "iterations: " << iterations << std::endl;

        csv << numThreads << ',' << seconds << std::endl;

        std::ofstream output("../src/testerPy/OutvecFromEach.bin", std::ios::binary);
        output.write(reinterpret_cast<const char*>(foundVecX.data()), sizeof(float) * foundVecX.size());
        output.close();
//#pragma omp barrier
    }
}
