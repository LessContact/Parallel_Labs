#include <mpi.h>
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <vector>

inline float NormOf2(float *vec, size_t size)
{
    float normVal = 0;
    for (size_t i = 0; i < size; i++)
    {
        normVal += vec[i] * vec[i];
    }
    return normVal;
}

float DotProduct(const float *a, const float *b, int N)
{
    float result = 0;
    for (int i = 0; i < N; ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

static const size_t VecSize = 50 * 50;
static size_t MatSize = VecSize * VecSize;
static int32_t processNum = 0;
static int32_t clusterSize = 0;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &clusterSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &processNum);

    char processorName[MPI_MAX_PROCESSOR_NAME] = {0};
    int32_t nameLength = 0;
    MPI_Get_processor_name(processorName, &nameLength);

    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processorName, processNum, clusterSize);

    const float tau = -0.05;
    const float epsilon = 0.001;
    //     const float tau = -0.05;
    // const float epsilon = 0.00001;

    std::vector<float> matA(MatSize);
    if (processNum == 0)
    {
        std::ifstream input("bin/matA.bin", std::ios::binary);
        if (!input.is_open())
        {
            // throw std::runtime_error("Failed to load data file.");
            exit(1);
        }
        input.read(reinterpret_cast<char *>(matA.data()), sizeof(float) * MatSize);
        input.close();
        if (matA.size() != MatSize)
        {
            printf("dangit!\n");
            return 1;
        }
    }

    std::vector<int32_t> adjustedVecSizes(clusterSize);
    std::vector<int32_t> adjustedVecOffsets(clusterSize);
    int32_t generalVecOffset = VecSize / clusterSize;
    int32_t additionalLines = VecSize % clusterSize;
    int32_t currentOffset = 0;
    for (int32_t i = 0; i < additionalLines; ++i)
    {
        int32_t givenSize = generalVecOffset + 1;
        adjustedVecSizes[i] = givenSize;
        adjustedVecOffsets[i] = currentOffset;
        currentOffset += givenSize;
    }

    for (int32_t i = additionalLines; i < clusterSize; ++i)
    {
        int32_t givenSize = generalVecOffset;
        adjustedVecSizes[i] = givenSize;
        adjustedVecOffsets[i] = currentOffset;
        currentOffset += givenSize;
    }

    std::vector<int32_t> adjustedMatOffsets(clusterSize);
    std::vector<int32_t> adjustedMatSizes(clusterSize);

    for (int32_t i = 0; i < adjustedMatSizes.size(); ++i)
    {
        adjustedMatSizes[i] = adjustedVecSizes[i] * VecSize;
    }
    for (int32_t i = 0; i < adjustedMatOffsets.size(); ++i)
    {
        adjustedMatOffsets[i] = adjustedVecOffsets[i] * VecSize;
    }

    std::vector<float> vecB(VecSize);

    std::ifstream inputVec("bin/vecB.bin", std::ios::binary);
    if (!inputVec.is_open())
    {
        // throw std::runtime_error("Failed to load data file.");
        exit(1);
    }
    inputVec.read(reinterpret_cast<char *>(vecB.data()), sizeof(float) * VecSize);
    inputVec.close();
    if (vecB.size() != VecSize)
    {
        printf("dangit!\n");
        return 0;
    }

    std::cout << "init finished" << std::endl;

    auto start = MPI_Wtime();
    std::vector<float> partMatA(adjustedMatSizes[processNum]);
    MPI_Scatterv(matA.data(), adjustedMatSizes.data(), adjustedMatOffsets.data(),
                 MPI_FLOAT, partMatA.data(), (int32_t)partMatA.size(),
                 MPI_FLOAT, 0, MPI_COMM_WORLD);

    std::cout << "scatter done" << std::endl;

    std::vector<float> vecX(VecSize, 0);

    float normVecB2 = NormOf2(vecB.data(), vecB.size());

    // std::vector<float> localAxb(adjustedVecOffsets[ProcessNum], 0);
    std::vector<float> localVecX(adjustedVecSizes[processNum], 0);

    size_t iterations = 0;
    // bool isDone = false;
    while (true)
    {
        float localNorm2 = 0.0f;
        for (int32_t i = 0; i < adjustedVecSizes[processNum]; ++i)
        {
            localVecX[i] = DotProduct(partMatA.data() + VecSize * i, vecX.data(), VecSize);
            localVecX[i] -= vecB[adjustedVecOffsets[processNum] + i];
            localNorm2 += localVecX[i] * localVecX[i];
            localVecX[i] *= -tau;
            localVecX[i] += vecX[adjustedVecOffsets[processNum] + i];
        }

        MPI_Allgatherv(localVecX.data(), adjustedVecSizes[processNum],
                       MPI_FLOAT, vecX.data(), adjustedVecSizes.data(),
                       adjustedVecOffsets.data(), MPI_FLOAT, MPI_COMM_WORLD);

        // float localepsilon = NormOf(localAxb.data(), localAxb.size()) * normVecB;

        float Norm2 = 0.0f;
        MPI_Allreduce(&localNorm2, &Norm2, 1, MPI_FLOAT, MPI_SUM,
                      MPI_COMM_WORLD);
        if (Norm2 < epsilon * epsilon * normVecB2)
            break;

        if (processNum == 0 && iterations % 50 == 0)
        {
            printf("Iterations: %zu\n", iterations + 1);
            printf("Norm2: %f\n", Norm2);
        }

        ++iterations;
    }

    if (processNum == 0)
    {
        printf("Done in %zu iterations\n", iterations);
    }

    double time = MPI_Wtime() - start;
    double maxTime = 0;
    MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (processNum == 0){
        printf("Time spent: %lf", maxTime);

        std::ofstream output("bin/out.bin", std::ios::binary);
        output.write(reinterpret_cast<const char *>(vecX.data()), sizeof(float) * vecX.size());
        output.close();
        std::ofstream timeF("time.txt");
        timeF << maxTime;
    }


    MPI_Finalize();

    return 0;
}
