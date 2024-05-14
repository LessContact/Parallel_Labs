#include <iostream>
#include <limits>
#include <mpi.h>
#include <vector>
#include <cstring>

#define VALUE_AT(i, j, k) ((i) * N * N + (j) * N + (k))

static const uint32_t N = 768;
static int32_t ProcessNum = 0;
static int32_t ClusterSize = 0;
const float paramA = 10e5;
const float epsilon = 10e-8;

const int32_t Dx = 2, Dy = 2, Dz = 2; // cube sizes
const int32_t X0 = -1, Y0 = -1, Z0 = -1; // cube origin coordinates

constexpr float Hx = Dx / static_cast<float>(N - 1); // "stride" of the grid
constexpr float Hy = Dy / static_cast<float>(N - 1); // "stride" of the grid
constexpr float Hz = Dz / static_cast<float>(N - 1); // "stride" of the grid

inline float Phi(const float x, const float y, const float z) {
    return x * x + y * y + z * z;
}

inline float Ro(const float x, const float y, const float z) {
    return 6 - paramA * Phi(x, y, z);
}

// these calculate the coords of the discrete "nodes" we use for calculations
inline float GetCoordX(const uint32_t i) {
    return X0 + static_cast<float>(i) * Hx;
}

// these calculate the coords of the discrete "nodes" we use for calculations
inline float GetCoordY(const uint32_t j) {
    return Y0 + static_cast<float>(j) * Hy;
}

// these calculate the coords of the discrete "nodes" we use for calculations
inline float GetCoordZ(const uint32_t k) {
    return Z0 + static_cast<float>(k) * Hz;
}

float CalculateSolutionDelta(const std::vector<float> &grid) {
    float maxDelta = std::numeric_limits<float>::min();

    for (uint32_t i = 0; i < N; ++i) {
        const float x = GetCoordX(i);

        for (uint32_t j = 0; j < N; ++j) {
            const float y = GetCoordY(j);

            for (uint32_t k = 0; k < N; ++k) {
                const float z = GetCoordZ(k);

                const float currentDelta = std::abs(grid[VALUE_AT(i, j, k)] - Phi(x, y, z));
                maxDelta = std::max(maxDelta, currentDelta);
            }
        }
    }
    return maxDelta;
}

// в области Ω с краевыми условиями 1-го рода (т.е. на границе G известны значения искомой функции φ)
// это те самые значения
void InitializeGrid(std::vector<float> &grid, const uint32_t layerSize, const int32_t layerCoord) {
    for (uint32_t i = 0; i < layerSize + 2; ++i) {
        const float z = GetCoordZ(layerCoord + i);

        for (uint32_t j = 0; j < N; ++j) {
            const float x = GetCoordX(j);

            for (uint32_t k = 0; k < N; ++k) {
                const float y = GetCoordY(k);

                if (k != 0 && k != N - 1 && j != 0 && j != N - 1 && z != Z0 && z != Z0 + Dz) { // if not on the edge
                    grid[VALUE_AT(i, j, k)] = 0;  // body
                } else {
                    grid[VALUE_AT(i, j, k)] = Phi(x, y, z);  // edges
                }
            }
        }
    }
}

float UpdateLayerAndGetDelta(std::vector<float> &layer, std::vector<float> &layerNew, uint32_t layerIdx,
                             int32_t layerCoord) {
    const int32_t globalLayerCoord = layerCoord + static_cast<int32_t>(layerIdx);

    float maxDelta = std::numeric_limits<float>::min();

    if (globalLayerCoord == 0 || globalLayerCoord == N - 1)  // lower/upper border -- and that means we have to preserve it
    {
        const uint32_t layerDataOffset = layerIdx * N * N;
        memcpy(layerNew.data() + layerDataOffset, layer.data() + layerDataOffset, N * N * sizeof(layer[0]));
        return 0.0f;
    }

    const float z = GetCoordZ(globalLayerCoord);
    for (uint32_t i = 0; i < N; ++i) {
        const float x = GetCoordX(i);
        for (uint32_t j = 0; j < N; ++j) {
            const float y = GetCoordY(j);

            if (i == 0 || i == N - 1 || j == 0 || j == N - 1)  // side parts and front we have to preserve those
            {
                layerNew[VALUE_AT(layerIdx, i, j)] = layer[VALUE_AT(layerIdx, i, j)];
                continue;
            }

            const float invHx = 1.0f / Hx, invHy = 1.f / Hy, invHz = 1.f / Hz;
            const float invHx2 = invHx * invHx, invHy2 = invHy * invHy, invHz2 = invHz * invHz;
            const float multiplier = 1.0f / (2 * invHx2 + 2 * invHy2 + 2 * invHz2 + paramA);

            const float partX = invHx2 * (layer[VALUE_AT(layerIdx, i + 1, j)] + layer[VALUE_AT(layerIdx, i - 1, j)]);
            const float partY = invHy2 * (layer[VALUE_AT(layerIdx, i, j + 1)] + layer[VALUE_AT(layerIdx, i, j - 1)]);
            const float partZ = invHz2 * (layer[VALUE_AT(layerIdx + 1, i, j)] + layer[VALUE_AT(layerIdx - 1, i, j)]);
            layerNew[VALUE_AT(layerIdx, i, j)] = multiplier * (partX + partY + partZ - Ro(x, y, z));

            maxDelta = std::max(maxDelta,
                                std::abs(layerNew[VALUE_AT(layerIdx, i, j)] - layer[VALUE_AT(layerIdx, i, j)]));
        }
    }

    return maxDelta;
}

std::vector<float> Solve() {
    const uint32_t layerSize = N / ClusterSize;
    const int32_t layerCoord = ProcessNum * static_cast<int32_t>(layerSize) - 1;

    std::vector<float> layerBuffer(N * N * (layerSize + 2), 0);  // storing next state
    std::vector<float> layer(N * N * (layerSize + 2), 0);        // current state
    InitializeGrid(layer, layerSize, layerCoord);

    MPI_Request lowerSend = 0, lowerRecv = 0;
    MPI_Request upperSend = 0, upperRecv = 0;
    float globalDelta = std::numeric_limits<float>::max();
    uint32_t iterationCount = 0;
    while (globalDelta > epsilon) {
        ++iterationCount;

        float localDelta = 0.0f;
        float localMaxDelta = std::numeric_limits<float>::min();

        if (ProcessNum > 0)  // send lower stripe, receive upper
        {
            MPI_Isend(layer.data() + N * N, N * N, MPI_FLOAT, ProcessNum - 1,
                      80085, MPI_COMM_WORLD, &lowerSend); // a
            MPI_Irecv(layer.data(), N * N, MPI_FLOAT, ProcessNum - 1,
                      MPI_ANY_TAG, MPI_COMM_WORLD, &lowerRecv); // b
        }

        if (ProcessNum < ClusterSize - 1)  // send upper stripe, receive lower
        {
            MPI_Isend(layer.data() + N * N * layerSize, N * N, MPI_FLOAT, ProcessNum + 1,
                      80085, MPI_COMM_WORLD, &upperSend); // b
            MPI_Irecv(layer.data() + N * N * (layerSize + 1), N * N, MPI_FLOAT, ProcessNum + 1,
                      MPI_ANY_TAG, MPI_COMM_WORLD, &upperRecv); // a
        }

        // update insides w/o overlap
        for (uint32_t layerIdx = 2; layerIdx < layerSize; ++layerIdx) {
            localDelta = UpdateLayerAndGetDelta(layer, layerBuffer, layerIdx, layerCoord);
            localMaxDelta = std::max(localMaxDelta, localDelta);
        }

        MPI_Status stat = {};

        if (ProcessNum > 0) {
            MPI_Wait(&lowerRecv, &stat);  // Wait receive finish
        }

        localDelta = UpdateLayerAndGetDelta(layer, layerBuffer, 1, layerCoord);  // upper overlap
        localMaxDelta = std::max(localMaxDelta, localDelta);

        if (ProcessNum < ClusterSize - 1) {
            MPI_Wait(&upperRecv, &stat);  // Wait receive finish
        }

        localDelta = UpdateLayerAndGetDelta(layer, layerBuffer, layerSize, layerCoord);  //  lower overlap
        localMaxDelta = std::max(localMaxDelta, localDelta);

        if (ProcessNum > 0) MPI_Wait(&lowerSend, &stat); // Wait send finish
        if (ProcessNum < ClusterSize - 1) MPI_Wait(&upperSend, &stat);  // Wait send finish

        MPI_Allreduce(&localMaxDelta, &globalDelta, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

        if (ProcessNum == 0) {
            std::cout << "Allreduce delta: " << globalDelta << std::endl;
            std::cout << "Epsilon: " << epsilon << std::endl;
        }

        memcpy(layer.data(), layerBuffer.data(), layer.size() * sizeof(layer[0]));
    }
    if (ProcessNum == 0) {
        std::cout << "Iterations: " << iterationCount << std::endl;
    }

    std::vector<float> result = {};
    if (ProcessNum == 0) {
        result.resize(N * N * N, 0);
    }
    const auto recvCount = static_cast<int32_t>(layerSize * N * N);
    MPI_Gather(layer.data() + N * N, recvCount, MPI_FLOAT, result.data(), recvCount, MPI_FLOAT, 0, MPI_COMM_WORLD);

    return result;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ClusterSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcessNum);

    if (ProcessNum == 0) {
        std::cout << "ClusterSize:" << ClusterSize << std::endl;
    }

    if (N % ClusterSize != 0) {
        if (ProcessNum == 0) std::cout << "N % clusterSize != 0\naka i dont want to bother with non evenly divisible stuff\n";
        return 0;
    }

    const double startTime = MPI_Wtime();
    const std::vector<float> solution = Solve();
    const double diff = MPI_Wtime() - startTime;

    double maxTime = 0.0;
    MPI_Allreduce(&diff, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (ProcessNum == 0) {
        std::cout << "Max time: " << maxTime << " seconds.\n";
        std::cout << "Max Delta: " << CalculateSolutionDelta(solution) << std::endl;
    }

    MPI_Finalize();

    return 0;
}