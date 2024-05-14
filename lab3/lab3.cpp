#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <array>
#include <eigen3/Eigen/Core>
#include <numeric>

#define X 0
//#define X 1
#define Y 1
//#define Y 0
#define VALUETYPE float

static int32_t ProcessNum = 0;
static int32_t ClusterSize = 0;

static const int32_t DimCount = 2;
static std::array<int32_t, DimCount> ProcessesInDim;
static std::array<int32_t, DimCount> ProcessCoordinates = {0};

static MPI_Comm GridCommunicator = {};
static MPI_Comm RowCommunicator = {};
static MPI_Comm ColCommunicator = {};

static const size_t n1 = 1600;
static const size_t n2 = 2500;
static const size_t n3 = 2800;
static const size_t p1 = 4;
static const size_t p2 = 4;

class Matrix{
public:
    Matrix() = default;
    Matrix(const int32_t rows, const int32_t columns) : Rows(rows), Columns(columns){
        Data.resize(rows * columns, 0);
    }
    ~Matrix() = default;

    void Dump(const std::string &outputFilePath) const {
        std::ofstream out(outputFilePath, std::ios::out);
        if (!out.is_open()) return;
        for (uint32_t row{}; row < Rows; ++row) {
            for (uint32_t col{}; col < Columns; ++col) {
                out << Data[row * Columns + col] << " ";
            }
            out << std::endl;
        }
        out.close();
    }

    Matrix MultiplyMatrixbyMatrix(Matrix &mat){
        Matrix matC(Rows, mat.Columns);

//        for (uint32_t rowA = 0; rowA < Rows; ++rowA) {
//            for (uint32_t colB = 0; colB < mat.Columns; ++colB) {
//                for (uint32_t colA = 0; colA < Columns; ++colA) {
//                    matC.Data[rowA * mat.Columns + colB] +=
//                                Data[rowA * Columns + colA] *
//                                        mat.Data[colA * mat.Columns + colB];
//                }
//            }
//        }

//        for (uint32_t i = 0; i < Rows; i++) {
//            for (uint32_t j = 0; j < mat.Rows; j++) {
//                for (uint32_t k = 0; k < mat.Columns; k++) {
//                    matC.Data[i * mat.Columns + k] += Data[i * Columns + j] * mat.Data[j * mat.Columns + k];
//                }
//            }
//        }

        for (uint32_t i = 0; i < Rows; i++) {
            for (uint32_t j = 0; j < mat.Columns; j++) {
                float sum = 0.0;
                for (uint32_t k = 0; k < Columns; k++)
                    sum = sum + Data[i * Columns + k] * mat.Data[k * mat.Columns + j];
                matC.Data[i * mat.Columns + j] = sum;
            }
        }

        return matC;
    }

    size_t Rows = 0;
    size_t Columns = 0;
    std::vector<VALUETYPE> Data = {};
};

bool AreMatricesEqual(const Matrix &lhs, const Matrix &rhs) {
    if (lhs.Rows != rhs.Rows || lhs.Columns != rhs.Columns) return false;
    for (uint32_t row{}; row < lhs.Rows; ++row) {
        for (uint32_t col{}; col < lhs.Columns; ++col)
            if (lhs.Data[row * lhs.Columns + col] != rhs.Data[row * lhs.Columns + col])return false;

    }
    return true;
}

void CreateCommunicators(){
    ProcessesInDim = {p1, p2};

    const int32_t periods[DimCount] = {0};
    const int32_t reorder = 1;  // mpi may reorder process ranks in grid.
                                // though it also says that they currently
                                // ignore reorder info????
    // global grid communicator
    MPI_Cart_create(MPI_COMM_WORLD, DimCount, ProcessesInDim.data(), periods, reorder, &GridCommunicator);
    // sub comm for X axis
    const int32_t remainDimsX[DimCount] = {0, 1}; // keep the connections along the Y
    MPI_Cart_sub(GridCommunicator, remainDimsX, &RowCommunicator);
    // sub comm for Y axis
    const int32_t remainDimsY[DimCount] = {1, 0}; // keep along X
    MPI_Cart_sub(GridCommunicator, remainDimsY, &ColCommunicator);
}
void DestroyCommunicators() {
    MPI_Comm_free(&RowCommunicator);
    MPI_Comm_free(&ColCommunicator);
    MPI_Comm_free(&GridCommunicator);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ClusterSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcessNum);

    char processorName[MPI_MAX_PROCESSOR_NAME] = {0};
    int32_t nameLength = 0;
    MPI_Get_processor_name(processorName, &nameLength);

    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processorName, ProcessNum, ClusterSize);

    CreateCommunicators();

    MPI_Cart_coords(GridCommunicator, ProcessNum, DimCount, ProcessCoordinates.data());
    printf("Proc %d has pos: (%d, %d)\n", ProcessNum, ProcessCoordinates[X], ProcessCoordinates[Y]);

    Matrix MatA(n1, n2);
    Matrix MatB(n2, n3);
    if (ProcessNum == 0) {
        std::ifstream input1("bin/MatA.bin", std::ios::binary);
        if (!input1.is_open()) {
            exit(1);
        }
        input1.read(reinterpret_cast<char *>(MatA.Data.data()), sizeof(VALUETYPE) * n1 * n2);
        input1.close();
//        printf("Mat A has been read.\n%p\n", MatA.Data.data());

        std::ifstream input2("bin/MatB.bin", std::ios::binary);
        if (!input2.is_open()) {
            exit(1);
        }
        input2.read(reinterpret_cast<char *>(MatB.Data.data()), sizeof(VALUETYPE) * n2 * n3);
        input2.close();
    } // mats init
    if(ProcessNum == 0) printf("mats are loaded\n");

//    std::vector<int32_t> rowsPerProcess(ProcessesInDim[X]);
//    for(size_t i = 0; i < rowsPerProcess.size(); ++i){
//        rowsPerProcess[i] = n1 / ProcessesInDim[X];
////        if(i < n1 % ProcessesInDim[X]) ++rowsPerProcess[i];
//    }
    int32_t rowsPerProcess = n1 / ProcessesInDim[X];
    printf("Proc %d has pos: (%d, %d)\n", ProcessNum, ProcessCoordinates[X], ProcessCoordinates[Y]);
    Matrix partMatA(rowsPerProcess, n2);
    { // sending mat A
        int32_t partMatASize = n2 * rowsPerProcess;
        if (ProcessCoordinates[Y] == 0) {
            MPI_Scatter(MatA.Data.data(), partMatASize, MPI_FLOAT,
                        partMatA.Data.data(), partMatASize, MPI_FLOAT, 0, ColCommunicator);
        }
        MPI_Bcast(partMatA.Data.data(), partMatASize, MPI_FLOAT, 0, RowCommunicator);
    }
    if(ProcessNum == 0) printf("mat A is sent\n");
    // mat A has been sent
    int32_t  colsPerProcess = n3 / ProcessesInDim[Y];
    Matrix partMatB (n2, colsPerProcess);
    { // sending mat B aka the "fun" part
        int32_t partMatBsize = n2 * colsPerProcess;
        if(ProcessCoordinates[X] == 0){
            MPI_Datatype columnType = 0, resizedColumnType = 0;
            MPI_Type_vector(n2, colsPerProcess, n3, MPI_FLOAT, &columnType);
//            MPI_Type_commit(&columnType); // is this actually needed?
            MPI_Type_create_resized(columnType, 0, partMatB.Columns * sizeof(VALUETYPE), &resizedColumnType);
            MPI_Type_commit(&resizedColumnType);

            MPI_Scatter(MatB.Data.data(), 1, resizedColumnType, partMatB.Data.data(),
                        partMatBsize, MPI_FLOAT, 0, RowCommunicator);

            MPI_Type_free(&columnType);
            MPI_Type_free(&resizedColumnType);
        }
        MPI_Bcast(partMatB.Data.data(), partMatBsize, MPI_FLOAT, 0, ColCommunicator);
    }
    if(ProcessNum == 0) printf("mat B is sent\n");
    //mat B has been sent
    auto start = MPI_Wtime();
    Matrix localMatC = partMatA.MultiplyMatrixbyMatrix(partMatB);
    if(ProcessNum == 0) printf("multiplication done\n");

    Matrix GatheredMatC(n1, n3);
    {
        MPI_Datatype blockType = 0, resizedBlockType = 0;
        auto partMatCSize = localMatC.Rows * localMatC.Columns;

        MPI_Type_vector(localMatC.Rows, localMatC.Columns, n3, MPI_FLOAT, &blockType);
//        MPI_Type_commit(&blockType);

        MPI_Type_create_resized(blockType, 0,
                                /*localMatC.Columns * */sizeof(VALUETYPE),
                                &resizedBlockType);
        MPI_Type_commit(&resizedBlockType);

//        std::vector<int32_t> displacements(ProcessesInDim[X] * ProcessesInDim[Y], 0);
//        int32_t i = 0;
//        for(int32_t row = 0; row < ProcessesInDim[X]; ++row){
//            for(int32_t col = 0; col < ProcessesInDim[Y]; ++col){
////                displacements[row * ProcessesInDim[X] + col] = row * ProcessesInDim[Y] * localMatC.Rows + col;
//                int32_t blockRow = i / ProcessesInDim[Y];
//                int32_t blockCol = i % ProcessesInDim[Y];
//                displacements[row * ProcessesInDim[X] + col] = blockRow * n3 * localMatC.Rows + blockCol * localMatC.Columns;
//                ++i;
//            }
//        }
        std::vector<int> offsets(ProcessesInDim[X] * ProcessesInDim[Y]);

        std::iota(offsets.begin(), offsets.end(), 0);
        std::transform(offsets.begin(), offsets.end(), offsets.begin(), [&](int num) {
            int blockRow = num / ProcessesInDim[Y];
            int blockCol = num % ProcessesInDim[Y];
            return blockRow * n3 * localMatC.Rows + blockCol * localMatC.Columns;
        });

        std::vector<int32_t> recvcounts(ProcessesInDim[X] * ProcessesInDim[Y], 1);
        MPI_Gatherv(localMatC.Data.data(), partMatCSize, MPI_FLOAT, GatheredMatC.Data.data(),
                    recvcounts.data(), offsets.data(), resizedBlockType, 0, MPI_COMM_WORLD);
        MPI_Type_free(&blockType);
        MPI_Type_free(&resizedBlockType);
    }
    // mat C gathered
    const double time = MPI_Wtime() - start;
    double maxTime = 0.0;
    MPI_Reduce(&time, &maxTime,
               1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (ProcessNum == 0) {
        printf("Time spent: %lf\n", maxTime);
        
//        Matrix MatCReal(n1, n3);
//        std::ifstream input3("bin/MatCReal.bin", std::ios::binary);
//        if (!input3.is_open()) {
//            exit(1);
//        }
//        input3.read(reinterpret_cast<char *>(MatCReal.Data.data()), sizeof(float) * n1 * n3);
//        input3.close();

//        using Mat = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
//        Mat matCorrect = Eigen::Map<const Mat>(MatCReal.Data.data(), n1, n3);
//        Mat matResultEigen = Eigen::Map<const Mat>(GatheredMatC.Data.data(), n1, n3);
        Matrix test = MatA.MultiplyMatrixbyMatrix(MatB);

        if(!AreMatricesEqual(GatheredMatC, test)){
//        if(!matResultEigen.isApprox(matCorrect)){
            printf("you've made a big mistake;\nyour Matrices are not equal!\n");
            std::ofstream output1("bin/MatCMine.bin", std::ios::binary);
            output1.write(reinterpret_cast<const char *>(GatheredMatC.Data.data()), sizeof(float) * GatheredMatC.Data.size());
            output1.close();
//            MatCReal.Dump("MatReal.txt");
//            GatheredMatC.Dump("CalcMat.txt");
        }
    }
    DestroyCommunicators();
    MPI_Finalize();
    return 0;
}