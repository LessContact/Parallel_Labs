#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <eigen3/Eigen/Core>

static const size_t n1 = 1600;
static const size_t n2 = 2500;
static const size_t n3 = 2800;

static std::random_device RandomDevice = {};
static std::mt19937 Gen(RandomDevice());

static const int32_t MatrixValueRange = 100;
static std::uniform_real_distribution<float> floatDistribution(
        -MatrixValueRange, MatrixValueRange);

#define VALUETYPE float

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

int main(){
    size_t rows = n1;
    size_t cols = n2;
    Matrix MatA(n1, n2);
//    MatA.resize(rows * cols);

    for (size_t row = 0; row < rows; ++row) {
        for (size_t col = 0; col < cols; ++col)
            MatA.Data[row * cols + col] = floatDistribution(Gen);
    }

    std::ofstream output1("bin/MatA.bin", std::ios::binary);
    if(!output1.is_open()){
        return 1;
    }
    output1.write(reinterpret_cast<const char *>(MatA.Data.data()), sizeof(float) * MatA.Data.size());
    output1.close();

    rows = n2;
    cols = n3;
//    std::vector<float> MatB;
//    MatB.resize(rows * cols);
Matrix MatB(n2, n3);
    for (size_t row = 0; row < rows; ++row) {
        for (size_t col = 0; col < cols; ++col)
            MatB.Data[row * cols + col] = floatDistribution(Gen);
    }
    std::ofstream output2("bin/MatB.bin", std::ios::binary);
    if(!output2.is_open()){
        return 1;
    }
    output2.write(reinterpret_cast<const char *>(MatB.Data.data()), sizeof(float) * MatA.Data.size());
    output2.close();

//    std::vector<float> matC(n1 * n3, 0);
Matrix MatC(n1, n3);
//    for (uint32_t rowA = 0; rowA < n1; ++rowA) {
//        for (uint32_t colB = 0; colB < n3; ++colB) {
//            for (uint32_t colA = 0; colA < n2; ++colA) {
//                matC[rowA * n3 + colB] +=
//                        MatA[rowA * n2 + colA] *
//                        MatB[colA * n3 + colB];
//            }
//        }
//    }

//    for (uint32_t i = 0; i < n1; i++) {
//        for (uint32_t j = 0; j < n3; j++) {
//            float sum = 0.0;
//            for (uint32_t k = 0; k < n2; k++)
//                sum = sum + Data[i * Columns + k] * mat.Data[k * mat.Columns + j];
//            matC.Data[i * mat.Columns + j] = sum;
//        }
//    }
    MatC = MatA.MultiplyMatrixbyMatrix(MatB);
//    using Mat = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
//
//    Mat matAEigen = Eigen::Map<const Mat>(MatA.data(), n1, n2);
//    Mat matBEigen = Eigen::Map<const Mat>(MatB.data(), n2, n3);
////    Mat matResultEigen = Eigen::Map<const Mat>(matC.data(), n1, n3);
//
//    Mat trustedResult = matAEigen * matBEigen;
//    trustedResult.data()

    std::ofstream output3("bin/MatCReal.bin", std::ios::binary);
    if(!output3.is_open()){
        return 1;
    }
    output3.write(reinterpret_cast<const char *>(MatC.Data.data()), sizeof(float) * MatC.Data.size());
//    output3.write(reinterpret_cast<const char *>(trustedResult.data()), sizeof(float) * trustedResult.size());
    output3.close();


}