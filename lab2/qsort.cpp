#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
//#include <omp.h>
#define SIZE 33554432

void quickSort(int64_t */*std::vector<int64_t>& */array, int64_t low, int64_t high)
{
    int64_t left = low;
    int64_t right = high;
    int64_t pivot = array[(left + right) / 2];

    while (left <= right)
    {
        while (array[left] < pivot)
            left++;
        while (array[right] > pivot)
            right--;
        if (left <= right)
        {
            std::swap(array[left], array[right]);
            left++;
            right--;
        }
    }
    if (right > low)
//#pragma omp task if (right - low >= 1000)
        quickSort(array, low, right);
    if (left < high)
//#pragma omp task if (high - left >= 1000)
        quickSort(array, left, high);
}

int main() {
//    omp_set_num_threads(16);
    std::ifstream inputFile("/home/less/Documents/OPP/src/Biggestdata.txt");
    if (!inputFile.is_open()) {
        std::cout << "Unable to open file." << std::endl;
        return 1;
    }

    std::vector<int64_t> numbers;
    numbers.reserve(SIZE * sizeof(int64_t));
//    numbers.reserve(2097152 * sizeof(int64_t));
    int64_t number;
    while (inputFile >> number) {
        numbers.push_back(number);
    }

    inputFile.close();
    std::chrono::time_point<std::chrono::high_resolution_clock> start;

        quickSort(numbers.data(), 0, numbers.size() - 1);

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout<< std::fixed << duration.count() << std::endl;
    // Output sorted numbers
    std::ofstream outfile("/home/less/Documents/OPP/src/sorted.txt");
//    for (int64_t num : numbers) {
//        outfile << num << std::endl;
//    }

    if(!std::is_sorted(numbers.begin(), numbers.end())){
        std::cout << "zamn" << std::endl;
        return 1;
    }

    return 0;
}
