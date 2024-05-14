#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <chrono>

#define SIZE 268435456
#define THOLD 1000

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
#pragma omp task if (right - low >= THOLD)
        quickSort(array, low, right);
    if (left < high)
#pragma omp task if (high - left >= THOLD)
        quickSort(array, left, high);
}

int main() {
    for(int i = 1; i < 20; i++) {
        omp_set_num_threads(i);
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
        auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel
        {
#pragma omp single
            quickSort(numbers.data(), 0, numbers.size() - 1);
        }
        auto end = std::chrono::high_resolution_clock::now();
//        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        auto elapsed = std::chrono::duration<double, std::milli> (end - start);
//        std::cout << elapsed.count() << std::endl;
        printf("%lf\n", elapsed.count()/1000);
//        printf("%ld\n", duration.count());
//        printf("%lf", elapsed.count()/1000);
//        std::cout << duration.count() << std::endl;
        // Output sorted numbers
        std::ofstream outfile("/home/less/Documents/OPP/src/sorted.txt");
//    for (int64_t num : numbers) {
//        outfile << num << std::endl;
//    }

        if (!std::is_sorted(numbers.begin(), numbers.end())) {
            std::cout << "zamn" << std::endl;
            return 1;
        }
    }
    return 0;
}
