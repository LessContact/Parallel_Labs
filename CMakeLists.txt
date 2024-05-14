cmake_minimum_required(VERSION 3.27)
project(OPP)

set(CMAKE_CXX_STANDARD 20)

include_directories(lab2)

find_package(OpenMP REQUIRED)


add_executable(OPP2_base
        lab2/src/BasicIteration.cpp)
add_executable(OPP2_Each
        lab2/src/EachParallel.cpp)
add_executable(OPP2_onceAndMore
        lab2/src/PalrallelOnceAndSchedule.cpp)
add_executable(OPP2_Once
        lab2/src/PalrallelOnce.cpp)
add_executable(QSort lab2/ParalellQSort.cpp)
add_executable(qS lab2/qsort.cpp)

# uncomment to create random number file
# add_custom_command(TARGET OPP2_onceAndMore POST_BUILD WORKING_DIRECTORY "lab2"
#        COMMAND bash -c "shuf -i 0-268435456 -r -n 268435456 > Biggestdata.txt")

target_link_libraries(OPP2_base PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(OPP2_base PRIVATE ${OpenMP_CXX_FLAGS} -O3)
target_link_libraries(OPP2_Each PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(OPP2_Each PRIVATE ${OpenMP_CXX_FLAGS} -O3)
target_link_libraries(OPP2_onceAndMore PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(OPP2_onceAndMore PRIVATE ${OpenMP_CXX_FLAGS} -O3)
target_link_libraries(OPP2_Once PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(OPP2_Once PRIVATE ${OpenMP_CXX_FLAGS} -O3)
target_link_libraries(QSort PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(QSort PRIVATE ${OpenMP_CXX_FLAGS} -O3)
target_compile_options(qS PRIVATE -O3)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_DIRS})

add_executable(lab1 lab1/src/main.cpp)
target_link_libraries(lab1 ${MPI_CXX_LIBRARIES})

#set(CUSTOM_BUILD_DIR ${CMAKE_BINARY_DIR}/../lab3)
##set(EXECUTABLE_OUTPUT_PATH ./lab3)
#set(EXECUTABLE_OUTPUT_PATH ${CUSTOM_BUILD_DIR})

add_executable(lab3 lab3/lab3.cpp)
target_link_libraries(lab3 ${MPI_CXX_LIBRARIES})
target_compile_options(lab3 PRIVATE -O3)

add_executable(lab4 lab4/src/main.cpp)
target_link_libraries(lab4 ${MPI_CXX_LIBRARIES})
target_compile_options(lab4 PRIVATE -O3)
#target_link_libraries(lab1 ${MPI_C_LIBRARIES})

add_executable(lab5 lab5/src/main.cpp)
target_link_libraries(lab5 ${MPI_CXX_LIBRARIES})
target_compile_options(lab5 PRIVATE -O3)
