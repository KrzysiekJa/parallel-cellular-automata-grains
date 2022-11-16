#include <map>
#include <iostream>
#include <string>
#include <chrono>
#include "classes/grainsgrowth.cpp"
#include "utils/utils_functions.cpp"


int main(void)
{
    INT_TYPE size = 9;
    INT_TYPE nucleusNum = 6;
    INT_TYPE iterMC = 5;


    auto t_start = std::chrono::high_resolution_clock::now();
    simula::GrainsGrowth grainsGrowth(size, "moore", "periodic", "3D", "monte_carlo");
    auto t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double, std::milli > init_time = t_finish - t_start;
    std::cout << "Initialization time: " << init_time.count() << " ms\n";


    t_start = std::chrono::high_resolution_clock::now();
    grainsGrowth.makeSimulation(nucleusNum, iterMC);
    t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double, std::milli > simula_time = t_finish - t_start;
    std::cout << "Simulation time:     " << simula_time.count() << " ms\n";


    std::string matrixStr = grainsGrowth.spaceToDisplay();
    std::cout << matrixStr;


    INT_TYPE *** matrix = grainsGrowth.getSpace();
    SIZE_TYPE *** energySpace = grainsGrowth.getEnergySpace();
    std::map< std::string, SIZE_TYPE > dimMap { {"first", size}, {"second", size}, {"third", size} };


    t_start = std::chrono::high_resolution_clock::now();
    utils::saveSpaceToCSV( "test.csv", matrix, energySpace, dimMap , nucleusNum, simula_time.count() );
    t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double, std::milli > saving_time = t_finish - t_start;
    std::cout << "Time for saving:     " << saving_time.count() << " ms\n";
}
