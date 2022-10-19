#include <iostream>
#include <string>
#include "classes/grainsgrowth.cpp"


int main(void)
{   
    unsigned size = 12;
    unsigned seedsNum = 7;

    simula::GrainsGrowth grainsGrowth(size);
    grainsGrowth.makeSimulation(seedsNum);
    std::string matrixStr = grainsGrowth.gridToString();
    std::cout << matrixStr;
    
}