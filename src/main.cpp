#include <iostream>
#include <string>
#include "classes/grainsgrowth.cpp"


int main(void)
{   
    INT_TYPE size = 15;
    INT_TYPE seedsNum = 9;

    simula::GrainsGrowth grainsGrowth(size, "moore", "periodic");
    grainsGrowth.makeSimulation(seedsNum);
    std::string matrixStr = grainsGrowth.gridToDisplay();
    std::cout << matrixStr;
}