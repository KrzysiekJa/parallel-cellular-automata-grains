#include <iostream>
#include <string>
#include "classes/grainsgrowth.cpp"


int main(void)
{   
    INT_TYPE size = 15;
    INT_TYPE nucleusNum = 9;
    INT_TYPE iterMC = 7;

    simula::GrainsGrowth grainsGrowth(size, "moore", "periodic", "monte_carlo");
    grainsGrowth.makeSimulation(nucleusNum, iterMC);
    std::string matrixStr = grainsGrowth.gridToDisplay();
    std::cout << matrixStr;
}
