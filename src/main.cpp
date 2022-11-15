#include <iostream>
#include <string>
#include "classes/grainsgrowth.cpp"


int main(void)
{   
    INT_TYPE size = 9;
    INT_TYPE nucleusNum = 6;
    INT_TYPE iterMC = 5;

    simula::GrainsGrowth grainsGrowth(size, "moore", "periodic", "3D", "monte_carlo");
    grainsGrowth.makeSimulation(nucleusNum, iterMC);
    std::string matrixStr = grainsGrowth.spaceToDisplay();
    std::cout << matrixStr;
}
