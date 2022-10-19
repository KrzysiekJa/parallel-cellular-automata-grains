#include <iostream>
#include <string>
#include "classes/grainsgrowth.cpp"


int main(void)
{   
    unsigned size = 12;

    simula::GrainsGrowth grainsGrowth(size);
    std::string matrixStr = grainsGrowth.gridToString();
    std::cout << matrixStr;
    
}