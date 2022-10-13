#include <iostream>
#include "utils/grainsgrowth.h"


int main(void)
{   
    int var = 12;
    //std::cout << "Input: ";
    //std::cin >> var;
    //std::cout << "var: " << var << " \n";

    GrainsGrowth gsg(var, var);
    
    return 0;
}