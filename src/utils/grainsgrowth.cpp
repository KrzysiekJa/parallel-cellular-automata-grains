#include "grainsgrowth.h"


GrainsGrowth::GrainsGrowth(unsigned long inputX, unsigned long inputY)
{
    X = inputX;
    Y = inputY;

    grid = new unsigned long * [inputY];
    
    for (int y = 0; y < inputY; ++y)
    {
        grid[y] = new unsigned long [inputX];
        
        for (int x = 0; x < inputX; ++x)
        {
            grid[y][x] = 0;
        }
    }
}

GrainsGrowth::~GrainsGrowth()
{
}

int GrainsGrowth::makeSimulation()
{

}
