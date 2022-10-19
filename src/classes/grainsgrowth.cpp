#include "../headers/grainsgrowth.hpp"


namespace simula{

    GrainsGrowth::GrainsGrowth(unsigned long inputSize)
    {
        size = inputSize;
        initGrid();
    }

    GrainsGrowth::~GrainsGrowth()
    {
    }

    void GrainsGrowth::initGrid()
    {
        grid = new unsigned long * [size];
        
        for (int y = 0; y < size; ++y)
        {
            grid[y] = new unsigned long [size];
            
            for (int x = 0; x < size; ++x)
            {
                grid[y][x] = 0;
            }
        }
    }

    int GrainsGrowth::makeSimulation()
    {
        return 0;
    }

    std::string GrainsGrowth::gridToString()
    {
        std::string outputStr = "";

        for (int y = 0; y < size; ++y)
        {            
            for (int x = 0; x < size; ++x)
            {
                outputStr += std::to_string(grid[y][x]);
                outputStr += " ";
            }
            outputStr += "\n";
        }
        return outputStr;
    }

    unsigned long GrainsGrowth::getSize()
    {
        return size;
    }

    void GrainsGrowth::setSize(unsigned long newSize)
    {
        size = newSize;
    }

    unsigned long ** GrainsGrowth::getGrid()
    {
        return grid;
    }

}