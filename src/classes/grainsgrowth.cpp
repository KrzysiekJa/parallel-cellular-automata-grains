#include <ctime>
#include "../headers/grainsgrowth.hpp"


namespace simula{

    /************************************************
     *      Initial actions
     ***********************************************/

    GrainsGrowth::GrainsGrowth(unsigned long inputSize)
    {
        dimSize = inputSize;
        initGrid();
    }

    GrainsGrowth::~GrainsGrowth()
    {
    }

    void GrainsGrowth::initGrid()
    {
        grid = new unsigned long * [dimSize];
        
        for (int y = 0; y < dimSize; ++y)
        {
            grid[y] = new unsigned long [dimSize];
            
            for (int x = 0; x < dimSize; ++x)
            {
                grid[y][x] = 0;
            }
        }
    }

    /************************************************
     *      Simulation operations
     ***********************************************/

    unsigned long ** GrainsGrowth::makeSimulation(int numberOfSeeds)
    {
        putSeedsInGrid(numberOfSeeds);

        return grid;
    }

    void GrainsGrowth::putSeedsInGrid(int numberOfSeeds)
    {
        srand (time( NULL ));
        unsigned long x, y;

        for (int i = 1; i < numberOfSeeds + 1; ++i)
        {
            do
            {
                x = rand() % dimSize;
                y = rand() % dimSize;
            } while (grid[y][x] != 0);

            grid[y][x] = i;
        }
    }

    /************************************************
     *      Outbounde actions
     ***********************************************/

    std::string GrainsGrowth::gridToString()
    {
        std::string outputStr = "";

        for (int y = 0; y < dimSize; ++y)
        {            
            for (int x = 0; x < dimSize; ++x)
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
        return dimSize;
    }

    void GrainsGrowth::setSize(unsigned long newDimSize)
    {
        dimSize = newDimSize;
    }

    unsigned long ** GrainsGrowth::getGrid()
    {
        return grid;
    }

}