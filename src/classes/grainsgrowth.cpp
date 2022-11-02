#include <ctime>
#include <algorithm>
#include "../headers/grainsgrowth.hpp"


namespace simula{

    /************************************************
     *      Initial actions
     ***********************************************/

    GrainsGrowth::GrainsGrowth(
        unsigned inputSize, 
        std::string inCA_NeighborhoodType, 
        std::string inBoundaryConditionsType
        )
    {
        dimSize = inputSize;
        checkInput_NeighborhoodType(inCA_NeighborhoodType);
        checkInput_BoundaryConditionsType(inBoundaryConditionsType);
        setNeighborhoodVector();
        grid = initGrid(grid);
        nextGrid = initGrid(nextGrid);
    }

    GrainsGrowth::~GrainsGrowth()
    {
    }

    void GrainsGrowth::checkInput_NeighborhoodType( std::string inCA_NeighborhoodType )
    {
        if (inCA_NeighborhoodType == "von_neumann" || 
                inCA_NeighborhoodType == "moore"
            )
        {
            CA_NeighborhoodType = inCA_NeighborhoodType;
        }
        else
        {
            CA_NeighborhoodType = "von_neumann";
        }
    }

    void GrainsGrowth::checkInput_BoundaryConditionsType( std::string inBoundaryConditionsType )
    {
        if (inBoundaryConditionsType == "absorbing" || 
                inBoundaryConditionsType == "periodic"
            )
        {
            boundaryConditionsType = inBoundaryConditionsType;
        }
        else
        {
            boundaryConditionsType = "absorbing";
        }
    }

    void GrainsGrowth::setNeighborhoodVector()
    {
        if(CA_NeighborhoodType == "von_neumann")
        {
            neighborhoodVector = {{0, -1}, {-1, 0}, {0, 1}, {1, 0}};
        }
        if(CA_NeighborhoodType == "moore")
        {
            neighborhoodVector = {{0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}};
        }
    }

    INT_TYPE ** GrainsGrowth::initGrid( INT_TYPE ** gridObject )
    {
        gridObject = new INT_TYPE * [dimSize];
        
        for (SIZE_TYPE x = 0; x < dimSize; ++x)
        {
            gridObject[x] = new INT_TYPE [dimSize];
            memset( gridObject[x], 0, dimSize * sizeof(INT_TYPE) );
        }
        return gridObject;
    }

    /************************************************
     *      Simulation operations
     ***********************************************/

    INT_TYPE ** GrainsGrowth::makeSimulation( unsigned numberOfSeeds )
    {
        bool change_happened = true;
        
        putSeedsInGrid(numberOfSeeds);
        nextGrid = copyGrid(grid);

        while (change_happened)
        {
            change_happened = makeIterationStep();
        }
        return nextGrid;
    }

    void GrainsGrowth::putSeedsInGrid( unsigned numberOfSeeds )
    {
        srand ( time( NULL ));
        SIZE_TYPE x, y;

        for (unsigned i = 1; i < numberOfSeeds + 1; ++i)
        {
            do
            {
                x = rand() % dimSize;
                y = rand() % dimSize;
            } while (grid[x][y] != 0);

            grid[x][y] = i;
        }
    }

    INT_TYPE ** GrainsGrowth::copyGrid( INT_TYPE ** modelGrid )
    {
        INT_TYPE ** targetGrid = new INT_TYPE * [dimSize];

        for (SIZE_TYPE x = 0; x < dimSize; ++x){
            targetGrid[x] = new INT_TYPE[dimSize];
            std::copy( &modelGrid[x][0],  &modelGrid[x][dimSize],  targetGrid[x] );
        }

        return targetGrid;
    }

    bool GrainsGrowth::makeIterationStep()
    {   
        bool change_happened = false;

        change_happened = iterateOverGrid();

        grid = copyGrid(nextGrid);

        return change_happened;
    }

    bool GrainsGrowth::iterateOverGrid()
    {
        bool change_happened = false;
        for (int x = 0; x < dimSize; ++x)
        {
            for (int y = 0; y < dimSize; ++y)
            {
                change_happened = checkIfCell_IdEqual_0(x, y) || change_happened;
                // if change happened variable keeps the state of truth
                std::cout << grid[x][y] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        return change_happened;
    }

    /************************************************
     *      Simulation operations -> detailed
     ***********************************************/

    bool GrainsGrowth::checkIfCell_IdEqual_0( INT_TYPE cell_x, INT_TYPE cell_y )
    {
        if (grid[cell_x][cell_y] == 0)
        {
            checkCell_NeighborhoodIds(cell_x, cell_y);
            return true;
        }
        else
        {
            return false;
        }
    }

    void GrainsGrowth::checkCell_NeighborhoodIds( INT_TYPE cell_x, INT_TYPE cell_y )
    {
        bool skip_process = false;
        for (auto iter = neighborhoodVector.begin(); iter != neighborhoodVector.end(); ++iter)
        {
            skip_process = inspectNeighborhoodCell(cell_x, cell_y, * iter);
            if ( skip_process )
            {
                break;
            }
        }
    }

    bool GrainsGrowth::inspectNeighborhoodCell( 
        INT_TYPE cell_x, INT_TYPE cell_y, std::tuple< short, short > iter
        )
    {
        long neigh_cell_x = cell_x + std::get<0>(iter);
        long neigh_cell_y = cell_y + std::get<1>(iter);

        if( neigh_cell_x >= 0 && neigh_cell_x < dimSize && 
            neigh_cell_y >= 0 && neigh_cell_y < dimSize )
        {
            if ( grid[neigh_cell_x][neigh_cell_y] != 0 )
            {
                // found cell in neighborhood with id =/= 0
                performCell_IdChange(cell_x, cell_y);
                return true;
            }
        }
        return false;
    }

    void GrainsGrowth::performCell_IdChange( INT_TYPE cell_x, INT_TYPE cell_y )
    {
        std::map< INT_TYPE, unsigned > neighborhoodMap = countNeighborhood(cell_x, cell_y);
        INT_TYPE newId = getNeighborhood_MaxFranction(neighborhoodMap);
        nextGrid[cell_x][cell_y] = newId;
    }

    std::map< INT_TYPE, unsigned > GrainsGrowth::countNeighborhood(
        INT_TYPE cell_x, INT_TYPE cell_y
        )
    {
        std::map<INT_TYPE, unsigned> rMap;
    
        for (auto iter = neighborhoodVector.begin(); iter != neighborhoodVector.end(); ++iter)
        {
            long neigh_cell_x = cell_x + std::get<0>( * iter);
            long neigh_cell_y = cell_y + std::get<1>( * iter);

            rMap = checkCellDuringCounting(neigh_cell_x, neigh_cell_y, rMap);
        }
        return rMap;
    }

    std::map< INT_TYPE, unsigned > GrainsGrowth::checkCellDuringCounting(
        long neigh_cell_x, long neigh_cell_y, std::map< INT_TYPE, unsigned > valMap
        )
    {
        if( neigh_cell_x >= 0 && neigh_cell_x < dimSize && 
            neigh_cell_y >= 0 && neigh_cell_y < dimSize )
        {
            INT_TYPE value = grid[neigh_cell_x][neigh_cell_y];
            if (value != 0)
            {
                valMap[ value ]++;
            }
        }
        else 
        {   // specially for periodic boundary conditions case
            if (boundaryConditionsType == "periodic" && 
                grid[ mapIfPeriodic(neigh_cell_x) ][ mapIfPeriodic(neigh_cell_y) ] != 0 )
            {
                INT_TYPE value = grid[ mapIfPeriodic(neigh_cell_x) ][ mapIfPeriodic(neigh_cell_y) ];
                valMap[ value ]++;
            }
        }
        return valMap;
    }

    INT_TYPE GrainsGrowth::getNeighborhood_MaxFranction( std::map< INT_TYPE, unsigned > neighborhoodMap )
    {
        std::map< INT_TYPE, unsigned >::iterator  dominant = std::max_element(
            std::begin(neighborhoodMap), std::end(neighborhoodMap),
            [] (const std::pair<INT_TYPE, unsigned> & pair1, const std::pair<INT_TYPE, unsigned> & pair2) 
            {
                return pair1.second < pair2.second;
            }
        );
        return dominant -> first;
    }

    INT_TYPE GrainsGrowth::mapIfPeriodic( long location )
    {
        if (location < 0) {
            return dimSize - 1;
        }
        if (location >= dimSize) {
            return 0;
        }
        else {
            return location;
        }
    }

    /************************************************
     *      Outbounde actions
     ***********************************************/

    std::string GrainsGrowth::gridToDisplay()
    {
        std::string outputStr = "";

        for (int x = 0; x < dimSize; ++x)
        {
            for (int y = 0; y < dimSize; ++y)
            {
                outputStr += std::to_string( nextGrid[x][y] );
                outputStr += " ";
            }
            outputStr += "\n";
        }
        return outputStr;
    }

    /************************************************
     *      Getters, setters
     ***********************************************/

    SIZE_TYPE GrainsGrowth::getSize()
    {
        return dimSize;
    }

    void GrainsGrowth::setSize( SIZE_TYPE newDimSize )
    {
        dimSize = newDimSize;
        grid = initGrid(grid);
        nextGrid = initGrid(nextGrid);
    }

    INT_TYPE ** GrainsGrowth::getGrid()
    {
        return nextGrid;
    }

}