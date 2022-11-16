#include <ctime>
#include <random>
#include <algorithm>
#include "../headers/grainsgrowth.hpp"


namespace simula{

    /************************************************
     *      Initial actions
     ***********************************************/

    GrainsGrowth::GrainsGrowth(
        unsigned inputSize, 
        std::string inCA_NeighborhoodType, 
        std::string inBoundaryConditionsType,
        std::string inDimensionalityType,
        std::string inSimulationType
        )
    {
        dimSize = inputSize;
        checkInput_NeighborhoodType(inCA_NeighborhoodType);
        checkInput_BoundaryConditionsType(inBoundaryConditionsType);
        checkInput_DimensionalityType(inDimensionalityType);
        checkInput_SimulationType(inSimulationType);
        setThirdDim();
        setNeighborhoodVector();
        setMonteCarloVector();
        space = initSpace(space);
        nextSpace = initSpace(nextSpace);
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

    void GrainsGrowth::checkInput_DimensionalityType( std::string inDimensionalityType )
    {
        if (inDimensionalityType == "2D" || 
                inDimensionalityType == "3D"
            )
        {
            dimensionalityType = inDimensionalityType;
        }
        else
        {
            dimensionalityType = "2D";
        }
    }

    void GrainsGrowth::checkInput_SimulationType( std::string inSimulationType )
    {
        if (inSimulationType == "none" || 
            inSimulationType == "monte_carlo"
            )
        {
            simulationType = inSimulationType;
        }
        else
        {
            CA_NeighborhoodType = "none";
        }
    }

    void GrainsGrowth::setThirdDim()
    {
        if (dimensionalityType == "2D")
        {
            thirdDim = 1;
        }
        if (dimensionalityType == "3D")
        {
            thirdDim = dimSize;
        }
        else
        {
            dimensionalityType = 1;
        }
    }

    void GrainsGrowth::setNeighborhoodVector()
    {
        if(CA_NeighborhoodType == "von_neumann")
        {
            if (dimensionalityType == "2D")
            {
                neighborhoodVector = {{0, -1, 0}, {-1, 0, 0}, {0, 1, 0}, {1, 0, 0}};
            }
            if (dimensionalityType == "3D")
            {
                neighborhoodVector = {{0, -1, 0}, {-1, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}, {0, 0, -1}};
            }
        }
        if(CA_NeighborhoodType == "moore")
        {
            if (dimensionalityType == "2D")
            {
                neighborhoodVector = {{0, -1, 0}, {-1, -1, 0}, {-1, 0, 0}, {-1, 1, 0}, 
                                        {0, 1, 0}, {1, 1, 0}, {1, 0, 0}, {1, -1, 0}};
            }
            if (dimensionalityType == "3D")
            {
                neighborhoodVector = {{-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, 
                                        {-1, 0, 1}, {-1, 1, -1}, {-1, 1, 0}, {-1, 1, 1}, {0, -1, -1}, 
                                        {0, -1, 0}, {0, -1, 1}, {0, 0, -1}, {0, 0, 1}, {0, 1, -1}, 
                                        {0, 1, 0}, {0, 1, 1}, {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, 
                                        {1, 0, -1}, {1, 0, 0}, {1, 0, 1}, {1, 1, -1}, {1, 1, 0}, {1, 1, 1}};
            }
        }
    }

    void GrainsGrowth::setMonteCarloVector()
    {
        for (SIZE_TYPE x = 0; x < dimSize; ++x)
        {
            for (SIZE_TYPE y = 0; y < dimSize; ++y)
            {
                for (SIZE_TYPE z = 0; z < thirdDim; ++z)
                {
                    monteCarloVector.push_back( { x, y, z} );
                }
            }
        }
    }

    INT_TYPE *** GrainsGrowth::initSpace( INT_TYPE *** spaceObject )
    {
        spaceObject = new INT_TYPE ** [dimSize];
        
        for (SIZE_TYPE x = 0; x < dimSize; ++x)
        {
            spaceObject[x] = new INT_TYPE * [dimSize];
            for (SIZE_TYPE y = 0; y < dimSize; ++y)
            {
                spaceObject[x][y] = new INT_TYPE [thirdDim];
                for (SIZE_TYPE z = 0; z < thirdDim; ++z)
                {
                    spaceObject[x][y][z] = 0;
                }
                // didn't worked: memset( spaceObject[x][y], 0, thirdDim * sizeof(INT_TYPE) );
            }
        }
        return spaceObject;
    }

    /************************************************
     *      Simulation preparation
     ***********************************************/

    INT_TYPE *** GrainsGrowth::makeSimulation( unsigned numberOfNucleus, unsigned numberOfMCExecutions )
    {        
        createBasicStructure(numberOfNucleus); // Cellular automata
        //std::cout << spaceToDisplay() << "\n";

        if (simulationType == "monte_carlo") 
        {
            numberOfMCExecutions = numberOfMCExecutions > 0 ? numberOfMCExecutions : numberOfNucleus;
            makeMonteCarloSimulation(numberOfMCExecutions);
        }
        return nextSpace;
    }

    INT_TYPE *** GrainsGrowth::copySpace( INT_TYPE *** modelSpace )
    {
        INT_TYPE *** targetSpace = new INT_TYPE ** [dimSize];

        for (SIZE_TYPE x = 0; x < dimSize; ++x)
        {
            targetSpace[x] = new INT_TYPE * [dimSize];

            for (SIZE_TYPE y = 0; y < dimSize; ++y)
            {
                targetSpace[x][y] = new INT_TYPE [dimSize];
                std::copy( &modelSpace[x][y][0],  &modelSpace[x][y][thirdDim],  targetSpace[x][y] );
            }
        }
        return targetSpace;
    }

    void GrainsGrowth::createBasicStructure( unsigned numberOfNucleus)
    {
        bool change_happened_iter = true;

        putNucleusInSpace(numberOfNucleus);

        while (change_happened_iter)
        {
            change_happened_iter = iterateOverSpace();
            space = copySpace(nextSpace);
        }
    }

    void GrainsGrowth::putNucleusInSpace( unsigned numberOfNucleus )
    {
        srand ( time( NULL ));
        SIZE_TYPE x, y, z;

        for (unsigned i = 1; i < numberOfNucleus + 1; ++i)
        {
            do
            {
                x = rand() % dimSize;
                y = rand() % dimSize;
                z = rand() % thirdDim;
            } while (space[x][y][z] != 0);

            space[x][y][z] = i;
        }
    }

    bool GrainsGrowth::iterateOverSpace()
    {
        bool zero_hit = false;
        for (SIZE_TYPE x = 0; x < dimSize; ++x)
        {
            for (SIZE_TYPE y = 0; y < dimSize; ++y)
            {
                for (SIZE_TYPE z = 0; z < thirdDim; ++z)
                {
                    zero_hit = checkIfCell_IdEqual_0(x, y, z) || zero_hit;
                // if change happened variable keeps the state of truth
                }
            }
        }
        return zero_hit;
    }

    /************************************************
     *      Cellular automata operations
     ***********************************************/

    bool GrainsGrowth::checkIfCell_IdEqual_0( INT_TYPE cell_x, INT_TYPE cell_y , INT_TYPE cell_z)
    {
        if (space[cell_x][cell_y][cell_z] == 0)
        {
            checkCell_NeighborhoodIds(cell_x, cell_y, cell_z);
            return true;
        }
        else
        {
            return false;
        }
    }

    void GrainsGrowth::checkCell_NeighborhoodIds( INT_TYPE cell_x, INT_TYPE cell_y, INT_TYPE cell_z )
    {
        bool skip_process = false;
        for (auto iter = neighborhoodVector.begin(); iter != neighborhoodVector.end(); ++iter)
        {
            skip_process = inspectNeighborhoodCell(cell_x, cell_y, cell_z, * iter);
            if ( skip_process )
            {
                break;
                // to reduce repetition of operation if in neighborhoob cell =/= 0
            }
        }
    }

    bool GrainsGrowth::inspectNeighborhoodCell( 
        INT_TYPE cell_x, INT_TYPE cell_y, INT_TYPE cell_z, std::tuple< short, short, short > iter
        )
    {
        long neigh_cell_x = cell_x + std::get<0>(iter);
        long neigh_cell_y = cell_y + std::get<1>(iter);
        long neigh_cell_z = cell_z + std::get<2>(iter);

        if( neigh_cell_x >= 0 && neigh_cell_x < dimSize && 
            neigh_cell_y >= 0 && neigh_cell_y < dimSize && 
            neigh_cell_z >= 0 && neigh_cell_z < thirdDim )
        {
            if ( space[neigh_cell_x][neigh_cell_y][neigh_cell_z] != 0 )
            {
                // found cell in neighborhood with id =/= 0
                performCell_IdChange(cell_x, cell_y, cell_z);
                return true;
            }
        }
        return false;
    }

    void GrainsGrowth::performCell_IdChange( INT_TYPE cell_x, INT_TYPE cell_y, INT_TYPE cell_z )
    {
        std::map< INT_TYPE, unsigned > neighborhoodMap = countNeighborhood(cell_x, cell_y, cell_z);
        INT_TYPE newId = getNeighborhood_MaxFranction(neighborhoodMap);
        nextSpace[cell_x][cell_y][cell_z] = newId;
    }

    std::map< INT_TYPE, unsigned > GrainsGrowth::countNeighborhood(
        INT_TYPE cell_x, INT_TYPE cell_y, INT_TYPE cell_z
        )
    {
        std::map<INT_TYPE, unsigned> rMap;
    
        for (auto iter = neighborhoodVector.begin(); iter != neighborhoodVector.end(); ++iter)
        {
            long neigh_cell_x = cell_x + std::get<0>( * iter);
            long neigh_cell_y = cell_y + std::get<1>( * iter);
            long neigh_cell_z = cell_z + std::get<2>( * iter);

            rMap = checkCellDuringCounting(neigh_cell_x, neigh_cell_y, neigh_cell_z, rMap);
        }
        return rMap;
    }

    std::map< INT_TYPE, unsigned > GrainsGrowth::checkCellDuringCounting(
        long neigh_cell_x, long neigh_cell_y, long neigh_cell_z, std::map< INT_TYPE, unsigned > valMap
        )
    {
        if( neigh_cell_x >= 0 && neigh_cell_x < dimSize && 
            neigh_cell_y >= 0 && neigh_cell_y < dimSize && 
            neigh_cell_z >= 0 && neigh_cell_z < thirdDim )
        {
            INT_TYPE value = space[neigh_cell_x][neigh_cell_y][neigh_cell_z];
            if (value != 0)
            {
                valMap[ value ]++;
            }
        }
        else 
        {   // specially for periodic boundary conditions case
            if (boundaryConditionsType == "periodic" && 
                space[ mapIfPeriodic(neigh_cell_x) ][ mapIfPeriodic(neigh_cell_y) ][ mapIfPeriodic(neigh_cell_z) ] != 0 )
            {
                INT_TYPE value = space[ mapIfPeriodic(neigh_cell_x) ][ mapIfPeriodic(neigh_cell_y) ][ mapIfPeriodic(neigh_cell_z) ];
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
     *      Monte Carlo Simulation
     ***********************************************/

    void GrainsGrowth::makeMonteCarloSimulation( unsigned numberOfIter )
    {
        while (numberOfIter)
        {
            shuffleMonteCarloVector();
            iterateOverMCVector();
            space = copySpace(nextSpace);

            numberOfIter--;
        }
    }

    void GrainsGrowth::shuffleMonteCarloVector()
    {
        std::default_random_engine random_engine = std::default_random_engine {};
        std::shuffle(std::begin(monteCarloVector), std::end(monteCarloVector), random_engine);
    }

    void GrainsGrowth::iterateOverMCVector()
    {
        std::vector< std::tuple< SIZE_TYPE, SIZE_TYPE, SIZE_TYPE >> iterVector = monteCarloVector;

        for (auto iter = iterVector.begin(); iter != iterVector.end(); )
        {
            SIZE_TYPE x = std::get<0>( * iter);
            SIZE_TYPE y = std::get<1>( * iter);
            SIZE_TYPE z = std::get<2>( * iter);

            inspectCellEnergy(x, y, z);

            iter = iterVector.erase(iter); // instrad of iter++
        }
    }

    void GrainsGrowth::inspectCellEnergy( INT_TYPE cell_x, INT_TYPE cell_y, INT_TYPE cell_z )
    {
        std::map< INT_TYPE, unsigned > neighborhoodMap = countNeighborhood(cell_x, cell_y, cell_z);
        unsigned initEnergy = calculateCellEnergy( neighborhoodMap , space[cell_x][cell_y][cell_z] );
        
        if (initEnergy > 0)
        {
            checkCellEnergeticStatus(neighborhoodMap, initEnergy, cell_x, cell_y, cell_z);
        }
    }

    void GrainsGrowth::checkCellEnergeticStatus( 
        std::map< INT_TYPE, unsigned > inMap, unsigned initEnergy, INT_TYPE cell_x, INT_TYPE cell_y, INT_TYPE cell_z 
        )
    {
        INT_TYPE maxFranctionId = getNeighborhood_MaxFranction(inMap);
        unsigned maxFranctionEnergy = calculateCellEnergy( inMap , maxFranctionId );
        if (initEnergy > maxFranctionEnergy)
        {
            nextSpace[cell_x][cell_y][cell_z] = maxFranctionId;
        }
    }

    unsigned GrainsGrowth::calculateCellEnergy( std::map< INT_TYPE, unsigned > inMAp , INT_TYPE id )
    {
        inMAp.erase( inMAp.find( id ) );
        // count all without indicated id
        unsigned sum = std::accumulate ( 
            inMAp.begin(), inMAp.end(), 0, 
            [] ( unsigned acc, std::pair< INT_TYPE, unsigned > p ) 
            {
                return ( acc + p.second ); 
            } 
        );
        return sum;
    }

    /************************************************
     *      Outbounde actions
     ***********************************************/

    std::string GrainsGrowth::spaceToDisplay()
    {
        std::string outputStr = "";

        for (int x = 0; x < dimSize; ++x)
        {
            for (int y = 0; y < dimSize; ++y)
            {
                for (int z = 0; z < thirdDim; ++z)
                {
                    outputStr += std::to_string( nextSpace[x][y][z] );
                    outputStr += " ";
                }
                outputStr += "\n";
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
        setThirdDim();
        space = initSpace(space);
        nextSpace = initSpace(nextSpace);
        setMonteCarloVector();
    }

    INT_TYPE *** GrainsGrowth::getSpace()
    {
        return nextSpace;
    }

    SIZE_TYPE *** GrainsGrowth::createEnergySpace() // private method
    {
        SIZE_TYPE *** spaceObject = new SIZE_TYPE ** [dimSize];
        
        for (SIZE_TYPE x = 0; x < dimSize; ++x)
        {
            spaceObject[x] = new INT_TYPE * [dimSize];
            for (SIZE_TYPE y = 0; y < dimSize; ++y)
            {
                spaceObject[x][y] = new INT_TYPE [thirdDim];
                for (SIZE_TYPE z = 0; z < thirdDim; ++z)
                {
                    spaceObject[x][y][z] = 0;
                }
            }
        }
        return spaceObject;
    }

    SIZE_TYPE *** GrainsGrowth::getEnergySpace()
    {
        SIZE_TYPE *** energySpace = createEnergySpace();

        for (SIZE_TYPE x = 0; x < dimSize; ++x)
        {
            for (SIZE_TYPE y = 0; y < dimSize; ++y)
            {
                for (SIZE_TYPE z = 0; z < thirdDim; ++z)
                {
                    std::map< INT_TYPE, unsigned > neighborhoodMap = countNeighborhood(x, y, z);
                    energySpace[x][y][z] = calculateCellEnergy( neighborhoodMap , space[x][y][z] );
                }
            }
        }
        return energySpace;
    }

    void GrainsGrowth::setNeighborhoodType( std::string inCA_NeighborhoodType )
    {
        checkInput_NeighborhoodType(inCA_NeighborhoodType);
        setNeighborhoodVector();
    }

    void GrainsGrowth::setBoundaryConditionsType( std::string inBoundaryConditionsType )
    {
        checkInput_BoundaryConditionsType(inBoundaryConditionsType);
    }

    void GrainsGrowth::setDimensionalityType( std::string inDimensionalityType )
    {
        checkInput_DimensionalityType(inDimensionalityType);
    }

    void GrainsGrowth::setSimulationType( std::string inSimulationType )
    {
        checkInput_SimulationType(inSimulationType);
    }

}