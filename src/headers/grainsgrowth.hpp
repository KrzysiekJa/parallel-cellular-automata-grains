#include <map>
#include <tuple>
#include <vector>
#include <string>
#include "variables.hpp"


namespace simula{

    class GrainsGrowth
    {
    private:
        SIZE_TYPE dimSize;
        INT_TYPE ** grid;
        INT_TYPE ** nextGrid;
        std::string CA_NeighborhoodType = "von_neumann"; // alternative: "moore"
        std::string boundaryConditionsType = "absorbing"; // alternative: "periodic"
        std::vector< std::tuple< short, short >> neighborhoodVector;
        void checkInput_NeighborhoodType( std::string );
        void checkInput_BoundaryConditionsType( std::string );
        void setNeighborhoodVector();
    protected:
        INT_TYPE ** initGrid( INT_TYPE ** );
        INT_TYPE ** copyGrid( INT_TYPE ** );
        void putSeedsInGrid( unsigned );
        bool makeIterationStep();
        bool iterateOverGrid();
        bool checkIfCell_IdEqual_0( INT_TYPE, INT_TYPE );
        void checkCell_NeighborhoodIds( INT_TYPE, INT_TYPE );
        bool inspectNeighborhoodCell( INT_TYPE, INT_TYPE, std::tuple< short, short > );
        void performCell_IdChange( INT_TYPE, INT_TYPE );
        std::map< INT_TYPE, unsigned > countNeighborhood( INT_TYPE, INT_TYPE );
        std::map< INT_TYPE, unsigned > checkCellDuringCounting( long, long, std::map< INT_TYPE, unsigned > );
        INT_TYPE mapIfPeriodic( long );
        INT_TYPE getNeighborhood_MaxFranction( std::map< INT_TYPE, unsigned > );
    public:
        GrainsGrowth( unsigned , std::string = "von_neumann", std::string = "absorbing");
        ~GrainsGrowth();
        INT_TYPE ** makeSimulation( unsigned );
        std::string gridToDisplay();
        SIZE_TYPE getSize();
        void setSize( SIZE_TYPE );
        INT_TYPE ** getGrid();
    };

}