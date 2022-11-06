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
        std::string CA_NeighborhoodType = "von_neumann"; // alternative: "moore", "monte_carlo"
        std::string boundaryConditionsType = "absorbing"; // alternative: "periodic"
        std::string simulationType = "standard"; // alternative: "monte_carlo"
        std::vector< std::tuple< short, short >> neighborhoodVector;
        std::vector< std::tuple< SIZE_TYPE, SIZE_TYPE >> monteCarloVector;
        void checkInput_NeighborhoodType( std::string );
        void checkInput_BoundaryConditionsType( std::string );
        void checkInput_SimulationType( std::string );
        void setNeighborhoodVector();
        void setMonteCarloVector();
    protected:
        INT_TYPE ** initGrid( INT_TYPE ** );
        INT_TYPE ** copyGrid( INT_TYPE ** );
        void putNucleusInGrid( unsigned );
        void createBasicStructure( unsigned );
        bool iterateOverGrid();
        bool checkIfCell_IdEqual_0( INT_TYPE, INT_TYPE );
        void checkCell_NeighborhoodIds( INT_TYPE, INT_TYPE );
        bool inspectNeighborhoodCell( INT_TYPE, INT_TYPE, std::tuple< short, short > );
        void performCell_IdChange( INT_TYPE, INT_TYPE );
        std::map< INT_TYPE, unsigned > countNeighborhood( INT_TYPE, INT_TYPE );
        std::map< INT_TYPE, unsigned > checkCellDuringCounting( long, long, std::map< INT_TYPE, unsigned > );
        INT_TYPE mapIfPeriodic( long );
        INT_TYPE getNeighborhood_MaxFranction( std::map< INT_TYPE, unsigned > );
        void makeMonteCarloSimulation( unsigned );
        void shuffleMonteCarloVector();
        void iterateOverMCVector( std::vector< std::tuple< INT_TYPE, INT_TYPE > > );
        void inspectCellEnergy( INT_TYPE, INT_TYPE );
        unsigned calculateCellEnergy( std::map< INT_TYPE, unsigned >, INT_TYPE );
        void checkCellEnergeticStatus( std::map< INT_TYPE, unsigned >, unsigned , INT_TYPE, INT_TYPE );
    public:
        GrainsGrowth( unsigned , std::string = "von_neumann", std::string = "absorbing", std::string = "standard");
        ~GrainsGrowth();
        INT_TYPE ** makeSimulation( unsigned, unsigned = 0 );
        std::string gridToDisplay();
        SIZE_TYPE getSize();
        void setSize( SIZE_TYPE );
        INT_TYPE ** getGrid();
        void setNeighborhoodType( std::string );
        void setBoundaryConditionsType( std::string );
        void setSimulationType( std::string );
    };

}