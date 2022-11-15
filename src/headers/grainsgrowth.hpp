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
        SIZE_TYPE thirdDim;
        INT_TYPE *** space;
        INT_TYPE *** nextSpace;
        std::string CA_NeighborhoodType = "von_neumann"; // alternative: "moore"
        std::string boundaryConditionsType = "absorbing"; // alternative: "periodic"
        std::string dimensionalityType = "2D"; // alternative: "3D"
        std::string simulationType = "none"; // alternative: "monte_carlo"
        std::vector< std::tuple< short, short, short >> neighborhoodVector;
        std::vector< std::tuple< SIZE_TYPE, SIZE_TYPE, SIZE_TYPE >> monteCarloVector;

        void checkInput_NeighborhoodType( std::string );
        void checkInput_BoundaryConditionsType( std::string );
        void checkInput_DimensionalityType( std::string );
        void checkInput_SimulationType( std::string );
        void setThirdDim();
        void setNeighborhoodVector();
        void setMonteCarloVector();
        INT_TYPE *** initSpace( INT_TYPE *** );

        INT_TYPE *** copySpace( INT_TYPE *** );
        void putNucleusInSpace( unsigned );
        void createBasicStructure( unsigned );
        bool iterateOverSpace();
        bool checkIfCell_IdEqual_0( INT_TYPE, INT_TYPE, INT_TYPE );
        void checkCell_NeighborhoodIds( INT_TYPE, INT_TYPE, INT_TYPE );
        bool inspectNeighborhoodCell( INT_TYPE, INT_TYPE, INT_TYPE, std::tuple< short, short, short > );
        void performCell_IdChange( INT_TYPE, INT_TYPE, INT_TYPE );
        std::map< INT_TYPE, unsigned > countNeighborhood( INT_TYPE, INT_TYPE, INT_TYPE );
        std::map< INT_TYPE, unsigned > checkCellDuringCounting( long, long, long, std::map< INT_TYPE, unsigned > );
        INT_TYPE mapIfPeriodic( long );
        INT_TYPE getNeighborhood_MaxFranction( std::map< INT_TYPE, unsigned > );

        void makeMonteCarloSimulation( unsigned );
        void shuffleMonteCarloVector();
        void iterateOverMCVector( std::vector< std::tuple< INT_TYPE, INT_TYPE, INT_TYPE > > );
        void inspectCellEnergy( INT_TYPE, INT_TYPE, INT_TYPE );
        void checkCellEnergeticStatus( std::map< INT_TYPE, unsigned >, unsigned , INT_TYPE, INT_TYPE, INT_TYPE );
        unsigned calculateCellEnergy( std::map< INT_TYPE, unsigned >, INT_TYPE );
    public:
        GrainsGrowth( 
            unsigned , std::string = "von_neumann", std::string = "absorbing", 
            std::string = "2D", std::string = "none" 
            );
        ~GrainsGrowth();
        INT_TYPE *** makeSimulation( unsigned, unsigned = 0 );
        std::string spaceToDisplay();
        SIZE_TYPE getSize();
        void setSize( SIZE_TYPE );
        INT_TYPE *** getSpace();
        void setNeighborhoodType( std::string );
        void setBoundaryConditionsType( std::string );
        void setDimensionalityType( std::string );
        void setSimulationType( std::string );
    };

}