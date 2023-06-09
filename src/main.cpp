#include <map>
#include <fstream>
#include <chrono>
#include "headers/cell.hpp"
#include "headers/grainsgrowth.hpp"
#include "utils/utils_functions.hpp"


int main( int argc, char * argv[] )
{
    /************************************************
     *      Reading from xml config file
     ***********************************************/
    
    // std::string inFileName( argv[1] );
    std::string inFileName = "data/2022-12-15-00:46:57:63_moore_periodic_monte_carlo_3D_10.xml";
    std::ifstream inputFile( inFileName );
    std::map< std::string, std::string > configDataMap = utils::readDataFromXML( inputFile );
    
    /************************************************
     *      Initializations of structure
     ***********************************************/

    auto t_start = std::chrono::high_resolution_clock::now();

    simula::GrainsGrowth grainsGrowth( configDataMap );

    auto t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double> init_time = t_finish - t_start;

    /************************************************
     *      Simulation CA & MC
     ***********************************************/

    t_start = std::chrono::high_resolution_clock::now();

    grainsGrowth.createBasicStructure( std::stoul( configDataMap["numberOfNucleus"] ) );

    t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double> simula_CA_time = t_finish - t_start;
    //    ---   **************   ---    //
    t_start = std::chrono::high_resolution_clock::now();

    grainsGrowth.makeMonteCarloSimulation( std::stoul( configDataMap["iterationsMC"] ) );
    
    t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double> simula_MC_time = t_finish - t_start;

    //std::string matrixStr = grainsGrowth.spaceToDisplay();

    /************************************************
     *      Saving sturcture to csv file
     ***********************************************/

    simula::CELL *** space_matrix = grainsGrowth.getSpace();
    SIZE_TYPE size = std::stoul( configDataMap["dimSize"] );
    std::map< std::string, SIZE_TYPE > dimMap { {"first", size}, {"second", size} };
    configDataMap["dimensionalityType"] == "3D" ? 
                dimMap.insert( {"third", size} ) : dimMap.insert( {"third", 1} );

    t_start = std::chrono::high_resolution_clock::now();

    utils::saveSpaceToCSV( 
        configDataMap["outputFile"], space_matrix, dimMap , 
        std::stoi( configDataMap["numberOfNucleus"] ) 
    );
    
    t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double> saving_time = t_finish - t_start;

    /************************************************
     *      Saving measured time
     ***********************************************/
    
    std::map< std::string, double > timeMeasurementsMap = { 
        {"initTime", init_time.count()}, 
        {"simulCATime", simula_CA_time.count()}, 
        {"simulMCTime", simula_MC_time.count()}, 
        {"savingTime", saving_time.count()} 
    };
    utils::saveTimeMeasurements( configDataMap["measurementsFile"], timeMeasurementsMap );

    return 0;
}
