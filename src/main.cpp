#include <map>
#include <fstream>
#include <chrono>
#include "classes/grainsgrowth.cpp"
#include "utils/utils_functions.cpp"


int main( int argc, char * argv[] )
{
    /************************************************
     *      Reading from xml config file
     ***********************************************/
    
    std::string inFileName( argv[1] );
    std::ifstream inputFile( "data/data.xml" );
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
    // -- ************** -- //
    t_start = std::chrono::high_resolution_clock::now();

    grainsGrowth.makeMonteCarloSimulation( std::stoul( configDataMap["iterationsMC"] ) );
    
    t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double> simula_MC_time = t_finish - t_start;

    //std::string matrixStr = grainsGrowth.spaceToDisplay();

    /************************************************
     *      Saving sturcture to csv file
     ***********************************************/

    INT_TYPE *** matrix = grainsGrowth.getSpace();
    SIZE_TYPE *** energySpace = grainsGrowth.getEnergySpace();
    SIZE_TYPE size = std::stoul( configDataMap["dimSize"] );
    std::map< std::string, SIZE_TYPE > dimMap { {"first", size}, {"second", size} };
    configDataMap["dimensionalityType"] == "3D" ? 
                dimMap.insert( {"third", size} ) : dimMap.insert( {"third", 1} );

    t_start = std::chrono::high_resolution_clock::now();

    utils::saveSpaceToCSV( 
        configDataMap["outputFile"], matrix, energySpace, dimMap , 
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
}
