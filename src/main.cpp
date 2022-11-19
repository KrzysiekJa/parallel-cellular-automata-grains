#include <map>
#include <fstream>
#include <chrono>
#include "classes/grainsgrowth.cpp"
#include "utils/utils_functions.cpp"


int main(void)
{
    /************************************************
     *      Reading from xml config file
     ***********************************************/

    std::ifstream inputFile( "data_file/data.xml" );
    std::map< std::string, std::string > configDataMap = utils::readDataFromXML( inputFile );

    /************************************************
     *      Initializations of structure
     ***********************************************/

    auto t_start = std::chrono::high_resolution_clock::now();

    simula::GrainsGrowth grainsGrowth( configDataMap );

    auto t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double> init_time = t_finish - t_start;

    /************************************************
     *      Simulation
     ***********************************************/

    t_start = std::chrono::high_resolution_clock::now();

    grainsGrowth.makeSimulation(
        std::stoul( configDataMap["numberOfNucleus"] ), 
        std::stoul( configDataMap["iterationsMC"] )
        );
    
    t_finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration< double> simula_time = t_finish - t_start;


    std::string matrixStr = grainsGrowth.spaceToDisplay();
    //std::cout << matrixStr;

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
    
    std::string measurementsFile = "data_file/time.txt";
    std::map< std::string, double > timeMeasurementsMap = { 
        {"initTime", init_time.count()}, 
        {"simulTime", simula_time.count()}, 
        {"savingTime", saving_time.count()} 
        };
    utils::saveTimeMeasurements( measurementsFile, timeMeasurementsMap);
}
