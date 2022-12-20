#ifndef UTILS_FUNCTIONS_H
#define UTILS_FUNCTIONS_H

#include <map>
#include <fstream>
#include <string>
#include "../headers/variables.hpp"
#include "../headers/cell.hpp"


namespace utils{

    std::map< std::string, std::string > readDataFromXML( std::istream & );

    void toRGBColor( SIZE_TYPE, SIZE_TYPE, float &, float &, float & );

    void saveSpaceToCSV( 
        std::string, simula::CELL ***, std::map< std::string, SIZE_TYPE >, SIZE_TYPE 
    );

    void saveTimeMeasurements( std::string, std::map< std::string, double > );

}
#endif /* UTILS_FUNCTIONS_H */