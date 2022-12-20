#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include "utils_functions.hpp"


namespace utils{

    std::map< std::string, std::string > readDataFromXML( std::istream & inStream )
    {
        std::map< std::string, std::string > outputMap;

        using boost::property_tree::ptree;
        ptree pt;
        read_xml( inStream, pt );
        ptree subtree = pt.get_child("GG_config");

        BOOST_FOREACH( ptree::value_type const& node, subtree.get_child( "" ) ) 
        {
            std::string label = node.first;
            std::string value = subtree.get<std::string>( label );
            outputMap.insert( {label, value} );
            
        }
        return outputMap;
    }


    void toRGBColor( SIZE_TYPE p, SIZE_TYPE np, float & r, float & g, float & b )
    {
        // returns: r,g,b in range <0,1>
        float inc = 6.0 / np;
        float x = p * inc;
        r = 0.0f; g = 0.0f; b = 0.0f;

        if ((0 <= x && x <= 1) || (5 <= x && x <= 6)) r = 1.0f;
        else if (4 <= x && x <= 5) r = x - 4;
        else if (1 <= x && x <= 2) r = 1.0f - (x - 1);
        if (1 <= x && x <= 3) g = 1.0f;
        else if (0 <= x && x <= 1) g = x - 0;
        else if (3 <= x && x <= 4) g = 1.0f - (x - 3);
        if (3 <= x && x <= 5) b = 1.0f;
        else if (2 <= x && x <= 3) b = x - 2;
        else if (5 <= x && x <= 6) b = 1.0f - (x - 5);
    }


    // ugly function
    void saveSpaceToCSV( 
        std::string outputFileName, 
        simula::CELL *** spaceMatrix,
        std::map< std::string, SIZE_TYPE > dimMap,
        SIZE_TYPE numberOfNucleus
        )
    {
        long unsigned id_count = 0;
        float r, g, b;
        std::string outputString = "";

        outputString += "id,x,y,z,grainId,colorRGB,energy\n";


        for (SIZE_TYPE x = 0; x < dimMap["first"]; ++x)
        {
            for (SIZE_TYPE y = 0; y < dimMap["second"]; ++y)
            {
                for (SIZE_TYPE z = 0; z < dimMap["third"]; ++z)
                {
                    outputString += std::to_string(id_count) + ",";
                    outputString += std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) + ",";
                    outputString += std::to_string( spaceMatrix[x][y][z].id ) + ","; // grain id
                    toRGBColor( spaceMatrix[x][y][z].id, numberOfNucleus, r, g, b );
                    outputString += std::to_string( int(r*255) ) + " " + std::to_string( int(g*255) ) + " ";
                    outputString += std::to_string( int(b*255) ) + ",";
                    outputString += std::to_string( spaceMatrix[x][y][z].energy );
                    outputString += "\n";
                    ++id_count;
                }
            }
        }
        outputString.pop_back(); // to remove last line -> last "\n"

        std::ofstream outputFile( outputFileName, std::ios::trunc );
        outputFile << outputString;
        outputFile.close();
    }


    void saveTimeMeasurements(
        std::string outputFileName, std::map< std::string, double > timeMeasurementsMap
        )
    {
        std::ofstream outputFile( outputFileName, std::ios::trunc );

        outputFile << timeMeasurementsMap["initTime"] << "\n";
        outputFile << timeMeasurementsMap["simulCATime"] << "\n"; // Cellular automata
        outputFile << timeMeasurementsMap["simulMCTime"] << "\n"; // Monte Carlo
        outputFile << timeMeasurementsMap["savingTime"];
        outputFile.close();
    }
}
