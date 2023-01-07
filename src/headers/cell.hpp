#ifndef CELL_H
#define CELL_H

#include <vector>
#include <string>
#include "variables.hpp"


namespace simula{

    class CELL
    {
    // simple container
    public:
        INT_TYPE id;
        unsigned energy;
        std::vector< CELL * > neighborhood;

        CELL();
        CELL( std::string, std::string );
        ~CELL();
        void addNeighbor( CELL * );
    };
}
#endif /* CELL_H */