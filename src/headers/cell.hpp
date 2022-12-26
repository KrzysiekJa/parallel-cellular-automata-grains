#ifndef CELL_H
#define CELL_H

#include <vector>
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
        CELL( INT_TYPE, unsigned );
        ~CELL();
        void populateNeighborhoodVec( std::vector< CELL * > );
    };
}
#endif /* CELL_H */