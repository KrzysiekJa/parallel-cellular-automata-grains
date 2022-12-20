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
        CELL( INT_TYPE id );
        ~CELL();
    };
}
#endif /* CELL_H */