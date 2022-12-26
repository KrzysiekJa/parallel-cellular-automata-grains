#include "../headers/cell.hpp"


namespace simula{
    
    // default
    CELL::CELL()
    {
        this->id = 0;
        this->energy = 0;
    }

    CELL::CELL( INT_TYPE id, unsigned energy )
    {
        this->id = id;
        this->energy = energy;
    }
    
    CELL::~CELL(){}

    void CELL::populateNeighborhoodVec( std::vector< CELL * > inVec )
    {
        this->neighborhood = inVec;
    }
}