#include "../headers/cell.hpp"


namespace simula{
    
    // default
    CELL::CELL()
    {
        this->id = 0;
    }

    CELL::CELL( INT_TYPE id )
    {
        this->id = id;
    }
    
    CELL::~CELL(){}
}