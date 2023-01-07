#include "../headers/cell.hpp"


namespace simula{
    
    // default
    CELL::CELL()
    {
        this->id = 0;
        this->energy = 0;
    }

    CELL::CELL( std::string CA_NeighborhoodType, std::string dimensionalityType )
    {
        this->id = 0;
        this->energy = 0;
        
        if (CA_NeighborhoodType == "von_neumann" && dimensionalityType == "2D")
        {
            this->neighborhood.reserve(4);
        }
        if (CA_NeighborhoodType == "moore" && dimensionalityType == "2D")
        {
            this->neighborhood.reserve(8);
        }
        if (CA_NeighborhoodType == "von_neumann" && dimensionalityType == "3D")
        {
            this->neighborhood.reserve(6);
        }
        if (CA_NeighborhoodType == "moore" && dimensionalityType == "3D")
        {
            this->neighborhood.reserve(26);
        }
    }
    
    CELL::~CELL(){}

    void CELL::addNeighbor( CELL * neigh )
    {
        this->neighborhood.emplace_back( neigh );
    }
}