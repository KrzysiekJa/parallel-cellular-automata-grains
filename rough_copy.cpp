#include <map>
#include <tuple>
#include <vector>
#include <algorithm>
#include <iostream>

#define dimSize 3

int ** grid = new int * [dimSize];


std::map<int, unsigned int> counter(
        const std::vector< std::tuple< int, int >> & positions
    )
{
    std::map<int, unsigned int> rv;
    
    for (auto pos = positions.begin(); pos != positions.end(); ++pos)
    {
        int value = grid[ std::get<0>( * pos) ][ std::get<1>( * pos) ];
        if (value != 0)
        {
            rv[ value ]++;
        }
    }

    return rv;
}

void display( const std::map<int, unsigned int> & counts )
{
    for (auto count = counts.begin(); count != counts.end(); ++count)
    {
        std::cout << "Value " << count->first << " has count "
                  << count->second << std::endl;
    }

    auto max = std::max_element(
        std::begin(counts), std::end(counts),
        [] (const std::pair<int, unsigned int> & p1, const std::pair<int, unsigned int> & p2) 
        {
            return p1.second < p2.second;
        }
    );

    std::cout << "Max caount for:" << max->first << std::endl;
}



int main(void)
{        
    for (int i = 0; i < dimSize; ++i)
    {
        grid[i] = new int [dimSize];
        memset( grid[i], 0, dimSize * sizeof(int));
    }

    grid[1][0] = 1;
    grid[0][1] = 2;
    grid[2][1] = 1;
    grid[1][2] = 0;
    //{<0,2,0>, <1,0,0>, <0,1,0>};

    std::vector< std::tuple< int, int >> vector {{1, 0}, {0, 1}, {2, 1}, {1, 2}};

    display( counter( vector ) );

    for (int i = 0; i < dimSize; i++)
    {
        for (size_t j = 0; j < dimSize; j++)
        {
            std::cout << grid[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    std::string CA_NeighborhoodType = "von_neumann";
    std::cout << CA_NeighborhoodType.compare("vn_neumann") << "\n";

    return 0;
}