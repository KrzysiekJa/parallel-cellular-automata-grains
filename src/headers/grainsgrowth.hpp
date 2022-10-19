#include <string>

namespace simula{

    class GrainsGrowth
    {
    private:
        unsigned long dimSize;
        unsigned long ** grid;
    protected:
        void initGrid();
        void putSeedsInGrid(int numberOfSeeds);
    public:
        GrainsGrowth(unsigned long size);
        ~GrainsGrowth();
        unsigned long ** makeSimulation(int numberOfSeeds);
        std::string gridToString();
        unsigned long getSize();
        void setSize(unsigned long newDimSize);
        unsigned long ** getGrid();
    };

}