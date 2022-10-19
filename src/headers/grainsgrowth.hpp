#include <string>

namespace simula{

    class GrainsGrowth
    {
    private:
        unsigned long size;
        unsigned long ** grid;
    public:
        GrainsGrowth(unsigned long size);
        ~GrainsGrowth();
        void initGrid();
        int makeSimulation();
        std::string gridToString();
        unsigned long getSize();
        void setSize(unsigned long newSize);
        unsigned long ** getGrid();
    };

}