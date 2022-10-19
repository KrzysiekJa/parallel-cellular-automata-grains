namespace simula{

    class GrainsGrowth
    {
    private:
        unsigned long X;
        unsigned long Y;
        unsigned long ** grid;
    public:
        GrainsGrowth(unsigned long inputX, unsigned long inputY);
        ~GrainsGrowth();
        int makeSimulation();
    };

}