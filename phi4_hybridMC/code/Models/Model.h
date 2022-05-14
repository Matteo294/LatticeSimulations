#pragma once
#include "../Lattices/Lattice.h"
#include "../Other/other.h"

class Model{
    public:
        Model(class Lattice* l);
        ~Model();
        virtual double evaluateAction(){return 0.;}
        virtual void ModelInfo(){};
        virtual double evaluateMDdrift(int nt, int nx){return 0.;}
    protected:
        class Lattice* lattice;
        
};