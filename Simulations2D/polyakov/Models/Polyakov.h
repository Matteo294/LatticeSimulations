#include "Model.h"
#include "../Other/other.h"

using namespace std;

class Polyakov : public Model{
    public:
        Polyakov(class Lattice* latt, int Nlinks, double beta);
        ~Polyakov(){;}
        void modelInfo();
        vector<vector<double>> copyConfiguration(){return Model::copyConfiguration();}
        void writeConfiguration(vector<vector<vector<double>>> U);
        double evaluateAction();
        double evaluateMDdrift(int nt, int nx, int ny, int nz); // drift term for Molecular Dynamics evolution
        vector<vector<vector<complex<double>>>> U; // Polyakov chain
    protected:
        int Nlinks;
        double beta;
};