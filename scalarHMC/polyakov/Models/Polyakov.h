#include "Model.h"
#include "../Other/other.h"

using namespace std;

class Polyakov : public Model{
    public:
        Polyakov(class Lattice* latt, int Nlinks, double beta);
        ~Polyakov(){;}
        void modelInfo();
        double evaluateAction();
        void saveConfiguration();
        void writeSavedConfiguration();
        double evaluateMDdrift(int nt, int nx, int ny, int nz); // drift term for Molecular Dynamics evolution
        vector<vector<vector<complex<double>>>> U; // Polyakov chain
        void printU(int idx);
        void printUcopy(int idx);
        void newConf(int idx, int copyConf=1);
    protected:
        vector<vector<vector<complex<double>>>> Ucopy; // Polyakov chain
        int Nlinks;
        double beta;
};