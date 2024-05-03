#include "tools.hh"
#include "ME2_Analytic.hh"
#include "cuts.hh"
#include "integration.hh"
#include "dsigma_nu.hh"

using namespace std;
using namespace Variables;

// Processes for neutrino scattering on a neutrino target


// Function for partonic cross-section for pp > Q Qb at LO
double dsigma_nunu_LO( KinematicData &Kin ){
        // nu nubar > f fbar process
        // Immediately check to see if kinematic passed
        if( !Kin.get_cuts() ){
                return 0.0;
        }
        // Partonic cross-section = |M|^2 / ( 2 shat ) * dphi_2
        double ME2 = ME2_Analytic::nunux_ffx(1,2,3,4,Kin,12,11);

        // dphi_2 factor (and integration jacobian stored within Kin.weight)    
        return ME2 * Kin.weight() / Kin.flux();
}

