#include "tools.hh"
#include "ME2_Analytic.hh"
#include "var.hh"
#include "cuts.hh"
#include "integration.hh"
#include "dsigma.hh"

using namespace std;
using namespace Variables;



// Master function for all 2 to 2 neutrino scattering processes
double dsigma_nu( KinematicData &Kin, const int pdg_codes[] ){
        // nu + nu scattering

        // Immediately check to see if kinematic passed
        if( !Kin.get_cuts() ){
                return 0.0;
        }

        for( int i = 0; i < 2; i++){
                if( !is_neutrino(pdg_codes[i]) ){
                        cerr << "dsigma_nu: pdg particle  = " << pdg_codes[i] << endl;
                        cerr << "this is not a neutrino\n";
                        abort();
                }
        }

        // Initialise to zero
        double ME2(0.);

        // The special nu nubar > f fbar case
        if( !is_neutrino(pdg_codes[3]) ){

                // Check that same flavour fermions are requested
                if( abs(pdg_codes[2]) != abs(pdg_codes[3]) ){
                        cerr << "dsigma_nu: unequal outgoing fermion flavours  = " << pdg_codes[2] << "\t" << pdg_codes[3] << endl;
                        abort();
                }

                // Set up correct flavours
                int f1 = (pdg_codes[0] > 0)? 1: 2;
                int f2 = (pdg_codes[1] > 0)? 1: 2;
                // 
                int f3 = (pdg_codes[2] > 0)? 3: 4;
                int f4 = (pdg_codes[3] > 0)? 3: 4;                
                // nu nubar > f fbar
                ME2 = ME2_Analytic::nunux_ffx(f1,f2,f3,f4,Kin,abs(pdg_codes[0]),abs(pdg_codes[2]));
        }
        // Otherwise all neutrino scattering
        else{
                const int array_pdg[4] = {pdg_codes[0],pdg_codes[1],pdg_codes[2],pdg_codes[3]};
                // Determine whether we are in the same or different flavour scattering configuration
                ME2 = ME2_Analytic::nunu_nunu(Kin,array_pdg);
        }    

        // Some tests of the analytic continuation
        // cross the process and take ratio
        // cout << ME2_Analytic::nu1nu1_nu1nu1(1,2,3,4,Kin,12,12) / ME2_Analytic::nu1nu1x_nu1nu1x(1,3,4,2,Kin,12,12) << endl;


        // dphi_2 factor (and integration jacobian stored within Kin.weight)    
        return ME2 * Kin.weight() / Kin.flux();
}



// Master function for all 2 to 2 neutrino scattering processes
double dsigma_gamma( KinematicData &Kin, const int pdg_codes[] ){
        // (anti)nu + photon scattering

        // Immediately check to see if kinematic passed
        if( !Kin.get_cuts() ){
                return 0.0;
        }

        if( !is_neutrino(pdg_codes[0]) ){
                cerr << "dsigma_gamma: pdg particle 1 = " << pdg_codes[0] << endl;
                abort();
        }

        if( pdg_codes[1] != 22 ){
                cerr << "dsigma_gamma: pdg particle 2 = " << pdg_codes[1] << endl;
                abort();             
        }


        // Initialise to zero
        double ME2(0.);

        // For now only consider the simple case  nu + gamma > nu + f + fbar

        if( !is_neutrino(pdg_codes[2]) ){
                cerr << "dsigma_gamma: pdg particle 3 = " << pdg_codes[2] << endl;
                abort();                 
        }

        ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, abs(pdg_codes[2]), abs(pdg_codes[3]));

        // dphi_3 factor (and integration jacobian stored within Kin.weight)    
        return ME2 * Kin.weight() / Kin.flux();
}
