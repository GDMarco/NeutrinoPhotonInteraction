#include "tools.hh"
#include "ME2_Analytic.hh"
#include "var.hh"
#include "cuts.hh"
#include "integration.hh"
#include "dsigma.hh"

using namespace std;
using namespace Variables;



// The previous routines were becoming cumbersome, so they are now updated with a channel dependent integer
double dsigma_channels( KinematicData &Kin, int channel_id ){

        // The summary of all channels is stored via the global map<int,str> process_map
        // The general form of the differential cross-section is always:
        // d sigma = \sum |M|^2 / flux * dphi_n

        // The phase space dphi_n is generated according to whether we are in 2to2 or 2to3, accessed via Kin.weight() which contains all jacobian factors
        // The flux is always 1 / (4p1.p2) for massless particles, accessed via Kin.flux()

        // Thus, here we evaluate |M|^2 given a channel id
        double ME2(0.);

        // If computing virtual corrections for channels 1-8, require use of recola
        if( active_virtual ){
                active_recola = true;
                if( channel_id > 8 ){
                        cerr << "dsigma_channels: requested virtual correction for channel " << channel_id << endl;
                        abort();
                }
        }

        // Option 1: use only |M|^2 results from Recola
        if( active_recola ){
                ME2 = ME2_Analytic::compute_process_recola(Kin,channel,active_virtual);
        }
        // Option 2: use analytic results from |M|^2 which I have computed
        // For 2to2 they are approximately 20 times faster
        else{
                // All equal flavour results obtained from |M|^2 for nu1 + nu1 > nu1 + nu1
                // 1) nu1 + nu1 > nu1 + nu1
                if( channel_id == 1 )       ME2 = ME2_Analytic::nu1nu1_nu1nu1(1,2,3,4, Kin, pdg_projectile, 12);
                // 2) nu1 + nu1bar > nu1 + nu1bar
                if( channel_id == 2  ) ME2 = ME2_Analytic::nu1nu1_nu1nu1(1,4,3,2, Kin, pdg_projectile, 12);
                // 3) nu1bar + nu1bar > nu1bar + nu1bar
                if( channel_id == 3  ) ME2 = ME2_Analytic::nu1nu1_nu1nu1(3,4,1,2, Kin, pdg_projectile, 12);
                // 4) nu1bar + nu1 > nu1bar + nu1
                if( channel_id == 4  ) ME2 = ME2_Analytic::nu1nu1_nu1nu1(3,2,1,4, Kin, pdg_projectile, 12);
                // Mixed flavour results obtained from |M|^2 for nu1 + nu2 > nu1 + nu2
                // 5) nu1 + nu2 > nu1 + nu2        
                if( channel_id == 5  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(1,2,3,4, Kin, pdg_projectile, 12);
                // 6) nu1 + nu2bar > nu1 + nu2bar        
                if( channel_id == 6  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(1,4,3,2, Kin, pdg_projectile, 12);
                // 7) nu1bar + nu2 > nu1bar + nu2        
                if( channel_id == 7  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(3,2,1,4, Kin, pdg_projectile, 12);
                // 8) nu1bar + nu2bar > nu1bar + nu2bar        
                if( channel_id == 8  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(3,4,1,2, Kin, pdg_projectile, 12);
                // neutrino annihilation channels (3 channels)        
                // 9) nu1 + nu1bar > nu2 + nu2bar
                if( channel_id == 9  ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 14);
                // 10-12) nu1 + nu1bar > l1 + l1bar (provide outgoing fermion PDG)
                if( channel_id == 10 ) ME2 = ME2_Analytic::nu1nu1bar_l1l1bar(1,2,3,4, Kin, 11);
                if( channel_id == 11 ) ME2 = ME2_Analytic::nu1nu1bar_l1l1bar(1,2,3,4, Kin, 13);
                if( channel_id == 12 ) ME2 = ME2_Analytic::nu1nu1bar_l1l1bar(1,2,3,4, Kin, 15);                
                // 13-18) nu1 + nu2bar > l1 + l2bar [CC process, no Z-coupling info needed]
                if( channel_id > 12 and channel_id < 19 )  ME2 = ME2_Analytic::nu1nu2bar_l1l2bar(1,2,3,4, Kin );
                // 19-27 nu nubar > f fbar
                if( channel_id == 19 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 1);
                if( channel_id == 20 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 2);
                if( channel_id == 21 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 3);
                if( channel_id == 22 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 4);
                if( channel_id == 23 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 5);
                if( channel_id == 24 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 6);
                if( channel_id == 25 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 11);
                if( channel_id == 26 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 13);
                if( channel_id == 27 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 15);
                // 28-33
                if( channel_id > 27 and channel_id < 34 ) ME2 = ME2_Analytic::nugumma_Wl(1,2,3,4, Kin);

                // testing channel
                if( channel_id == 99 ) ME2 = ME2_Analytic::eeB0g0NCM(1,2,3,4, Kin, pdg_projectile, pdg_fermion);


                // Channels 101+
                // 101) nu_1 gamma -> nu_1 l_2 l_2~
                if( channel_id == 101 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 102) nu_1 gamma -> nu_1 q q~
                if( channel_id == 102 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion); 
                // 103) nu_1bar gamma -> nu_1bar l_2~ l_2
                if( channel_id == 103 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(3, 2, 1, 5, 4, Kin, pdg_projectile, pdg_fermion);
                // 104) nu_1bar gamma -> nu_1bar q~ q
                if( channel_id == 104 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(3, 2, 1, 5, 4, Kin, pdg_projectile, pdg_fermion);
                // 105) nu_1 gamma -> l_1 nu_2 l_2~
                if( channel_id == 105 ) ME2 = ME2_Analytic::nu1gamma_l1nu2l2x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 106) nu_1 gamma -> l_1 q q~
                if( channel_id == 106 ) ME2 = ME2_Analytic::nu1gamma_l1qqbar(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 107) nu_1x gamma -> l_1x nu_2x l_2
                if( channel_id == 107 ) ME2 = ME2_Analytic::nu1gamma_l1nu2l2x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);                
                // 108) nu_1x gamma -> l_1x qx q [c.c. of the above squared amplitude so should be the same]
                if( channel_id == 108 ) ME2 = ME2_Analytic::nu1gamma_l1qqbar(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 109) nu_1 gamma -> l_1 nu_1 l_1~
                if( channel_id == 109 ) ME2 = ME2_Analytic::nu1gamma_l1nu1l1x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 110) nu_1~ gamma -> l_1~ nu_1~ l_1
                if( channel_id == 110 ) ME2 = ME2_Analytic::nu1gamma_l1nu1l1x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);

                if( channel_id > 108 ){
                        cout << "Computation still in progress, use recola for now\n";
                        abort();
                }

        }

        // For performing phase-space point tests of the squared amplitudes
        bool check_recola = false;
        // Phase-space point test against recola
        if( check_recola and !active_recola ){
                cout << "channel = " << channel << endl;
                cout << "process      = " << process_map.at(channel) << endl;
                cout << "process(rcl) = " << process_map_rcl.at(channel) << endl;                
                double ME2_rcl = ME2_Analytic::compute_process_recola(Kin,channel,active_virtual);
                // RECOLA averages over neutrino spins... so must correct this factor of 2 per incoming neutrino
                double nu_spin_correction = ( channel < 100 )? 4.: 2.;
                if( channel > 27 and channel < 34 ) nu_spin_correction = 2.;
                cout << "Apply correction to Recola |M|^2 for neutrino spin averaging of " << nu_spin_correction << endl;
                cout << setprecision(15) << "ME2 / ME2_recola = " << ME2 / ( ME2_rcl * nu_spin_correction ) << endl;
                // cout << "ME2 = " << ME2 << endl;
                // cout << "ME2_rcl = " << ME2_rcl << endl;                
                cout << endl;
        }

        return ME2 * Kin.weight() / Kin.flux();
}
