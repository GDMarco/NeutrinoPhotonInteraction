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

        // Return |M|^2 and the phase-space point
        // cout << "Test of recola\n";
        // Kin.p(1).print();
        // Kin.p(2).print();
        // Kin.p(3).print();
        // Kin.p(4).print();
        // cout << "ME2 = " << ME2 << endl;
        abort();



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


        cout << "Ratio is here\n";
        cout << ME2_Analytic::compute_process_recola(Kin,1,active_virtual) / ME2 << endl;      

        cout << "ME2 = " << ME2 << endl;
        abort();

        // dphi_3 factor (and integration jacobian stored within Kin.weight)    
        return ME2 * Kin.weight() / Kin.flux();
}






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
                if( channel_id == 1 )       ME2 = ME2_Analytic::nu1nu1_nu1nu1(1,2,3,4, Kin, pdg_projectile, pdg_fermion);
                // 2) nu1 + nu1bar > nu1 + nu1bar
                else if( channel_id == 2  ) ME2 = ME2_Analytic::nu1nu1_nu1nu1(1,4,3,2, Kin, pdg_projectile, pdg_fermion);
                // 3) nu1bar + nu1bar > nu1bar + nu1bar
                else if( channel_id == 3  ) ME2 = ME2_Analytic::nu1nu1_nu1nu1(3,4,1,2, Kin, pdg_projectile, pdg_fermion);
                // 4) nu1bar + nu1 > nu1bar + nu1
                else if( channel_id == 4  ) ME2 = ME2_Analytic::nu1nu1_nu1nu1(3,2,1,4, Kin, pdg_projectile, pdg_fermion);
                // Mixed flavour results obtained from |M|^2 for nu1 + nu2 > nu1 + nu2
                // 5) nu1 + nu2 > nu1 + nu2        
                else if( channel_id == 5  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(1,2,3,4, Kin, pdg_projectile, pdg_fermion);
                // 6) nu1 + nu2bar > nu1 + nu2bar        
                else if( channel_id == 6  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(1,4,3,2, Kin, pdg_projectile, pdg_fermion);
                // 7) nu1bar + nu2 > nu1bar + nu2        
                else if( channel_id == 7  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(3,2,1,4, Kin, pdg_projectile, pdg_fermion);
                // 8) nu1bar + nu2bar > nu1bar + nu2bar        
                else if( channel_id == 8  ) ME2 = ME2_Analytic::nu1nu2_nu1nu2(3,4,1,2, Kin, pdg_projectile, pdg_fermion);
                // neutrino annihilation channels (3 channels)        
                // 9) nu1 + nu1bar > nu2 + nu2bar
                else if( channel_id == 9  ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, 12, 14);
                // 10) nu1 + nu1bar > l1 + l1bar (provide outgoing fermion PDG)
                else if( channel_id == 10 ) ME2 = ME2_Analytic::nu1nu1bar_l1l1bar(1,2,3,4, Kin, pdg_projectile, pdg_fermion);
                // 11) nu1 + nu1bar > f + fbar (provide outoing fermion PDG)
                else if( channel_id == 11 ) ME2 = ME2_Analytic::nu1nu2bar_l1l2bar(1,2,3,4, Kin, pdg_projectile, pdg_fermion);
                //
                else if( channel_id == 12 ) ME2 = ME2_Analytic::nunux_ffx(1,2,3,4, Kin, pdg_projectile, pdg_fermion);
                //
                else if( channel_id == 13 ) ME2 = ME2_Analytic::nugumma_Wl(1,2,3,4, Kin);
                //
                else if( channel_id == 14 ) ME2 = ME2_Analytic::nugumma_Wl(1,2,3,4, Kin);                                                                


                //
                // Channels 101+
                // 101) nu_1 gamma -> nu_1 l_2 l_2~
                else if( channel_id == 101 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 102) nu_1 gamma -> nu_1 q q~
                else if( channel_id == 102 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion); 
                // 103) nu_1bar gamma -> nu_1bar l_2~ l_2
                else if( channel_id == 103 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(3, 2, 1, 5, 4, Kin, pdg_projectile, pdg_fermion);
                // 104) nu_1bar gamma -> nu_1bar q~ q
                else if( channel_id == 104 ) ME2 = ME2_Analytic::nu1gamma_nu1f2f2x(3, 2, 1, 5, 4, Kin, pdg_projectile, pdg_fermion);
                // 105) nu_1 gamma -> l_1 nu_2 l_2~
                else if( channel_id == 105 ) ME2 = ME2_Analytic::nu1gamma_l1nu2l2x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 106) nu_1 gamma -> l_1 q q~
                else if( channel_id == 106 ) ME2 = ME2_Analytic::nu1gamma_l1qqbar(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 107) nu_1x gamma -> l_1x nu_2x l_2
                else if( channel_id == 107 ) ME2 = ME2_Analytic::nu1gamma_l1nu2l2x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);                
                // 108) nu_1x gamma -> l_1x qx q [c.c. of the above squared amplitude so should be the same]
                else if( channel_id == 108 ) ME2 = ME2_Analytic::nu1gamma_l1qqbar(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 109) nu_1 gamma -> l_1 nu_1 l_1~
                else if( channel_id == 109 ) ME2 = ME2_Analytic::nu1gamma_l1nu1l1x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);
                // 110) nu_1~ gamma -> l_1~ nu_1~ l_1
                else if( channel_id == 110 ) ME2 = ME2_Analytic::nu1gamma_l1nu1l1x(1, 2, 3, 4, 5, Kin, pdg_projectile, pdg_fermion);

                if( channel_id > 108 ){
                        cout << "Computation still in progress, use recola for now\n";
                        abort();
                }

        }

        // For performing phase-space point tests of the squared amplitudes
        bool check_recola = true;
        // Phase-space point test against recola
        if( check_recola and !active_recola ){
                cout << "channel = " << channel << endl;
                cout << "process      = " << process_map.at(channel) << endl;
                cout << "process(rcl) = " << process_map_rcl.at(channel) << endl;                
                double ME2_rcl = ME2_Analytic::compute_process_recola(Kin,channel,active_virtual);
                cout << setprecision(15) << "ME2 / ME2_recola = " << ME2 / ME2_rcl << endl;
                // cout << "ME2 = " << ME2 << endl;
                // cout << "ME2_rcl = " << ME2_rcl << endl;                
                cout << endl;
        }

        return ME2 * Kin.weight() / Kin.flux();
}
