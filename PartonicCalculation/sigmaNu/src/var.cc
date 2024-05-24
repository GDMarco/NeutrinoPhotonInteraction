#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <array>
#include "var.hh"
#include "recola.hpp"

using namespace std;
using namespace Variables;

// Declaration of all the external variables from var.hh
// std::string Variables::s_process, Variables::s_target, Variables::s_projectile, Variables::s_correction;

// Process map for internal processes and recola interface
std::map<int,std::string> Variables::process_map, Variables::process_map_rcl;
// Channel selection, then PDG of the neutrino projectile and outgoing fermion line (via particle 4)
int Variables::channel, Variables::pdg_projectile, Variables::pdg_fermion;
// Recola interface
bool Variables::active_recola;
// Some integration options (cuba dimensions, treatment of resonances)
int Variables::cuba_dimensions;

// Global variables/couplings etc.
std::complex<double> Variables::ALPHA, Variables::SW2;
std::complex<double> Variables::gLl, Variables::gLd, Variables::gLu, Variables::gLnu;
std::complex<double> Variables::gRl, Variables::gRd, Variables::gRu;
double Variables::ALPHAS;
// Global scales
double Variables::muf_var, Variables::mur_var, Variables::mu0, Variables::mu_loop;
// Virtual
bool Variables::active_virtual;
// Recola interface
// Some flavour info for QCD
int Variables::nf_as, Variables::nf_pdf;
double Variables::nf_var;
// Global CoM and beam definition
double Variables::Ecms, Variables::Ecms2;

int Variables::scale_opt;

//// Neutrino cuts/selections

// pT
bool Variables::active_pTnu_min,	Variables::active_pTnu_max;
double Variables::pTnu_min, 		Variables::pTnu_max;



// Function to update the cuba integration dimenions for every implemented process
void Variables::update_process_dimensions(){
	// 2to2 processes registered with < 100
	// 2to3 processes registered with > 100
	// Fixes the phase-space integration dimensions
	cuba_dimensions = (channel < 100)? 1: 4;
}

// Initialise the Electroweak scheme
void Variables::init_scheme(int ischeme){
	// Initiate various EW schemes

	// alpha_0 CMS scheme
	if( ischeme == 0 ){
		ALPHA = ALPHA_ZERO;
		SW2 = SW2_OS;
		cout << "init_scheme: Initiated alpha_0 scheme" << endl;		
	}
    // alpha_GF CMS
	else if( ischeme == 1 ){
		// complex<double> agf = sqrt(2.0) * gf * MW2C * SW2_OS / M_PI;
		complex<double> agf = sqrt(2) / M_PI * gf * fabs( MW2C * SW2_OS );
		ALPHA = agf;
		SW2 = SW2_OS;
		cout << "init_scheme: Initiated GF scheme" << endl;
	}	
	else{
		cout << "init_scheme: Unsupported EW input scheme" << endl;
		abort();
	}
	cout << setprecision(15);
	cout << "ALPHA: " << ALPHA << endl;
	cout << "SW2: " << SW2 << endl;
    // // left
    gLl  = (-1.0/2.0 - Ql * SW2)*sqrt(1.0/(SW2*(1.0-SW2)) );
    gLnu = ( 1.0/2.0           )*sqrt(1.0/(SW2*(1.0-SW2)) );
    gLd  = (-1.0/2.0 - Qd * SW2)*sqrt(1.0/(SW2*(1.0-SW2)) );
    gLu  = ( 1.0/2.0 - Qu * SW2)*sqrt(1.0/(SW2*(1.0-SW2)) );    
    // // right
    gRl = -sqrt( SW2/(1.0-SW2) ) * Ql;
    gRd = -sqrt( SW2/(1.0-SW2) ) * Qd;
    gRu = -sqrt( SW2/(1.0-SW2) ) * Qu;

}

// Processes for neutrino scattering on a neutrino target
bool Variables::is_neutrino( int pdg ){
        if( abs(pdg) == 12 or abs(pdg) == 14 or abs(pdg) == 16 )
                return true;
        else{
                return false;
        }
}

double Variables::get_mass_fermion_pdg( int pdg ){
	// Check if neutrino
	if( is_neutrino(pdg) ) return 0.0;
	// quarks
	if( abs(pdg) == 1 ) return md;
	if( abs(pdg) == 2 ) return mu;
	if( abs(pdg) == 3 ) return ms;
	if( abs(pdg) == 4 ) return mc;
	if( abs(pdg) == 5 ) return mb;
	if( abs(pdg) == 6 ) return mt;
	// leptons
	if( abs(pdg) == 11 ) return me;
	if( abs(pdg) == 13 ) return mm;
	if( abs(pdg) == 15 ) return mtau;

	cerr << "get_fermion_mass_pdg: no value for pdg = " << pdg << endl;
	abort();
}

// Set-up of default settings and cuts
void Variables::init_default(){

	// Default values
	scale_opt = 1;
	muf_var = 1.0;
	mur_var = 1.0;
	mu0 = 100.0; // Fixed-scale if not using dynamic scale variation
	mu_loop = 100.0; // Fixed-scale entering Loops, wave function renormalisation constants etc.

	// Whether or not recola and virtual corrections are active be default
	active_virtual = false;
	active_recola = false;

	// Initialise
	nf_var = 5;
	nf_pdf = 5;
	nf_as = 5;

	// Cuba dimensions
	cuba_dimensions = 1; // number of dimensions of the integrand

	// Cuts on outgoing particles: if we wish to be differential
	active_pTnu_min = false; pTnu_min = 0.;
	active_pTnu_max = false; pTnu_max = 0.;

	// Effective EW scheme
	int ischeme_default = 1;
	init_scheme(ischeme_default);
}


void Variables::print_settings(){
	cout << "\033[0;31m Summary of global settings\033[0m" << endl;
	cout << "channel = " << channel << endl;
	cout << "cuba_dimensions = " << cuba_dimensions << endl;
}

// A function to initialise Recola and set up parameters/scheme
void Variables::init_recola(){
	// Particle masses
	Recola::set_pole_mass_z_rcl (mz,gz);
	Recola::set_pole_mass_w_rcl (mw,gw);
	// charged leptons
	Recola::set_pole_mass_electron_rcl(me);
	Recola::set_pole_mass_muon_rcl(mm,0.);
	Recola::set_pole_mass_tau_rcl(mtau,0.);
	// quarks
	Recola::set_pole_mass_top_rcl(mt,0.0);
	Recola::set_pole_mass_bottom_rcl(mb,0.0);
	Recola::set_pole_mass_charm_rcl(mc,0.0);
	Recola::set_pole_mass_strange_rcl(ms);
	Recola::set_pole_mass_strange_rcl(0.0);
	Recola::set_pole_mass_strange_rcl(0.0);
	// When to ignore fermion masses
	Recola::set_light_fermions_rcl(0.5*me);

	Recola::set_complex_mass_scheme_rcl();
	Recola::use_gfermi_scheme_and_set_alpha_rcl( ALPHA.real() );

	// Recola::set_resonant_particle_rcl ("Z");
	// Recola::set_resonant_particle_rcl ("W+");
	// Recola::set_resonant_particle_rcl ("W-");
}

// Function to initialise all scattering channels into a map
void Variables::init_channels(){

	////////////////////
	// 2to2 processes //
	////////////////////

	// All neutrino scattering (8 channels)
	// same flavour channels
	// 1) nu1 + nu1 > nu1 + nu1
	process_map.emplace( 1, "nu_1 nu_1 -> nu_1 nu_1");
	// 2) nu1 + nu1bar > nu1 + nu1bar
	process_map.emplace( 2, "nu_1 nu_1~ -> nu_1 nu_1~");
	// 3) nu1bar + nu1bar > nu1bar + nu1bar
	process_map.emplace( 3, "nu_1~ nu_1~ -> nu_1~ nu_1~");
	// 4) nu1bar + nu1 > nu1bar + nu1
	process_map.emplace( 4, "nu_1~ nu_1 -> nu_1~ nu_1");
	// mixed-flavour channels
	// 5) nu1 + nu2 > nu1 + nu2
	process_map.emplace( 5, "nu_1 nu_2 -> nu_1 nu_2");
	// 6) nu1 + nu2bar > nu1 + nu2bar
	process_map.emplace( 6, "nu_1 nu_2~ -> nu_1 nu_2~");
	// 7) nu1bar + nu2 > nu1bar + nu2
	process_map.emplace( 7, "nu_1~ nu_2 -> nu_1~ nu_2");
	// 8) nu1bar + nu2bar > nu1bar + nu2bar
	process_map.emplace( 8, "nu_1~ nu_2 -> nu_1~ nu_2");	
	// neutrino annihilation channels (3 channels)
	// 9) nu1 + nu1bar > nu2 + nu2bar
	process_map.emplace( 9, "nu_1 nu_1~ -> nu_2 nu_2~");
	// 10) nu1 + nu1bar > l1 + l1bar (provide outgoing fermion PDG)
	process_map.emplace( 10, "nu_1 nu_1~ -> l_1 l_1~");	
	// 11) nu1 + nu2bar > l1 + l2bar (provide PDG1 and PDG4 [l2bar flavour])
	process_map.emplace( 11, "nu_1 nu_2~ -> l1 l2bar");
	// 11) nu1 + nu1bar > f + fbar (provide outoing fermion PDG)
	process_map.emplace( 12, "nu_1 nu_1~ -> f f~");		
	// Onshell W processes
	// 13) nu1 + gamma > W+ l-
	process_map.emplace( 13, "nu_1 gamma -> W+ l_1");
	// 14) nu1bar + gamma > W- l+
	process_map.emplace( 14, "nu_1~ gamma -> W- l_1~");	

	////////////////////
	// 2to3 processes //
	////////////////////

	// neutrino  + photon scattering NC exchanges (2 channels + c.c.)
	// 101) nu1  + gamma > nu1 + l2 + l2bar
	process_map.emplace( 101, "nu_1 gamma -> nu_1 l_2 l_2~");
	// 102) nu1  + gamma > nu1 + q + qbar
	process_map.emplace( 102, "nu_1 gamma -> nu_1 q q~");
	// 103) nu1~ + gamma > nu1 + l2 + l2bar
	process_map.emplace( 103, "nu_1~ gamma -> nu_1~ l~ l_2");
	// 104) nu1~ + gamma > nu1 + q + qbar
	process_map.emplace( 104, "nu_1~ gamma -> nu_1~ q~ q");

	// neutrino  + photon scattering CC exchanges (2 channels + c.c.)
	// 105) nu1  + gamma > l1 + nu2 + l2~
	process_map.emplace( 105, "nu_1 gamma -> l_1 nu_2 l_2~");
	// 106) nu1  + gamma > l1 + u + d~ (massless quarks)
	process_map.emplace( 106, "nu_1 gamma -> l_1 u d~ [massless quarks]");
	// 107) nu1~ + gamma > l1~ + nu2~ + l2
	process_map.emplace( 107, "nu_1~ gamma -> l_1~ nu_2~ l_2");
	// 108) nu1~ + gamma > l1~ + u~ + d (massless quarks)
	process_map.emplace( 108, "nu_1~ gamma -> l_1~ u~ d [massless quarks]");

	// neutrino  + photon same lepton flavour cases (1 channel + c.c.)
	// 109) nu1  + gamma > l1 + nu1 + l1~
	process_map.emplace( 109, "nu_1 gamma -> l_1 nu_1 l_1~");
	// 110) nu1~ + gamma > l1~ + nu1~ + l1
	process_map.emplace( 110, "nu_1~ gamma -> l_1~ nu_1~ l_1");

	// Implementation of the Recola channels as above

	// The recola map is the same as above, but with some explicit fermion names
	// 1) nu1 + nu1 > nu1 + nu1
	process_map_rcl.emplace( 1, "nu_e nu_e -> nu_e nu_e");
	// 2) nu1 + nu1bar > nu1 + nu1bar
	process_map_rcl.emplace( 2, "nu_e nu_e~ -> nu_e nu_e~");
	// 3) nu1bar + nu1bar > nu1bar + nu1bar
	process_map_rcl.emplace( 3, "nu_e~ nu_e~ -> nu_e~ nu_e~");
	// 4) nu1bar + nu1 > nu1bar + nu1
	process_map_rcl.emplace( 4, "nu_e~ nu_e -> nu_e~ nu_e");
	// mixed-flavour channels
	// 5) nu1 + nu2 > nu1 + nu2
	process_map_rcl.emplace( 5, "nu_e nu_mu -> nu_e nu_mu");
	// 6) nu1 + nu2bar > nu1 + nu2bar
	process_map_rcl.emplace( 6, "nu_e nu_mu~ -> nu_e nu_mu~");
	// 7) nu1bar + nu2 > nu1bar + nu2
	process_map_rcl.emplace( 7, "nu_e~ nu_mu -> nu_e~ nu_mu");
	// 8) nu1bar + nu2bar > nu1bar + nu2bar
	process_map_rcl.emplace( 8, "nu_e~ nu_mu~ -> nu_e~ nu_mu~");
	// neutrino annihilation channels (3 channels)
	// 9) nu1 + nu1bar > nu2 + nu2bar
	process_map_rcl.emplace( 9, "nu_e nu_e~ -> nu_mu nu_mu~");
	// 10) nu1 + nu1bar > l1 + l1bar (provide outgoing fermion PDG)
	process_map_rcl.emplace( 10, "nu_tau nu_tau~ -> tau- tau+");	
	// 11) nu1 + nu1bar > f + fbar (provide outoing fermion PDG)
	process_map_rcl.emplace( 11, "nu_e nu_mu~ -> e- mu+");
	// 12) nu1 + nu2bar > l1 + l2bar (provide outoing fermion PDG)
	process_map_rcl.emplace( 12, "nu_e nu_e~ -> c c~");
	// 13) nu1 + gamma > W+ l-
	process_map_rcl.emplace( 13, "nu_mu gamma -> W+ mu-");
	// 14) nu1bar + gamma > W- l+
	process_map_rcl.emplace( 14, "nu_mu~ gamma -> W- mu+");	

	// neutrino  + photon scattering NC exchanges (2 channels + c.c.)
	// 101) nu1  + gamma > nu1 + l2 + l2bar
	process_map_rcl.emplace( 101, "nu_e gamma -> nu_e mu- mu+");
	// 102) nu1  + gamma > nu1 + q + qbar
	process_map_rcl.emplace( 102, "nu_e gamma -> nu_e c c~");
	// 103) nu1~ + gamma > nu1 + l2 + l2bar
	process_map_rcl.emplace( 103, "nu_e~ gamma -> nu_e~ mu+ mu-");
	// 104) nu1~ + gamma > nu1 + q + qbar
	process_map_rcl.emplace( 104, "nu_e~ gamma -> nu_e~ c~ c");

	// The following channels were computed only in a resonance approximation
	// So explicitly check them vs Recola for specific fermion choices
	// neutrino  + photon scattering CC exchanges (2 channels + c.c.)
	// 105) nu1  + gamma > l1 + nu2 + l2~
	process_map_rcl.emplace( 105, "nu_e gamma -> e- nu_mu mu+");
	// 106) nu1  + gamma > l1 + u + d~ (massless quarks)
	process_map_rcl.emplace( 106, "nu_e gamma -> e- u d~");
	// 107) nu1~ + gamma > l1~ + nu2~ + l2
	process_map_rcl.emplace( 107, "nu_e~ gamma -> e+ nu_mu~ mu-");
	// 108) nu1~ + gamma > l1~ + u~ + d (massless quarks)
	process_map_rcl.emplace( 108, "nu_e~ gamma -> e+ u~ d");

	// neutrino  + photon same lepton flavour cases (1 channel + c.c.)
	// 109) nu1  + gamma > l1 + nu1 + l1~
	process_map_rcl.emplace( 109, "nu_mu gamma -> mu- nu_mu mu+");
	// 110) nu1~ + gamma > l1~ + nu1~ + l1
	process_map_rcl.emplace( 110, "nu_mu~ gamma -> mu+ nu_mu~ mu-");

}

// Function for printing available channels
void Variables::print_channels(){
	// Process registration instructions
	cout << "\033[0;31m Channel selection\033[0m" << endl;
	cout << "provide the program with an integer to select the required channel\n";
	cout << "-c i\n\n";
	cout << "additionally provide the PDG code of the projectile (which identifies neutrino flavour)\n";
	cout << "-p PDG (where PDG = +/- 12, 14, 16)\n\n";
	cout << "additionally provide the PDG code of outgoing fermion 4 (that will fix the masses etc.)\n";
	cout << "-f PDG\n\n";
	cout << "This information gives the program all infomration it needs for channel selection from the following list:\n";
	cout << "All fermion masses included unless stated explicitly\n";
	cout << "\nFor channels 1-8, NLO corrections are available via Recola interface\n";
	cout << "This is controlled via global 'bool active_virtual = true/false'\n";

	// Print available channels
	for( auto element: process_map){
		cout << element.first << " " << element.second << endl;
	}
}

// Register and generate processes with recola
void Variables::init_recola_processes(){
	// Define/register
	for( auto &element: process_map_rcl ){

		// For all 2to2 processes, register as NLO (to provide access to the virtual corr.)		
		string order = "LO";
		if( element.first < 100 ) order = "NLO";
		else{
			order = "LO";
		}

		// Speed up by only defining requested channel?
		// if( element.first != channel ) continue;

		Recola::define_process_rcl(element.first,element.second,order);
		// First element is the integer, second is the string process declaration
	}
	// Generate
	Recola::generate_processes_rcl();
}


// A function that writes the settings of the program to a text file
// i.e. the settings used for the computation

void Variables::write_settings(ofstream &infile, string process ){

	// Save information on the electroweak scheme
	infile << "# Program settings for SigmaNu\n";

	// Electroweak inputs
	infile << "# Electroweak scheme and input parameters\n";
	infile << "# Scheme = Complex mass scheme\n";
	infile << "# alpha_GF = " << setprecision(15) << ALPHA.real() << endl;
	infile << "# Pole masses/widths for particles [GeV]\n";
	infile << "# mw = " << mw << endl;
	infile << "# gw = " << gw << endl;	
	infile << "# mz = " << mz << endl;
	infile << "# gz = " << gz << endl;
	infile << "# s^2_thetaw (derived) = " << SW2 << endl;
	// Fermion masses
	infile << "# charged leptons\n";
	infile << "# m_e = " << me << endl;
	infile << "# m_m = " << mm << endl;	
	infile << "# m_t = " << mtau << endl;
	infile << "# quarks\n";
	infile << "# m_d = " << md << endl;
	infile << "# m_u = " << mu << endl;	
	infile << "# m_s = " << ms << endl;
	infile << "# m_c = " << mc << endl;
	infile << "# m_b = " << mb << endl;	
	infile << "# m_t = " << mt << endl;	

	// Ecms
	// infile << "# Ecms = " << Ecms << endl;


	return;
}

