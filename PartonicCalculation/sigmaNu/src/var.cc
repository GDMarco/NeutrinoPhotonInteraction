#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <array>
#include "var.hh"

using namespace std;
using namespace Variables;

// Declaration of all the external variables from var.hh
std::string Variables::s_process, Variables::s_target, Variables::s_projectile, Variables::s_correction;

// Some integration options (cuba dimensions, treatment of resonances)
int Variables::cuba_dimensions;

// Global variables/couplings etc.
std::complex<double> Variables::ALPHA, Variables::SW2;
std::complex<double> Variables::gLl, Variables::gLd, Variables::gLu, Variables::gLnu;
std::complex<double> Variables::gRl, Variables::gRd, Variables::gRu;
double Variables::ALPHAS;
// Global scales
double Variables::muf_var, Variables::mur_var, Variables::mu0, Variables::mu_loop;
// Some flavour info for QCD
int Variables::nf_as, Variables::nf_pdf;
double Variables::nf_var;
// Global switch to include Ioperator IR-subtraction in virtuals
bool Variables::active_Iop;
// Global CoM and beam definition
double Variables::Ecms, Variables::Ecms2;

int Variables::scale_opt;



//// Neutrino cuts/selections

// pT
bool Variables::active_pTnu_min,	Variables::active_pTnu_max;
double Variables::pTnu_min, 		Variables::pTnu_max;




// Function to update the cuba integration dimenions for every implemented process
void Variables::update_process_dimensions(){

	if( s_target == "nu" or s_target == "nux" )
		cuba_dimensions = 1;
	else if( s_target == "gamma" ){
		cuba_dimensions = 4;
	}
	else{
		cerr << "process s_target = " << s_target << " which is unsupported\n";
		abort();
	}

}

void Variables::update_process( int pdg_codes[] ){

	// Check the projecile p1
	if( is_neutrino(pdg_codes[0]) )
		s_projectile = (pdg_codes[0]>0)? "nu": "nux";

	// Check the target p2
	if( is_neutrino(pdg_codes[1]) ){
		s_target = (pdg_codes[1]>0)? "nu": "nux";
	}
	else if( pdg_codes[1] == 22 ){
		s_target = "gamma";
	}
	else{
		cout << "update_process: unspported target pdg code " << pdg_codes[1] << endl;
		abort();
	}

	// Update the cuba dimensions accordingly
	cuba_dimensions = (s_target=="gamma")? 4: 1;
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

// Set-up of default settings and cuts
void Variables::init_default(){

	// proc = nunu, nugamma
	s_process = "";
	s_correction = "";

	// Default values
	scale_opt = 1;
	muf_var = 1.0;
	mur_var = 1.0;
	mu0 = 100.0; // Fixed-scale if not using dynamic scale variation
	mu_loop = 100.0; // Fixed-scale entering Loops, wave function renormalisation constants etc.

	// Automatically include IR subtraction in Virtuals (used in Catani Seymour)
	active_Iop = true;

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
	cout << "projectile = " << s_projectile << endl;
	cout << "target = " << s_target << endl;	
	// cout << "proc = " << s_process << endl;
	// cout << "correction = " << s_correction << endl;
	//
	cout << "scale options:\n";
	cout << "scale_opt = " << scale_opt << endl;
	cout << "mu0 = " << mu0 << endl;
	cout << "muf_var = " << muf_var << endl;
	cout << "mur_var = " << mur_var << endl;
	//
	cout << "cuba_dimensions = " << cuba_dimensions << endl;
}

// A function that writes the settings of the program to a text file
// i.e. the settings used for the computation

void Variables::write_settings(ofstream &infile, const string pdfset, string process ){

	// Ecms
	infile << "# Ecms = " << Ecms << endl;
	// Scales
	infile << "# Scale option = " << scale_opt << " (-1=Fixed scale, 1,2,3,...=Dynamic)" << endl;
	infile << "# mu0 = " << mu0 << " (used if Scale option = -1 above)" << endl;
	infile << "# mur_var = " << mur_var << " (muR multiplicative factor)" << endl;
	infile << "# muf_var = " << muf_var << " (muF multiplicative factor)" << endl;
	// Process options
	infile << "# proc = " << s_process << endl;
	infile << "# correction = " << s_correction << endl;
	// Flavour scheme options
	infile << "# nf_max pdf = " << nf_pdf << " (max number of flavours considered in pdfs)" << endl;
	infile << "# nf_max as = " << nf_as << " (max number of flavours considered in as)" << endl;	
	infile << "# nf_var = " << nf_var << " (nf variable used to dynamically compute b0[nf], b1[nf] terms)" << endl;	

	// If process not specified, print all program settings
	if( process == ""){
		// Write out all cuts
	}

	return;

}

