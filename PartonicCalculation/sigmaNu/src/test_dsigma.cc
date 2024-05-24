#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <array>
#include <cstring>
#include <getopt.h>
// General files
#include <cuba.h>
#include "integration.hh"
#include "var.hh"
#include "tools.hh"
// Process files
#include "dsigma.hh"
//Timing
#include <sys/time.h>

using namespace std;
using namespace Variables;

// Define some globals, to ease the set-up
// pdf_info pdf_cache;
int seed_cache;
int grid_cache = 2;

// Write more user-friendly command line reader
void print_usage() {
	cout << endl
		<< "Usage, dSigma:" << endl
		<< " -s --seed			<s>		seed number for integration\n"
		<< " -o --iobs			<o>		potentially if we consider some specific observeables or channel?\n"
		// Then the PDG of particles 1,2,3,4
		<< " -p --pdg_1			<p>		pdg of (anti)neutrino projecile\n" 
		<< " -f --pdg_4			<f>		pdg of fermion 4, fixes masses, couplings etc.\n"
		<< " -c --chan			<c>		channel selection\n"
		<< endl;
	// Supported channels
	print_channels();
	exit(0);
	return;
}

// Introduce a more user friendly interface
void read_arguments(int argc, char* argv[], int &seed, int &iobs, int &ichan, int &flav_1, int &flav_4) {
	// provide it length 4 integer array for PDG codes (these must all be entered)
	const char* const short_options = "s:o:p:f:c:";
	const struct option long_options[] = { { "help", 0, NULL, 'h' },
		   { "seed", 1, NULL,  's' },
		   { "iobs", 1, NULL,  'o' },
		   { "pdg_1",1, NULL,  'p' },
		   { "pdg_4",1, NULL,  'f' },
		   { "chan", 1, NULL,  'c' },
		   { NULL, 0, NULL, 0 } };
	int next_option;
	do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
		switch (next_option) {
			case 's':
				seed = stoi(optarg, NULL);
				break;						
			case 'o':
				iobs = stoi(optarg, NULL);
				break;
			case 'c':
				ichan = stoi(optarg, NULL);
				break;				
			case 'p':
				flav_1 = stoi(optarg, NULL);
				break;
			case 'f':
				flav_4 = stoi(optarg, NULL);
				break;															
			case '?':
				print_usage();
			case -1: break;
			default: abort();
		}
	}
	while (next_option != -1);
	return;
}


// Use a general CUBA/Vegas interface for all processes
// Control it with some globally defined variables (correction type:LO,R, ... ; calculation type:M,0,L; ...)
int Vegas_Interface(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
	// The integrand value (thing to be integrated), initialise to zero

	double dsigma_summed(0.);

	// Define some variables to control the phase-space and random variable setup
	int nrandom(0); // No additional integrations required for parton level computations

	// 2to3 if target is a photon, otherwise 2to2
	int nfinal = ( channel > 100 )? 3: 2;

	// Create array of masses of final state particles
	double masses[nfinal] = {0.};

	// For two-to-two scattering set the outgoing fermion mass if relevant
	if( nfinal == 2 ){
		double m_fermion = get_mass_fermion_pdg( pdg_fermion );
		masses[0] = m_fermion;
		masses[1] = m_fermion;

		// Special cases
		//channel 12: nu1 + nu2bar > l1 l2bar
		if( channel == 11 ){
			// Find mass of l1 from correlation with the neutrino projectile
			masses[0] = get_mass_fermion_pdg( abs(pdg_projectile)-1 );
		}
		//channels 13,14: nu + gamma > W l
		else if( channel == 13 or channel == 14 ){
			masses[0] = mw;
			masses[1] = m_fermion;
		}
	}

	// Assign masses based channel selection
	// For now fix to muon mass
	if( nfinal == 3 ){
		// Collect the mass information on the outgoing fermion
		double m_fermion = get_mass_fermion_pdg( pdg_fermion );

		// The channel break down is:
		// 101-105: nu gamma > neutrino f fbar [NC type]
		if( channel < 105 ){
			masses[0] = 0.0;
			masses[1] = m_fermion;
			masses[2] = m_fermion;
		}
		// 105+: nu gamma > lepton f f'~ [where f,f' may be quarks or leptons]
		// if f f'~ are quarks, ignore masses
		// if f f~ are leptons, include charged lepton mass effects
		else{
			// Get the mass of lepton particle 3 (it is correlated with the incoming neutrino flavour)
			// Particle 3 = the charged lepton mass correlated with initial neutrino flavour
			double m_3 = get_mass_fermion_pdg( abs(pdg_projectile)-1 );
			// If particle 4,5 are quarks ignore masses
			double m_4 = 0.0;
			double m_5 = 0.0;
			// If particle 4 is a neutrino, include charged lepton mass
			if( is_neutrino( pdg_fermion) ){
				// e.g. electron neutrino (pdg = 12), electron (pdg = 11)
				m_5 = get_mass_fermion_pdg( abs(pdg_fermion)-1 );
			}
			masses[0] = m_3;
			masses[1] = m_4;
			masses[2] = m_5;
		}
	}

	// Create phase-space
	KinematicData Kin = Generate_Phase_Space( xx, nfinal, masses, nrandom, "ee" );
	// Manually set as(muR) currently not used
	Kin.set_as_mur( 0.118 );


	// New function ordered by channels (integers)
	dsigma_summed = dsigma_channels( Kin, channel );

	// Return the (integrand) differential cross-section in pb
	ff[0] = dsigma_summed * hbarc2;
	return 0;
}



// This is the main program, i.e. the one which is run by the executable
// we can supply it with some inputs at command line (lets just use integers)
int main(int argc, char *argv[])
{
	// Initialise the channels
	init_channels();
	if( argc < 2 ) print_usage();
  // Initialise program defaults
	init_default();
	// Some default values
	channel = 1;
	pdg_projectile = 12;
	pdg_fermion = 12;

  // Initialise some process specific scales/variables
	scale_opt = -1;
	// Scale options
	mur_var = 1.;
	muf_var = 1.;	
  mu0 = mz;
	mu_loop = 100.;	// No results depend on mu_reg value

	// Set the collision environment (pp collisions at LHC 13 TeV)
	Ecms = 100.0;
	Ecms2 = pow(Ecms,2);

	int isetup(0); // If we perhaps want to select some specific processes

	// Print available channels
	print_channels();

	// Now read the specific set of command line arguments
	read_arguments(argc,argv,seed_cache,isetup,channel,pdg_projectile,pdg_fermion);	

	// Update cuba dimensions according to the process
	update_process_dimensions();

	// Print main program settings
	print_settings();

	// Initialise recola settings
	init_recola();

	// Initialises the process registration and constructs <int,str> map
	init_recola_processes();

	////////////
	// Setups //
	////////////
	// Setup 0, most simple inclusive LHC setup
	if( isetup == 0 ){
		cout << "Inclusive cross-section at fixed Ecms = " << Ecms << endl;
	}
	else if( isetup == 1 ){
		cout << "Inclusive cross-section for an energy scan in Ecms\n";
	}


	///////////////////////////
	// Information for Vegas //
	///////////////////////////	
	double warmup_precision = 4e-2;
	double integration_precision = 1e-3;
	int grid_no = 3;
	grid_cache = grid_no;
	int max_evaluations = 1e7;

	// When evaluating virtual corrections, work with reduced precision for NLO coefficient
	if( active_virtual ){
		warmup_precision = 1e-1;
		integration_precision = 1e-2;
	}

	//////////////////////////////
	// Store cross-section data //
	//////////////////////////////
	const std::string s_analysis[2] = {"SigmaIncl","SigmaIncl_Ecms"};
  string outfile = s_analysis[isetup]+"_channel"+to_string(channel);
  ofstream ofile_results;
  // open file
  ofile_results.open(outfile+"_s"+to_string(seed_cache)+".txt");
  // save general program settings
	write_settings(ofile_results,"");

	/////////////////
	// Start timer //
	/////////////////	
  struct timeval t0, t1;
  gettimeofday(&t0,NULL);

  // Fiducial results for isetup < 3
  if( isetup == 0 ){
		// Run a test integral, returns an array with two entries < integral, error >, accessed via test_setup[0] and test_setup[1] respectively.
		array<double,2> sigma_fiducial = {0.};
		//			
		sigma_fiducial = integration::vegasC4(Vegas_Interface, warmup_precision, integration_precision, cuba_dimensions, NULL, seed_cache, grid_cache, max_evaluations );

		// Output the results
		cout << "Integral = " << sigma_fiducial[0] << endl;
		cout << "Error = " << sigma_fiducial[1] << endl;	
		// Save this information to the file
		ofile_results << "# Sigma Fiducial" << endl;
		// Include a function which writes all relevant information to the text file
		ofile_results << sigma_fiducial[0] << "\t" << sigma_fiducial[1] << endl;
	}

	// Set up the energy scan
	if( isetup == 1 ){

		// Either derive Ecms for a varying Enu, or directly fix Ecms
		double Ecms_low = 10.0;
		double Ecms_high = 1e5;
		// Number of bins to consider
		int n_bins = 30;	
		vector<double> Ecms_values = linspace( log(Ecms_low),log(Ecms_high), n_bins);
		// ^^^^^^^^^^^^^^^^^^^^^^^
		// Above to be adjusted to match values required by Gaetano	

		// Vector to store the computed cross-section
		vector< array<double,2> > sigma;
		for( int i(0); i < n_bins; i++ ){
			// Select Ecms
			Ecms = exp(Ecms_values[i]);
			Ecms2 = pow(Ecms,2);
			// Compute the cross-section
			array<double,2> sigma_incl = integration::vegasC4(Vegas_Interface, warmup_precision, integration_precision, cuba_dimensions, NULL, seed_cache, grid_cache, max_evaluations );
			// Save the results in the vector
			sigma.push_back( sigma_incl );
		}

		ofile_results << "# Ecms\tsigma[pb]\tsigma_error[pb]\n";
		// Write the results to the file
		for( unsigned int i=0; i < sigma.size(); i++ ){
		// for( auto i: sigma ){
			array<double,2> sig = sigma[i];
			ofile_results << exp(Ecms_values[i]) << "\t"	<< sig[0] << "\t" << sig[1] << endl;
		}		
	}


	ofile_results.close();

  gettimeofday(&t1,NULL);
  double ta=t1.tv_sec-t0.tv_sec+(t1.tv_usec-t0.tv_usec)*0.000001;
  cout << "Total time: " << ta << "s" << endl;		

	return 0;
}






	// // Setup 2 is a mass scan, in an energy slice
	// if( isetup == 2 ){
	// 	double m_low = 0.01;
	// 	double m_max = 3.0;
	// 	int n_bins = 15;	
	// 	vector<double> m_bins = linspace( log(m_low),log(m_max), n_bins);
	// 	// Vector to store the computed cross-section
	// 	vector< array<double,2> > sigma;	

	// 	for( int i(0); i < n_bins; i++ ){
	// 		// Select mQ
	// 		mQ = exp(m_bins[i]);
	// 		mH = mQ;
	// 		// Update hadronic info
	// 		hadr_info.mHadron = mQ;

	// 		// Compute the differential cross-section
	// 		array<double,2> sigma_fiducial = {0.};
	// 		if( s_process == "eeQQb" )
	// 			sigma_fiducial = integration::vegasC4(Vegas_Interface_eeQQb, warmup_precision, integration_precision, cuba_dimensions, &hadr_info, seed_cache, grid_cache, max_evaluations );
	// 		// Store the cross-section and MC error in a vector
	// 		sigma.push_back( sigma_fiducial );
	// 	}
	// 	// Save result
	// 	for( unsigned int i=0; i < sigma.size(); i++ ){
	// 	// for( auto i: sigma ){
	// 		array<double,2> sig = sigma[i];
	// 		ofile_results << exp(m_bins[i]) << "\t"	<< sig[0] << "\t" << sig[1] << endl;
	// 	}
	// }

