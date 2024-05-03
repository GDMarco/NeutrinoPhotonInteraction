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
#include "dsigma_nu.hh"
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
		<< " -s --seed			<s>		seed number for integration"  << endl
		<< " -t --targ			<t>		target particle (nu,nux,gamma)" << endl		
		<< " -p --proc			<p>   process (nunu or nugamma)" << endl
		<< " -c --corr			<c>		correction (LO supported)" << endl		
		<< " -o --iobs			<o>		potentially if we consider some specific observeables or channel?" << endl
		<< endl;
	exit(0);
	return;
}

// Introduce a more user friendly interface
void read_arguments(int argc, char* argv[], int &seed, string &targ, string &proc, string &corr, int &iobs ) {
	const char* const short_options = "s:t:p:c:o:";
	const struct option long_options[] = { { "help", 0, NULL, 'h' },
		   { "seed", 1, NULL,  's' },
		   { "targ", 1, NULL,  't' },
		   { "proc", 1, NULL,  'p' },		      
		   { "corr", 1, NULL,  'c' },
		   { "iobs", 1, NULL,  'o' },
		   { NULL, 0, NULL, 0 } };
	int next_option;
	ostringstream s_corr, s_proc, s_targ;	
	do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
		switch (next_option) {
			case 's':
				seed = stoi(optarg, NULL);
				break;
			case 't':
				s_targ << optarg;
				corr = s_targ.str();
				break;				
			case 'p':
				s_proc << optarg;
				proc = s_proc.str();
				break;
			case 'c':
				s_corr << optarg;
				corr = s_corr.str();
				break;							
			case 'o':
				iobs = stoi(optarg, NULL);
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

	if( s_process != "nunu" and s_process != "nugamma" ){
		cerr << "Vegas_Interface:: unknown process selection s_process = " << s_process << endl;
		abort();
	}

	// Define some variables to control the phase-space and random variable setup
	int nrandom(0); // No additional integrations required for parton level computations
	// Final-state particle multiplicity (based on correction type)
	int nfinal = nfinal_2to2_Born( s_correction );

	// Create array of masses
	double masses[nfinal] = {0.};

	// Assign these particle masses according to the selected subprocess
	// For now take them to be massless

	// Create phase-space
	KinematicData Kin = Generate_Phase_Space( xx, nfinal, masses, nrandom, "ee" );
	// Manually set as(muR) probably not needed
	Kin.set_as_mur( 0.118 );

	if( s_process == "nunu" ){
		//
		dsigma_summed = dsigma_nunu_LO( Kin );
	}

	// ...

	// Return the (integrand) differential cross-section in pb
	ff[0] = dsigma_summed * hbarc2;
	return 0;
}



// This is the main program, i.e. the one which is run by the executable
// we can supply it with some inputs at command line (lets just use integers)
int main(int argc, char *argv[])
{
	if( argc < 2 ) print_usage();
  // Initialise program defaults
	init_default();


	// Set up specific values for some global definitions
	s_process = "nunu"; 	// heavy-quark pair production
	s_target  = "nu"; // nu, nux, gamma
  // Initialise some process specific scales/variables
	scale_opt = -1;
	// Scale options
	mur_var = 1.;
	muf_var = 1.;	
  mu0 = mz;
	mu_loop = 100.;	// No results depend on mu_reg value

	// Set the collision environment (pp collisions at LHC 13 TeV)
	Ecms = 20.0;
	Ecms2 = pow(Ecms,2);

	int isetup(0); // If we perhaps want to select some specific processes

	// Now read the specific set of command line arguments
	read_arguments(argc,argv,seed_cache,s_target,s_process,s_correction,isetup);

	// Update the number of CUBA dimensions based on the assigned type of correction, process etc.
	update_process_dimensions();

	// Print main program settings
	print_settings();

	////////////
	// Setups //
	////////////
	// Setup 0, most simple inclusive LHC setup
	if( isetup == 0 ){
		cout << "Inclusive cross-section at fixed Ecms = " << Ecms << endl;
	}
	else if( isetup == 1 ){
		cout << "Inclusive cross-section for a scan in Ecms\n";
	}


	///////////////////////////
	// Information for Vegas //
	///////////////////////////	
	double warmup_precision = 4e-2;
	double integration_precision = 1e-3;
	int grid_no = 3;
	grid_cache = grid_no;
	int max_evaluations = 1e7;

	//////////////////////////////
	// Store cross-section data //
	//////////////////////////////
	const std::string s_analysis[2] = {"sigma","sigma_Ecms"};
  string outfile = s_analysis[isetup]+"_" + s_process + "_" + s_correction;
  ofstream ofile_results;
  ofile_results.open(outfile+"_s"+to_string(seed_cache)+".txt");
	write_settings(ofile_results,"no pdfs",s_process);

	/////////////////
	// Start timer //
	/////////////////	
  struct timeval t0, t1;
  gettimeofday(&t0,NULL);

  // Fiducial results for isetup < 3
  if( isetup < 2 ){
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

