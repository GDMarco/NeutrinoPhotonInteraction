#pragma once
#include <cmath>
#include <complex>
#include <array>
#include <map>

// Introduce custom typedef for switching precision
typedef double real_type;

// Declaration of all variable types
// Using units of GeV for masses/energies
namespace Variables {

	// Real/Imaginary projections
	const std::complex<double> im(0.0,1.0);
	const std::complex<double> re(1.0,0.0);	
	const std::complex<double> zero(0.0,0.0);
	const double Finite[3] = {1.,0.,0.};

	// Constants
	const double hbarc2 = 3.8937966e8;
	const double fm2GeV = 5.067731161976211;
	const double pb2GeVm2 = 2.56819e-9;
	const double pi = M_PI;
	const double pisq = M_PI*M_PI;
	const double zeta2 = pisq/6.0;
	const double zeta3 = 1.2020569031595942854;
	const double gammaE = 0.5772156649015328606;
	// Technical cut, a small cut for numerical dealing with numerical instablities
	const double tech_cut = 1e-10;
	const double m2_threshold = 1e-6; // Numerical threshold for considering particles massive/massless
	// Gauge-boson
	const double mw = 80.35797361; // Following pole values to be updated before next publication
	const double gw = 2.084298894;
	const double mz = 91.15348062;
	const double gz = 2.494266379;
	// Fermion-masses
	// leptons
	const double me = 0.5109989461e-3;
	const double mm = 105.65837e-3;
	const double mtau = 1.77686;
	// quarks
	const double md = 0.0;
	const double mu = 0.0;
	const double mc = 1.51;
	const double mb = 4.92;
	const double mt = 173.0;
	const double ms = 0.2;
	//Nucleons
	const double mp = 0.9382720813;
	const double mn = 0.9395654133;
	//Fermis constant
	const double gf = 1.16638e-5;//1.1663787e-5;
	//Couplings
	const double alpha_zero = 0.0072973525664;
	// const double alpha_zero = 0.00729927;
	const std::complex<double> ALPHA_ZERO(alpha_zero,0.0);
	//Declare alpha as static, and assign it within Main program
	const std::complex<double> MW2C(mw*mw,-gw*mw);
	const std::complex<double> MZ2C(mz*mz,-gz*mz);
	const std::complex<double> SW2_OS(1.0-real(MW2C/MZ2C),-imag(MW2C/MZ2C));
	//Define charges
	const double Qd = -1.0/3.0;
	const double Qu =  2.0/3.0;
	const double Ql = -1.0;    
	//Define complex couplings (global variables)
	extern std::complex<double> ALPHA, SW2;
	extern std::complex<double> gLl, gLnu, gLd, gLu;
	extern std::complex<double> gRl,       gRd, gRu;

	// Global variables
	extern double Ecms, Ecms2, Enu;

	// Function to get fermion mass based on pdg number
	double get_mass_fermion_pdg( int );

	// Function to print the global settings
	void print_settings();

	// Function to initalise all default settings
	void init_default();

	// Function to initiate different EW schemes
	void init_scheme(int);

	// Function to initiate Recola settings
	void init_recola();

	// Function that will write the program settings to a text file (passed by reference)
	void write_settings(std::ofstream &, std::string extra_info);

	// Boolean for neutrino pdg code
	bool is_neutrino(int);

	// Global variable to control integration dimensions of integrands	
	extern int cuba_dimensions;
	// Function to assign the global variable based on channel selection
	void update_process_dimensions();

	// -1 = fixed scale, 1 = et_V
	extern int scale_opt;
	// mu0 scale choice
	extern double mu0;
	// mu for loop integrals
	extern double mu_loop;
	// When computing viritual corrections for specific channels
	extern bool active_virtual;
	extern bool active_recola;
	// Multiplication for mu_0 scale
	extern double muf_var, mur_var;


	// Strong coupling
	extern double ALPHAS;
	// QCD constants
	const double nc = 3.0;
	const double cf = 4.0/3.0;
	const double ca = 3.0;//Nc
	const double tr = 1.0/2.0;
	// Define variable nf
	extern double nf_var;
	extern int nf_pdf, nf_as; // Number of active flavours in pdf/alphas running

	// Some process-dependent and channel dependent options
	// extern std::string s_process, s_target, s_projectile;
	// extern std::string s_correction; // LO, R, V, VV, RV, RR

	//// Any cuts on final-state?
	// pT of neutrino
	extern bool active_pTnu_min,	active_pTnu_max;
	extern double pTnu_min, 		pTnu_max;

	// Function to initialise Recola processes
	void init_channels();
	// Prints channel information
	void print_channels();

	// Function to setup particle masses per channel
	void init_channel_information(int, int *pdgs, double *masses);


	extern int channel, pdg_projectile;
	void init_recola_processes();
	// Construct a map between process strings and integers
	extern std::map<int,std::string> process_map;
	extern std::map<int,std::string> process_map_rcl;

}  



