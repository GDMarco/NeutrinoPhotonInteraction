#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <cuba.h>
#include <vector>
#include "var.hh"

using namespace std;
using namespace Variables;

// four-momentum class
class p4vec{
private:
	std::array<double, 4> p;
public:
	p4vec ();
	p4vec (double, double, double, double);
	~p4vec (){};
	void boost(double, double, double); // boost by {b_x,b_y,b_z}
    void rot(double, int); // rotation of cos[alpha] around axis i
	// Print four-vector (E,px,py,pz)
    void print(){
        std::cout << std::setprecision(10)      
            << p[0] << "\t" 
            << p[1] << "\t" 
            << p[2] << "\t" 
            << p[3] << std::endl;
	}
    double rap (void);
    double eta (void);
    double pT  (void);
    double ET  (void);
    double m2  (void);
    double m   (void);
    double modp(void);
    // For accessing individual four-momentum
    double E  (void);    
    double px (void);
    double py (void);
    double pz (void);
    // Extra useful functions
    double plus (void); // light-cone + = (E+pz)/sqrt(2)
    double minus (void);// light-cone - = (E-pz)/sqrt(2)
    double pi (int);    
    // Addition/Scalar multiplication operator
    p4vec operator+(const p4vec& b);        
    p4vec operator-(const p4vec& b);
    p4vec operator+=(const p4vec& b);
    p4vec operator-=(const p4vec& b);               
    p4vec operator*(double x);
    p4vec operator/(double x);    
};

// Useful kinematic functions
double kallen(const double, const double, const double);
double dotp4(p4vec &, p4vec &);
double dy_ij(p4vec &, p4vec &);
double deta_ij(p4vec &, p4vec &);
double dphi_ij(p4vec &, p4vec &);
// delta R = sqrt(dphi_ij^2 + dy_ij^2)
double dR_ij(p4vec &, p4vec &);
// Lepton specific observables
// Collins-soper angles
double costheta_CS(p4vec &, p4vec &);
double phi_CS(p4vec &, p4vec &);
// Phi-star variable
double phistar(p4vec &, p4vec &);
// Transverse mass
double MT(p4vec &, p4vec &); // For W boson production

// KinematicData class (contains momenta of a phase space point + other kinematic/useful info.)
class KinematicData {
private:
	// nparticles in the scattering process
	int npar;

	// Some generic booleans
	bool pass_cuts = false;
	// Define variable length information via pointers
	// (they will be assigned length when class instance initialised)
	// Set of particle momenta (partonic CoM)
	std::unique_ptr<p4vec[]> pset;
	// Set of particle pdg codes
	std::unique_ptr<int[]> pdgs;
	// Store necessary kinematic information on the phase-space point
	double mu_f = 0.0, mu_r = 0.0, as_mu_r = 0.0;
	double psp_weight = 0.0, Flux = 0.0;
    // Storage of additional random [0,1] variables (for fragmentation)
    std::unique_ptr<double[]> r_i;
    int n_ran = 0; // number of random variables stored

public:
	// Constructors
	KinematicData(int npar);

	// Functions to set private information (rather than directly altering)
	void set_pi(p4vec pin, int i) { pset[i - 1] = pin; } // The partonic momenta   
    // Allocates memory for storing r_i doubles (see tools.cc)
    void activate_ran( int ); // Allocate storage for n random variables
    void set_nran( int i ){ n_ran = i; }
    void set_ri(double ri, int i){ r_i[i-1] = ri; }; // Storage of random variables

	void set_muf(double muf_in) { mu_f = muf_in; }
	void set_mur(double mur_in) { mu_r = mur_in; }
    void set_as_mur(double as_mur_in) { as_mu_r = as_mur_in; }    
	void set_weight(double w_in) { psp_weight = w_in; }
	void set_flux(double flux_in) { Flux = flux_in; }
	// void set_pdg(){} A function which will set pdgs equal to array of integers
	void set_pdg_i(int i, int a) { pdgs[i - 1] = a; }
	void set_cuts(bool accept) { pass_cuts = accept; }

	// Access Kinenamtic information
	double sij(int, int); // access stored SH info
	double pij(int, int); // pij = sij / 2
    double pijext(int, p4vec); // dot production of pi with an external four vector

	// invariants
	// void fill_sij();    // Fill array with invariants for numerical efficiency?

	// Simple return functions to access private information
	p4vec p(int);
    p4vec p_Lab(int, double, double, double);
    int length() { return this->npar; }
    int pdg(int i) { return this->pdgs[i - 1]; }
    int nran() { return this->n_ran; } 
    bool get_cuts() { return this->pass_cuts; }
    // Access the random numbers (if stored)
    double r(int i);

    double muf() { return this->mu_f; }
    double mur() { return this->mu_r; }
    double as_mur() { return this->as_mu_r; }
    double weight() { return this->psp_weight; }
	double flux() { return this->Flux; }
};


// Some additional tools for setting up differential binning
template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

// A function for debugging a phase-space point
void debug_PS_point( KinematicData & );

// General function to generate kinematic data for all collision environments
// Random numbers [], number of final state particles, masses of final state particles, n random variables to be be stored in KinematicData, collision environment
KinematicData Generate_Phase_Space( const cubareal* rand, const int npar_out, const double *masses, const int n_random, const std::string environment );

// This function is designed for massless particles:
// inputs: random numbers, npar_in, npar_out, phasespace setup option (e.g. how the phase-space is setup)

// Same as above but for fixed Ecms2
KinematicData Gen_2to2_Massive( const cubareal*, const int, const double* );
KinematicData Gen_2to3_Massive( const cubareal*, const int, const double* );
KinematicData Gen_2to4_Massive( const cubareal*, const int, const double* );

//// The general kinematic scale routine, sets the scale based on kinematics
double KinematicScales( KinematicData &, int);