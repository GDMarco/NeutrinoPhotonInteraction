#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>
#include <sstream>
#include <array>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <cuba.h>
//glashow++ header files
// #include "var.hh"
#include "integration.hh"

using namespace std;

////////////////////////////////////////////////////////////
// INTEGRALS                                              //
////////////////////////////////////////////////////////////


////////////////////////////////
// MULTIDIMENSIONAL INTEGRALS //
////////////////////////////////

// Definitions for CUBA interface
// These definintions are required for CUBA
#define NCOMP 1
#define NVEC 1
//#define USERDATA NULL
#define EPSABS 1e-16
#define VERBOSE 1
#define SMOOTH 0
#define LAST 4
#define SEED 0
// Lower evaluations for warmup
#define MINEVAL_WU 10000
#define MINEVAL 50000
//#define MAXEVAL 4000000000
#define MAXEVAL 5e8//1800000000
long long int maxeval = 2e8; // for ll version
#define GRIDNO 5
#define STATEFILE ""
#define SPIN NULL
#define NBATCH 1000
#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.
#define KEY 0
#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
//#define LDXGIVEN NDIM
#define NEXTRA 0
#define RETAINSTATEFILE FALSE

// General routine for VEGAS, takes integrand of ndim and returns result
// Performs a set-up of the grid first, before re-evaluating the integrand with refined grid
std::array<double,2> integration::vegasC4(int f(const int *ndim, const cubareal xx[], const int *ncomp, 
  cubareal ff[], void *userdata), double wu_prec, double prec, int NDIM, void* data, int iseed, int gridnum, int max_eval ){
  #undef NSTART
  #undef NINCREASE
  #undef GRIDNO
  #undef STATEFILE
  #undef MAXEVAL
  #define NSTART 1e4
  #define NINCREASE 0//100000
  #define GRIDNO gridnum
  #define STATEFILE ""
  #define MAXEVAL max_eval
  // long long int neval;
  int neval;
  int fail;  
  cubareal wu_integral[1], wu_error[1], wu_prob[1];
  // First perform a warmup integration, storing the statefile
  // #define RETAINSTATEFILE TRUE
  Vegas(NDIM, NCOMP, f, data, NVEC,
    wu_prec, EPSABS, VERBOSE+SMOOTH, iseed,
    MINEVAL_WU, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, wu_integral, wu_error, wu_prob);
  cout << endl << "Warm-up run comlete, now performing final integration" << endl;
  #undef NSTART
  #undef NINCREASE
  #define NSTART 2e4
  #define NINCREASE 100000//1000000
  cubareal integral[1], error[1], prob[1]; 
  Vegas(NDIM, NCOMP, f, data, NVEC,
    prec, EPSABS, VERBOSE+SMOOTH, iseed,
    MINEVAL_WU, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, prob);   
  cout << endl << "Production-up run comlete" << endl;
  std::array<double,2> result = {integral[0],error[0]};
  return result;
}

/////////////////////////////////////
// Treatment of plus distributions //
/////////////////////////////////////

// Inputs are:
// A) 5 component array
// 0 f[z] g[z]
// 1 k[z]
// 2 f[1] g[z]
// 3 f[1] G[x]
// 4 h[z] delta(1-z)

// B) PDF[z]
// C) PDF[1.]
// D) z
// E) x

double integration::plus_distr( array<double,5> funcs, double pdf_1, double pdf_z, double x, double z ){
    double jacobian = (1. - x);
    double plus = 
        + ( funcs[0] + funcs[1] ) * pdf_z * jacobian / z    // (f[z]g[z] + k[z]) * pdf(x/z) / z
        - funcs[2] * pdf_1 * jacobian                       // -f[1]g[z] * pdf(x)
        - funcs[3] * pdf_1                                  // -f[1]G[x] * pdf(x)
        + funcs[4] * pdf_1;                                 // +h[z]delta(1-z) * pdf(x)
    return plus;
}

array<double,2> integration::combine2( array<double,2> A, array<double,2> B){
  array<double,2> output = {0.};
  output[0] = A[0]+B[0];
  output[1] = pow(pow(A[1],2) + pow(B[1],2),0.5);
  if( output[1] < 1e-16 ) output[1] = 0;
  return output;
}

array<double,2> integration::combine3( array<double,2> A, array<double,2> B, array<double,2> C){
  array<double,2> output = {0.};
  output[0] = A[0]+B[0]+C[0];
  output[1] = pow(pow(A[1],2) + pow(B[1],2) + pow(C[1],2),0.5);
  if( output[1] < 1e-16 ) output[1] = 0;  
  return output;
}


// Some integration maps from unit hypercube [0,1] -> [A,B]
// Output is mapped variable, and corresponding jacobian
array<double,2> integration::int_map( const cubareal x, int option ){
    array<double,2> out = {0.};
    // map [0,1] to [0,infinity]
    if( option == 1 ){
        out[0] = x / ( 1. - pow(x,2) );
        out[1] = ( 1. + pow(x,2) ) / pow( 1. - pow(x,2), 2 );
    }
    else if( option == -1 ){
        out[0] = - x / ( 1. - pow(x,2) );
        out[1] = ( 1. + pow(x,2) ) / pow( 1. - pow(x,2), 2 );       
    }
    else{
        cout << "currently only implemented the x/(1+x^2) mapping" << endl;
        abort();
        // 1/x - 1, jacob = 1/x^2
        // tanh[x], j = 1 - tanh[x]^2
        // x = atanh(xr);
        // xjacob = 1. / ( 1. - pow(xr,2) );        
    }
    // Return the mapped variable and jacobian
    return out;
}


////////////////////////////////
// 1-D Integrals with GSL Lib //
////////////////////////////////
// Some 1-d integrals with gauss integration

// double integration::gauss(double f(double,void*), double xmin, double xmax, double prec, double *error, void *par){
//   double result, err;
//   gsl_function F;
//   gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
//   F.function = f;
//   F.params=par;
//   gsl_integration_qags (&F, xmin, xmax, 0, prec, 1000,    w, &result, &err);
//   gsl_integration_workspace_free (w);
//   if(error) *error = err;
//   return result;
// }

// double integration::gauss_inf(double f(double,void*), double xmin, double prec, double *error, void *par){
//   double result, err;
//   gsl_function F;
//   gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
//   F.function = f;
//   F.params=par;
//   gsl_integration_qagiu (&F, xmin, 0, prec, 100000, w, &result, &err);
//   gsl_integration_workspace_free (w);
//   if(error) *error = err;
//   return result;
// }

