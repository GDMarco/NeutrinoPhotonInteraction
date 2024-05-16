#pragma once
#include <cuba.h>
#include <string>
#include <array>
////////////////////////////////////////////////////////////
// INTEGRALS                                              //
////////////////////////////////////////////////////////////

namespace integration {

	// Multi-D integrals

	// The function is written such that it first does a run of the integration (storing the importance sampling grid), then evaluates the integral
	extern std::array<double,2> vegasC4(int f(const int*, const cubareal*, const int*, cubareal*, void*),double wu_prec, double prec, int NDIM, void*, int iseed, int gridno, int max_eval);	

	// These are helper functions, for combining the results of integrals and dealing with plus distributions
 	extern double plus_distr( std::array<double,5> funcs, double pdf_1, double pdf_z, double z, double x );
 	std::array<double,2> combine2( std::array<double,2> A, std::array<double,2> B);
	std::array<double,2> combine3( std::array<double,2> A, std::array<double,2> B, std::array<double,2> C);

	// Some integration maps, which return mapped variable and jacobian
	std::array<double,2> int_map( const cubareal, int );

	// For now I comment out the 1D integrals, but leave them there if we want to use them, or something similar
	// 1-D gauss integrals
 	// extern double gauss(double f(double,void*), double xmin, double xmax, double prec, double *error=NULL, void *par=NULL);
 	// extern double gauss_inf(double f(double,void*), double xmin, double prec, double *error=NULL, void *par=NULL);
};

