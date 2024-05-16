#pragma once
#include <cmath>
#include <complex>
#include "tools.hh"

// The general routine for scattering on a neutrino target: provide it with array of pdg codes of the four external particles
extern double dsigma_nu( KinematicData &, const int* );

extern double dsigma_gamma( KinematicData &, const int* );