#pragma once
#include <cmath>
#include <complex>
#include "tools.hh"
#include "recola.hpp"

// Declaration of all ME2 functions
namespace ME2_Analytic {

    //////////////////////////////////////////////////////////
    //// (anti)neutrino + [anti]neutrino scattering results //
    //////////////////////////////////////////////////////////

    // Massive fermion production in nu nubar annihilation
    extern double nunux_ffx(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2);

    // Same flavour neutrino scattering
    // Master function which takes the pdg codes of the 4 neutrinos, and returns the appropriate |M|^2
    extern double nunu_nunu(KinematicData &Kin, const int pdg_ids[4] );    
    // nu1 + nu1 > nu1 + nu1
    extern double nu1nu1_nu1nu1(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2);
    // nu1 + nu2 > nu1 + nu2 
    extern double nu1nu2_nu1nu2(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2);
    // All neutrino scattering configurations obtained from these two functions via crossings


    ////////////////////////////////////////////////
    // (anti)neutrino + photon scattering results //
    ////////////////////////////////////////////////

    // Neutral Current exchange
    // nu + photon > nu + f fbar [fermion line can be quarks or charged leptons]
    extern double nu1gamma_nu1f2f2x(int i1, int i2, int i3, int i4, int i5, KinematicData &Kin, int f1, int f2);

    // Charged Current exchange

    // nu + photon > nu + f fbar [fermion line can be quarks or charged leptons]
    extern double nu1gamma_l1qqbar(int i1, int i2, int i3, int i4, int i5, KinematicData &Kin, int f1, int f2);

    // nu1 + photon > l1 + nul + lx [different flavour leptons]
    extern double nu1gamma_l1nu2l2x(int i1, int i2, int i3, int i4, int i5, KinematicData &Kin, int f1, int f2);


    // nu1 + photon > l1 + nul1 + l1x [same flavour leptons]
    extern double nu1gamma_l1nu1l1x(int i1, int i2, int i3, int i4, int i5, KinematicData &Kin, int f1, int f2);


    // Recola interface
    extern double compute_process_recola(KinematicData &Kin, int proc, bool virt);

    // Perhaps update to provide it with fixed length array?


    // The below functions were used only for testing the crossing relations
    // nu1 + nu1bar > nu1 + nu1bar [only used for testing the crossing relations]
    // extern double nu1nu1x_nu1nu1x(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2);
    // Different flavour neutrino scattering
    // nu1x + nu2 > nu1x + nu2 [only used for testing the crossing relations]
    // extern double nu1xnu2_nu1xnu2(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2);



}