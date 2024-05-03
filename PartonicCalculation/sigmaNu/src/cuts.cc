#include <cmath>
#include <iostream>
#include "var.hh"
#include "tools.hh"
#include "cuts.hh"

using namespace Variables;

// Some particle specific cuts
bool pass_cuts_neutrino( p4vec &pnu ){
        return true;

        // Currently not set up cuts
}



// General function for applying cuts (set true if all cuts are passed)
void apply_cuts( KinematicData &Kin ){

        // How many particles in the Kinematic data structure
        // int nparticles = Kin.length();

        Kin.set_cuts(true);
        return;
}