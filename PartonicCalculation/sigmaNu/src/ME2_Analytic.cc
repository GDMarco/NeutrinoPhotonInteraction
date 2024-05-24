#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include "var.hh"
#include "tools.hh"
#include "ME2_Analytic.hh"
#include "recola.hpp"
// Scalar one-loop integrals from OneLOop external libraries
// #include "cavh_olo.h" // Link to OneLOop

using namespace std;
using namespace Variables;
using namespace ME2_Analytic;

// Expressions for heavy-quark pair production in e-(p1) + e+(p2) > gamma/Z > Q(p3) + Qbar(p4) [ + g(p5) ]
void assign_EW_charges(double &Q, complex<double> &gL, complex<double> &gR, int pdg){
	if (pdg==11 or pdg==13 or pdg==15){
		gL = gLl;
		gR = gRl;
		Q = Ql;
	}
	else if ( pdg==12 or pdg==14 or pdg==16){
		gL = gLnu;
		gR = 0;
		Q = 0;
	}	
	else if (pdg==1 or pdg==3 or pdg==5){   //quarks as incoming particles to check q+qb->l+lb process
		gL = gLd;
		gR = gRd;
		Q = Qd;
	}
	else if (pdg==2 or pdg==4 or pdg==6){
		gL = gLu;
		gR = gRu;
		Q = Qu;
	}
	else{
		cerr << "assign_EW_charges: unsupported pdf code " << pdg << endl;
		abort();
	}
	return;
}


//////////////////////////////////////////////////////////
//// (anti)neutrino + [anti]neutrino scattering results //
//////////////////////////////////////////////////////////


// nu nubar > Q Qbar (applicable to all fermion lines, and accounts for nc factors)
double ME2_Analytic::nunux_ffx(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2 ){

	if( !is_neutrino(f1) ){
		cerr << "nunu_ffx: fermion line 1 is not a neutrino\n";
		abort();
	}
	// f1 line should always be the neutrino line, check this

	// Check spin sum for special case of different flavour neutrinos

	// Kinematics
	double s12 = 2.0*Kin.pij(i1,i2);
	double s13 = 2.0*Kin.pij(i1,i3);
	double mfsq = Kin.p(i3).m2();
	// Propagators
	complex<double> chiz = 1.0 / ( s12 - MZ2C );
	// f3, f4 a generic fermion line
	complex<double> gLf, gRf; double Qf;
	assign_EW_charges(Qf,gLf,gRf,abs(f2));

	complex<double> ME2 = 48.* ALPHA * conj(ALPHA) * pisq *chiz*gLnu*conj(chiz)*conj(gLnu) * 
	( conj(gLf)*(gRf*mfsq*s12 + gLf*pow(s12 - s13,2)) + conj(gRf)*(gLf*mfsq*s12 + gRf*pow(s13,2)) );

	// Remove colour factor
	if( abs(f2) > 6 ) ME2 /= 3.0;

	return ME2.real();
}


// Master function to access all the various scattering configurations in 4 neutrino combinations
double ME2_Analytic::nunu_nunu(KinematicData &Kin, const int pdg_ids[4] ){

	// equal flavour 4 neutrino cases
	if( abs(pdg_ids[0]) == abs(pdg_ids[1]) && abs(pdg_ids[0]) == abs(pdg_ids[2]) ){

		// nu nu nu nu
		if( pdg_ids[0] > 0 && pdg_ids[1] > 0 ){
			return nu1nu1_nu1nu1(1,2,3,4, Kin, pdg_ids[0], pdg_ids[0]);
		}
		// nubar nubar > nubar nubar
		else if( pdg_ids[0] < 0 && pdg_ids[1] < 0 ){
			return nu1nu1_nu1nu1(2,1,4,3, Kin, pdg_ids[0], pdg_ids[0]);
		}
		// else the mixed cases: nu nubar > X + c.c.
		else{
			// nu nubar > X
			if( pdg_ids[0] > 0 ){
				// nu nubar > nu nbar [using line reversal]
				if( pdg_ids[2] > 0 ) return nu1nu1_nu1nu1(1,4,3,2, Kin, pdg_ids[0], pdg_ids[0]);
				// nu nubar > nubar nu
				else{
					return nu1nu1_nu1nu1(1,3,2,4, Kin, pdg_ids[0], pdg_ids[0]);
				}
			}
			// nubar nu > X
			else{
				// nubar nu > nu nubar
				if( pdg_ids[2] > 0 ) return nu1nu1_nu1nu1(4,2,3,1, Kin, pdg_ids[0], pdg_ids[0]);
				// nubar nu > nubar nu
				else{
					return nu1nu1_nu1nu1(3,2,1,4, Kin, pdg_ids[0], pdg_ids[0]);
				}
			}
		}
	}
	// Flavour annihilation channel: nu1 nu1x > nu2 nu2x
	else if( abs(pdg_ids[0]) == abs(pdg_ids[1]) and abs(pdg_ids[2]) == abs(pdg_ids[3]) ){
        // Set up correct flavours
        int f1 = (pdg_ids[0] > 0)? 1: 2;
        int f2 = (pdg_ids[1] > 0)? 1: 2;
        // 
        int f3 = (pdg_ids[2] > 0)? 3: 4;
        int f4 = (pdg_ids[3] > 0)? 3: 4; 
        // nu nubar > f fbar
        return ME2_Analytic::nunux_ffx(f1,f2,f3,f4,Kin,abs(pdg_ids[0]),abs(pdg_ids[2]));
	}
	// Otherwise we are in the mixed flavour cases
	else{

		// nu1 + nu2
		if( pdg_ids[0] > 0  and pdg_ids[1] > 0 ){
			// nu1 + nu2 > nu1 + nu2 
			if( pdg_ids[0] == pdg_ids[2] )
				return nu1nu2_nu1nu2(1,2,3,4,Kin,pdg_ids[0],pdg_ids[1]);
			// nu1 + nu2 > nu2 + nu1
			else if( pdg_ids[0] == pdg_ids[3] )
				return nu1nu2_nu1nu2(1,2,4,3,Kin,pdg_ids[0],pdg_ids[1]);
		}
		// nu1 + nu2x
		else if( pdg_ids[0] > 0  and pdg_ids[1] < 0 ){
			// nu1 + nu2x > nu1 + nu2x
			if( pdg_ids[0] == pdg_ids[2] )
				return nu1nu2_nu1nu2(1,4,3,2,Kin,pdg_ids[0],pdg_ids[1]);
			// nu1 + nu2x > nu2x + nu1
			else if( pdg_ids[0] == pdg_ids[3] )
				return nu1nu2_nu1nu2(1,2,4,3,Kin,pdg_ids[0],pdg_ids[1]);
		}
		// nu1x + nu2
		else if( pdg_ids[0] < 0  and pdg_ids[1] > 0 ){
			// nu1x + nu2 > nu1x + nu2
			if( pdg_ids[0] == pdg_ids[2] )
				return nu1nu2_nu1nu2(3,2,1,4,Kin,pdg_ids[0],pdg_ids[1]);
			// nu1x + nu2 > nu2 + nu1x
			else if( pdg_ids[0] == pdg_ids[3] )
				return nu1nu2_nu1nu2(4,2,1,3,Kin,pdg_ids[0],pdg_ids[1]);
		}
		// nu1x nu2x
		else if( pdg_ids[0] < 0  and pdg_ids[1] < 0 ){
			// nu1x + nu2x > nu1x + nu2x
			if( pdg_ids[0] == pdg_ids[2] )
				return nu1nu2_nu1nu2(3,4,1,2,Kin,pdg_ids[0],pdg_ids[1]);
			// nu1x + nu2x > nu2x + nu1x
			else if( pdg_ids[0] == pdg_ids[3] )
				return nu1nu2_nu1nu2(4,3,1,2,Kin,pdg_ids[0],pdg_ids[1]);
		}

		// One should never reach here if I accounted for all configurations
		cout << pdg_ids[0] << endl;
		cout << pdg_ids[1] << endl;
		cout << pdg_ids[2] << endl;
		cout << pdg_ids[3] << endl;
		cerr << "nunu_nunu: should not reach here\n";
		abort();
		return 0.0;
	}
}

// Same flavour neutrino scattering

// nu1 + nu1 > nu1 + nu1
double ME2_Analytic::nu1nu1_nu1nu1(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2 ){

	// f1 != f2

	// Analytic continuation of momenta
	double ac_i1 = (i1 < 3)?  +1.: -1.;
	double ac_i2 = (i2 < 3)?  +1.: -1.;
	double ac_i3 = (i3 < 3)?  -1.: +1.;
	double ac_i4 = (i4 < 3)?  -1.: +1.;
	// Kinematics
	double s12 = ac_i1 * ac_i2 * 2.0*Kin.pij(i1,i2);
	double s13 = ac_i1 * ac_i3 * 2.0*Kin.pij(i1,i3);

	complex<double> ME2 = 16.*ALPHA*(2.*MZ2C + s12)*conj(ALPHA)*(s12 + 2.*conj(MZ2C))*pow(gLnu,2)*pow(pi,2)*pow(s12,2)*
	pow(MZ2C + s12 - s13,-1)*pow(MZ2C + s13,-1)*pow(conj(gLnu),2)*pow(s12 - s13 + conj(MZ2C),-1)*
	pow(s13 + conj(MZ2C),-1);

	// Must take care when equal flavour neutrinos in final state
	double averaging(1.);
	if( ac_i3 * ac_i4 > 0 ) averaging /= 2.;
	return ME2.real() * averaging;
}


// Different flavour neutrino scattering

// nu1 + nu2 > nu1 + nu2
double ME2_Analytic::nu1nu2_nu1nu2(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2 ){

	// Analytic continuation of momenta
	double ac_i1 = (i1 < 3)?  +1.: -1.;
	double ac_i2 = (i2 < 3)?  +1.: -1.;
	double ac_i3 = (i3 < 3)?  -1.: +1.;
	// Kinematics
	double s12 = ac_i1 * ac_i2 * 2.0*Kin.pij(i1,i2);
	double s13 = ac_i1 * ac_i3 * 2.0*Kin.pij(i1,i3);

	complex<double> ME2 = 16.*ALPHA*conj(ALPHA)*pow(gLnu,2)*pow(pi,2)*pow(s12,2)*pow(MZ2C + s13,-1)*pow(conj(gLnu),2)*
   pow(s13 + conj(MZ2C),-1);

   return ME2.real();
}

// nu1 + nu1bar > l1 + l1bar
double ME2_Analytic::nu1nu1bar_l1l1bar(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2 ){

	// Kinematics
	double s12 = 2.0*Kin.pij(i1,i2);
	double s13 = 2.0*Kin.pij(i1,i3);
	double mfsq = Kin.p(i3).m2();

	// f3, f4 a generic fermion line
	complex<double> gLf, gRf; double Qf;
	assign_EW_charges(Qf,gLf,gRf,abs(f2));	

	complex<double> ME2 = -2.*ALPHA*conj(ALPHA)*pow(MW2C,-1)*pow(pi,2)*pow(MZ2C - s12,-1)*(conj(MW2C)*(s12 - conj(MZ2C) + 2.*conj(gLf)*conj(gLnu)*(mfsq - s13 - conj(MW2C))*conj(SW2))*
      (mfsq*s12*(mfsq*(MZ2C - s12) + 4.*gLnu*gRf*MW2C*(-mfsq + MW2C + s13)*SW2) + 2.*MW2C*(MZ2C - s12 + 2.*gLf*gLnu*(-mfsq + MW2C + s13)*SW2)*pow(s12 - s13,2)) + 
     (mfsq*s12 - mfsq*conj(MZ2C) + 4.*conj(gLnu)*conj(gRf)*(mfsq - s13 - conj(MW2C))*conj(MW2C)*conj(SW2))*
      (mfsq*MW2C*s12*(MZ2C - s12 + 2.*gLf*gLnu*(-mfsq + MW2C + s13)*SW2) + 
        ((4.*gLnu*gRf*MW2C*(MW2C + s13)*SW2 + mfsq*(MZ2C - s12 - 4.*gLnu*gRf*MW2C*SW2))*pow(s13,2))/2.))*pow(-mfsq + MW2C + s13,-1)*pow(SW2,-1)*
   pow(mfsq - s13 - conj(MW2C),-1)*pow(conj(MW2C),-1)*pow(s12 - conj(MZ2C),-1)*pow(conj(SW2),-1);
	return ME2.real();
}

// nu1 + nu2bar > l1 + l2bar
double ME2_Analytic::nu1nu2bar_l1l2bar(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2 ){
	// Kinematics
	double s12 = 2.0*Kin.pij(i1,i2);
	double s13 = 2.0*Kin.pij(i1,i3);
	double ml1sq = Kin.p(i3).m2();
	double ml2sq = Kin.p(i4).m2();
	complex<double> ME2 = ALPHA*conj(ALPHA)*(-2.*ml1sq*ml2sq*MW2C*s12 + ml1sq*ml2sq*(ml1sq - ml2sq - s13)*s13 - 2.*ml1sq*(ml2sq*s12 + 2.*MW2C*(s12 - s13))*conj(MW2C) + 
     4.*MW2C*(s12 - s13)*(ml2sq - s12 + s13)*conj(MW2C))*pow(MW2C,-1)*pow(pi,2)*pow(-ml1sq + MW2C + s13,-1)*pow(SW2,-1)*pow(ml1sq - s13 - conj(MW2C),-1)*
   pow(conj(MW2C),-1)*pow(conj(SW2),-1);
	return ME2.real();
}

// nu1 + gamma > W+ + l1
double ME2_Analytic::nugumma_Wl(int i1, int i2, int i3, int i4, KinematicData &Kin){

	// Kinematics
	double s12 = 2.0*Kin.pij(i1,i2);
	double s13 = 2.0*Kin.pij(i1,i3);
	double mw2r = pow(mw,2);
	double mlsq = Kin.p(i4).m2();
	complex<double> sw2 = SW2.real();
	complex<double> ME2 = 8.*ALPHA*s12*conj(ALPHA)*conj(pow(sw2,-1/2.))*pow(mw2r,-1)*pow(pi,2)*
   ((-mw2r - 2*s12 + s13)*pow(mlsq,3) + pow(mlsq,4) + pow(mlsq,2)*(2*mw2r*s12 - 3*pow(mw2r,2) + pow(s12 - s13,2)) - 
     2*mw2r*(mw2r - s13)*(-2*mw2r*s12 + pow(mw2r,2) + pow(s12,2) + pow(s13,2)) + 
     mlsq*(-((4*s12 + 3*s13)*pow(mw2r,2)) + 5*pow(mw2r,3) + s13*pow(s12 - s13,2) + mw2r*(6*s12*s13 + pow(s12,2) + pow(s13,2))))*pow(mlsq - mw2r + s13,-2)*
   pow(mlsq - mw2r - s12 + s13,-2)*pow(sw2,-1/2.);
// 	complex<double> ME2 = 4.*ALPHA*s12*conj(ALPHA)*conj(pow(SW2,-1/2.))*pow(MW2C,-1)*pow(pi,2)*pow(-mlsq + MW2C - s13,-1)*pow(-mlsq + MW2C + s12 - s13,-1)*
   // ((-4*s12 + 2*s13)*pow(mlsq,3) + 2*pow(mlsq,4) + pow(mlsq,2)*(MW2C*s13 - 6.*pow(MW2C,2) + 2*pow(s12 - s13,2)) + 
   //   2.*MW2C*((2*s12 + s13)*pow(MW2C,2) + 2*s13*(pow(s12,2) + pow(s13,2)) - MW2C*(3*s12*s13 + 2*pow(s12,2) + pow(s13,2))) - 
   //   conj(MW2C)*((-4*s12 + s13)*pow(mlsq,2) + 2*pow(mlsq,3) + 2.*MW2C*(s13*(s12 + s13) - MW2C*(2*s12 + s13) + 2.*pow(MW2C,2)) + 
   //      mlsq*(-(s12*s13) + MW2C*(2*s12 + s13) - 6.*pow(MW2C,2) + 2*pow(s12,2) + pow(s13,2))) + 
   //   mlsq*(-((6*s12 + 5*s13)*pow(MW2C,2)) + 4.*pow(MW2C,3) + 2*s13*pow(s12 - s13,2) + MW2C*(11*s12*s13 + 4*pow(s12,2) + 3*pow(s13,2))))*pow(SW2,-1/2.)*
   // pow(mlsq + s13 - conj(MW2C),-1)*pow(mlsq - s12 + s13 - conj(MW2C),-1);
   // Photon averaging
   return ME2.real() / 2.;
}


////////////////////////////////////////////////
// (anti)neutrino + photon scattering results //
////////////////////////////////////////////////

double ME2_Analytic::nu1gamma_nu1f2f2x(int i1, int i2, int i3, int i4, int i5, KinematicData &Kin, int f1, int f2){

	if( !is_neutrino(f1) ){
		cerr << "nu1gamma_nu1f2f2x: fermion line 1 is not a neutrino\n";
		abort();
	}

	// Allow the analytic continuation for momenta 1 and 3 (i.e. for the massless neutrino line)
	// Analytic continuation of momenta
	double ac_i1 = (i1 < 3)?  +1.: -1.;
	double ac_i3 = (i3 < 3)?  -1.: +1.;
	// Take fermion mass from kinematic structure
	double mfsq = Kin.p(i4).m2();
	// Other kinematics
	double s12 = ac_i1 * 2.0*Kin.pij(i1,i2);
	double s13 = ac_i1 * ac_i3 * 2.0*Kin.pij(i1,i3);
	double s15 = ac_i1 * 2.0*Kin.pij(i1,i5);
	double s24 = 2.0*Kin.pij(i2,i4);
	double s25 = 2.0*Kin.pij(i2,i5);

	// f3, f4 a generic fermion line
	complex<double> gLf, gRf; double Qf;
	assign_EW_charges(Qf,gLf,gRf,abs(f2));

	// Include mass for quarks. The f2 f2x quark system should also be at least as massive as pion
	complex<double> chiz13 = 1. / ( MZ2C + s13 );
	// alpha alpha*  (Re[alpha])
	complex<double> prefactor = 128.*chiz13*gLnu*conj(chiz13)*conj(gLnu)*pow(pi,3)*pow(Qf,2)*pow(s24,-2)*pow(s25,-2)
	 * ALPHA * conj(ALPHA) * ALPHA.real();

	complex<double> ME2 = (conj(gRf)*(2.*gLf*mfsq*s24*s25*(s12*(s24 + s25) + s13*(-s13 + s24 + s25) - pow(s12,2)) + 
        2.*gRf*mfsq*(s24 + s25)*(-(s12*s24*(2*s13 + 2*s15 + s25)) + s24*pow(s12,2) + 
           (s24 + s25)*pow(s13 + s15,2)) + 
        gRf*s13*s24*s25*(4*s13*s15 - 2*s12*(s13 + s15) - 2*s13*s24 - 2*s15*s24 + pow(s12,2) + 
           2*pow(s13,2) + 2*pow(s15,2) + pow(s24,2)) - 2.*gLf*s13*pow(mfsq,2)*pow(s24 + s25,2)) + 
     conj(gLf)*(2.*gRf*mfsq*s24*s25*(s12*(s24 + s25) + s13*(-s13 + s24 + s25) - pow(s12,2)) + 
        2.*gLf*mfsq*(s24 + s25)*(-(s12*s24*(2*s15 + s25)) + s24*pow(s12,2) + (s24 + s25)*pow(s15,2)) + 
        gLf*s13*s24*s25*(2*s15*s25 - 2*s12*(s15 + s25) + pow(s12,2) + 2*pow(s15,2) + pow(s25,2)) - 
        2.*gRf*s13*pow(mfsq,2)*pow(s24 + s25,2)));

	// Averaging over photon polarisations
	prefactor /= 2.;
	// Account for colour sum
	if( abs(f2) < 7 ) prefactor *= 3.;
	// Include the prefactor and averaging/sum factors
	ME2 *= prefactor;
	return ME2.real();
	// Validated against Recola up to a factor of two
}



double ME2_Analytic::nu1gamma_l1qqbar(int i1, int i2, int i3, int i4, int i5, KinematicData &Kin, int f1, int f2){

	if( !is_neutrino(f1) ){
		cerr << "nu1gamma_nu1f2f2x: fermion line 1 is not a neutrino\n";
		abort();
	}

	// Allow the analytic continuation for momenta 1 and 3 (i.e. for the massless neutrino line)
	// Analytic continuation of momenta
	// double ac_i1 = (i1 < 3)?  +1.: -1.;
	// double ac_i3 = (i3 < 3)?  -1.: +1.;
	// Take fermion mass from kinematic structure
	// double mfsq = Kin.p(i4).m2();
	// Other kinematics
	// double s13 = ac_i1 * ac_i3 * 2.0*Kin.pij(i1,i3);
	// double s15 = ac_i1 * 2.0*Kin.pij(i1,i5);
	// double s24 = 2.0*Kin.pij(i2,i4);
	// double s25 = 2.0*Kin.pij(i2,i5);
	// double s34 = 2.0*Kin.pij(i3,i4);

	// Using another basis of kinematics
	double s12 = 2.0*Kin.pij(i1,i2);
	double s14 = 2.0*Kin.pij(i1,i4);	
	double s23 = 2.0*Kin.pij(i2,i3);
	double s35 = 2.0*Kin.pij(i3,i5);		
	double s45 = 2.0*Kin.pij(i4,i5);
	// fermion masses
	double mlsq = Kin.p(3).m2();

	// f3, f4 a generic fermion line
	// complex<double> gLf, gRf; double Qf;
	// assign_EW_charges(Qf,gLf,gRf,abs(f2));

	// Resonant propagator
	complex<double> chiw45 = 1. / (s45 - MW2C);
	complex<double> mw2c = MW2C;
	// The following option to drop width effects in numerator (a sort of pole approximation)
	// mw2c = pow(mw,2);
	// Non-resonant propagator
	complex<double> chiw13 = 1. / (mw2c + s12 - s23 - s45 );
	//	
	complex<double> prefactor = 4.*chiw13*chiw45*conj(ALPHA)*conj(chiw13)*conj(chiw45)*pow(ALPHA,2)*pow(pi,3);

	complex<double> ME2 = (16*pow(s23,-2)*pow(mlsq - s12 + s14 + s35,-1)*pow(mlsq + s14 - s23 + s35,-1)*pow(SW2,-1)*
     (9*pow(mlsq,4)*(2*(s14 - s45)*pow(s12,2) - 
          s12*(s23*s35 - s23*s45 + s14*(s23 + 4*s45) - 4*pow(s45,2)) + 
          s45*(s23*s35 - s23*s45 + s14*(s23 + 2*s45) - pow(s23,2) - 2*pow(s45,2))) - 
       3*pow(mlsq,3)*(6*(s14 - s45)*pow(s12,3) + 
          pow(s12,2)*(6*s45*(2*s35 + s45) + 3*s14*(3*s23 - 4*s35 + 4*s45) - s23*(6*s35 + 11*s45) - 
             18*pow(s14,2)) - s45*(9*(s23 + 2*s45)*pow(s14,2) + 
             s14*(6*s23*(2*s35 - 5*s45) + 12*(s35 - 2*s45)*s45 - 17*pow(s23,2)) + 
             (-11*s35 + 14*s45)*pow(s23,2) + 8*pow(s23,3) + 6*(-2*s35 + s45)*pow(s45,2) + 
             s23*(-(s35*s45) + 9*pow(s35,2) + 23*pow(s45,2))) + 
          s12*(9*(s23 + 4*s45)*pow(s14,2) + (-3*s35 + 7*s45)*pow(s23,2) + 
             6*(-4*s35 + s45)*pow(s45,2) - 
             3*s14*(-4*s23*s35 + 13*s23*s45 - 8*s35*s45 + 2*pow(s23,2) + 14*pow(s45,2)) + 
             s23*(5*s35*s45 + 9*pow(s35,2) + 34*pow(s45,2)))) + 
       s23*(s12 - s23 - s45)*s45*((3*s14*(7*s35 + 6*s45) - (2*s23 - 3*s35)*(9*s35 + 8*s45))*
           pow(s12,2) + (-6*s14 + 6*s23 - 9*s35)*pow(s12,3) + 
          (3*s14 - 2*s23 + 3*s35)*(3*s35 + 2*s45)*
           (2*s23*s45 + 2*s35*s45 - 2*s14*(s23 + s45) + pow(s14,2) + pow(s23,2) + pow(s35,2) + 
             2*pow(s45,2)) + s12*(3*(6*s23 - 3*s35 + 4*s45)*pow(s14,2) - 6*pow(s14,3) + 
             (2*s23 - 3*s35)*(6*s23*s45 + 16*s35*s45 + 3*pow(s23,2) + 9*pow(s35,2) + 10*pow(s45,2)) - 
             6*s14*(s23*(-3*s35 + 4*s45) + 3*pow(s23,2) + 4*(s35*s45 + pow(s35,2) + pow(s45,2))))) + 
       mw2c*conj(mw2c)*(-3*(s12*(6*s14 - s23 - 6*s45) - 
             (3*s14 - s23 - 3*s45)*(6*s14 - 5*s23 + 4*s35 - 2*s45))*pow(mlsq,3) + 
          9*(2*s14 - s23 - 2*s45)*pow(mlsq,4) + 
          pow(mlsq,2)*(84*s23*s35*s45 + 9*s23*pow(s12,2) - 18*(5*s23 - 4*s35 + 5*s45)*pow(s14,2) + 
             54*pow(s14,3) + 14*s35*pow(s23,2) - 26*s45*pow(s23,2) + 3*s23*pow(s35,2) - 
             18*s45*pow(s35,2) - 12*s23*pow(s45,2) + 36*s35*pow(s45,2) - 
             3*s12*(4*s23*s35 + 17*s23*s45 - 6*s35*s45 - 6*s14*(2*s23 - s35 + 3*s45) + 12*pow(s14,2) + 
                2*pow(s23,2) + 6*pow(s45,2)) + 
             s14*(s23*(-75*s35 + 108*s45) + 37*pow(s23,2) + 
                18*(-6*s35*s45 + pow(s35,2) + 2*pow(s45,2)))) + 
          s23*((-3*s14 + 2*s23)*pow(s12,3) + 
             pow(s12,2)*(-12*s14*s23 + 15*s14*s35 - 10*s23*s35 - 6*s35*s45 + 9*pow(s14,2) + 
                4*pow(s23,2)) + (3*s14 - 2*s23 + 3*s35)*(3*s14 - 2*(s23 + s45))*
              (2*s23*s45 + 2*s35*s45 - 2*s14*(s23 + s45) + pow(s14,2) + pow(s23,2) + pow(s35,2) + 
                2*pow(s45,2)) + s12*(2*(4*s23 - 9*s35 - 6*s45)*pow(s14,2) - 3*pow(s14,3) + 
                2*(6*s35*s45*(s35 + s45) - 2*(2*s35 + s45)*pow(s23,2) + pow(s23,3) + 
                   s23*(4*s35*s45 + 7*pow(s35,2) - 2*pow(s45,2))) + 
                s14*(24*s23*s35 + 14*s23*s45 - 12*s35*s45 - 7*pow(s23,2) - 21*pow(s35,2) + 
                   6*pow(s45,2)))) - mlsq*
           (s23*(-18*s14 + 11*s23 - 12*s35)*pow(s12,2) + 3*s23*pow(s12,3) + 
             9*(3*s23 - 4*s35 + 4*s45)*pow(s14,3) - 18*pow(s14,4) + 34*s35*s45*pow(s23,2) - 
             4*s35*pow(s23,3) + 18*s45*pow(s23,3) + 10*pow(s23,4) - 27*s23*s45*pow(s35,2) + 
             7*pow(s23,2)*pow(s35,2) - 6*s23*pow(s35,3) - 
             s14*(36*s35*s45*(-s35 + s45) + 5*(3*s35 + 4*s45)*pow(s23,2) + 27*pow(s23,3) - 
                3*s23*(-36*s35*s45 + pow(s35,2) - 6*pow(s45,2))) + 24*s23*s35*pow(s45,2) + 
             24*pow(s23,2)*pow(s45,2) - 18*pow(s35,2)*pow(s45,2) + 
             pow(s14,2)*(s23*(54*s35 - 33*s45) + 8*pow(s23,2) - 
                18*(-4*s35*s45 + pow(s35,2) + pow(s45,2))) + 
             s12*(-6*(5*s23 - 3*s35 + 6*s45)*pow(s14,2) + 18*pow(s14,3) + 
                s14*(18*s45*(-2*s35 + s45) + 3*s23*(5*s35 + 27*s45) + 13*pow(s23,2)) - 
                3*(6*(s35 + 2*s45)*pow(s23,2) - 6*s35*pow(s45,2) + 
                   s23*(-9*s35*s45 - 5*pow(s35,2) + 8*pow(s45,2)))) + 12*s23*pow(s45,3))) - 
       pow(mlsq,2)*(9*pow(s12,3)*(s14*(-4*s23 + 2*s35 - 6*s45) + 2*s45*(-s35 + s45) + 
             s23*(s35 + 4*s45) + 4*pow(s14,2)) - 
          s45*(27*(s23 + 2*s45)*pow(s14,3) - 
             3*pow(s14,2)*(6*s45*(-4*s35 + 5*s45) + s23*(-15*s35 + 51*s45) + 25*pow(s23,2)) + 
             (48*s35 - 45*s45)*pow(s23,3) - 21*pow(s23,4) + 
             pow(s23,2)*(44*s35*s45 - 48*pow(s35,2) - 74*pow(s45,2)) - 
             18*s35*(s35 - 2*s45)*pow(s45,2) + 
             s14*((-93*s35 + 145*s45)*pow(s23,2) + 69*pow(s23,3) + 
                18*s45*(-6*s35*s45 + pow(s35,2) + 2*pow(s45,2)) + 
                3*s23*(-37*s35*s45 + 15*pow(s35,2) + 54*pow(s45,2))) + 
             3*s23*(4*s45*pow(s35,2) + 9*pow(s35,3) + 34*s35*pow(s45,2) - 10*pow(s45,3))) - 
          3*pow(s12,2)*(22*s23*s35*s45 - 6*(3*s23 - 4*s35 + s45)*pow(s14,2) + 18*pow(s14,3) - 
             6*s35*pow(s23,2) + 12*s23*pow(s35,2) - 6*s45*pow(s35,2) + 23*s23*pow(s45,2) - 
             3*s14*(8*s35*s45 + s23*(2*s35 + s45) - 2*pow(s35,2) + 8*pow(s45,2)) + 12*pow(s45,3)) + 
          3*s12*(9*(s23 + 4*s45)*pow(s14,3) - 
             3*pow(s14,2)*(16*s45*(-s35 + s45) + s23*(-5*s35 + 23*s45) + 4*pow(s23,2)) - 
             s45*pow(s23,3) + 6*pow(s45,2)*(3*s35*s45 - 2*pow(s35,2) + pow(s45,2)) - 
             pow(s23,2)*(3*s35*s45 + 6*pow(s35,2) + 19*pow(s45,2)) + 
             s14*((-12*s35 + 31*s45)*pow(s23,2) + 3*pow(s23,3) + 
                6*s45*(-11*s35*s45 + 2*pow(s35,2) + pow(s45,2)) + 
                s23*(-43*s35*s45 + 15*pow(s35,2) + 63*pow(s45,2))) + 
             s23*(16*s45*pow(s35,2) + 9*pow(s35,3) + 53*s35*pow(s45,2) + pow(s45,3)))) - 
       mlsq*(6*s23*s45*pow(s12,4) + 3*pow(s12,3)*
           (-3*s14*s23*(s35 - 6*s45) + 6*s14*s45*(-2*s35 + s45) - 6*(2*s23 - s35 + 2*s45)*pow(s14,2) + 
             6*pow(s14,3) + 6*s14*pow(s23,2) - (3*s35 + 8*s45)*pow(s23,2) + 
             s23*(-2*s35*s45 + 3*pow(s35,2) - 14*pow(s45,2)) + 6*s35*pow(s45,2)) + 
          pow(s12,2)*(9*(3*s23 - 4*s35)*pow(s14,3) - 18*pow(s14,4) + 21*s45*pow(s23,3) - 
             18*s35*(s35 + 2*s45)*pow(s45,2) + 
             6*pow(s14,2)*(6*s23*s35 + 8*s23*s45 + 6*s35*s45 - 3*pow(s35,2) + 9*pow(s45,2)) + 
             pow(s23,2)*(42*s35*s45 + 18*pow(s35,2) + 83*pow(s45,2)) - 
             3*s14*(23*s45*pow(s23,2) + 3*pow(s23,3) - 12*s45*(s35*s45 + pow(s35,2) - pow(s45,2)) + 
                s23*(8*s35*s45 + 3*pow(s35,2) + 39*pow(s45,2))) + 
             s23*(-6*s45*pow(s35,2) - 18*pow(s35,3) + 33*s35*pow(s45,2) + 78*pow(s45,3))) + 
          s45*(-9*(s23 + 2*s45)*pow(s14,4) + 
             3*pow(s14,3)*(-6*s23*(s35 - 4*s45) + 12*s45*(-s35 + s45) + 11*pow(s23,2)) + 
             (30*s35 - 8*s45)*pow(s23,4) - 6*pow(s23,5) + 
             pow(s23,3)*(56*s35*s45 - 30*pow(s35,2) - 6*pow(s45,2)) - 
             pow(s14,2)*((-69*s35 + 97*s45)*pow(s23,2) + 45*pow(s23,3) + 
                18*s45*(-4*s35*s45 + pow(s35,2) + pow(s45,2)) + 
                3*s23*(-39*s35*s45 + 6*pow(s35,2) + 29*pow(s45,2))) + 
             s14*((-81*s35 + 51*s45)*pow(s23,3) + 27*pow(s23,4) + 36*s35*(s35 - s45)*pow(s45,2) + 
                pow(s23,2)*(-135*s35*s45 + 51*pow(s35,2) + 58*pow(s45,2)) - 
                6*s23*(-5*s45*pow(s35,2) + 3*pow(s35,3) + 24*s35*pow(s45,2) - 3*pow(s45,3))) - 
             18*pow(s35,2)*pow(s45,3) + pow(s23,2)*
              (16*s45*pow(s35,2) + 33*pow(s35,3) + 106*s35*pow(s45,2) + 12*pow(s45,3)) + 
             s23*(3*s45*pow(s35,3) - 9*pow(s35,4) - 9*pow(s35,2)*pow(s45,2) + 60*s35*pow(s45,3) + 
                12*pow(s45,4))) + s12*(9*(s23 + 4*s45)*pow(s14,4) - 
             9*pow(s14,3)*(-2*s23*s35 + 11*s23*s45 - 8*s35*s45 + 2*pow(s23,2) + 6*pow(s45,2)) + 
             3*pow(s14,2)*(6*s35*(2*s35 - 7*s45)*s45 + (-9*s35 + 22*s45)*pow(s23,2) + 3*pow(s23,3) + 
                s23*(-51*s35*s45 + 6*pow(s35,2) + 25*pow(s45,2))) + 
             s14*(3*(3*s35 + 4*s45)*pow(s23,3) + 
                18*pow(s45,2)*(2*s35*s45 - 4*pow(s35,2) + pow(s45,2)) + 
                pow(s23,2)*(69*s35*s45 - 18*pow(s35,2) + 10*pow(s45,2)) + 
                3*s23*(-7*s45*pow(s35,2) + 6*pow(s35,3) + 59*s35*pow(s45,2) + 15*pow(s45,3))) - 
             3*(5*s45*pow(s23,4) + 11*pow(s23,3)*pow(s45,2) - 6*s35*(2*s35 + s45)*pow(s45,3) + 
                pow(s23,2)*(17*s45*pow(s35,2) + 3*pow(s35,3) + 37*s35*pow(s45,2) + 24*pow(s45,3)) + 
                s23*(-5*s45*pow(s35,3) - 3*pow(s35,4) - 2*pow(s35,2)*pow(s45,2) + 29*s35*pow(s45,3) + 
                   18*pow(s45,4))))) - (9*pow(mlsq,4)*
           (s23*s35 - 2*s23*s45 + s14*(s23 + 4*s45) + s12*(-4*s14 + s23 + 4*s45) - pow(s23,2) - 
             4*pow(s45,2)) + 3*pow(mlsq,3)*
           (-5*s23*s35*s45 + (12*s14 - s23 - 12*s45)*pow(s12,2) + 9*(s23 + 4*s45)*pow(s14,2) + 
             s14*(3*s23*(4*s35 - 17*s45) + 24*(s35 - 2*s45)*s45 - 17*pow(s23,2)) - 11*s35*pow(s23,2) + 
             19*s45*pow(s23,2) - s12*(2*s23*s35 + 27*s23*s45 - 24*s35*s45 - 
                6*s14*(5*s23 - 4*s35 + 6*s45) + 36*pow(s14,2) + pow(s23,2)) + 8*pow(s23,3) + 
             9*s23*pow(s35,2) + 40*s23*pow(s45,2) - 24*s35*pow(s45,2) + 12*pow(s45,3)) + 
          pow(mlsq,2)*(-9*s23*pow(s12,3) + 
             3*pow(s12,2)*(12*s45*(-s35 + s45) - 12*s14*(2*s23 - s35 + 3*s45) + s23*(7*s35 + 32*s45) + 
                24*pow(s14,2)) + 27*(s23 + 4*s45)*pow(s14,3) + 58*s35*s45*pow(s23,2) + 
             48*s35*pow(s23,3) - 45*s45*pow(s23,3) - 21*pow(s23,4) + 15*s23*s45*pow(s35,2) - 
             48*pow(s23,2)*pow(s35,2) + 27*s23*pow(s35,3) + 186*s23*s35*pow(s45,2) - 
             100*pow(s23,2)*pow(s45,2) - 36*pow(s35,2)*pow(s45,2) - 
             3*pow(s14,2)*(-15*s23*s35 + 81*s23*s45 - 48*s35*s45 + 25*pow(s23,2) + 60*pow(s45,2)) - 
             3*s12*(-12*(4*s23 - 4*s35 + 3*s45)*pow(s14,2) + 36*pow(s14,3) - 
                (12*s35 + 5*s45)*pow(s23,2) + 8*pow(s23,3) + 
                s14*(s23*(-31*s35 + 33*s45) + 4*pow(s23,2) + 
                   12*(-5*s35*s45 + pow(s35,2) - pow(s45,2))) + 
                12*s45*(s35*s45 - pow(s35,2) + pow(s45,2)) + 
                s23*(51*s35*s45 + 13*pow(s35,2) + 24*pow(s45,2))) + 
             s14*((-93*s35 + 182*s45)*pow(s23,2) + 69*pow(s23,3) + 
                36*s45*(-6*s35*s45 + pow(s35,2) + 2*pow(s45,2)) + 
                3*s23*(-62*s35*s45 + 15*pow(s35,2) + 90*pow(s45,2))) - 42*s23*pow(s45,3) + 
             72*s35*pow(s45,3)) - s23*(-3*(s14 - s23)*pow(s12,4) - 9*s45*pow(s14,4) + 
             pow(s12,3)*(-15*s35*s45 + 2*s23*(-9*s35 + s45) - 3*s14*(5*s23 - 5*s35 + s45) + 
                9*pow(s14,2) + 6*pow(s23,2)) + 
             pow(s12,2)*(3*s35*s45*(13*s35 + 14*s45) + 3*(3*s23 - 6*s35 - 7*s45)*pow(s14,2) - 
                3*pow(s14,3) - 2*(9*s35 + 7*s45)*pow(s23,2) + 3*pow(s23,3) + 
                2*s23*(10*s35*s45 + 18*pow(s35,2) - 9*pow(s45,2)) - 
                3*s14*(-13*s23*s35 - 12*s23*s45 + 2*s35*s45 + 3*pow(s23,2) + 7*pow(s35,2) - 
                   8*pow(s45,2))) - (2*s23 - 3*s35)*(s23 + s45)*(3*s35 + 4*s45)*
              (2*s23*s45 + 2*s35*s45 + pow(s23,2) + pow(s35,2) + 2*pow(s45,2)) + 
             pow(s14,3)*(9*s23*(s35 + 4*s45) + 30*pow(s45,2)) - 
             pow(s14,2)*((24*s35 + 53*s45)*pow(s23,2) + 6*(s35 + 7*s45)*pow(s45,2) + 
                s23*(12*s35*s45 - 9*pow(s35,2) + 86*pow(s45,2))) + 
             s12*(-3*(11*s23 - 3*s35 + 9*s45)*pow(s14,3) + 9*pow(s14,4) - 
                4*(3*s35 - 5*s45)*pow(s23,3) + 6*pow(s23,4) - 
                3*s35*s45*(24*s35*s45 + 11*pow(s35,2) + 18*pow(s45,2)) + 
                pow(s23,2)*(s35*s45 + 18*pow(s35,2) + 44*pow(s45,2)) + 
                pow(s14,2)*(-30*s23*s35 + 74*s23*s45 + 3*s35*s45 + 45*pow(s23,2) + 9*pow(s35,2) + 
                   54*pow(s45,2)) - 2*s23*
                 (23*s45*pow(s35,2) + 15*pow(s35,3) + 17*s35*pow(s45,2) - 16*pow(s45,3)) - 
                s14*((-33*s35 + 67*s45)*pow(s23,2) + 27*pow(s23,3) - 9*s45*pow(s35,2) - 9*pow(s35,3) - 
                   6*s35*pow(s45,2) + s23*(12*s35*s45 + 33*pow(s35,2) + 100*pow(s45,2)) + 42*pow(s45,3)
                   )) + s14*((21*s35 + 34*s45)*pow(s23,3) - 6*(pow(s35,2) - 4*pow(s45,2))*pow(s45,2) + 
                pow(s23,2)*(18*s35*s45 - 18*pow(s35,2) + 80*pow(s45,2)) + 
                s23*(9*pow(s35,3) + 18*s35*pow(s45,2) + 76*pow(s45,3)))) + 
          mlsq*(3*s23*(-6*s14 + 5*s23 - 4*s35 + s45)*pow(s12,3) + 3*s23*pow(s12,4) + 
             9*(s23 + 4*s45)*pow(s14,4) - 
             3*pow(s14,3)*(24*s45*(-s35 + s45) + s23*(-6*s35 + 33*s45) + 11*pow(s23,2)) - 
             52*s35*s45*pow(s23,3) - 30*s35*pow(s23,4) - 2*s45*pow(s23,4) + 6*pow(s23,5) - 
             23*s45*pow(s23,2)*pow(s35,2) + 30*pow(s23,3)*pow(s35,2) + 3*s23*s45*pow(s35,3) - 
             33*pow(s23,2)*pow(s35,3) + 9*s23*pow(s35,4) - 140*s35*pow(s23,2)*pow(s45,2) - 
             12*pow(s23,3)*pow(s45,2) + 36*s23*pow(s35,2)*pow(s45,2) + 
             pow(s12,2)*(-6*(11*s23 - 6*s35 + 12*s45)*pow(s14,2) + 36*pow(s14,3) - 
                (51*s35 + 76*s45)*pow(s23,2) + 
                3*s14*(12*s45*(-2*s35 + s45) + s23*(2*s35 + 51*s45) + 8*pow(s23,2)) + 6*pow(s23,3) + 
                3*s23*(11*s35*s45 + 8*pow(s35,2) - 20*pow(s45,2)) + 36*s35*pow(s45,2)) + 
             pow(s14,2)*((-69*s35 + 89*s45)*pow(s23,2) + 45*pow(s23,3) + 
                36*s45*(-4*s35*s45 + pow(s35,2) + pow(s45,2)) + 
                3*s23*(-57*s35*s45 + 6*pow(s35,2) + 40*pow(s45,2))) + 
             s14*(3*(27*s35 - 8*s45)*pow(s23,3) - 27*pow(s23,4) + 
                pow(s23,2)*(150*s35*s45 - 51*pow(s35,2) - 38*pow(s45,2)) + 
                72*s35*(-s35 + s45)*pow(s45,2) + 3*s23*s35*(-11*s35*s45 + 6*pow(s35,2) + 84*pow(s45,2))
                ) - 84*s23*s35*pow(s45,3) - 36*pow(s23,2)*pow(s45,3) + 36*pow(s35,2)*pow(s45,3) + 
             s12*(18*(3*s23 - 4*s35 + 2*s45)*pow(s14,3) - 36*pow(s14,4) + 
                3*pow(s14,2)*(3*s23*(10*s35 + s45) + 8*pow(s23,2) + 
                   12*(3*s35*s45 - pow(s35,2) + pow(s45,2))) - 
                s14*((-21*s35 + 95*s45)*pow(s23,2) + 66*pow(s23,3) + 
                   36*s45*(-2*pow(s35,2) + pow(s45,2)) + 
                   6*s23*(26*s35*s45 + pow(s35,2) + 27*pow(s45,2))) + 
                3*((-12*s35 + 17*s45)*pow(s23,3) + 8*pow(s23,4) - 12*s35*(s35 + s45)*pow(s45,2) + 
                   pow(s23,2)*(31*s35*s45 + 23*pow(s35,2) + 42*pow(s45,2)) + 
                   s23*(-13*s45*pow(s35,2) - 8*pow(s35,3) + 8*s35*pow(s45,2) + 26*pow(s45,3)))) - 
             24*s23*pow(s45,4)))*pow(mw2c*conj(mw2c),1/2.))*pow(conj(pow(SW2,1/2.)),-2))/3.;

	// Averaging over photon polarisations
	prefactor /= 2.;
	// Include the prefactor and averaging/sum factors
	ME2 *= prefactor;
	return ME2.real();
	// Validated against Recola up to a factor of two
}


double ME2_Analytic::nu1gamma_l1nu2l2x(int i1, int i2, int i3, int i4, int i5, KinematicData &Kin, int f1, int f2){

	if( !is_neutrino(f1) ){
		cerr << "nu1gamma_nu1f2f2x: fermion line 1 is not a neutrino\n";
		abort();
	}

	// Using another basis of kinematics
	double s12 = 2.0*Kin.pij(i1,i2);
	double s14 = 2.0*Kin.pij(i1,i4);	
	double s23 = 2.0*Kin.pij(i2,i3);
	double s25 = 2.0*Kin.pij(i2,i5);	
	double s35 = 2.0*Kin.pij(i3,i5);		
	double s45 = 2.0*Kin.pij(i4,i5);
	// fermion masses
	double ml1sq = Kin.p(3).m2();
	double ml2sq = Kin.p(5).m2();	
	// f3, f4 a generic fermion line
	// complex<double> gLf, gRf; double Qf;
	// assign_EW_charges(Qf,gLf,gRf,abs(f2));

	complex<double> chiw45 = 1. / (ml2sq + s45 - MW2C);
	complex<double> mw2c = MW2C;
	// The following option to drop width effects in numerator (a sort of pole approximation)
	// mw2c = pow(mw,2);
	complex<double> chiw13 = 1. / (mw2c + s12 - s23 - s45 );
	//
	complex<double> prefactor = 4.*chiw13*chiw45*conj(ALPHA)*conj(chiw13)*conj(chiw45)*pow(ALPHA,2)*pow(pi,3);

	complex<double> ME2 = 4.*pow(s23,-2)*pow(s25,-2)*pow(SW2,-2)*(-(ml2sq*(5*s23 + 2*s45)*pow(ml1sq,5)) + 
     pow(ml1sq,4)*(2*(12*s14 - 3*s23 - 5*s45)*pow(ml2sq,2) + 8*(s14 - s45)*pow(s12,2) + 
        ml2sq*(-2*s14*s23 - 10*s23*s35 + 28*s14*s45 + 9*s23*s45 - 4*s35*s45 + 
           4*s12*(-8*s14 + s23 + 4*s45) + 6*pow(s23,2) - 14*pow(s45,2)) - 
        4*s12*(s23*s35 - s23*s45 + s14*(s23 + 4*s45) - 4*pow(s45,2)) - 
        4*s45*(-(s23*s35) + s23*s45 - s14*(s23 + 2*s45) + pow(s23,2) + 2*pow(s45,2))) + 
     4.*mw2c*conj(mw2c)*((2*ml2sq*(s14 - s23 - s45) + 
           (2*s14 - s23 - 2*s45)*(3*s14 - 2*s23 + 2*s35 - s45))*pow(ml1sq,3) + 
        (2*s14 - s23 - 2*s45)*pow(ml1sq,4) - 
        pow(ml1sq,2)*(s12*s23*s35 + 2*s12*s23*s45 - 8*s23*s35*s45 + 2*(s14 - s45)*pow(ml2sq,2) - 
           s23*pow(s12,2) + 2*(5*s23 - 4*s35 + 5*s45)*pow(s14,2) - 6*pow(s14,3) + s12*pow(s23,2) - 
           3*s35*pow(s23,2) + 3*s45*pow(s23,2) + pow(s23,3) + s23*pow(s35,2) + 2*s45*pow(s35,2) + 
           ml2sq*(2*s12*s23 + 5*s14*s23 + 2*s23*s35 + 8*s14*s45 - 5*s23*s45 - 4*pow(s14,2) - 
              2*pow(s23,2) - 4*pow(s45,2)) + 2*s23*pow(s45,2) - 4*s35*pow(s45,2) - 
           s14*(-9*s23*s35 + 12*s23*s45 - 12*s35*s45 + 5*pow(s23,2) + 2*pow(s35,2) + 4*pow(s45,2))) + 
        ml1sq*(2*(-s14 + s23 + s45)*pow(ml2sq,3) + (-3*s23 + 4*s35 - 4*s45)*pow(s14,3) + 
           2*pow(s14,4) + s14*(-(s23*s35*(s35 - 12*s45)) + 4*s35*s45*(-s35 + s45) - 
              s12*s23*(s23 + 3*s35 + 4*s45) + 2*s23*pow(s12,2) + (3*s35 - s45)*pow(s23,2) + pow(s23,3))
             + pow(ml2sq,2)*(-4*s12*s23 + 3*s14*s23 - 4*s14*s35 + 2*s23*s35 + 3*s23*s45 + 4*s35*s45 - 
              2*pow(s14,2) + 2*pow(s23,2) + 2*pow(s45,2)) + 
           pow(s14,2)*(s23*(-6*s35 + 5*s45) + 2*(-4*s35*s45 + pow(s35,2) + pow(s45,2))) - 
           (s23 - s35)*(-(s12*s23*(s23 + s35 + 2*s45)) + s23*pow(s12,2) + 
              s35*(3*s23*s45 + pow(s23,2) + 2*pow(s45,2))) + 
           2*ml2sq*(-2*s12*s23*(s35 + s45) + s23*pow(s12,2) - (s23 + 3*s45)*pow(s14,2) + pow(s14,3) + 
              s35*(2*s23*s45 + s45*(s35 + 2*s45) + pow(s23,2)) + 
              s14*(-2*s12*s23 + s23*s35 + 4*s23*s45 - 2*s35*s45 - pow(s35,2) + 2*pow(s45,2)))) + 
        s23*((-2*s12 + s14 + 2*s23 + 2*s35 + 3*s45)*pow(ml2sq,3) + pow(ml2sq,4) + 
           pow(ml2sq,2)*(3*s23*s35 + 3*s23*s45 + 4*s35*s45 - s12*(4*s14 + 3*s23 + 3*s35 + 2*s45) + 
              s14*(-s23 + 3*s35 + 4*s45) + pow(s12,2) + pow(s23,2) + pow(s35,2) + 2*pow(s45,2)) + 
           s14*(s14 - s23 + s35)*(2*s23*s45 + 2*s35*s45 - 2*s14*(s23 + s45) - 2*s12*(s35 + s45) + 
              pow(s12,2) + pow(s14,2) + pow(s23,2) + pow(s35,2) + 2*pow(s45,2)) + 
           ml2sq*(s35*(s23 + s45)*(s23 + s35 + 2*s45) + (2*s14 + s23 + s35)*pow(s12,2) - 
              (2*s23 - 2*s35 + s45)*pow(s14,2) - 
              s12*((s23 + s35)*(s23 + s35 + 2*s45) + s14*(-3*s23 + 5*s35 + 4*s45) + 2*pow(s14,2)) + 
              pow(s14,3) + s14*(-3*s23*s35 + s23*s45 + 4*s35*s45 + pow(s23,2) + 3*pow(s35,2) + 
                 4*pow(s45,2))))) + pow(ml1sq,3)*
      ((40*s14 + 8*s23 - 6*s45)*pow(ml2sq,3) + 
        pow(ml2sq,2)*(10*s23*s35 - 2*s12*(24*s14 + 5*s23 - 8*s45) + 45*s23*s45 - 16*s35*s45 + 
           s14*(-50*s23 + 48*s35 + 16*s45) + 56*pow(s14,2) + 10*pow(s23,2) - 4*pow(s45,2)) + 
        4*(pow(s12,2)*(s14*(-6*s23 + 4*s35 - 8*s45) + 2*s45*(-2*s35 + s45) + s23*(s35 + 6*s45) + 
              6*pow(s14,2)) + s45*(3*(s23 + 2*s45)*pow(s14,2) + (-4*s35 + 5*s45)*pow(s23,2) + 
              3*pow(s23,3) + 2*(-2*s35 + s45)*pow(s45,2) - 
              2*s14*(-2*s23*s35 + 5*s23*s45 - 2*s35*s45 + 3*pow(s23,2) + 4*pow(s45,2)) + 
              s23*(-(s35*s45) + 3*pow(s35,2) + 7*pow(s45,2))) - 
           s12*(3*(s23 + 4*s45)*pow(s14,2) + 
              s14*(4*s23*(s35 - 4*s45) + 8*(s35 - 2*s45)*s45 - 3*pow(s23,2)) + 
              (-2*s35 + 3*s45)*pow(s23,2) + 4*(-2*s35 + s45)*pow(s45,2) + 
              s23*(3*pow(s35,2) + 13*pow(s45,2)))) + 
        ml2sq*(40*s23*s35*s45 + (8*s14 + 5*s23 - 8*s45)*pow(s12,2) + 3*(5*s23 + 26*s45)*pow(s14,2) + 
           7*s35*pow(s23,2) + 5*s45*pow(s23,2) - 
           s12*(9*s23*s35 + s14*(-92*s23 + 64*s35 - 48*s45) + 62*s23*s45 - 32*s35*s45 + 
              80*pow(s14,2) + 9*pow(s23,2)) + 7*pow(s23,3) - s23*pow(s35,2) - 2*s45*pow(s35,2) + 
           s14*(6*s23*s35 - 74*s23*s45 + 60*s35*s45 - 22*pow(s23,2) - 52*pow(s45,2)) + 
           46*s23*pow(s45,2) - 28*s35*pow(s45,2) + 8*pow(s45,3))) + 
     ml1sq*((-8*s14 - 3*s23 + 8*s45)*pow(ml2sq,5) + 
        pow(ml2sq,4)*(-10*s23*s35 + 2*s12*(8*s14 + 3*s23 - 8*s45) - 21*s23*s45 + 16*s35*s45 + 
           2*s14*(s23 - 8*(s35 + s45)) - 8*pow(s14,2) - 6*pow(s23,2) + 24*pow(s45,2)) + 
        pow(ml2sq,3)*(-52*s23*s35*s45 + (-8*s14 - 3*s23 + 8*s45)*pow(s12,2) - 
           (23*s23 + 40*s45)*pow(s14,2) + 8*pow(s14,3) - 33*s35*pow(s23,2) - 29*s45*pow(s23,2) + 
           s12*(-16*s14*(s23 - 2*s35 - s45) - 32*s45*(s35 + s45) + s23*(29*s35 + 42*s45) + 
              16*pow(s14,2) + 5*pow(s23,2)) + pow(s23,3) + s23*pow(s35,2) + 8*s45*pow(s35,2) - 
           46*s23*pow(s45,2) + 48*s35*pow(s45,2) + 
           2*s14*(-7*s23*s35 + 15*s23*s45 - 24*s35*s45 + 14*pow(s23,2) - 4*pow(s35,2) + 
              4*pow(s45,2)) + 24*pow(s45,3)) + 
        4*(s23*s35*s45*pow(s12,3) + pow(s12,2)*
            ((-6*s23 + 4*s35 - 4*s45)*pow(s14,3) + 2*pow(s14,4) + 
              s14*(s23*(11*s35 - 4*s45)*s45 + 4*s35*s45*(-s35 + s45) + 2*(s35 - 4*s45)*pow(s23,2) - 
                 2*pow(s23,3)) + (s35 + 2*s45)*pow(s23,3) + 
              s23*s35*(-(s35*s45) + pow(s35,2) - 7*pow(s45,2)) - 
              2*pow(s23,2)*(2*s35*s45 + pow(s35,2) - pow(s45,2)) + 2*pow(s35,2)*pow(s45,2) + 
              pow(s14,2)*(s23*(-7*s35 + 10*s45) + 6*pow(s23,2) + 
                 2*(-4*s35*s45 + pow(s35,2) + pow(s45,2)))) + 
           s45*((s23 + 2*s45)*pow(s14,4) + 
              pow(s14,3)*(2*s23*(s35 - 4*s45) + 4*(s35 - s45)*s45 - 4*pow(s23,2)) + 
              (-4*s35 + 3*s45)*pow(s23,4) + pow(s23,5) + 4*pow(s23,3)*pow(s35 - s45,2) + 
              pow(s14,2)*((-8*s35 + 13*s45)*pow(s23,2) + 6*pow(s23,3) + 
                 2*s45*(-4*s35*s45 + pow(s35,2) + pow(s45,2)) + 
                 s23*(-13*s35*s45 + 2*pow(s35,2) + 11*pow(s45,2))) - 
              s14*(-10*(s35 - s45)*pow(s23,3) + 4*pow(s23,4) + 4*s35*(s35 - s45)*pow(s45,2) + 
                 pow(s23,2)*(-17*s35*s45 + 6*pow(s35,2) + 11*pow(s45,2)) - 
                 2*s23*(-2*s45*pow(s35,2) + pow(s35,3) + 8*s35*pow(s45,2) - 2*pow(s45,3))) + 
              s23*s35*(-(s45*pow(s35,2)) + pow(s35,3) + s35*pow(s45,2) - 6*pow(s45,3)) + 
              2*pow(s35,2)*pow(s45,3) + pow(s23,2)*(-4*pow(s35,3) - 11*s35*pow(s45,2) + 2*pow(s45,3)))\
            - s12*(-((3*s23 - 2*s35 + 2*s45)*(s23 + 4*s45)*pow(s14,3)) + (s23 + 4*s45)*pow(s14,4) + 
              s45*pow(s23,4) + pow(s23,3)*(-(s35*s45) + pow(s35,2) + 5*pow(s45,2)) + 
              pow(s14,2)*((-4*s35 + 17*s45)*pow(s23,2) + 3*pow(s23,3) + 
                 4*s45*(-4*s35*s45 + pow(s35,2) + pow(s45,2)) + 
                 s23*(-20*s35*s45 + 2*pow(s35,2) + 21*pow(s45,2))) - 
              pow(s23,2)*(4*s45*pow(s35,2) + 2*pow(s35,3) + 13*s35*pow(s45,2) - 4*pow(s45,3)) + 
              4*pow(s35,2)*pow(s45,3) + s23*(pow(s35,4) - 12*s35*pow(s45,3)) - 
              s14*(-2*(s35 - 4*s45)*pow(s23,3) + pow(s23,4) + 8*s35*(s35 - s45)*pow(s45,2) + 
                 pow(s23,2)*(-13*s35*s45 + 3*pow(s35,2) + 18*pow(s45,2)) + 
                 s23*(4*s45*pow(s35,2) - 2*pow(s35,3) - 27*s35*pow(s45,2) + 8*pow(s45,3))))) + 
        ml2sq*(4*s23*s35*pow(s12,3) + 4*(s23 + 4*s45)*pow(s14,4) - 
           8*pow(s14,3)*(-(s23*s35) + 7*s23*s45 - 4*s35*s45 + 2*pow(s23,2) + 3*pow(s45,2)) + 
           pow(s12,2)*(8*s35*s45*(s35 + 2*s45) - (11*s23 + 24*s45)*pow(s14,2) + 8*pow(s14,3) + 
              4*(s35 - 3*s45)*pow(s23,2) - 8*pow(s23,3) + 
              s14*(s23*(5*s35 + 24*s45) + 11*pow(s23,2) - 8*(2*s35*s45 + pow(s35,2) - 2*pow(s45,2))) - 
              8*s23*(5*s35*s45 + pow(s35,2) + 2*pow(s45,2))) + 
           pow(s14,2)*((-36*s35 + 66*s45)*pow(s23,2) + 25*pow(s23,3) - 
              8*s45*(8*s35*s45 - 2*pow(s35,2) + pow(s45,2)) + 
              2*s23*(-52*s35*s45 + 4*pow(s35,2) + 29*pow(s45,2))) + 
           s14*(15*(3*s35 - 2*s45)*pow(s23,3) - 17*pow(s23,4) + 
              8*pow(s45,2)*(2*s35*s45 - 5*pow(s35,2) + 2*pow(s45,2)) - 
              6*pow(s23,2)*(-15*s35*s45 + 4*pow(s35,2) + 7*pow(s45,2)) + 
              2*s23*s35*(-12*s35*s45 + 4*pow(s35,2) + 53*pow(s45,2))) - 
           2*s12*(-2*(13*s23 - 8*s35 + 4*s45)*pow(s14,3) + 8*pow(s14,4) + 
              pow(s14,2)*(s23*(-44*s35 + 27*s45) + 25*pow(s23,2) + 
                 8*(-4*s35*s45 + pow(s35,2) - 2*pow(s45,2))) - 
              s14*((-25*s35 + 19*s45)*pow(s23,2) + pow(s23,3) + 
                 s23*(-57*s35*s45 + 10*pow(s35,2) - 8*pow(s45,2)) + 8*s45*(3*pow(s35,2) - 2*pow(s45,2))
                 ) - 2*(-3*s35*pow(s23,3) + 3*pow(s23,4) - 8*s35*(s35 + s45)*pow(s45,2) + 
                 pow(s23,2)*(19*s35*s45 + 6*pow(s35,2) + 2*pow(s45,2)) + 
                 s23*(5*s45*pow(s35,2) - 2*pow(s35,3) + 26*s35*pow(s45,2) + 8*pow(s45,3)))) + 
           4*((-4*s35 + s45)*pow(s23,4) + pow(s23,5) + 
              pow(s23,3)*(-5*s35*s45 + 3*pow(s35,2) + 2*pow(s45,2)) + 
              2*s35*(3*s35 + 2*s45)*pow(s45,3) - 
              pow(s23,2)*(12*s45*pow(s35,2) + 3*pow(s35,3) + 22*s35*pow(s45,2) + 2*pow(s45,3)) + 
              s23*(2*s45*pow(s35,3) + pow(s35,4) - 4*pow(s35,2)*pow(s45,2) - 16*s35*pow(s45,3) - 
                 4*pow(s45,4)))) + pow(ml2sq,2)*
         (-8*(3*s23 - 2*s35)*pow(s14,3) + 8*pow(s14,4) - 89*s35*s45*pow(s23,2) + 
           pow(s12,2)*(2*s14*(5*s23 - 8*s35) + 8*s45*(2*s35 + s45) - s23*(23*s35 + 16*s45) - 
              8*pow(s14,2) + pow(s23,2)) + 5*s35*pow(s23,3) - 4*s45*pow(s23,3) - 8*pow(s23,4) - 
           19*s23*s45*pow(s35,2) - 35*pow(s23,2)*pow(s35,2) + 12*s23*pow(s35,3) + 
           pow(s14,2)*(-(s23*(52*s35 + 9*s45)) + 10*pow(s23,2) + 
              8*(-4*s35*s45 + pow(s35,2) - 6*pow(s45,2))) - 82*s23*s35*pow(s45,2) - 
           36*pow(s23,2)*pow(s45,2) + 24*pow(s35,2)*pow(s45,2) + 
           s14*((23*s35 + 31*s45)*pow(s23,2) + 17*pow(s23,3) - 
              32*s45*(s35*s45 + pow(s35,2) - pow(s45,2)) + 
              s23*(28*s35*s45 - 8*pow(s35,2) + 44*pow(s45,2))) + 
           s12*((34*s23 + 64*s45)*pow(s14,2) - 16*pow(s14,3) + 2*(9*s35 + 17*s45)*pow(s23,2) + 
              7*pow(s23,3) + s14*(s23*(9*s35 - 68*s45) - 35*pow(s23,2) + 
                 16*(4*s35*s45 + pow(s35,2) - 2*pow(s45,2))) - 
              16*s45*(4*s35*s45 + pow(s35,2) + pow(s45,2)) + 
              s23*(86*s35*s45 + 15*pow(s35,2) + 68*pow(s45,2))) - 44*s23*pow(s45,3) + 
           48*s35*pow(s45,3) + 8*pow(s45,4))) + 
     pow(ml1sq,2)*(2*(4*s14 + 5*(s23 + s45))*pow(ml2sq,4) + 
        pow(ml2sq,3)*(30*s23*s35 + s14*(-38*s23 + 32*s35 - 28*s45) + 27*s23*s45 + 4*s35*s45 - 
           4*s12*(5*s23 + 4*s45) + 48*pow(s14,2) - 6*pow(s23,2) + 34*pow(s45,2)) + 
        4*(pow(s12,2)*(-2*s35*(s35 - 2*s45)*s45 - 2*(6*s23 - 4*s35 + 5*s45)*pow(s14,2) + 
              6*pow(s14,3) - 2*(s35 + 3*s45)*pow(s23,2) + 
              s23*(7*s35*s45 + 2*pow(s35,2) - 4*pow(s45,2)) + 
              2*s14*(-3*s23*s35 + 8*s23*s45 - 6*s35*s45 + 3*pow(s23,2) + pow(s35,2) + 2*pow(s45,2))) - 
           s12*(3*(s23 + 4*s45)*pow(s14,3) + 
              pow(s14,2)*(s23*(5*s35 - 29*s45) + 4*(4*s35 - 5*s45)*s45 - 6*pow(s23,2)) + 
              (s35 - 3*s45)*pow(s23,3) - 4*s35*(s35 - 2*s45)*pow(s45,2) - 
              2*pow(s23,2)*(2*pow(s35,2) + 7*pow(s45,2)) + 
              s14*((-6*s35 + 20*s45)*pow(s23,2) + 3*pow(s23,3) + 
                 4*s45*(-6*s35*s45 + pow(s35,2) + 2*pow(s45,2)) + 
                 s23*(-19*s35*s45 + 5*pow(s35,2) + 34*pow(s45,2))) + 
              s23*(2*s45*pow(s35,2) + 3*pow(s35,3) + 17*s35*pow(s45,2) - 8*pow(s45,3))) + 
           s45*(3*(s23 + 2*s45)*pow(s14,3) + (6*s35 - 7*s45)*pow(s23,3) - 3*pow(s23,4) + 
              pow(s14,2)*(5*s23*s35 - 17*s23*s45 + 8*s35*s45 - 9*pow(s23,2) - 10*pow(s45,2)) + 
              pow(s23,2)*(6*s35*s45 - 6*pow(s35,2) - 9*pow(s45,2)) - 2*s35*(s35 - 2*s45)*pow(s45,2) + 
              s14*((-11*s35 + 18*s45)*pow(s23,2) + 9*pow(s23,3) + 
                 2*s45*(-6*s35*s45 + pow(s35,2) + 2*pow(s45,2)) + 
                 s23*(-13*s35*s45 + 5*pow(s35,2) + 18*pow(s45,2))) + 
              s23*(3*pow(s35,3) + 10*s35*pow(s45,2) - 4*pow(s45,3)))) - 
        2*pow(ml2sq,2)*(-36*s23*s35*s45 + (4*s14 - 7*s23 - 4*s45)*pow(s12,2) + 
           (34*s23 - 32*s35 - 13*s45)*pow(s14,2) - 20*pow(s14,3) + 5*s35*pow(s23,2) + 
           20*s45*pow(s23,2) - 16*s23*pow(s35,2) + 3*s45*pow(s35,2) - 8*s23*pow(s45,2) - 
           10*s35*pow(s45,2) + 2*s12*(7*s23*s35 + 12*s23*s45 - 2*s14*(7*s23 - 4*s35 + 10*s45) + 
              16*pow(s14,2) + 12*pow(s45,2)) + 
           s14*(26*s23*s35 + 6*s35*s45 - 11*pow(s23,2) - 12*pow(s35,2) + 26*pow(s45,2)) - 20*pow(s45,3)
           ) + ml2sq*(-4*s23*pow(s12,3) + 16*(s23 + 4*s45)*pow(s14,3) - 23*s35*s45*pow(s23,2) + 
           15*s35*pow(s23,3) - 20*s45*pow(s23,3) - 12*pow(s23,4) + 47*s23*s45*pow(s35,2) - 
           7*pow(s23,2)*pow(s35,2) + 8*s23*pow(s35,3) + 
           pow(s14,2)*(24*s23*s35 - 139*s23*s45 + 96*s35*s45 - 44*pow(s23,2) - 62*pow(s45,2)) + 
           66*s23*s35*pow(s45,2) - 52*pow(s23,2)*pow(s45,2) - 14*pow(s35,2)*pow(s45,2) + 
           pow(s12,2)*(-6*s14*s23 + 9*s23*s35 - 32*s14*s45 + 20*s23*s45 + 16*pow(s14,2) + 
              3*pow(s23,2) + 16*pow(s45,2)) + 
           s14*(4*s35*(8*s35 - 23*s45)*s45 + (-43*s35 + 93*s45)*pow(s23,2) + 41*pow(s23,3) + 
              4*s23*(-23*s35*s45 + 4*pow(s35,2) + 26*pow(s45,2))) + 
           s12*(4*(35*s23 + 12*(-2*s35 + s45))*pow(s14,2) - 64*pow(s14,3) + 
              2*(6*s35 + 29*s45)*pow(s23,2) - 3*pow(s23,3) + 
              16*s45*(-2*s35*s45 + pow(s35,2) - 2*pow(s45,2)) + 
              s23*(-90*s35*s45 - 25*pow(s35,2) + 4*pow(s45,2)) + 
              s14*(s23*(107*s35 - 116*s45) - 79*pow(s23,2) + 32*(3*s35*s45 - pow(s35,2) + pow(s45,2))))
             - 12*s23*pow(s45,3) + 32*s35*pow(s45,3) + 16*pow(s45,4))) - 
     4*s23*((-3*s12 + 2*s14 - s23 + 5*s35 + 3*s45)*pow(ml2sq,5) + pow(ml2sq,6) + 
        s35*(s14 - s23 + s35)*s45*(-s12 + s23 + s45)*
         (2*s23*s45 + 2*s35*s45 - 2*s14*(s23 + s45) - 2*s12*(s35 + s45) + pow(s12,2) + pow(s14,2) + 
           pow(s23,2) + pow(s35,2) + 2*pow(s45,2)) + 
        pow(ml2sq,4)*(s23*s35 - 5*s23*s45 + 16*s35*s45 + s14*(s23 + 5*s35 + 6*s45) - 
           s12*(6*s14 - 4*s23 + 13*s35 + 6*s45) + 3*pow(s12,2) + pow(s14,2) - 2*pow(s23,2) + 
           8*pow(s35,2) + 3*pow(s45,2)) + 
        pow(ml2sq,3)*((6*s14 - 5*s23 + 11*s35 + 3*s45)*pow(s12,2) - pow(s12,3) + 
           (-s23 + s35 + 3*s45)*pow(s14,2) - 3*s35*pow(s23,2) - 7*s45*pow(s23,2) - pow(s23,3) + 
           5*s23*pow(s35,2) + 23*s45*pow(s35,2) + 5*pow(s35,3) - 9*s23*pow(s45,2) + 
           19*s35*pow(s45,2) + 2*s14*(7*s35*s45 + s23*(s35 + 2*s45) + pow(s23,2) + 2*pow(s35,2) + 
              3*pow(s45,2)) - s12*(28*s35*s45 - 4*s23*(s35 + 3*s45) + 2*s14*(7*s35 + 6*s45) + 
              3*pow(s14,2) - 3*pow(s23,2) + 17*pow(s35,2) + 3*pow(s45,2)) + pow(s45,3)) + 
        pow(ml2sq,2)*((-2*s14 + 2*s23 - 3*s35)*pow(s12,3) + (s23 + s35)*pow(s14,3) - 
           7*s35*s45*pow(s23,2) - 4*s45*pow(s23,3) - pow(s23,4) + 9*s23*s45*pow(s35,2) - 
           2*pow(s23,2)*pow(s35,2) + pow(s12,2)*
            (-3*s14*s23 - 8*s23*s35 - 7*s23*s45 + 13*s35*s45 + 6*s14*(2*s35 + s45) + 3*pow(s14,2) + 
              11*pow(s35,2)) + 4*s23*pow(s35,3) + 12*s45*pow(s35,3) + pow(s35,4) - 
           4*s23*s35*pow(s45,2) - 9*pow(s23,2)*pow(s45,2) + 24*pow(s35,2)*pow(s45,2) + 
           pow(s14,2)*(-2*s23*(s35 + 2*s45) - 3*pow(s23,2) + pow(s35,2) + 3*pow(s45,2)) + 
           s12*(s23*s45*(11*s35 + 12*s45) + (s23 - 2*(s35 + 3*s45))*pow(s14,2) + 
              (3*s35 + 7*s45)*pow(s23,2) + pow(s23,3) - 
              s14*(27*s35*s45 + s23*(s35 + s45) + 2*pow(s23,2) + 8*pow(s35,2) + 6*pow(s45,2)) - 
              s35*(33*s35*s45 + 8*pow(s35,2) + 19*pow(s45,2))) - 7*s23*pow(s45,3) + 
           10*s35*pow(s45,3) + s14*((s35 + 8*s45)*pow(s23,2) + 3*pow(s23,3) + 8*s45*pow(s35,2) + 
              pow(s35,3) + 15*s35*pow(s45,2) + s23*(7*s35*s45 + pow(s35,2) + 6*pow(s45,2)) + 
              2*pow(s45,3))) - ml2sq*(-((2*s35*s45 + s23*(s35 + s45))*pow(s14,3)) + 
           pow(s12,3)*(-3*s23*s35 + s14*(-2*s23 + 3*s35) + s35*(2*s35 + s45) + pow(s14,2) + 
              pow(s23,2)) + 3*s35*s45*pow(s23,3) + s35*pow(s23,4) + s45*pow(s23,4) + 
           s45*pow(s23,2)*pow(s35,2) - pow(s23,3)*pow(s35,2) - 5*s23*s45*pow(s35,3) + 
           pow(s23,2)*pow(s35,3) - s23*pow(s35,4) - 2*s45*pow(s35,4) + 8*s35*pow(s23,2)*pow(s45,2) + 
           3*pow(s23,3)*pow(s45,2) - 6*s23*pow(s35,2)*pow(s45,2) - 9*pow(s35,3)*pow(s45,2) - 
           pow(s12,2)*(-(s14*s23*(s35 + 3*s45)) + s14*s35*(4*s35 + 13*s45) + 
              (s35 + 3*s45)*pow(s14,2) + 3*s35*(4*s35*s45 + pow(s35,2) + pow(s45,2)) - 
              s23*(9*s35*s45 + 3*pow(s35,2) + 2*pow(s45,2))) + 
           pow(s14,2)*(3*(s35 + s45)*pow(s23,2) - s45*(-3*s35*s45 + 2*pow(s35,2) + pow(s45,2)) + 
              s23*(7*s35*s45 - pow(s35,2) + 3*pow(s45,2))) + 5*s23*s35*pow(s45,3) + 
           4*pow(s23,2)*pow(s45,3) - 11*pow(s35,2)*pow(s45,3) - 
           s14*(3*(s35 + s45)*pow(s23,3) + 2*s35*s45*(2*s35*s45 + pow(s35,2) + 4*pow(s45,2)) + 
              pow(s23,2)*(8*s35*s45 - 2*pow(s35,2) + 6*pow(s45,2)) + 
              s23*(-(s45*pow(s35,2)) + pow(s35,3) + 11*s35*pow(s45,2) + 3*pow(s45,3))) + 
           s12*((s23 + s35)*pow(s14,3) + (s35 - 3*s45)*pow(s23,3) - pow(s23,4) + 
              pow(s14,2)*(-(s23*(s35 + 3*s45)) - 3*pow(s23,2) + pow(s35,2) + 3*pow(s45,2)) - 
              pow(s23,2)*(5*s35*s45 + 2*pow(s35,2) + 5*pow(s45,2)) + 
              s14*(-((s35 - 6*s45)*pow(s23,2)) + 3*pow(s23,3) + 
                 s23*(5*s35*s45 + pow(s35,2) + 2*pow(s45,2)) + 
                 s35*(10*s35*s45 + pow(s35,2) + 17*pow(s45,2))) + 
              s23*(pow(s35,3) - 11*s35*pow(s45,2) - 4*pow(s45,3)) + 
              s35*(12*s45*pow(s35,2) + pow(s35,3) + 20*s35*pow(s45,2) + 4*pow(s45,3))) + 
           2*s23*pow(s45,4) - 2*s35*pow(s45,4))) - 
     4*(pow(ml1sq,4)*(s14*s23 + s23*s35 + ml2sq*(6*s14 - 4*s45) + 4*s14*s45 - 2*s23*s45 + 
           s12*(-4*s14 + s23 + 4*s45) - pow(s23,2) - 4*pow(s45,2)) + 
        pow(ml1sq,3)*(4*s14*s23*s35 - 17*s14*s23*s45 + 8*s14*s35*s45 - 3*s23*s35*s45 + 
           (8*s14 - 4*s45)*pow(ml2sq,2) + 3*s23*pow(s14,2) + 12*s45*pow(s14,2) + 
           2*ml2sq*(-2*s12*s14 - 7*s14*s23 + 6*s14*s35 + 2*s23*s35 + 2*s12*s45 - 6*s14*s45 + 
              4*s23*s45 - 4*s35*s45 + 8*pow(s14,2)) - 6*s14*pow(s23,2) - 4*s35*pow(s23,2) + 
           7*s45*pow(s23,2) - s12*(-13*s14*s23 - s23*(s35 - 11*s45) + 8*s14*(s35 - 2*s45) + 
              4*s45*(-2*s35 + s45) + 12*pow(s14,2) + pow(s23,2)) + 3*pow(s23,3) + 3*s23*pow(s35,2) - 
           16*s14*pow(s45,2) + 12*s23*pow(s45,2) - 8*s35*pow(s45,2) + 4*pow(s45,3)) + 
        pow(ml1sq,2)*(-22*s14*s23*s35*s45 - 2*(s14 - 2*s45)*pow(ml2sq,3) + 
           s23*(s35 + 3*s45)*pow(s12,2) - s23*pow(s12,3) + 5*s23*s35*pow(s14,2) - 
           27*s23*s45*pow(s14,2) + 16*s35*s45*pow(s14,2) + 3*s23*pow(s14,3) + 12*s45*pow(s14,3) - 
           11*s14*s35*pow(s23,2) + 23*s14*s45*pow(s23,2) + 9*s35*s45*pow(s23,2) - 
           9*pow(s14,2)*pow(s23,2) + 9*s14*pow(s23,3) + 6*s35*pow(s23,3) - 8*s45*pow(s23,3) - 
           3*pow(s23,4) + 5*s14*s23*pow(s35,2) + 4*s14*s45*pow(s35,2) - s23*s45*pow(s35,2) - 
           6*pow(s23,2)*pow(s35,2) + 3*s23*pow(s35,3) + 30*s14*s23*pow(s45,2) - 
           24*s14*s35*pow(s45,2) + 18*s23*s35*pow(s45,2) - 20*pow(s14,2)*pow(s45,2) - 
           12*pow(s23,2)*pow(s45,2) - 4*pow(s35,2)*pow(s45,2) + 
           2*pow(ml2sq,2)*(s23*s35 + 2*s12*(s14 - s23 - s45) + 4*s23*s45 - 
              2*s14*(3*s23 - s35 + 5*s45) + 6*pow(s14,2) - pow(s23,2) + 6*pow(s45,2)) - 
           s12*(-4*s35*(s35 - 2*s45)*s45 - 2*(11*s23 - 8*s35 + 10*s45)*pow(s14,2) + 12*pow(s14,3) - 
              (2*s35 + 7*s45)*pow(s23,2) + 2*pow(s23,3) + 
              s23*(16*s35*s45 + pow(s35,2) - 4*pow(s45,2)) + 
              s14*(-15*s23*s35 + 28*s23*s45 - 24*s35*s45 + 8*pow(s23,2) + 4*pow(s35,2) + 8*pow(s45,2)))
             + 8*s14*pow(s45,3) - 6*s23*pow(s45,3) + 8*s35*pow(s45,3) + 
           ml2sq*(13*s23*s35*s45 + 4*s23*pow(s12,2) + (-23*s23 + 20*s35 - 12*s45)*pow(s14,2) + 
              14*pow(s14,3) - 6*s35*pow(s23,2) - 5*s45*pow(s23,2) + pow(s23,3) + 7*s23*pow(s35,2) - 
              4*s45*pow(s35,2) + s12*(5*s14*s23 - 3*s23*s35 + 16*s14*s45 - 11*s23*s45 - 8*pow(s14,2) + 
                 pow(s23,2) - 8*pow(s45,2)) + 
              s14*(-16*s23*s35 + 15*s23*s45 - 24*s35*s45 + 8*pow(s23,2) + 6*pow(s35,2) - 
                 8*pow(s45,2)) + 4*s23*pow(s45,2) + 8*s35*pow(s45,2) + 8*pow(s45,3))) + 
        ml1sq*(-4*(s14 - s45)*pow(ml2sq,4) - s23*(2*s14 - 2*s23 + s35)*pow(s12,3) - 
           19*s23*s35*s45*pow(s14,2) + 2*s23*s35*pow(s14,3) - 11*s23*s45*pow(s14,3) + 
           8*s35*s45*pow(s14,3) + s23*pow(s14,4) + 4*s45*pow(s14,4) + 20*s14*s35*s45*pow(s23,2) - 
           8*s35*pow(s14,2)*pow(s23,2) + 13*s45*pow(s14,2)*pow(s23,2) - 4*pow(s14,3)*pow(s23,2) + 
           10*s14*s35*pow(s23,3) - 9*s14*s45*pow(s23,3) - 9*s35*s45*pow(s23,3) + 
           6*pow(s14,2)*pow(s23,3) - 4*s14*pow(s23,4) - 4*s35*pow(s23,4) + 3*s45*pow(s23,4) + 
           pow(s23,5) - 5*s14*s23*s45*pow(s35,2) + 2*s23*pow(s14,2)*pow(s35,2) + 
           4*s45*pow(s14,2)*pow(s35,2) - 6*s14*pow(s23,2)*pow(s35,2) + s45*pow(s23,2)*pow(s35,2) + 
           4*pow(s23,3)*pow(s35,2) + s23*pow(s12,2)*
            (-4*s23*s35 - 5*s23*s45 + 2*s35*s45 + 3*s14*(s35 + 2*s45) + pow(s35,2)) + 
           2*s14*s23*pow(s35,3) - s23*s45*pow(s35,3) - 4*pow(s23,2)*pow(s35,3) + s23*pow(s35,4) + 
           28*s14*s23*s35*pow(s45,2) + 16*s23*pow(s14,2)*pow(s45,2) - 16*s35*pow(s14,2)*pow(s45,2) - 
           8*pow(s14,3)*pow(s45,2) - 12*s14*pow(s23,2)*pow(s45,2) - 14*s35*pow(s23,2)*pow(s45,2) + 
           4*pow(s23,3)*pow(s45,2) - 8*s14*pow(s35,2)*pow(s45,2) + 4*s23*pow(s35,2)*pow(s45,2) + 
           2*pow(ml2sq,3)*(-2*s23*s35 + s14*(s23 - 4*s35 - 2*s45) + 2*s12*(s14 - s23 - s45) + 
              4*s35*s45 - 2*pow(s14,2) + 4*pow(s45,2)) - 
           s12*((-9*s23 + 8*s35 - 8*s45)*pow(s14,3) + 4*pow(s14,4) + (4*s35 - 2*s45)*pow(s23,3) + 
              s14*(-(s23*s35*(s35 - 26*s45)) + 8*s35*s45*(-s35 + s45) + (s35 - 5*s45)*pow(s23,2) + 
                 5*pow(s23,3)) - 3*pow(s23,4) + s23*s35*(3*s35*s45 + pow(s35,2) - 6*pow(s45,2)) + 
              4*pow(s35,2)*pow(s45,2) - pow(s23,2)*(9*s35*s45 + 6*pow(s35,2) + 2*pow(s45,2)) + 
              pow(s14,2)*(s23*(-13*s35 + 15*s45) + 3*pow(s23,2) + 
                 4*(-4*s35*s45 + pow(s35,2) + pow(s45,2)))) - 4*s14*s23*pow(s45,3) + 
           8*s14*s35*pow(s45,3) - 8*s23*s35*pow(s45,3) + 4*pow(s14,2)*pow(s45,3) + 
           2*pow(s23,2)*pow(s45,3) + 4*pow(s35,2)*pow(s45,3) + 
           pow(ml2sq,2)*(-5*s23*s35*s45 + 6*s23*pow(s12,2) - (7*s23 + 16*s45)*pow(s14,2) + 
              4*pow(s14,3) - 8*s35*pow(s23,2) - 3*s45*pow(s23,2) + 
              s12*(s14*(-9*s23 + 8*s35) - 4*s45*(2*s35 + s45) - s23*(s35 + 5*s45) + 4*pow(s14,2) + 
                 7*pow(s23,2)) + pow(s23,3) - 3*s23*pow(s35,2) + 4*s45*pow(s35,2) - 4*s23*pow(s45,2) + 
              16*s35*pow(s45,2) + s14*(-2*s23*s35 + 17*s23*s45 - 16*s35*s45 + 4*pow(s23,2) - 
                 4*pow(s35,2) + 8*pow(s45,2)) + 4*pow(s45,3)) - 
           2*ml2sq*(s23*(-3*s14 + 4*s23 - 2*s35 - 3*s45)*pow(s12,2) + s23*pow(s12,3) + 
              2*(2*s23 - 2*s35 + s45)*pow(s14,3) - 2*pow(s14,4) + 6*s35*s45*pow(s23,2) - 
              s35*pow(s23,3) + s45*pow(s23,3) + pow(s23,4) + s23*s45*pow(s35,2) + 
              5*pow(s23,2)*pow(s35,2) - s23*pow(s35,3) - 
              pow(s14,2)*(-8*s23*s35 + 3*s23*s45 - 8*s35*s45 + pow(s23,2) + 2*pow(s35,2) - 
                 4*pow(s45,2)) + 4*s23*s35*pow(s45,2) + 2*pow(s23,2)*pow(s45,2) - 
              4*pow(s35,2)*pow(s45,2) + 2*s12*
               (s35*s45*(s35 + 2*s45) - (s23 + 3*s45)*pow(s14,2) + pow(s14,3) - 
                 (3*s35 + 2*s45)*pow(s23,2) + 
                 s14*(2*s23*s35 + 5*s23*s45 - 2*s35*s45 - pow(s35,2) + 2*pow(s45,2))) + 
              2*s23*pow(s45,3) - 4*s35*pow(s45,3) - 
              s14*(3*s35*pow(s23,2) + 2*pow(s23,3) - 6*s45*pow(s35,2) + 
                 s23*(12*s35*s45 - pow(s35,2) + 6*pow(s45,2)) + 4*pow(s45,3)))) + 
        s23*(-((s12 + s14 - 3*s23 + 3*s35 - 2*s45)*pow(ml2sq,4)) + 
           pow(ml2sq,3)*(2*s23*s35 + 9*s23*s45 - 5*s35*s45 + s14*(-2*s23 - 2*s35 + s45) - 
              s12*(s14 + 3*s23 - 3*s35 + 5*s45) + 2*pow(s12,2) - pow(s14,2) + 3*pow(s23,2) - 
              7*pow(s35,2) + 4*pow(s45,2)) + 
           pow(ml2sq,2)*(7*s23*s35*s45 + (4*s14 - 2*s23 + s35 + 3*s45)*pow(s12,2) - pow(s12,3) + 
              (-s23 + s35 - 3*s45)*pow(s14,2) + pow(s14,3) + 4*s35*pow(s23,2) + 6*s45*pow(s23,2) + 
              pow(s23,3) - 4*s23*pow(s35,2) - 13*s45*pow(s35,2) - 5*pow(s35,3) - 
              s14*(5*s23*s35 + 3*s23*s45 + 2*s35*s45 + pow(s23,2) + pow(s35,2) - 6*pow(s45,2)) + 
              s12*(4*s23*s35 + s14*(6*s23 + s35 - 8*s45) - 7*s23*s45 - 6*pow(s23,2) + 7*pow(s35,2) - 
                 4*pow(s45,2)) + 8*s23*pow(s45,2) - 2*s35*pow(s45,2) + 2*pow(s45,3)) - 
           ml2sq*((2*s14 - 2*s23 + s35)*pow(s12,3) + (4*s23 + s45)*pow(s14,3) - pow(s14,4) + 
              pow(s12,2)*(6*s14*s23 + s35*(s35 - 2*s45) + s23*(6*s35 + s45) - 3*s14*(s35 + 2*s45) - 
                 2*pow(s14,2) - 4*pow(s23,2)) - 5*s35*s45*pow(s23,2) - 
              pow(s14,2)*(5*s23*s45 + 3*s35*s45 + 6*pow(s23,2)) - 3*s45*pow(s23,3) - pow(s23,4) + 
              3*s23*s45*pow(s35,2) - 2*pow(s23,2)*pow(s35,2) + 4*s23*pow(s35,3) + 7*s45*pow(s35,3) + 
              pow(s35,4) - 6*s23*s35*pow(s45,2) - 4*pow(s23,2)*pow(s45,2) + 8*pow(s35,2)*pow(s45,2) + 
              s14*(7*s45*pow(s23,2) + 4*pow(s23,3) + 2*s23*(4*s35*s45 + pow(s35,2) + 2*pow(s45,2)) - 
                 s45*(pow(s35,2) + 4*pow(s45,2))) + 
              s12*((-s23 + 3*s35 + s45)*pow(s14,2) + pow(s14,3) + 6*(s35 + s45)*pow(s23,2) + 
                 pow(s23,3) + s23*(s35*s45 - 8*pow(s35,2) + 2*pow(s45,2)) + 
                 s35*(-7*s35*s45 - 3*pow(s35,2) + 2*pow(s45,2)) + 
                 s14*(-9*s23*s35 - 7*s23*s45 + 2*s35*s45 - pow(s23,2) + pow(s35,2) + 8*pow(s45,2))) - 
              2*s23*pow(s45,3)) - (s14 - s23 + s35)*
            ((-2*s14*s35 + 3*s23*s35 - 3*s14*s45 + 2*s23*s45 + s35*s45)*pow(s12,2) + 
              (s14 - s23)*pow(s12,3) + (s23 + s45)*(s35 + 2*s45)*pow(s14,2) - s45*pow(s14,3) + 
              s35*(s23 + s45)*(2*s23*s45 + 2*s35*s45 + pow(s23,2) + pow(s35,2) + 2*pow(s45,2)) - 
              s14*(2*s23*s45*(2*s35 + s45) + (2*s35 + s45)*pow(s23,2) + 
                 s45*(4*s35*s45 + pow(s35,2) + 2*pow(s45,2))) + 
              s12*(-2*s35*s45*(s35 + s45) - (3*s23 + 2*s45)*pow(s14,2) + pow(s14,3) - 
                 2*s45*pow(s23,2) - pow(s23,3) - s23*(4*s35*s45 + 3*pow(s35,2) + 2*pow(s45,2)) + 
                 s14*(4*s23*s45 + 3*pow(s23,2) + pow(s35 + 2*s45,2))))))*pow(mw2c*conj(mw2c),1/2.));

	// Averaging over photon polarisations
	prefactor /= 2.;
	// Include the prefactor and averaging/sum factors
	ME2 *= prefactor;
	return ME2.real();
	// Validated against Recola up to a factor of two
}


double ME2_Analytic::nu1gamma_l1nu1l1x(int i1, int i2, int i3, int i4, int i5, KinematicData &Kin, int f1, int f2){

if( !is_neutrino(f1) ){
		cerr << "nu1gamma_nu1f2f2x: fermion line 1 is not a neutrino\n";
		abort();
	}

	// Using another basis of kinematics
	double s12 = 2.0*Kin.pij(i1,i2);
	double s14 = 2.0*Kin.pij(i1,i4);	
	double s23 = 2.0*Kin.pij(i2,i3);
	// double s25 = 2.0*Kin.pij(i2,i5);	
	double s35 = 2.0*Kin.pij(i3,i5);		
	double s45 = 2.0*Kin.pij(i4,i5);
	// fermion masses
	double mlsq = Kin.p(3).m2();

	// Implemented as
	// |M|^2 = |M_W|^2 + |M_Z|^2 + M_Z M_W* + M_W M_Z*
	// where the intereference term is: M_Z M_W* + M_W M_Z*

	cout << abs(f1)-1 << endl;
	// f3, f4 a generic fermion line
	complex<double> gLf, gRf; double Qf;
	assign_EW_charges(Qf,gLf,gRf,abs(f1)-1);


	complex<double> chiw45 = 1. / (mlsq + s45 - MW2C);
	complex<double> mw2c = MW2C;
	// The following option to drop width effects in numerator (a sort of pole approximation)
	complex<double> chiw13 = 1. / (mw2c + s12 - s23 - s45 );	
	// prefactor
	complex<double> prefac_interference = 8.*gLnu*conj(ALPHA)*pow(ALPHA,2)*pow(MW2C,-1)*pow(pi,3)*pow(MZ2C + s14,-1)*pow(s23,-2)*pow(2*mlsq + s14 - s23 + s35,-2)*pow(SW2,-1);
	// interference term
	complex<double> interference = -(chiw45*(2*mlsq + s14 - s23 + s35)*(gLf*((8*s14*(-1. + chiw13*s23) - 4.*chiw13*s23*(s23 + s45))*pow(mlsq,4) + 
           2*pow(mlsq,3)*(s12*s23*(2. + chiw13*(8.*mw2c + s23)) - 2*s14*s35 + s14*s23*(3. - 2.*chiw13*s12 + 2.*chiw13*s35) - 2.*chiw13*s23*(4.*mw2c + s35)*(s23 + s45) + 
              2.*(-1. + chiw13*s23)*pow(s14,2) - 2.*chiw13*s14*pow(s23,2)) + 
           pow(mlsq,2)*(-8.*mw2c*(-(s12*s23*(-1. + chiw13*(4*s23 + 2*s35 + 3*s45))) + chiw13*s23*pow(s12,2) + 2.*(-1. + chiw13*s23)*pow(s14,2) + 
                 (s23 + s45)*(-2*s45 + s23*(-1. + 3.*chiw13*s35 + 4.*chiw13*s45) + chiw13*pow(s23,2)) - 
                 s14*(-4*s45 + s23*(-2. + chiw13*s12 + 3.*chiw13*s45) + 2.*chiw13*pow(s23,2))) + 
              s23*(2*s12*s35 + 3*s14*s35 - chiw13*s12*(2*s14 - s23)*(s14 + s35) - 2*s14*s23*(1. + chiw13*s35) - 2*pow(s12,2) + pow(s14,2) + 
                 chiw13*(-s23 + s45)*pow(s14,2) + 2.*chiw13*s14*pow(s23,2) - chiw13*(s23 + s45)*pow(s35,2))) - 
           4.*mlsq*mw2c*(s23*(-1. + chiw13*(4*s23 + s35 + s45))*pow(s12,2) + 2.*(-1. + chiw13*s23)*pow(s14,3) - 
              pow(s14,2)*(2*(s35 - 2*s45) + s23*(-2. + chiw13*(s12 - 2*s35 + 3*s45)) + chiw13*pow(s23,2)) + 
              s14*(2*(2*s35 - s45)*s45 + s12*s23*(1. + chiw13*s23 - chiw13*(s35 + s45)) + s23*(s35 - 2.*chiw13*s35*s45 + s45*(-3. + 2.*chiw13*s45)) + 
                 chiw13*s23*pow(s12,2) - (-1. + chiw13*(s35 + s45))*pow(s23,2) - 3.*chiw13*pow(s23,3)) + 
              s12*s23*(s35 - 5.*chiw13*s35*s45 + s23*(1. - 7.*chiw13*s35 - 8.*chiw13*s45) + s45*(2. - 3.*chiw13*s45) - 2.*chiw13*pow(s23,2) - chiw13*pow(s35,2)) + 
              (s23 + s45)*(-2*s35*s45 + (-1. + chiw13*(s35 + 3*s45))*pow(s23,2) + 2.*chiw13*pow(s23,3) + 
                 s23*(s45*(-1. + 2.*chiw13*s45) + s35*(-1. + 6.*chiw13*s45) + 3.*chiw13*pow(s35,2)))) + 
           2.*mw2c*s23*(chiw13*s23*pow(s12,3) - (-2. + chiw13*(2*s23 + s45))*pow(s14,3) + 
              pow(s14,2)*(s35 - 2*s23*(2. + chiw13*(s35 - 3*s45)) - chiw13*s35*s45 + 2*s45*(-2. + chiw13*s45) + 4.*chiw13*pow(s23,2)) - 
              pow(s12,2)*(s14*(-1. + chiw13*s45) + chiw13*(s35*s45 + 3*s23*(s35 + s45) + pow(s23,2))) + 
              s14*((2. + 4.*chiw13*s35 - 5.*chiw13*s45)*pow(s23,2) - 2.*chiw13*pow(s23,3) + (1. - chiw13*s45)*pow(s35,2) - 
                 s23*(s35*(2. - 4.*chiw13*s45) + s45*(-4. + 5.*chiw13*s45) + chiw13*pow(s35,2)) + (3. - 2.*chiw13*s45)*pow(s45,2)) - 
              s35*(s23 + s45)*(-s45 + s23*(-1. + 3.*chiw13*s45) + 2.*chiw13*pow(s23,2) + chiw13*(2*s35*s45 + pow(s35,2) + 2*pow(s45,2))) + 
              s12*(-(s23*(s23 + s45)) + s14*(s23 - 2*(s35 + s45)) + 
                 chiw13*(s35*s45*(2*s35 + 3*s45) + (-s23 + s35 - 2*s45)*pow(s14,2) + pow(s14,3) + (s35 + 3*s45)*pow(s23,2) + 2*pow(s23,3) + 
                    s23*(6*s35*s45 + 3*pow(s35,2) + 2*pow(s45,2)) + s14*(2*s23*s45 - 2*pow(s23,2) + 3*pow(s45,2)))))) + 
        gRf*mlsq*(4*pow(mlsq,2)*(s12*s23*(-1. + chiw13*(s23 + s45)) - 2.*mw2c*(s14*(2. - 2.*chiw13*s23) + chiw13*s23*(2*s12 + s23 + s45)) - 
              (-1. + chiw13*s23)*(3*s23*s45 + pow(s23,2) + 2*pow(s45,2))) + 
           2.*mlsq*(-((s23 + s45)*(-(s23*s35) - s23*s45 - 2*s35*s45 + s14*(-1. + chiw13*s23)*(s23 + 2*s45) + chiw13*s23*(s23 + s35 + s45)*(s23 + 2*s45))) - 
              s23*(-1. + chiw13*(2*s23 + s45))*pow(s12,2) + 2.*mw2c*
               (s14*(s23 - 4.*chiw13*s12*s23 - 2*s35 + 2.*chiw13*s23*s35) + 
                 s23*(s12*(2. + chiw13*s23 - 4.*chiw13*s35) - 2.*chiw13*s35*(s23 + s45) + 2.*chiw13*pow(s12,2)) + 2.*(-1. + chiw13*s23)*pow(s14,2)) + 
              s12*s23*(-s35 - 2*s45 + chiw13*s35*s45 + s23*(-1. + chiw13*s35 + 7.*chiw13*s45) + s14*(-1. + chiw13*(s23 + s45)) + 3.*chiw13*pow(s23,2) + 
                 3.*chiw13*pow(s45,2))) + s23*(3*s14*s23*s45 + s23*s35*s45 + chiw13*s23*pow(s12,3) + s14*pow(s23,2) - 4.*chiw13*s14*s45*pow(s23,2) - 
              4.*chiw13*s35*s45*pow(s23,2) + s12*(-(s23*s45) + s14*(s23 + s45)*(-2. + 2.*chiw13*s23 + 3.*chiw13*s45) + 
                 chiw13*(s23 + s45)*(3*s35*s45 + 2*s23*(s35 + s45) + pow(s23,2))) - 
              pow(s12,2)*(s14*(-1. + chiw13*(s23 + s45)) + chiw13*(s35*s45 + s23*(s35 + 3*s45) + 2*pow(s23,2))) - chiw13*s14*pow(s23,3) - chiw13*s35*pow(s23,3) + 
              2.*mw2c*(2.*(-1. + chiw13*(s14 + s35))*pow(s12,2) + (-1. + chiw13*(s23 + s45))*pow(s14,2) + 
                 s12*(s14*(4. - 3.*chiw13*s23 - 4.*chiw13*s35) + s35*(2. + chiw13*s23 - 2.*chiw13*s35) - 2.*chiw13*pow(s14,2)) + 
                 s14*(-2*s23 + s35 + 2.*chiw13*pow(s23,2)) - chiw13*(s23 + s45)*pow(s35,2)) + 3*s14*pow(s45,2) - 5.*chiw13*s14*s23*pow(s45,2) + s35*pow(s45,2) - 
              5.*chiw13*s23*s35*pow(s45,2) - 2.*chiw13*s14*pow(s45,3) - 2.*chiw13*s35*pow(s45,3))))) - 
   chiw13*s23*(gRf*mlsq*(4.*(2.*mw2c*(2*s12 + s14) + (s12 - s23 - s45)*(s12 - 2*s14 - s23 - 2*s35 - s45))*pow(mlsq,2) + 8*(-s12 + s23 + s45)*pow(mlsq,3) + 
         2.*mlsq*(2*s14*s23*s35 + 3*s14*s23*s45 + 2*s14*s35*s45 + s23*s35*s45 + mw2c*(4*s12*(s14 - s23 + 2*s35) + 2*s14*(2*s14 - s23 + 4*s35) - 4*pow(s12,2)) + 
            (2*s14 + s23 + s35)*pow(s12,2) + s23*pow(s14,2) + s45*pow(s14,2) + 2*s14*pow(s23,2) + s35*pow(s23,2) + s45*pow(s23,2) + s23*pow(s35,2) + 
            s45*pow(s35,2) - s12*(2*s23*s35 + s23*s45 + 2*s35*s45 + 2*s14*(2*s23 + s35 + 2*s45) + pow(s14,2) + pow(s23,2) + pow(s35,2)) + 2*s14*pow(s45,2) + 
            s23*pow(s45,2)) + (s14 - s23 + s35)*(4.*mw2c*s12*s35 + 2.*mw2c*s14*(s14 + 2*s23 + 3*s35) + s12*s23*s45 + s14*s23*s45 - s23*s35*s45 - 
            2*s12*s14*(s23 + s45) - 4.*mw2c*pow(s12,2) + s14*pow(s12,2) + s14*pow(s23,2) + s14*pow(s45,2) - s35*pow(s45,2))) + 
      gLf*((8*s12 - 4*s14)*pow(mlsq,4) - 2*pow(mlsq,3)*(s14*(2*s14 - s23) + s12*(-6*s14 + 2*s23 - 4*s35) + 8.*mw2c*(s12 - s23 - s45) + 2*pow(s12,2)) + 
         pow(mlsq,2)*(8.*mw2c*(-2*s12*(s14 + s23 + s35 + s45) + (s23 + s45)*(2*s35 + s45) + s14*(s35 + 2*s45) + pow(s12,2)) - 
            (s14 - s23 + s35)*(s14*(s14 - 2*s23 - s35) - 2*s12*(2*s14 + s35) + 2*pow(s12,2))) + 
         4.*mlsq*mw2c*((2*s14 + s23 + s35)*pow(s12,2) + (s23 + s45)*pow(s14,2) - 
            s12*(-(s14*s23) + 4*s14*(s35 + s45) + s23*(2*s35 + s45) + s35*(s35 + 2*s45) + pow(s14,2)) + (s23 + s45)*(s23*(-s35 + s45) + pow(s23,2) + pow(s35,2)) + 
            s14*(s23*(s35 - 3*s45) - 2*pow(s23,2) + 2*(3*s35*s45 + pow(s35,2) + pow(s45,2)))) + 
         2.*mw2c*(s14 - s23 + s35)*(s12*s23*(s23 + s45) - s12*s14*(s23 + 2*(s35 + s45)) + s14*pow(s12,2) - s35*pow(s14,2) + 
            s14*(2*s23*s35 + 4*s35*s45 + pow(s35,2) + pow(s45,2)) - s35*pow(s23 + s45,2))));

    // Averaging over photon polarisations
	prefac_interference /= 2.;
	// Result for the interference
	complex<double> ME2 = interference * prefac_interference;

	// Add to it the results for |M_W^2|
	double Mw_piece = nu1gamma_l1nu2l2x(i1,i2,i3,i4,i5,Kin,f1,f2);
	// Then the |M_Z|^2, obtained from nu1gamma_nu1f2f2x (with leptons, and the momentum swap for neutrino position)
	// i4 neutrino in third entry i3 (where it is defined for nu1gamma_nu1f2f2x)
	double Mz_piece = nu1gamma_nu1f2f2x(i1,i2,i4,i3,i5,Kin,f1,abs(f2)-1);

	cout << "|M_w|^2 = " << Mw_piece << endl;
	cout << "|M_z|^2 = " << Mz_piece << endl;
	cout << "|M_inter|^2 = " << ME2 << endl;	

	double mfsq = mlsq;
	complex<double> chiz14 = 1. / ( MZ2C + s14 );
	// direct MZ calculation part
	complex<double> ME2_Z = -(conj(gRf)*(-(gRf*(2*mfsq + s23)*s45*(s23 + s45)) + 2.*gLf*s14*pow(mfsq,2))*pow(s23,-2)) - 
   conj(gLf)*(gLf*(2*mfsq + s23)*(s14 - s45)*(-s14 + s23 + s45) + 2.*gRf*s14*pow(mfsq,2))*pow(s23,-2) + 
   conj(gRf)*(2.*gLf*mfsq*s14*(mfsq + s14 - s23 + s35) + 
      gRf*(s12 - s23 - s45)*(2*mfsq*(s12 - s14 - s35) + (s14 - s23 + s35)*s45 - 4*pow(mfsq,2)))*
    pow(2*mfsq + s14 - s23 + s35,-2) + conj(gLf)*
    (2.*gRf*mfsq*s14*(mfsq + s14 - s23 + s35) + 
      gLf*(s14 - s23 + s35)*(s14 - s23 - s45)*(-s12 + s35 + s45) - 4.*gLf*s12*pow(mfsq,2) + 
      2.*gLf*mfsq*((s14 - s23 + s35)*(s14 - s23 - s45) - s12*(s35 + s45) + pow(s12,2)))*
    pow(2*mfsq + s14 - s23 + s35,-2) + conj(gRf)*pow(s23,-1)*pow(2*mfsq + s14 - s23 + s35,-1)*
    (2.*gLf*mfsq*(2*mfsq*s12 + s12*(s14 + s35) + s14*(-s23 + s35) - pow(s12,2)) + 
      gRf*(2*mfsq*(s12 - s23 - s45)*(s12 - s14 - s35 - s45) + s14*s23*s45 - s23*s35*s45 - 
         s12*(s23*s45 + 2*s14*(s23 + s45)) + 4*(-s12 + s23 + s45)*pow(mfsq,2) + s14*pow(s12,2) + 
         s14*pow(s23,2) + 2*s45*pow(s23,2) + s14*pow(s45,2) + 2*s23*pow(s45,2) - s35*pow(s45,2))) + 
   conj(gLf)*pow(s23,-1)*pow(2*mfsq + s14 - s23 + s35,-1)*
    (2.*gRf*mfsq*(2*mfsq*s12 + s12*(s14 + s35) + s14*(-s23 + s35) - pow(s12,2)) + 
      gLf*(-2*s14*s23*s45 + 4*s14*s35*s45 - s12*s23*(s23 + s45) + s12*s14*(s23 - 2*(s35 + s45)) + 
         4*(-s12 + s23 + s45)*pow(mfsq,2) + s14*pow(s12,2) + 
         2*mfsq*(s14*(-s23 + s35 + s45) + (s23 + s45)*(s23 + s35 + s45) - 
            s12*(s14 + s23 + s35 + 2*s45) + pow(s12,2)) - s35*pow(s14,2) + s35*pow(s23,2) + 
         2*s45*pow(s23,2) + s14*pow(s35,2) + s14*pow(s45,2) + 2*s23*pow(s45,2) - s35*pow(s45,2)));

    ME2_Z *= 32.*chiz14*gLnu*conj(ALPHA)*conj(chiz14)*conj(gLnu)*pow(ALPHA,2)*pow(pi,3);

    cout << "Z ratio part = " << Mz_piece / ME2_Z.real() << endl;

	return (2.*ME2.real() + Mw_piece + Mz_piece*2);
	// Validated against Recola up to a factor of two
}


// An interface to recola squared amplitudes
double ME2_Analytic::compute_process_recola(KinematicData &Kin, int process, bool virt){
	// First translate the Kinematic data
	int n = Kin.length();
	// Now construct the momentum set
	double pset[n][4];
	for( int i = 0; i < n; i++){
		for( int j = 0; j < 4; j++)
			pset[i][j] = Kin.p(i+1).pi(j);
	}

	// for( int i = 0; i < n; i++)
	// 	cout << pset[i][0] << " " 
	// 		 << pset[i][1] << " "
	// 		 << pset[i][2] << " "
	// 		 << pset[i][3] << "\n";
	string order = (virt)? "NLO":"LO";

	double temp[2];
	Recola::compute_process_rcl(process,pset,order,temp);
	// Return the sum of squared amplitude at LO and NLO
	if( order == "LO" ) return temp[0];
	// Return just 2 Re | M1 M0 |, obtain with lower stat. precision for efficiency
	else{
		return temp[1];
	}

	// |M0|^2 + 2 Re | M1 M0 |
	// if( order == "NLO" ) return temp[0]+temp[1];
}



