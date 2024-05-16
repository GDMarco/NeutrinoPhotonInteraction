#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include "var.hh"
#include "tools.hh"
#include "ME2_Analytic.hh"
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

	complex<double> ME2 = 192.* ALPHA * conj(ALPHA) * pisq *chiz*gLnu*conj(chiz)*conj(gLnu) * 
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
	// Kinematics
	double s12 = ac_i1 * ac_i2 * 2.0*Kin.pij(i1,i2);
	double s13 = ac_i1 * ac_i3 * 2.0*Kin.pij(i1,i3);

	complex<double> ME2 = 64.*ALPHA*(2.*MZ2C + s12)*conj(ALPHA)*(s12 + 2.*conj(MZ2C))*pow(gLnu,2)*pow(pi,2)*pow(s12,2)*
   pow(MZ2C + s12 - s13,-1)*pow(MZ2C + s13,-1)*pow(conj(gLnu),2)*pow(s12 - s13 + conj(MZ2C),-1)*
   pow(s13 + conj(MZ2C),-1);

   return ME2.real();
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

	complex<double> ME2 = 64.*ALPHA*conj(ALPHA)*pow(gLnu,2)*pow(pi,2)*pow(s12,2)*pow(MZ2C + s13,-1)*pow(conj(gLnu),2)*
   pow(s13 + conj(MZ2C),-1);

   return ME2.real();
}




////////////////////////////////////////////////
// (anti)neutrino + photon scattering results //
////////////////////////////////////////////////

double ME2_Analytic::nu1gamma_nu1f2f2x(int i1, int i2, int i3, int i4, int i5, KinematicData &Kin, int f1, int f2){

	if( !is_neutrino(f1) ){
		cerr << "nu1gamma_nu1f2f2x: fermion line 1 is not a neutrino\n";
		abort();
	}

	// Take fermion mass from kinematic structure
	double mfsq = Kin.p(i4).m2();
	// Other kinematics
	double s12 = 2.0*Kin.pij(i1,i2);
	double s13 = 2.0*Kin.pij(i1,i3);
	double s15 = 2.0*Kin.pij(i1,i5);
	double s24 = 2.0*Kin.pij(i2,i4);
	double s25 = 2.0*Kin.pij(i2,i5);

	// f3, f4 a generic fermion line
	complex<double> gLf, gRf; double Qf;
	assign_EW_charges(Qf,gLf,gRf,abs(f2));

	// Include mass for quarks. The f2 f2x quark system should also be at least as massive as pion
	complex<double> chiz13 = 1. / ( MZ2C + s13 );
	// alpha alpha*  (Re[alpha])
	complex<double> prefactor = 32.*chiz13*gLnu*conj(chiz13)*conj(gLnu)*pow(pi,3)*pow(Qf,2)*pow(s24,-2)*pow(s25,-2)
	 * ALPHA * conj(ALPHA) * ALPHA.real();

	complex<double> ME2 = 3.*(conj(gRf)*(2.*gLf*mfsq*s24*s25*(s12*(s24 + s25) + s13*(-s13 + s24 + s25) - pow(s12,2)) + 
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
	if( abs(f2) > 6 ) prefactor /= 3.;
	// Include the prefactor and averaging/sum factors
	ME2 *= prefactor;
	return ME2.real();
}







// Below functions were used for validating crossing relations

// nu1 + nu1bar > nu1 + nu1bar [should be related to the above result via crossing i2 <> i4]
// double ME2_Analytic::nu1nu1x_nu1nu1x(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2 ){

// 	// Analytic continuation of momenta
// 	double ac_i1 = (i1 < 3)?  +1.: -1.;
// 	double ac_i2 = (i2 < 3)?  +1.: -1.;
// 	double ac_i3 = (i3 < 3)?  -1.: +1.;
// 	// Kinematics
// 	double s12 = ac_i1 * ac_i2 * 2.0*Kin.pij(i1,i2);
// 	double s13 = ac_i1 * ac_i3 * 2.0*Kin.pij(i1,i3);

// 	complex<double> ME2 = 64.*ALPHA*(2.*MZ2C - s12 + s13)*conj(ALPHA)*(s12 - s13 - 2.*conj(MZ2C))*pow(gLnu,2)*pow(pi,2)*
//    pow(MZ2C - s12,-1)*pow(s12 - s13,2)*pow(MZ2C + s13,-1)*pow(conj(gLnu),2)*pow(s12 - conj(MZ2C),-1)*
//    pow(s13 + conj(MZ2C),-1);

//    return ME2.real();
// }

// // nu1x + nu2 > nu1x + nu2 [related to the above by crossing, not needed]
// double ME2_Analytic::nu1xnu2_nu1xnu2(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2 ){
// 	return 0.0;
// }



