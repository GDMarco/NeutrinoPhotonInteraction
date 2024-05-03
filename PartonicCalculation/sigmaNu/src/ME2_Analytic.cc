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

// Expressions for heavy-quark pair production in e-(p1) + e+(p2) > gamma/Z > Q(p3) + Qbar(p4) [ + g(p5) ]
void assign_EW_charges(double &Q, complex<double> &gL, complex<double> &gR, int pdg){
	if (pdg == 11){
		gL = gLl;
		gR = gRl;
		Q = Ql;
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


// // e- (p1) + e+(p2) > Q(p3) + Qbar(p4)
// double ME2_Analytic::eeB0g0NCM(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2){
// 	double s12 = 2.0*Kin.pij(i1,i2);
// 	double s13 = 2.0*Kin.pij(i1,i3);
// 	double msqQ = Kin.p(i4).m2();
// 	// Assign the charges to fermions f1, and fermions f2
// 	complex<double> gLf1, gRf1; double Q1;
// 	assign_EW_charges(Q1,gLf1,gRf1,f1);
// 	complex<double> gLf2, gRf2; double Q2;
// 	assign_EW_charges(Q2,gLf2,gRf2,f2);
// 	// the summed/averaged ME2
// 	complex<double> prefac = 48.0 * ALPHA * conj(ALPHA) * pisq;
// 	complex<double> chiz_s = 1.0 / ( s12 - MZ2C );
// 	complex<double> denom = 1.0 / (pow(s12,2) );
// 	complex<double> ME2 = ((msqQ*s12*(Q1*Q2 + gLf1*gRf2*chiz_s*s12) + (Q1*Q2 + gLf1*gLf2*chiz_s*s12)*(s12-s13)*(s12-s13))* (Q1*Q2+s12*conj(gLf1)*conj(gLf2)*conj(chiz_s))
// 		+(msqQ*s12*(Q1*Q2 + gRf1*gRf2*chiz_s*s12) + (Q1*Q2 + gLf2*gRf1*chiz_s*s12)*s13*s13)*(Q1*Q2+s12*conj(gLf2)*conj(gRf1)*conj(chiz_s))
// 		+(msqQ*s12*(Q1*Q2 + gLf1*gLf2*chiz_s*s12) + (Q1*Q2 + gLf1*gRf2*chiz_s*s12)*s13*s13)*(Q1*Q2+s12*conj(gLf1)*conj(gRf2)*conj(chiz_s))
// 		+(msqQ*s12*(Q1*Q2 + gLf2*gRf1*chiz_s*s12) + (Q1*Q2 + gRf1*gRf2*chiz_s*s12)*(s12-s13)*(s12-s13))*(Q1*Q2+s12*conj(gRf1)*conj(gRf2)*conj(chiz_s)));

// 	// Include prefactorsand denom.
// 	ME2*=prefac*denom;
// 	return ME2.real();
// }


double ME2_Analytic::nunux_ffx(int i1, int i2, int i3, int i4, KinematicData &Kin, int f1, int f2 ){

	// f1 line should always be the neutrino line, check this

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

	return ME2.real();

}

