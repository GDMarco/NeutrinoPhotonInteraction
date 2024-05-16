#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include "tools.hh"
#include "cuts.hh"

using namespace std;

p4vec::p4vec() {
	p = {0.,0.,0.,0.};
}

p4vec::p4vec(double E, double px, double py, double pz){
	p = {E, px, py, pz};
}

// General boost
void p4vec::boost(double bx, double by, double bz){
    // beta_v = {bx,by,pz}
    double b2 = bx*bx + by*by + bz*bz;
    double gam  = 1./sqrt(1.-b2);
    double bp = p[1]*bx + p[2]*by + p[3]*bz;
    double gam2 = (gam-1.)/b2;
    if( b2 <= 0 ) gam2 = 0.; // For numerical stability when beta = 0.
    p[1] = p[1] + gam2*bp*bx + gam*bx*p[0];
    p[2] = p[2] + gam2*bp*by + gam*by*p[0];
    p[3] = p[3] + gam2*bp*bz + gam*bz*p[0];
    p[0] = gam*( p[0] + bp );
}

// General rotation around axis = i, by alpha (given by cos[alpha])
void p4vec::rot(double c_a, int i){
    p[0] = p[0]; // Energy component un-touched
    array< double, 3 > temp = {p[1],p[2],p[3]};
    array< array<double, 3>, 3 > rotation;
    double s_a = sqrt(1.-pow(c_a,2));
    if( i == 1 ){
        rotation[0] = { { 1.0, 0.0, 0.0} };
        rotation[1] = { { 0.0, c_a, s_a} };
        rotation[2] = { { 0.0,-s_a, c_a} };        
    }
    else if( i == 2 ){
        rotation[0] = { { c_a, 0.0,-s_a} };
        rotation[1] = { { 0.0, 1.0, 0.0} };
        rotation[2] = { { s_a, 0.0, c_a} };         
    }
    else if( i == 3 ){
        rotation[0] = { { c_a, s_a, 0.0} };
        rotation[1] = { {-s_a, c_a, 0.0} };
        rotation[2] = { { 0.0, 0.0, 1.0} };           
    }
    else{
        cout << "p4vec::rot unsuported axis rotation " << i << endl;
        abort();
    }
    for( int j = 0; j < 3; j++){
        p[j+1] = 0.;
        for( int k = 0; k < 3; k++){
            p[j+1] += rotation[j][k] * temp[k];
        }
    }
}

double p4vec::rap(){
    return (0.5*log( (p[0]+p[3])/(p[0]-p[3]) ) );
    // For better numerical stability
}

double p4vec::modp(){
    return sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
}

double p4vec::eta(){
    double modp = sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
    return (0.5*log( (modp+p[3])/(modp-p[3]) ) );
}

double p4vec::pT(){
    return sqrt(p[1]*p[1]+p[2]*p[2]);
}

double p4vec::ET(){
    double ET2 = (p[0] + p[3]) * (p[0] - p[3]);
    // For numerical instabilities with extreme kinematics    
    if( ET2 < 0. ){ // Should never be encountered
        cout << "p4vec::ET repaired\n";
        cout << "new ET2 = " << abs(ET2) << endl;
        return sqrt( abs(ET2) );
    }
    return sqrt( ET2 );
}

double p4vec::m2(){
    return p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
}

double p4vec::m(){
    return sqrt(p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3]);
}

double p4vec::E(){
    return p[0];
}
double p4vec::px(){
    return p[1];
}
double p4vec::py(){
    return p[2];
}
double p4vec::pz(){
    return p[3];
}
double p4vec::plus(){// (E + Pz) / sqrt(2)
        return (p[0] + p[3])/sqrt(2.);
}
double p4vec::minus(){// (E - Pz) / sqrt(2)
        return (p[0] - p[3])/sqrt(2.);
}
// Access individual elements i = 0 (E), 1 (px), 2 (py), 3 (pz) of four-momentum
double p4vec::pi(int i){
    if( i==0 or i==4 ) return p[0];
    else if( i>0 and i<4 )
        return p[i];
    else{
        std::cout << "p4vec::pi trying to access element i = " << i << std::endl;
        return 0.0;        
    }
}

p4vec p4vec::operator+(const p4vec& vec){
    p4vec p4out;
    p4out.p[1] = this->p[1] + vec.p[1];
    p4out.p[2] = this->p[2] + vec.p[2];
    p4out.p[3] = this->p[3] + vec.p[3];
    p4out.p[0] = this->p[0] + vec.p[0];
    return p4out;
}

p4vec p4vec::operator+=(const p4vec& vec){
    p4vec p4out;
    this->p[1] += vec.p[1];
    this->p[2] += vec.p[2];
    this->p[3] += vec.p[3];
    this->p[0] += vec.p[0];
    return p4out;
}

p4vec p4vec::operator-(const p4vec& vec){
    p4vec p4out;
    p4out.p[1] = this->p[1] - vec.p[1];
    p4out.p[2] = this->p[2] - vec.p[2];
    p4out.p[3] = this->p[3] - vec.p[3];
    p4out.p[0] = this->p[0] - vec.p[0];
    return p4out;
}

p4vec p4vec::operator-=(const p4vec& vec){
    p4vec p4out;
    this->p[1] -= vec.p[1];
    this->p[2] -= vec.p[2];
    this->p[3] -= vec.p[3];
    this->p[0] -= vec.p[0];
    return p4out;
}

p4vec p4vec::operator*(double x){
    p4vec p4out;
    p4out.p[1] = this->p[1] * x;
    p4out.p[2] = this->p[2] * x;
    p4out.p[3] = this->p[3] * x;
    p4out.p[0] = this->p[0] * x;
    return p4out;
}

p4vec p4vec::operator/(double x){
    p4vec p4out;
    p4out.p[1] = this->p[1] / x;
    p4out.p[2] = this->p[2] / x;
    p4out.p[3] = this->p[3] / x;
    p4out.p[0] = this->p[0] / x;
    return p4out;
}

// General kinematic functions given below

// Kallen function
double kallen(const double x, const double y, const double z){
    return pow(x,2) + pow(y,2) + pow(z,2) - 2.*x*y - 2.*x*z - 2.*y*z;
}

// Dot product of four vectors, pi.pj
double dotp4(p4vec &pi, p4vec &pj){
    return pi.E()*pj.E()-pi.px()*pj.px()-pi.py()*pj.py()-pi.pz()*pj.pz();
}

// dy_ij = y_i = y_j
double dy_ij(p4vec &pi, p4vec &pj){
    return pi.rap()-pj.rap();
}

// deta_ij = eta_i - eta_j
double deta_ij(p4vec &pi, p4vec &pj){
    return pi.eta()-pj.eta();
}

// dphi_ij = phi_i - phi_j
double dphi_ij(p4vec &pi, p4vec &pj){
    double phi1 = atan2( pi.py(), pi.px() );
    // [0,2Pi]
    if( phi1 < 0 ) phi1 += 2*M_PI;
    double phi2 = atan2( pj.py(), pj.px() );
    if( phi2 < 0 ) phi2 += 2*M_PI;
    double dphi = fabs(phi1-phi2);
    // [0,Pi]
    if( dphi > M_PI) dphi = 2*M_PI-dphi;
    return dphi;
}

// dR_ij = sqrt( dphi_ij^2 + dy_ij^2 )
double dR_ij(p4vec &pi, p4vec &pj){
    double dy2 = pow(pi.rap()-pj.rap(),2);
    double phi1 = atan2( pi.py(), pi.px() );
    // [0,2Pi]
    if( phi1 < 0 ) phi1 += 2*M_PI;
    double phi2 = atan2( pj.py(), pj.px() );
    if( phi2 < 0 ) phi2 += 2*M_PI;
    double dphi = fabs(phi1-phi2);
    // [0,Pi]
    if( dphi > M_PI) dphi = 2*M_PI-dphi;
    return sqrt(dy2+dphi*dphi);
}

// costheta_CS, i.e. in the Collins-Soper frame
double costheta_CS(p4vec &lm, p4vec &lp){
    // Construct V
    p4vec V = lm + lp;
    double mll = V.m();
    double ptll = V.pT();
    // + and - components
    double p1p = lm.E()+lm.pz();
    double p1m = lm.E()-lm.pz();
    double p2p = lp.E()+lp.pz();
    double p2m = lp.E()-lp.pz();
    // Fix sign of pz
    double sign_pz = 1.0;
    if( V.pz() < 0 ) sign_pz *= -1.0;
    //
    return sign_pz*((p1p*p2m)-(p1m*p2p))/mll/sqrt(mll*mll+ptll*ptll);
}

// phi_CS, i.e. in the Collins-Soper frame
// Also can compute costheta_CS this way
double phi_CS(p4vec &lm, p4vec &lp){
    // Construct V
    p4vec V = lm + lp;
    double QV2 = V.m2();
    double QV = sqrt(QV2);
    double PTV = V.pT();
    double PTV2 = pow(PTV,2);
    double XTV = sqrt(QV2+PTV2);
    //
    double sign_pz = 1.0;
    if( V.pz() < 0 ) sign_pz *= -1.0;
    // 
    double pl_x = 
        +PTV/QV/XTV * ( -V.E()*lm.E() +V.pz()*lm.pz() )
        +XTV/QV/PTV * ( +V.px()*lm.px() +V.py()*lm.py() );
    double pl_y = sign_pz/PTV * ( -V.py()*lm.px() +V.px()*lm.py() );
    // double pl_z = sign_pz/XTV * (
    //      -V.pz()*lm.E()
    //      +V.E()*lm.pz() );
    // Now evaluate cos_theta_cs, phi_cs
    // double cos_th_cs = pl_z / sqrt( pow(pl_x,2) + pow(pl_y,2) + pow(pl_z,2) );
    // double th_cs = atan2( sqrt( pow(pl_x,2) + pow(pl_y,2)), pl_z );
    double phi_cs = atan2( pl_y, pl_x );
    if( phi_cs < 0 ) phi_cs += 2.*M_PI;
    return phi_cs;
}

// phistar variable
double phistar(p4vec &lm, p4vec &lp){
    double dphi = dphi_ij( lm, lp );
    double deta = lm.rap() - lp.rap(); // eta_lm - eta_lp, massless leptons
    double cos_th = tanh(deta/2.0); 
    double sin_th = sqrt( (1.+cos_th)*(1.-cos_th) ); // independent of sign of deta
    return tan( (M_PI-dphi)/2.0 ) * sin_th;
}

// mT variable, Charged Current (see 1701.07420)
double MT(p4vec &nu, p4vec &l){
    double dphi = dphi_ij( nu, l );
    return sqrt(2. * nu.pT() * l.pT() * ( 1. - cos(dphi) ));
}

// Could also pass an array by reference to a void function?
// Function for defining energies and momenta in a particle decay of the form mij -> mi + mj in the mij CoM frame
// Forms the basis for generating particle momenta using sequences of 2-body phase spaces
std::array<double,3> CoM_1to2_kinematics(double mij, double mi, double mj){
    // Define dimensionless variables
    double xi = max(0.,pow(mi/mij,2));
    double xj = max(0.,pow(mj/mij,2));
    // 3d array to contain Ei, Ej, |pi|
    std::array<double,3> output = {0.};
    // Ei = mij / 2. ( 1 + xi - xj )
    output[0] = mij / 2. * ( 1. + xi - xj );
    // Ej = mij / 2. ( 1 - xi + xj )
    output[1] = mij / 2. * ( 1. - xi + xj );
    // |pi| = |pj| = mij/2. * Kallen^{1/2}(1.,xi,xj)
    double vel = sqrt( kallen(1.,xi,xj) );    
    if( xi == 0.0 ){
        vel = ( 1. - xj );
    }
    if( xj == 0.0 ){
        vel = ( 1. - xi );
    }
    output[2] = mij / 2. * vel;    
    // output[2] = mij / 2. * sqrt( kallen(1.,xi,xj) );
    return output; 
}

std::array<double,3> CoM_1to2_kinematics_sq(double msq_ij, double msq_i, double msq_j){
    // Define dimensionless variables
    double xi = max(0.,msq_i/msq_ij);
    double xj = max(0.,msq_j/msq_ij);
    // 3d array to contain Ei, Ej, |pi|
    std::array<double,3> output = {0.};
    double mij = sqrt(msq_ij);
    // Ei = mij / 2. ( 1 + xi - xj )
    output[0] = mij / 2. * ( 1. + xi - xj );
    // Ej = mij / 2. ( 1 - xi + xj )
    output[1] = mij / 2. * ( 1. - xi + xj );
    // |pi| = |pj| = mij/2. * Kallen^{1/2}(1.,xi,xj)
    // sqrt( 1 -)
    double vel = sqrt( kallen(1.,xi,xj) );
    if( xi == 0.0 ){
        vel = ( 1. - xj );
    }
    if( xj == 0.0 ){
        vel = ( 1. - xi );
    }
    output[2] = mij / 2. * vel;
    // Perform expansion?
    return output; 
}



////////////////////////////////////////////////////////////
// Information concerning the kinematic information class //
////////////////////////////////////////////////////////////


// The constructor for kinematic data
// Create the variable length objects accordingly
KinematicData::KinematicData(int npar) : npar(npar), pset(new p4vec[npar]), pdgs(new int[npar]) {
}

// Kinematic functions
// sij (currently didn't implement the SH functions)

// pi.pj
double KinematicData::pij(int i, int j) {
	if( i > npar or j > npar ) {
		cout << "KinematicData::pij, Error, p(" << max(i, j) << ") not initalised" << endl;
		cout << npar << endl;
		abort();
	}
	return dotp4(pset[i - 1], pset[j - 1]);
}

// pijext(int, p4vec)
double KinematicData::pijext(int i, p4vec j){
    if( i > npar ) {
        cout << "KinematicData::pijext, Error, p(" << i << ") not initalised" << endl;
        cout << npar << endl;
        abort();     
    }
    return dotp4(pset[i-1], j);
}

double KinematicData::sij(int i, int j) {
	if( i > npar or j > npar or min(i, j) < 1 ) {
		cout << "KinematicData::Za, p_i(" << i << ") not initalised" << endl;
		cout << "KinematicData::Za, p_j(" << j << ") not initalised" << endl;
		abort();
	}
    return 2. * pij(i,j);
}

// Return function to access the privately stored momentum information
p4vec KinematicData::p(int i) {
	if( i > npar ) {
		cout << "KinematicData::p, Error, p(" << i << ") not initalised" << endl;
		abort();
	}
	return pset[i - 1];
}

// Access Lab momenta (for cuts)
p4vec KinematicData::p_Lab(int i, double b_x, double b_y, double b_z) {
    if( i > npar ) {
        cout << "KinematicData::p, Error, p(" << i << ") not initalised" << endl;
        abort();
    }
    // Perform the boost to the Lab frame
    p4vec p_lab = pset[i - 1]; p_lab.boost(b_x,b_y,b_z);
    return p_lab;
}

// Function to activate storage of n random doubles [0,1]
void KinematicData::activate_ran( int nran ) {
    // Check if it was already activated
    if( r_i ){
        cout << "r_i has previously been activated" << endl;
    }

    // Otherwise create the nran length array of double
    r_i.reset( new double[nran] );
    // Also adjust the n_ran variable
    n_ran = nran;
    return;
}

// Access random variables (if set)
double KinematicData::r(int i) {
    if( i > n_ran ) {
        cout << "KinematicData::r, Error, r(" << i << ") not initalised" << endl;
        abort();
    }
    return r_i[i - 1];
}



// General wrapper for phase-space generation
KinematicData Generate_Phase_Space( const cubareal rand[], const int npar_out, const double masses[], const int n_random, const std::string environment ){

    if( npar_out > 3 ){
        cerr << "Generate_Phase_Space: implement n-body final-state phase space for n = " << npar_out << endl;
        abort();
    }

    // Double check integration dimenions matches that required
    int integration_dimensions(0);
    // Number of phase-space dimensions
    int ps_dimensions = 1 + (npar_out-2) * 3;
    integration_dimensions += ps_dimensions;

    // The random numbers stored in r_i container
    // Those are provided via n_random as additional integration parameters

    // Random numbers requested (e.g. for collinear convolutions at V, RV, VV levels)
    if( (integration_dimensions+n_random) != cuba_dimensions ){
        cerr << "Generate_Phase_Space: integration dimensions = " << integration_dimensions << endl;
        cout << "n random = " << n_random << endl;
        cerr << "cuba_dimensions = " << cuba_dimensions << endl;
        cerr << "PS_dimensions = " << ps_dimensions << endl;
        abort();
    }

    int mij_integration = 2; // Log

    // 2 to 2 Massive
    if( npar_out == 2 ){
        KinematicData Kin = Gen_2to2_Massive(rand,mij_integration,masses);
        // Store random variables
        if( n_random > 0 ){
            Kin.activate_ran(n_random);
            for( int i(0); i < n_random; i++)
                Kin.set_ri(rand[i+ps_dimensions],i+1);                
        }
        return Kin;
    }
    else if( npar_out == 3 ){
        KinematicData Kin = Gen_2to3_Massive(rand,mij_integration,masses);
        // Store random variables
        if( n_random > 0 ){
            Kin.activate_ran(n_random);
            for( int i(0); i < n_random; i++)
                Kin.set_ri(rand[i+ps_dimensions],i+1);                
        }
        return Kin;       
    }
    else if( npar_out == 4 ){
        cerr << "PS: 2to4 currently not implemented yet\n";
        abort();
    }


    cerr << "Generate_Phase_Space: should not reach end of this function\n";
    abort();
    return KinematicData(npar_out+2);
}




// For 2 -> n scattering, use phase-space factorisation
// Always writing dphi_n(shat;p3,..,pn) in factorised form, a product of many two-body phase spaces
// dphi_n(shat;p3,..,pn) -> dphi_{n-1}(shat;pX,...,pn) dmX^2/(2Pi) dphi_2(mx2;p3,p4) ...
// Note there I can still integrate over the trivial dphi only of the 'first two-body phase space' piece


// A general version of the 2to2 scattering kinematics for arbitrary final-state particle masses, at fixed Ecms2
KinematicData Gen_2to2_Massive( const cubareal rand[], const int ps_opt, const double masses[]){
    double jacob(0.), costh(0.);
    double m3 = max(masses[0],0.0);
    double m4 = max(masses[1],0.0);
    // (p1 + p2)^2 = Ecms2 = s12
    double s12 = Ecms2;

    // Just a single variable to consider: dcostheta
    // Perform the mapping of int_{-1}^{1} dcosth onto the hypercube [0,1]
    costh = -1.0 + 2.0 * rand[0];
    // Factor of two for the dcostheta jacobian/mapping
    jacob = 2.;
    // Derive sin(theta)
    double sinth = sqrt( (1.0+costh)*(1.0-costh) );
    // Derive the particle kinematics
    // Declare p1, p2 directly in Ecms2 frame (E,x,y,z)
    p4vec p1 ( Ecms/2.0, 0.0, 0.0, Ecms/2.0 );
    p4vec p2 ( Ecms/2.0, 0.0, 0.0,-Ecms/2.0 );

    // Information for p3, p4 in CoM frame
    double m12 = sqrt(s12);
    std::array<double,3> CoM_Kins = CoM_1to2_kinematics( m12, m3, m4 );
    double E3 = CoM_Kins[0];
    double E4 = CoM_Kins[1];
    double mod_p3 = CoM_Kins[2];    
    // This now defines the particle four momenta in partonic CoM frame
    p4vec p3 ( E3, 0, mod_p3*sinth, mod_p3*costh );
    p4vec p4 ( E4, 0,-mod_p3*sinth,-mod_p3*costh );

    // Generate 4 particle Kinematic Data structure
    KinematicData KinOut(4);
    // Add momenta and x1, x2 to KinematicData structure
    KinOut.set_pi( p1, 1 );
    KinOut.set_pi( p2, 2 );
    KinOut.set_pi( p3, 3 );
    KinOut.set_pi( p4, 4 );

    double scale = KinematicScales( KinOut, scale_opt );
    KinOut.set_muf( scale * muf_var );
    KinOut.set_mur( scale * mur_var );

    // The differential cross-section
    // A single two-body phase-space factor, and include analytically the 2Pi from dphi integration
    // Have used that lambda^{1/2}(1,s_3,s_4) = 2 |p3| / sqrt(s12) in partonic CoM frame, m12=sqrt(s12)
    double ps_factor = (2.*M_PI) / ( 2. * pow(4.*M_PI,2) ) * ( 2. * mod_p3 / m12 );
    // Weight it the phase space factor, the jacobian mappings for each integration, and jacobians for any change of variables
    double weight = jacob * ps_factor;
    KinOut.set_weight( weight );
    KinOut.set_flux( 2.0 * s12 );

    // Apply analysis cuts (check if phase-space point passes cuts or not)
    apply_cuts( KinOut );
    return KinOut;
}

KinematicData Gen_2to3_Massive( const cubareal rand[], const int ps_opt, const double masses[]){
    // 6 integration variables variables and jacobian factor
    double jacob(1.), msq45(0.), cth_12(0.), cth_45(0.), phi_45(0.);
    // Mass information
    double m3 = masses[0];
    double m4 = masses[1];
    double m5 = masses[2];
    // ######################################################## //
    // ######################################################## //
    // Start the generation of the above differential variables //
    if( Ecms2 < pow(m3+m4+m5,2) ){
        cout << "Gen_2to3_Massive: not enough energy to produce massive particles\n";
        abort();
    }
    ///////////////////////////
    // [1] m45^2 integration //
    ///////////////////////////
    // 1) Generate m45^2: min, max
    double msq45_min = pow(m4+m5,2)+tech_cut;//max(tech_cut, pow(m4+m5,2)+tech_cut );
    // If m3 and m45 are produced at rest and s12 = Ecms2, can maximise m45 for the condition
    double m45_max = Ecms - (m3+tech_cut);
    double msq45_max = pow(m45_max,2);
    // Consider different sampling options for m45^2
    // ps_opt = 1, Linear sampling in m45^2
    if( ps_opt == 1 ){
        msq45 = msq45_min + ( msq45_max - msq45_min ) * rand[0];
        jacob *= ( msq45_max - msq45_min );       
    }
    // ps_opt = 2, Log sampling in m45^2
    if( ps_opt == 2 ){
        double lmsq45 = log(msq45_min) + ( log(msq45_max) - log(msq45_min) ) * rand[0];
        msq45 = exp(lmsq45);
        jacob *= msq45 * (log(msq45_max) - log(msq45_min));
    }
    // Integrate over m45 in various ways (optimised for gauge-boson resonance etc.)
    if( ps_opt == 0 or ps_opt > 2 ){
        cout << "Gen_2to3_Massive: ps_option " << ps_opt << endl;
        cout << "Gen_2to3_Massive: not implemented yet " << endl;
        abort();
    }
    // (p4+p5)^2 = m45^2 = m4^2 + m5^2 + 2 p4.p5, with 2 p4.p5 = s45
    // double s45 = msq45 - pow(m4,2) - pow(m5,2);
    double m45 = sqrt(msq45);

    ////////////////////////////////////////////
    // Energies and |p_out| in the (12) frame //
    ////////////////////////////////////////////
    double s12 = Ecms2; // (p1+p2)^2 = s_12 = Ecms2
    double m12 = sqrt(s12); // special case since m1^2 = m2^2 = 0., or s12 = m2_12
    // Use these variables to derive Energies and |p| of particle 3 and 45
    // Function which returns Ei, Ej, |pi| = |pj| in CoM frame (mij;mi,mj)
    // std::array<double,3> CoM_Kins_12 = CoM_1to2_kinematics( m12, m3, m45 );
    std::array<double,3> CoM_Kins_12 = CoM_1to2_kinematics_sq( s12, pow(m3,2), msq45 );    
    double E3_12  = CoM_Kins_12[0];
    double E45_12 = CoM_Kins_12[1];
    double p3_12  = CoM_Kins_12[2];

    /////////////////////////////////////////
    // [1] costh_12 integration (12-frame) //
    /////////////////////////////////////////
    // this is scattering angle between p1 and p3 in CoM frame (p1,p2 aligned on z)
    double cth_12_min = -1. + tech_cut; // tech_cut is a tiny quantity (1e-10) for numerical stability (for emission of unresolved particles)
    double cth_12_max = +1. - tech_cut;
    cth_12 = cth_12_min + (cth_12_max - cth_12_min) * rand[1];
    jacob *= (cth_12_max - cth_12_min) * ( 2. * M_PI ); // dOmega_12 ~ cth_12 dphi_12, [where dphi_12 integration just dropped, i.e. dPhi_12 -> 2 Pi]
    //////////////////////////////////////////////////////////////////////
    // [1] Alternatively, consider t_hat integration: t_hat = (p1-p3)^2 //
    //////////////////////////////////////////////////////////////////////    
    // t_hat = (p1-p3)^2 = m3^2 - 2 (E1 E3 - |p_in| |p_out| cos_theta ), with |p_in| = sqrt(s12)/2 = m12/2
    // t_hat_min [cos_theta -> -1]
    bool integration_t_chan = false; bool linear = false;
    if( integration_t_chan ){
        double t_hat_max = min( pow(m3,2) - m12 * ( E3_12 - p3_12 ), -m12*tech_cut );
        double t_hat_min = pow(m3,2) - m12 * ( E3_12 + p3_12 );
        double t_hat(0.);
        // Linear sampling
        if( linear ){
            t_hat = t_hat_min + (t_hat_max - t_hat_min) * rand[1];
            jacob *= (t_hat_max - t_hat_min) / (cth_12_max - cth_12_min) / ( m12 * p3_12);
        }
        // Log sampling
        if( !linear ){
            // How about integration in log[-that]
            double ml_t_hat = log(-t_hat_min) + (log(-t_hat_max) - log(-t_hat_min)) * rand[1];
            t_hat = - exp(ml_t_hat);
            jacob *= t_hat * (log(-t_hat_max) - log(-t_hat_min)) / (cth_12_max - cth_12_min) / ( m12 * p3_12 );
        }
        // dervie cos theta = ( t_hat - m3^2 )/|p3|
        cth_12 = ( m12 * E3_12 + t_hat - pow(m3,2) ) / ( m12 * p3_12 );
    }
    // Whatever approach, derive s^2_theta
    double sth_12 = sqrt( (1.-cth_12)*(1.+cth_12) );    

    /////////////////////////////////////////
    // [2] costh_45 integration (45-frame) //
    // [3] phi_45 integration (45-frame)   //    
    /////////////////////////////////////////
    cth_45 = -1.0 + 2.0 * rand[2];
    phi_45 = 2.0 * M_PI * rand[3];
    jacob *= 2.0 * ( 2. * M_PI ); // dOmega_45 ~ costh_45 * dphi_45, jacobian
    // derive sinth_45 in this frame
    double sth_45 = sqrt( (1.0+cth_45)*(1.-cth_45) );    

    // All differential variables have been generated //
    // ############################################## //
    // ############################################## //   

    ///////////////////////////////////////////////////////
    // Now generate the momenta of the phase-space point // 
    p4vec p1( m12/2.0, 0, 0,+m12/2.0 );
    p4vec p2( m12/2.0, 0, 0,-m12/2.0 );
    // Define p3, p45 in the (12) frame
    p4vec p3(   E3_12, 0.0, +p3_12*sth_12, +p3_12*cth_12);
    p4vec p45( E45_12, 0.0, -p3_12*sth_12, -p3_12*cth_12);
    // For extra stability could consider:
    p45 = p1 + p2 - p3;
    /////////////////////////////////////////
    // Momenta of p4, p5 in the (45) frame //
    /////////////////////////////////////////
    std::array<double,3> CoM_Kins_45 = CoM_1to2_kinematics( m45, m4, m5 );
    double E4_45  = CoM_Kins_45[0];
    double E5_45  = CoM_Kins_45[1];
    double p4_45  = CoM_Kins_45[2];
    // Define p4,p5 in the (45) frame
    p4vec p4( E4_45, +p4_45*sth_45*cos(phi_45), +p4_45*sth_45*sin(phi_45), +p4_45*cth_45 );
    p4vec p5( E5_45, -p4_45*sth_45*cos(phi_45), -p4_45*sth_45*sin(phi_45), -p4_45*cth_45 );
    // Boost these two vectors from (45) -> (12) frame, i.e. the inverse boost that puts p45 to rest
    p4.boost( p45.px()/p45.E(), p45.py()/p45.E(), p45.pz()/p45.E() );
    p5.boost( p45.px()/p45.E(), p45.py()/p45.E(), p45.pz()/p45.E() );
    // For further extra stability, could consider:
    p5 = p1 + p2 - p3 - p4;

    // Generate 5 particle Kinematic Data structure
    KinematicData KinOut(5);
    // Add partonic CoM momenta
    KinOut.set_pi( p1, 1 );
    KinOut.set_pi( p2, 2 );
    KinOut.set_pi( p3, 3 );
    KinOut.set_pi( p4, 4 );
    KinOut.set_pi( p5, 5 );
    // Also default pdg information to 0
    KinOut.set_pdg_i( 1, 0 );
    KinOut.set_pdg_i( 2, 0 );
    KinOut.set_pdg_i( 3, 0 );
    KinOut.set_pdg_i( 4, 0 );
    KinOut.set_pdg_i( 5, 0 );    

    double scale = KinematicScales( KinOut, scale_opt );
    KinOut.set_muf( scale * muf_var );
    KinOut.set_mur( scale * mur_var );
    // PS factors
    double ps2body_12 = 1. / ( 2. * pow(4.*M_PI,2) ) * (2. * p3_12 / m12 );
    double ps2body_45 = 1. / ( 2. * pow(4.*M_PI,2) ) * (2. * p4_45 / m45 );    
    double ps_factor = ps2body_12 * ps2body_45 / ( 2. * M_PI );
    // Weight it the phase space factor, the jacobian mappings for each integration, and jacobians for any change of variables
    double weight = jacob * ps_factor;
    KinOut.set_weight( weight );
    KinOut.set_flux( 2.0 * s12 );
  
    // Apply analysis cuts (check if phase-space point passes cuts or not)
    apply_cuts( KinOut );

    // Debug phase-space point
    if( scale < 0 or isnan(scale) ){
        debug_PS_point( KinOut );
        KinOut.set_cuts(false);
    }

    // Some de-bugging for phase-space points
    if( KinOut.get_cuts() ){
        double ET_sum = p3.ET() + p4.ET() + p5.ET();
        if( isnan(ET_sum) or ET_sum < 0 ){
            debug_PS_point( KinOut );
        }
        // Perform sanity check of momentum conservation for a well defined PS point
        p4vec mom_cons = p1 + p2 - p3 - p4 - p5;
        // Check of momentum conservation
        if( abs(mom_cons.E()) > 5e-6 ){
        // if( abs(p3.m2()) > 1e-8 ){
            cout << setprecision(15) << "Gen_2to3_Massive: Momentum conservation not satisfied so well" << endl;
            cout << "(p1+p2-p3-p4-p5) =";
            mom_cons.print();
            cout << "" << endl;
            debug_PS_point( KinOut );
        }
        // Could consider also defining p5 = p1 + p2 - p3 - p4 (for numerical stability?)
    }
    return KinOut;
}



// Function which determines kinematic scales from a general phase-space
double KinematicScales( KinematicData &Kin, int option ){

    // Option -1: a fixed-scale (the parameter mu0)
    if( option == -1 ){
        return mu0;
    }            
    // Option  0: sqrt(shat) = s12
    else if( option == 0 ){
        return 2. * Kin.pij(1,2); // s12
    }
    // Option  1: ET_(n,n-1) particle, transverse energy of the nth and nth-1 particle
    else if( option == 1 ){
        p4vec p_combined = Kin.p( Kin.length() ) + Kin.p( Kin.length()-1 );
        return p_combined.ET();
    }
    // Option  2: \sum ET_i / 2, some of all final-state particle transverse energies over 2
    else if( option == 2 ){
        double ET_sum(0.);
        for( int i = 3; i <= Kin.length(); i++ ){
            ET_sum += Kin.p(i).ET();
        }
        return ET_sum / 2.;
    }
    // Option  3: \sum ET_i / 4, some of all final-state particle transverse energies over 2
    else if( option == 3 ){
        double ET_sum(0.);
        for( int i = 3; i <= Kin.length(); i++ ){
            ET_sum += Kin.p(i).ET();
        }
        return ET_sum / 4.;
    }
    // Option  4: m_(n,n-1), mass of the combined system of the nth and nth-1 particle
    else if( option == 4 ){
        p4vec p_combined = Kin.p( Kin.length() ) + Kin.p( Kin.length()-1 );
        return p_combined.m();
    }
    // Option  5: ( ET(n) + ET(n-1) )/2
    else if( option == 5 ){
        p4vec p_combined = Kin.p( Kin.length() ) + Kin.p( Kin.length()-1 );
        return ( Kin.p( Kin.length() ).ET() + Kin.p( Kin.length()-1 ).ET() ) / 2.;
    }    
    else{
        cout << "KinematicScales: unspecifed option " << option << endl;
        abort();
    }
    cout << "KinematicScales: unspecifed option " << option << endl;
    abort();
    return 0.;
}


// Print some useful information
void debug_PS_point( KinematicData &Kin ){
    cout << "Debugging PS point\n";
    cout << setprecision(15);
    cout << "n particles = " << Kin.length() << endl;
    for( int i = 1; i <= Kin.length(); i++ ){
        cout << "i ";
        Kin.p(i).print();
    }

    // Compute the ET of the particles
    for( int i = 1; i <= Kin.length(); i++ )
        cout << "ET = " << Kin.p(i).ET() << endl;

    cout << "scale = " << Kin.mur() << endl;
    return;
}
