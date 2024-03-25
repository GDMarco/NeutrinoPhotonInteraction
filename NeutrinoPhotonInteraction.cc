#include "NeutrinoPhotonInteraction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

//const std::string myPropertyName = "counter";

// The parent's constructor need to be called on initialization!
NeutrinoPhotonInteraction::NeutrinoPhotonInteraction(ref_ptr<PhotonField> photonField, bool haveSecondaries,  double limit) : Module() { //double thinning,
    setPhotonField(photonField);
    setHaveSecondaries(haveSecondaries);
    setLimit(limit);
    //setThinning(thinning);
    
}

void NeutrinoPhotonInteraction::setPhotonField(ref_ptr<PhotonField> photonField) {
    this->photonField = photonField;
    std::string fname = photonField->getFieldName();
    setDescription("NeutrinoPhotonInteraction::Module" + fname);
    initRate(getDataPath("NeutrinoPhotonInteraction/rate_" + fname + ".txt"));
}

void NeutrinoPhotonInteraction::setHaveSecondaries(bool haveSecondaries) {
    this->haveSecondaries = haveSecondaries;
}

void NeutrinoPhotonInteraction::setLimit(double limit) {
    this->limit = limit;
}
/**
void NeutrinoPhotonInteraction::setThinning(double thinning) {
    this->thinning = thinning;
}
 */
void NeutrinoPhotonInteraction::initRate(std::string filename) {
    std::ifstream infile(filename.c_str());
    
    if (!infile.good())
        throw std::runtime_error("NeutrinoPhotonInteraction: could not open file" + filename);
    
    tabEnergy.clear();
    tabRate.clear();
    
    while (infile.good()) {
        if (infile.peek() != '#') {
            double a, b;
            infile >> a >> b;
            if (infile) {
                tabEnergy.push_back(pow(10, a) * eV);
                tabRate.push_back(b / Mpc);
            }
        }
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    }
    infile.close();
}

void NeutrinoPhotonInteraction::performInteraction(Candidate *candidate) const {
    
    candidate->setActive(false);

    if (not haveSecondaries)
        return;

    double mass_W = 1.; //value in kg //W(+) particle ID ==24
    double w = 1.;
    // Use assumption of Secke 98
    // W boson produced on shell
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    double Ee = (E - mass_W * c_squared);

    Random &random = Random::instance();
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    if (haveSecondaries)
        candidate->addSecondary(11, Ee / (1 + z), pos, w, interactionTag);
}

void NeutrinoPhotonInteraction::process(Candidate *candidate) const
{
    // To enable parallelization, the modules have to be stateless - the
    // process method should thus not modify internal variables!
    std::cout << "NeutrinoPhotonInteraction::Module::process() called\n";
    if (candidate->current.getId() != 12)
        return;
   
    // scale the electron energy instead of background photons
        double z = candidate->getRedshift();
        double E = (1 + z) * candidate->current.getEnergy();

        // check if in tabulated energy range
        if (E < tabEnergy.front() or (E > tabEnergy.back()))
            return;

        // interaction rate
        double rate = interpolate(E, tabEnergy, tabRate);
        rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);

        // check for interaction
        Random &random = Random::instance();
        double randDistance = -log(random.rand()) / rate;
        double step = candidate->getCurrentStep();
        if (step < randDistance) {
            candidate->limitNextStep(limit / rate);
            return;
        } else { // after performing interaction photon ceases to exist (hence return)
            performInteraction(candidate);
            return;
        }
}

void NeutrinoPhotonInteraction::setInteractionTag(std::string tag) {
    interactionTag = tag;
}

std::string NeutrinoPhotonInteraction::getInteractionTag() const {
    return interactionTag;
}

} // end namespace crpropa
