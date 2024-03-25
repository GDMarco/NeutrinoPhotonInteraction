/// Example plugin for CRPropa.
///
/// Please consider sharing the awesome plugin with you fellow researchers by
/// creating a eperate repository for your project. We maintain a list of
/// plugins to CRPropa on our webpage and are happy to add a link to your
/// project, just send us: (name of the plugin, short description, url)
#include <fstream>
#include <cmath>

#include <crpropa/Module.h>
#include "crpropa/PhotonBackground.h"

namespace crpropa {
/// A custom C++ module for Neutrino-Photon Interaction in astrophysical scenarios. On-shell production of the W boson.
///
class NeutrinoPhotonInteraction : public Module {
private:
    ref_ptr<PhotonField> photonField;
    bool haveSecondaries;
    double limit;
    //double thinning;
    
    std::string interactionTag = "NGI";
    
    std::vector<double> tabEnergy; //!< neutrino energy in [J]
    std::vector<double> tabRate; //!< interaction rate in [1/m]
    
    std::vector<double> tabE; //!< neutrino energy in [J]
    std::vector<double> tabs; //!< s_kin = s - m^2 in [J**2**]
    std::vector<std::vector<double>> tabCDF; //!< cumulative interaction rate
public:
    /// The parent's constructor need to be called on initialization!
    NeutrinoPhotonInteraction(ref_ptr<PhotonField> photonField, bool haveSecondaries = false, double limit = 0.1); //double thinning = 0,
    
    // set the target photon field
    void setPhotonField(ref_ptr<PhotonField> photonField);

    // decide if secondary electrons are added to the simulation
    void setHaveSecondaries(bool haveSecondaries);

    /** Limit the propagation step to a fraction of the mean free path
     * @param limit fraction of the mean free path
     */
    void setLimit(double limit);

    
    /** Apply thinning with a given thinning factor
     * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
     */
    /**
    void setThinning(double thinning);
    */
    
    /** set a custom interaction tag to trace back this interaction
     * @param tag string that will be added to the candidate and output
     */
    void setInteractionTag(std::string tag);
    std::string getInteractionTag() const;
    
    void initRate(std::string filename);
    
    void process(Candidate *candidate) const;
    void performInteraction(Candidate *candidate) const;
};

} // end namespace crpropa


