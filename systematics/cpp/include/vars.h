/**
 * @file vars.h
 * @brief This file contains the reconstructed variables and systematic
 * parameters for the systematics code.
 * @author justin.mueller@colostate.edu
*/

#ifndef VARS_H
#define VARS_H

#include "types.h"

struct RecoVar
{
    std::string name;
    uint32_t nbins;
    double xmin;
    double xmax;
};

std::vector<RecoVar> reco_vars = {
    {"reco_muon_ke", 25, 0, 2000},
    {"reco_proton_ke", 25, 0, 600},
    {"reco_visible_energy", 25, 0, 3000},
    {"reco_muon_pt", 25, 0, 750},
    {"reco_proton_pt", 25, 0, 750},
    {"reco_muon_polar_angle", 25, 0, 3.14159},
    {"reco_muon_azimuthal_angle", 25, 0, 3.14159},
    {"reco_opening_angle", 25, 0, 3.14159},
    {"reco_delta_pT", 25, 0, 1200},
    {"reco_delta_phiT", 25, 0, 3.14159},
    {"reco_delta_alphaT", 25, 0, 3.14159},
    {"muon_softmax", 25, 0, 1},
    {"proton_softmax", 25, 0.8, 1}};

systs_t systs({
    {"GENIEReWeight_ICARUS_v2_multisim_ZExpAVariationResponse", 52},
    {"GENIEReWeight_ICARUS_v2_multisim_RPA_CCQE", 53},
    {"GENIEReWeight_ICARUS_v2_multisim_CoulombCCQE", 54},
    {"GENIEReWeight_ICARUS_v2_multisim_NormCCMEC", 55},
    {"GENIEReWeight_ICARUS_v2_multisim_NormNCMEC", 56},
    {"GENIEReWeight_ICARUS_v2_multisim_NCELVariationResponse", 57},
    {"GENIEReWeight_ICARUS_v2_multisim_CCRESVariationResponse", 58},
    {"GENIEReWeight_ICARUS_v2_multisim_NCRESVariationResponse", 59},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvpCC1pi", 60},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvpCC2pi", 61},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvpNC1pi", 62},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvpNC2pi", 63},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvnCC1pi", 64},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvnCC2pi", 65},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvnNC1pi", 66},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvnNC2pi", 67},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvbarpCC1pi", 68},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvbarpCC2pi", 69},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvbarpNC1pi", 70},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvbarpNC2pi", 71},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvbarnCC1pi", 72},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvbarnCC2pi", 73},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvbarnNC1pi", 74},
    {"GENIEReWeight_ICARUS_v2_multisim_NonRESBGvbarnNC2pi", 75},
    {"GENIEReWeight_ICARUS_v2_multisim_RDecBR1gamma", 76},
    {"GENIEReWeight_ICARUS_v2_multisim_RDecBR1eta", 77},
    {"GENIEReWeight_ICARUS_v2_multisim_COHVariationResponse", 78},
    {"GENIEReWeight_ICARUS_v2_multisim_DISBYVariationResponse", 79},
    {"GENIEReWeight_ICARUS_v2_multisim_FSI_pi_VariationResponse", 80},
    {"GENIEReWeight_ICARUS_v2_multisim_FSI_N_VariationResponse", 81},
    {"expskin_Flux", 115},
    {"horncurrent_Flux", 116},
    {"kminus_Flux", 117},
    {"kplus_Flux", 118},
    {"kzero_Flux", 119},
    {"nucleoninexsec_Flux", 120},
    {"nucleonqexsec_Flux", 121},
    {"nucleontotxsec_Flux", 122},
    {"piminus_Flux", 123},
    {"pioninexsec_Flux", 124},
    {"pionqexsec_Flux", 125},
    {"piontotxsec_Flux", 126},
    {"piplus_Flux", 127}});

#endif