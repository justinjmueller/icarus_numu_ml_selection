/**
 * @file analysis.h
 * @brief Header file for definitions of complete variables with relevant cuts.
 * @author justin.mueller@colostate.edu
*/

#include "definitions.h"
#include "cuts.h"
#include "variables.h"
#include "numu_variables.h"

using namespace ana;

namespace ana
{
    // "Simple" variables that attach to a single interaction attribute.
    VARDLP_RECO(kCountParticles,vars::count_particles,cuts::no_cut);
    VARDLP_RECO(kCountPrimaries,vars::count_primaries,cuts::no_cut);

    VARDLP_TRUE(kCountParticlesTruth,vars::count_particles,cuts::no_cut);
    VARDLP_TRUE(kCountPrimariesTruth,vars::count_primaries,cuts::no_cut);

    VARDLP_RECO(kVisibleEnergy_1mu1p, vars::visible_energy, cuts::all_1mu1p_cut);
    VARDLP_RECO(kVisibleEnergy_1muNp, vars::visible_energy, cuts::all_1muNp_cut);
    VARDLP_RECO(kVisibleEnergy_1muX, vars::visible_energy, cuts::all_1muX_cut);

    // Define "variables" for binning interactions by some categorical class.
    DEFINECAT();

    // Define variables that are broadcasted across each selection cut.
    TCATVAR(kCountTTP,count);
    RCATVAR(kCountPTT,count);
    TCATVAR(kVisibleEnergyTTP,visible_energy);
    RCATVAR(kVisibleEnergyPTT,visible_energy);
    RCATVAR(kFlashTimePTT,flash_time);
    
    // Variables for confusion matrix.
    PVARDLP_TRUE(kPrimaryTruth,vars::primary,cuts::no_cut,cuts::matched);
    PVARDLP_TRUE(kPIDTruth,vars::pid,cuts::no_cut,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth,vars::primary_pid,cuts::no_cut,cuts::matched);
    PVAR_TTP(kPrimary,vars::primary,cuts::no_cut,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID,vars::pid,cuts::no_cut,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID,vars::primary_pid,cuts::no_cut,cuts::no_cut,cuts::no_cut);

    // Variables for confusion matrix (well-reconstructed).
    PVARDLP_TRUE(kPrimaryWellRecoTruth,vars::primary,cuts::wellreco,cuts::matched);
    PVARDLP_TRUE(kPIDWellRecoTruth,vars::pid,cuts::wellreco,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDWellRecoTruth,vars::primary_pid,cuts::wellreco,cuts::matched);
    PVAR_TTP(kPrimaryWellReco,vars::primary,cuts::wellreco,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPIDWellReco,vars::pid,cuts::wellreco,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPIDWellReco,vars::primary_pid,cuts::wellreco,cuts::no_cut,cuts::no_cut);

    // Variables for neutrino-only confusion matrix.
    PVARDLP_TRUE(kPrimaryTruth_Neutrino,vars::primary,cuts::neutrino,cuts::matched);
    PVARDLP_TRUE(kPIDTruth_Neutrino,vars::pid,cuts::neutrino,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth_Neutrino,vars::primary_pid,cuts::neutrino,cuts::matched);
    PVAR_TTP(kPrimary_Neutrino,vars::primary,cuts::neutrino,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID_Neutrino,vars::pid,cuts::neutrino,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID_Neutrino,vars::primary_pid,cuts::neutrino,cuts::no_cut,cuts::no_cut);

    // Variables for neutrino-only confusion matrix (well-reconstructed).
    PVARDLP_TRUE(kPrimaryWellRecoTruth_Neutrino,vars::primary,cuts::wellreco_neutrino,cuts::matched);
    PVARDLP_TRUE(kPIDWellRecoTruth_Neutrino,vars::pid,cuts::wellreco_neutrino,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDWellRecoTruth_Neutrino,vars::primary_pid,cuts::wellreco_neutrino,cuts::matched);
    PVAR_TTP(kPrimaryWellReco_Neutrino,vars::primary,cuts::wellreco_neutrino,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPIDWellReco_Neutrino,vars::pid,cuts::wellreco_neutrino,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPIDWellReco_Neutrino,vars::primary_pid,cuts::wellreco_neutrino,cuts::no_cut,cuts::no_cut);

    // Variables for cosmic-only confusion matrix.
    PVARDLP_TRUE(kPrimaryTruth_Cosmic,vars::primary,cuts::cosmic,cuts::matched);
    PVARDLP_TRUE(kPIDTruth_Cosmic,vars::pid,cuts::cosmic,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth_Cosmic,vars::primary_pid,cuts::cosmic,cuts::matched);
    PVAR_TTP(kPrimary_Cosmic,vars::primary,cuts::cosmic,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID_Cosmic,vars::pid,cuts::cosmic,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID_Cosmic,vars::primary_pid,cuts::cosmic,cuts::no_cut,cuts::no_cut);

    // Variables for 2D "true vs. reco" style plots.
    PVARDLP_TRUE(kCSDATruth_muon,vars::ke_init,cuts::neutrino,cuts::matched_muon);
    PVAR_TTP(kCSDA_muon,vars::csda_ke,cuts::neutrino,cuts::muon,cuts::no_cut);
    PVAR_TTP(kCSDA_muon2muon,vars::csda_ke_muon,cuts::neutrino,cuts::muon,cuts::no_cut);

    PVARDLP_BIAS(kCSDA_muon_bias,vars::ke_init,vars::csda_ke_muon,cuts::neutrino,cuts::muon,cuts::no_cut);
    PVARDLP_BIAS(kCSDA_noncc_muon_bias,vars::ke_init,vars::csda_ke_muon,cuts::neutrino,cuts::contained_tpc_muon,cuts::no_cut);
    PVARDLP_BIAS(kCSDA_wellreco_muon_bias,vars::ke_init,vars::csda_ke_muon,cuts::neutrino,cuts::wellreco_muon,cuts::no_cut);

    // Variables for 2D bias plots.
    VARDLP_BIAS(kEnergy_1mu1p_signal_bias,vars::visible_energy,vars::visible_energy,cuts::signal_1mu1p,cuts::all_1mu1p_cut);
    VARDLP_BIAS(kEnergy_1mu1p_othernu_bias,vars::visible_energy,vars::visible_energy,cuts::other_nu_1mu1p,cuts::all_1mu1p_cut);
    VARDLP_BIAS(kEnergy_1mu1p_cosmic_bias,vars::visible_energy,vars::visible_energy,cuts::cosmic,cuts::all_1mu1p_cut);
    VARDLP_BIAS(kEnergy_1muNp_1p_signal_bias,vars::visible_energy,vars::visible_energy,cuts::signal_1mu1p,cuts::all_1muNp_cut);
    VARDLP_BIAS(kEnergy_1muNp_Np_signal_bias,vars::visible_energy,vars::visible_energy,cuts::signal_1muNp_Nnot1,cuts::all_1muNp_cut);
    VARDLP_BIAS(kEnergy_1muNp_othernu_bias,vars::visible_energy,vars::visible_energy,cuts::other_nu_1muNp,cuts::all_1muNp_cut);
    VARDLP_BIAS(kEnergy_1muNp_cosmic_bias,vars::visible_energy,vars::visible_energy,cuts::cosmic,cuts::all_1muNp_cut);
    VARDLP_BIAS(kEnergy_1muX_1p_signal_bias,vars::visible_energy,vars::visible_energy,cuts::signal_1mu1p,cuts::all_1muX_cut);
    VARDLP_BIAS(kEnergy_1muX_Np_signal_bias,vars::visible_energy,vars::visible_energy,cuts::signal_1muNp_Nnot1,cuts::all_1muX_cut);
    VARDLP_BIAS(kEnergy_1muX_X_bias,vars::visible_energy,vars::visible_energy,cuts::signal_1muX_notNp,cuts::all_1muX_cut);
    VARDLP_BIAS(kEnergy_1muX_othernu_bias,vars::visible_energy,vars::visible_energy,cuts::other_nu_1muX,cuts::all_1muX_cut);

    VARDLP_BIAS(kNuEnergy_1mu1p_signal_bias,vars::neutrino_energy,vars::visible_energy,cuts::signal_1mu1p,cuts::all_1mu1p_cut);
    VARDLP_BIAS(kNuEnergy_1mu1p_othernu_bias,vars::neutrino_energy,vars::visible_energy,cuts::other_nu_1mu1p,cuts::all_1mu1p_cut);
    VARDLP_BIAS(kNuEnergy_1muNp_1p_signal_bias,vars::neutrino_energy,vars::visible_energy,cuts::signal_1mu1p,cuts::all_1muNp_cut);
    VARDLP_BIAS(kNuEnergy_1muNp_Np_signal_bias,vars::neutrino_energy,vars::visible_energy,cuts::signal_1muNp_Nnot1,cuts::all_1muNp_cut);
    VARDLP_BIAS(kNuEnergy_1muNp_othernu_bias,vars::neutrino_energy,vars::visible_energy,cuts::other_nu_1muNp,cuts::all_1muNp_cut);
    VARDLP_BIAS(kNuEnergy_1muX_1p_signal_bias,vars::neutrino_energy,vars::visible_energy,cuts::signal_1mu1p,cuts::all_1muX_cut);
    VARDLP_BIAS(kNuEnergy_1muX_Np_signal_bias,vars::neutrino_energy,vars::visible_energy,cuts::signal_1muNp_Nnot1,cuts::all_1muX_cut);
    VARDLP_BIAS(kNuEnergy_1muX_X_bias,vars::neutrino_energy,vars::visible_energy,cuts::signal_1muX_notNp,cuts::all_1muX_cut);
    VARDLP_BIAS(kNuEnergy_1muX_othernu_bias,vars::neutrino_energy,vars::visible_energy,cuts::other_nu_1muX,cuts::all_1muX_cut);

    // Variables for match validation.
    PVARDLP_TRUE(kLowXTruth,vars::lowx,cuts::no_cut,cuts::matched);
    PVAR_TTP(kLowX,vars::lowx,cuts::no_cut,cuts::no_cut,cuts::no_cut);

    // Simple particle variables.
    PVARDLP_TRUE(kCCOverlap,vars::overlap,cuts::no_cut,cuts::cathode_crossing_muon);
    PVARDLP_TRUE(kNonCCOverlap,vars::overlap,cuts::no_cut,cuts::non_cathode_crossing_muon);
    VARDLP_TRUE(kProtonScattering,vars::proton_scattering_cosine,cuts::no_cut);
    VARDLP_TRUE(kLeadingProtonOverlap,vars::leading_proton_overlap,cuts::no_cut);
}