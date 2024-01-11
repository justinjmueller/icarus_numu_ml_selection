/**
 * @file analysis.h
 * @brief Header file for definitions of complete variables with relevant cuts.
 * @author justin.mueller@colostate.edu
*/

#include "definitions.h"
#include "cuts.h"
#include "variables.h"

using namespace ana;

namespace ana
{
    // "Simple" variables that attach to a single interaction attribute.
    VARDLP_RECO(kCountParticles,vars::count_particles,cuts::no_cut);
    VARDLP_RECO(kCountPrimaries,vars::count_primaries,cuts::no_cut);

    VARDLP_TRUE(kCountParticlesTruth,vars::count_particles,cuts::no_cut);
    VARDLP_TRUE(kCountPrimariesTruth,vars::count_primaries,cuts::no_cut);

    // Define "variables" for binning interactions by some categorical class.
    DEFINECAT();

    // Define variables that are broadcasted across each selection cut.
    TCATVAR(kCountTTP,count);
    RCATVAR(kCountPTT,count);
    TCATVAR(kVisibleEnergyTTP,visible_energy);
    RCATVAR(kFlashTime,flash_time);
    
    // Variables for confusion matrix.
    PVARDLP_TRUE(kPrimaryTruth,vars::primary,cuts::no_cut,cuts::matched);
    PVARDLP_TRUE(kPIDTruth,vars::pid,cuts::no_cut,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth,vars::primary_pid,cuts::no_cut,cuts::matched);
    PVAR_TTP(kPrimary,vars::primary,cuts::no_cut,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID,vars::pid,cuts::no_cut,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID,vars::primary_pid,cuts::no_cut,cuts::no_cut,cuts::no_cut);

    // Variables for neutrino-only confusion matrix.
    PVARDLP_TRUE(kPrimaryTruth_Neutrino,vars::primary,cuts::neutrino,cuts::matched);
    PVARDLP_TRUE(kPIDTruth_Neutrino,vars::pid,cuts::neutrino,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth_Neutrino,vars::primary_pid,cuts::neutrino,cuts::matched);
    PVAR_TTP(kPrimary_Neutrino,vars::primary,cuts::neutrino,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID_Neutrino,vars::pid,cuts::neutrino,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID_Neutrino,vars::primary_pid,cuts::neutrino,cuts::no_cut,cuts::no_cut);

    // Variables for cosmic-only confusion matrix.
    PVARDLP_TRUE(kPrimaryTruth_Cosmic,vars::primary,cuts::cosmic,cuts::matched);
    PVARDLP_TRUE(kPIDTruth_Cosmic,vars::pid,cuts::cosmic,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth_Cosmic,vars::primary_pid,cuts::cosmic,cuts::matched);
    PVAR_TTP(kPrimary_Cosmic,vars::primary,cuts::cosmic,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID_Cosmic,vars::pid,cuts::cosmic,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID_Cosmic,vars::primary_pid,cuts::cosmic,cuts::no_cut,cuts::no_cut);

    // Variables for 2D "true vs. reco" style plots.
    PVARDLP_TRUE(kCSDATruth,vars::csda_ke,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kCSDA,vars::csda_ke,cuts::no_cut,cuts::no_cut,cuts::no_cut);

    PVARDLP_TRUE(kLowXTruth,vars::lowx,cuts::no_cut,cuts::matched);
    PVAR_TTP(kLowX,vars::lowx,cuts::no_cut,cuts::no_cut,cuts::no_cut);
}