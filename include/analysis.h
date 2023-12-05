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
    VARDLP_RECO(kCountParticles,vars::count_particles,cuts::no_cut);
    VARDLP_RECO(kCountPrimaries,vars::count_primaries,cuts::no_cut);

    VARDLP_TRUE(kCountParticlesTruth,vars::count_particles,cuts::no_cut);
    VARDLP_TRUE(kCountPrimariesTruth,vars::count_primaries,cuts::no_cut);

    DEFINECAT();
    TCATVAR(kCountTTP,count);
    RCATVAR(kCountPTT,count);
    TCATVAR(kVisibleEnergyTTP,visible_energy)

    PVARDLP_RECO(kCSDA,vars::csda_ke,cuts::no_cut);
    PVARDLP_TRUE(kCSDATruth,vars::csda_ke,cuts::no_cut);
    
    PVARDLP_TRUE(kPrimaryTruth,vars::primary,cuts::matched);
    PVARDLP_TRUE(kPIDTruth,vars::pid,cuts::matched);
    PVARDLP_TRUE(kPrimaryPIDTruth,vars::primary_pid,cuts::matched);
    PVAR_TTP(kPrimary,vars::primary,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPID,vars::pid,cuts::no_cut,cuts::no_cut);
    PVAR_TTP(kPrimaryPID,vars::primary_pid,cuts::no_cut,cuts::no_cut);

    PVARDLP_TRUE(kPrimaryTruth_Neutrino,vars::primary,cuts::matched_neutrino);
    PVARDLP_TRUE(kPIDTruth_Neutrino,vars::pid,cuts::matched_neutrino);
    PVARDLP_TRUE(kPrimaryPIDTruth_Neutrino,vars::primary_pid,cuts::matched_neutrino);
    PVAR_TTP(kPrimary_Neutrino,vars::primary,cuts::neutrino,cuts::no_cut);
    PVAR_TTP(kPID_Neutrino,vars::pid,cuts::neutrino,cuts::no_cut);
    PVAR_TTP(kPrimaryPID_Neutrino,vars::primary_pid,cuts::neutrino,cuts::no_cut);

    PVARDLP_TRUE(kPrimaryTruth_Neutrino,vars::primary,cuts::matched_cosmic);
    PVARDLP_TRUE(kPIDTruth_Neutrino,vars::pid,cuts::matched_cosmic);
    PVARDLP_TRUE(kPrimaryPIDTruth_Neutrino,vars::primary_pid,cuts::matched_cosmic);
    PVAR_TTP(kPrimary_Neutrino,vars::primary,cuts::cosmic,cuts::no_cut);
    PVAR_TTP(kPID_Neutrino,vars::pid,cuts::cosmic,cuts::no_cut);
    PVAR_TTP(kPrimaryPID_Neutrino,vars::primary_pid,cuts::cosmic,cuts::no_cut);
}