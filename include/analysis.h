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

    TCATVAR(kCountTTP,count);
    RCATVAR(kCountPTT,count);

    PVARDLP_RECO(kCSDA,vars::csda_ke,cuts::no_cut);
    PVARDLP_TRUE(kCSDATruth,vars::csda_ke,cuts::no_cut);
    PVARDLP_TRUE(kPIDTruth,vars::pid,cuts::matched);

    PVAR_TTP(kPID_muon,vars::pid,cuts::muon,cuts::no_cut);
    PVAR_TTP(kPID,vars::pid,cuts::no_cut,cuts::no_cut);
}