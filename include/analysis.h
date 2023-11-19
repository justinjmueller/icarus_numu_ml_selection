/**
 * @file analysis.h
 * @brief Header file for definitions of complete variables with relevant cuts.
 * @author justin.mueller@colostate.edu
*/

#include "definitions.h"
#include "cuts.h"
#include "variables.h"

namespace ana
{
    VARDLP_RECO(kCountParticles,vars::count_particles,cuts::no_cut);
    VARDLP_RECO(kCountPrimaries,vars::count_primaries,cuts::no_cut);

    VARDLP_TRUE(kCountParticlesTruth,vars::count_particles,cuts::no_cut);
    VARDLP_TRUE(kCountPrimariesTruth,vars::count_primaries,cuts::no_cut);

    VARDLP_TTP(kCount_1mu1p_NoCut,vars::count,cuts::signal_1mu1p,cuts::no_cut);
    VARDLP_TTP(kCount_OtherNu_NoCut,vars::count,cuts::other_nu,cuts::no_cut);
    VARDLP_TTP(kCount_Cosmic_NoCut,vars::count,cuts::cosmic,cuts::no_cut);

    VARDLP_TTP(kCount_1mu1p_FVCut,vars::count,cuts::signal_1mu1p,cuts::fiducial_cut);
    VARDLP_TTP(kCount_OtherNu_FVCut,vars::count,cuts::other_nu,cuts::fiducial_cut);
    VARDLP_TTP(kCount_Cosmic_FVCut,vars::count,cuts::cosmic,cuts::fiducial_cut);

    VARDLP_TTP(kCount_1mu1p_FVConCut,vars::count,cuts::signal_1mu1p,cuts::fiducial_containment_cut);
    VARDLP_TTP(kCount_OtherNu_FVConCut,vars::count,cuts::other_nu,cuts::fiducial_containment_cut);
    VARDLP_TTP(kCount_Cosmic_FVConCut,vars::count,cuts::cosmic,cuts::fiducial_containment_cut);

    VARDLP_TTP(kCount_1mu1p_FVConTopCut,vars::count,cuts::signal_1mu1p,cuts::fiducial_containment_topological_cut);
    VARDLP_TTP(kCount_OtherNu_FVConTopCut,vars::count,cuts::other_nu,cuts::fiducial_containment_topological_cut);
    VARDLP_TTP(kCount_Cosmic_FVConTopCut,vars::count,cuts::cosmic,cuts::fiducial_containment_topological_cut);

    VARDLP_TTP(kCount_1mu1p_AllCut,vars::count,cuts::signal_1mu1p,cuts::all_cut);
    VARDLP_TTP(kCount_OtherNu_AllCut,vars::count,cuts::other_nu,cuts::all_cut);
    VARDLP_TTP(kCount_Cosmic_AllCut,vars::count,cuts::cosmic,cuts::all_cut);
}