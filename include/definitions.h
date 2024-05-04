/**
 * @file definitions.h
 * @brief Header file for definitions of various helper preprocessor macros.
 * @author justin.mueller@colostate.edu
*/
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <vector>
#include <map>

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

/**
 * Preprocessor wrapper for looping over reco interactions. The SpillMultiVar
 * accepts a vector as a result of some function running over the top-level
 * StandardRecord. This wrapper broadcasts a function across all interactions
 * within the reco interaction.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each interaction passing
 * the cut SEL.
*/
#define VARDLP_RECO(NAME,VAR,SEL)                             \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)  \
    {                                                         \
        std::vector<double> var;                              \
        for(auto const& i : sr->dlp)                          \
        {                                                     \
            if(SEL(i))                                        \
                var.push_back(VAR(i));                        \
        }                                                     \
        return var;                                           \
    })

/**
 * Preprocessor wrapper for looping over true interactions. The SpillMultiVar
 * accepts a vector as a result of some function running over the top-level
 * StandardRecord. This wrapper broadcasts a function across all interactions
 * within the true interaction.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the interactions.
 * @param SEL function to select interactions.
 * @return a vector with the result of VAR called on each interaction passing
 * the cut SEL.
*/
#define VARDLP_TRUE(NAME,VAR,SEL)                             \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)  \
    {                                                         \
        std::vector<double> var;                              \
        for(auto const& i : sr->dlp_true)                     \
        {                                                     \
            if(SEL(i))                                        \
                var.push_back(VAR(i));                        \
        }                                                     \
        return var;                                           \
    })

/**
 * Preprocessor wrapper for looping over true interactions and broadcasting a
 * SpillMultiVar over the matched (truth->reco) interaction. The SpillMultiVar
 * accepts a vector as a result of some function running over the top-level
 * StandardRecord.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the reco interactions.
 * @param CAT function that defines the truth category.
 * @param SEL function to select reco interactions.
 * @return a vector with the result of VAR called on each reco interaction
 * passing SEL that is matched to by the true interaction passing category
 * cut CAT.
*/
#define VARDLP_TTP(NAME,VAR,CAT,SEL)                                     \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)             \
    {                                                                    \
        std::vector<double> var;                                         \
        for(auto const& i : sr->dlp_true)                                \
        {                                                                \
            if(CAT(i) && i.match.size() > 0 && SEL(sr->dlp[i.match[0]])) \
                var.push_back(VAR(sr->dlp[i.match[0]]));                 \
        }                                                                \
        return var;                                                      \
    })

/**
 * Preprocessor wrapper for looping over reco interactions which match to
 * truth interactions of the specified category and applying a SpillMultiVar.
 * The SpillMultiVar accepts a vector as a result of some function running over
 * the top-level StandardRecord.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the reco interactions.
 * @param CAT function that defines the truth category.
 * @param SEL function to select reco interactions.
 * @return a vector with the result of VAR called on each reco interaction
 * passing SEL that is matched to by the true interaction passing category
 * cut CAT.
*/
#define VARDLP_PTT(NAME,VAR,CAT,SEL)                                          \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)                  \
    {                                                                         \
        std::vector<double> var;                                              \
        for(auto const& i : sr->dlp)                                          \
        {                                                                     \
            if(SEL(i) && i.match.size() > 0 && CAT(sr->dlp_true[i.match[0]])) \
                var.push_back(VAR(i));                                        \
        }                                                                     \
        return var;                                                           \
    })

/**
 * Preprocessor wrapper for looping over true interactions and broadcasting a
 * SpillMultiVar over the matched (truth->reco) interaction. The SpillMultiVar
 * accepts a vector as a result of some function running over the top-level
 * StandardRecord. This wrapper will calculate the bias between two variables
 * between truth and reco.
 * @param NAME of the resulting SpillMultiVar.
 * @param TVAR function to broadcast over the true interactions.
 * @param RVAR function to broadcast over the reco interactions.
 * @param CAT function that defines the truth category.
 * @param SEL function to select reco interactions.
 * @return a vector with the bias between TVAR and RVAR as called on reco
 * interactions passing SEL that are matched to by the true interaction passing
 * category cut CAT.
*/
#define VARDLP_BIAS(NAME,TVAR,RVAR,CAT,SEL)                                     \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)                    \
    {                                                                           \
        std::vector<double> var;                                                \
        for(auto const& i : sr->dlp_true)                                       \
        {                                                                       \
            if(CAT(i) && i.match.size() > 0 && SEL(sr->dlp[i.match[0]]))        \
                var.push_back((RVAR(sr->dlp[i.match[0]]) - TVAR(i)) / TVAR(i)); \
        }                                                                       \
        return var;                                                             \
    })

/**
 * Preprocessor wrapper for looping over true particles and broadcasting a
 * SpillMultiVar over the matched (truth->reco) particle. The SpillMultiVar
 * accepts a vector as a result of some function running over the top-level
 * StandardRecord. This wrapper will calculate the bias between two variables
 * between truth and reco.
 * @param NAME of the resulting SpillMultiVar.
 * @param TVAR function to broadcast over the true particles.
 * @param RVAR function to broadcast over the reco particles.
 * @param ICAT function that defines the truth category (interactions).
 * @param PCAT function that defines the truth category (particles).
 * @param SEL function to select reco particles.
 * @return a vector with the bias between TVAR and RVAR as called on reco
 * particles passing SEL that are matched to by the true particle passing
 * category cut CAT.
*/
#define PVARDLP_BIAS(NAME,TVAR,RVAR,ICAT,PCAT,SEL)                                               \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)                                     \
    {                                                                                            \
        std::vector<double> var;                                                                 \
        std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleDLP> *> reco_particles;    \
        for(auto const& i : sr->dlp)                                                             \
        {                                                                                        \
            for(auto const& p : i.particles)                                                     \
                reco_particles.insert(std::make_pair(p.id, &p));                                 \
        }                                                                                        \
        for(auto const& i : sr->dlp_true)                                                        \
        {                                                                                        \
            for(auto const& p : i.particles)                                                     \
            {                                                                                    \
                if(ICAT(i) && PCAT(p) && p.match.size() > 0 && SEL(*reco_particles[p.match[0]])) \
                    var.push_back((RVAR(*reco_particles[p.match[0]]) - TVAR(p)) / TVAR(p));      \
            }                                                                                    \
        }                                                                                        \
        return var;                                                                              \
    })

/**
 * Preprocessor wrapper for looping over reco particles. The SpillMultiVar
 * accepts a vector as a result of some function running over the top-level
 * StandardRecord. This wrapper broadcasts a function across all particles
 * within the reco interactions.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the particles.
 * @param SEL function to select particles.
 * @return a vector with the result of VAR called on each particle passing
 * the cut SEL.
*/
#define PVARDLP_RECO(NAME,VAR,SEL)                            \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)  \
    {                                                         \
        std::vector<double> var;                              \
        for(auto const& i : sr->dlp)                          \
        {                                                     \
            for(auto const& p : i.particles)                  \
            {                                                 \
                if(SEL(p))                                    \
                    var.push_back(VAR(p));                    \
            }                                                 \
        }                                                     \
        return var;                                           \
    })

/**
 * Preprocessor wrapper for looping over true particles. The SpillMultiVar
 * accepts a vector as a result of some function running over the top-level
 * StandardRecord. This wrapper broadcasts a function across all particles
 * within the true interactions.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the particles.
 * @param ISEL function to select interactions.
 * @param PSEL function to select particles.
 * @return a vector with the result of VAR called on each particle passing
 * the cut SEL.
*/
#define PVARDLP_TRUE(NAME,VAR,ISEL,PSEL)                      \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)  \
    {                                                         \
        std::vector<double> var;                              \
        for(auto const& i : sr->dlp_true)                     \
        {                                                     \
            for(auto const& p : i.particles)                  \
            {                                                 \
                if(ISEL(i) && PSEL(p))                        \
                    var.push_back(VAR(p));                    \
            }                                                 \
        }                                                     \
        return var;                                           \
    })

/**
 * Preprocessor wrapper for looping over true particles and broadcasting a
 * SpillMultiVar over the matched (truth->reco) particle. The SpillMultiVar
 * accepts a vector as a result of some function running over the top-level
 * StandardRecord.
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the reco particles.
 * @param ICAT function that defines the truth category (interactions).
 * @param PCAT function that defines the truth category (particles).
 * @param SEL function to select reco particles.
 * @return a vector with the result of VAR called on each reco particle
 * passing SEL that is matched to by the true particle passing category
 * cut CAT.
*/
#define PVAR_TTP(NAME,VAR,ICAT,PCAT,SEL)                                                         \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)                                     \
    {                                                                                            \
        std::vector<double> var;                                                                 \
        std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleDLP> *> reco_particles;    \
        for(auto const& i : sr->dlp)                                                             \
        {                                                                                        \
            for(auto const& p : i.particles)                                                     \
                reco_particles.insert(std::make_pair(p.id, &p));                                 \
        }                                                                                        \
        for(auto const& i : sr->dlp_true)                                                        \
        {                                                                                        \
            for(auto const& p : i.particles)                                                     \
            {                                                                                    \
                if(ICAT(i) && PCAT(p) && p.match.size() > 0 && SEL(*reco_particles[p.match[0]])) \
                    var.push_back(VAR(*reco_particles[p.match[0]]));                             \
            }                                                                                    \
        }                                                                                        \
        return var;                                                                              \
    })

/**
 * Preprocessor wrapper for looping over true interactions which match to a
 * specified category of reco interaction and broadcasting a variable on the
 * true interaction. The SpillMultiVar accepts a vector as a result of some
 * function running over the top-level StandardRecord.
 * @param NAME of the SpillMultiVar.
 * @param VAR function to broadcast over the true interactions.
 * @param SEL function to select reco interactions.
 * @return a vector containing the results of VAR called on each true
 * interaction which is matched to a reco interaction of the specified category.
*/
#define VARDLP_TCAT(NAME,VAR,SEL)                              \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)   \
    {                                                          \
        std::vector<double> var;                               \
        for(auto const& i : sr->dlp_true)                      \
        {                                                      \
            if(i.match.size() > 0 && SEL(sr->dlp[i.match[0]])) \
                var.push_back(VAR(i));                         \
        }                                                      \
        return var;                                            \
    })

/**
 * Preprocessor wrapper for looping over reco interactions which match to a
 * true interaction and broadcasting a variable on the reco or true interaction
 * (reco if not categorical). The SpillMultiVar accepts a vector as a result of
 * some function running over the top-level StandardRecord.
 * @param NAME of the SpillMultiVar.
 * @param VAR function to broadcast over the true interactions.
 * @param SEL function to select reco interactions.
 * @return a vector containing the results of VAR called on each true
 * interaction which is matched to by a reco interaction of the specified category.
*/
#define VARDLP_RCAT(NAME,VAR,SEL)                                                  \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)                       \
    {                                                                              \
        std::vector<double> var;                                                   \
        for(auto const& i : sr->dlp)                                               \
        {                                                                          \
            if(SEL(i) && i.match.size() > 0)                                       \
                var.push_back(VAR(sr->dlp_true[i.match[0]]));                      \
        }                                                                          \
        return var;                                                                \
    })

/**
 * Preprocessor macro for defining categorical variables for all cuts.
*/
#define DEFINECAT()                                                                                             \
    VARDLP_TCAT(kCategoryTTP_NoCut,vars::category,cuts::no_cut);                                                \
    VARDLP_TCAT(kCategoryTTP_FVCut,vars::category,cuts::fiducial_cut);                                          \
    VARDLP_TCAT(kCategoryTTP_FVConCut,vars::category,cuts::fiducial_containment_cut);                           \
    VARDLP_TCAT(kCategoryTTP_FVConTop1mu1pCut,vars::category,cuts::fiducial_containment_topological_1mu1p_cut); \
    VARDLP_TCAT(kCategoryTTP_FVConTop1muNpCut,vars::category,cuts::fiducial_containment_topological_1muNp_cut); \
    VARDLP_TCAT(kCategoryTTP_FVConTop1muXCut,vars::category,cuts::fiducial_containment_topological_1muX_cut);   \
    VARDLP_TCAT(kCategoryTTP_All1mu1pCut,vars::category,cuts::all_1mu1p_cut);                                   \
    VARDLP_TCAT(kCategoryTTP_All1muNpCut,vars::category,cuts::all_1muNp_cut);                                   \
    VARDLP_TCAT(kCategoryTTP_All1muXCut,vars::category,cuts::all_1muX_cut);                                     \
    VARDLP_RCAT(kCategoryPTT_NoCut,vars::category,cuts::no_cut);                                                \
    VARDLP_RCAT(kCategoryPTT_FVCut,vars::category,cuts::fiducial_cut);                                          \
    VARDLP_RCAT(kCategoryPTT_FVConCut,vars::category,cuts::fiducial_containment_cut);                           \
    VARDLP_RCAT(kCategoryPTT_FVConTop1mu1pCut,vars::category,cuts::fiducial_containment_topological_1mu1p_cut); \
    VARDLP_RCAT(kCategoryPTT_FVConTop1muNpCut,vars::category,cuts::fiducial_containment_topological_1muNp_cut); \
    VARDLP_RCAT(kCategoryPTT_FVConTop1muXCut,vars::category,cuts::fiducial_containment_topological_1muX_cut);   \
    VARDLP_RCAT(kCategoryPTT_All1mu1pCut,vars::category,cuts::all_1mu1p_cut);                                   \
    VARDLP_RCAT(kCategoryPTT_All1muNpCut,vars::category,cuts::all_1muNp_cut);                                   \
    VARDLP_RCAT(kCategoryPTT_All1muXCut,vars::category,cuts::all_1muX_cut);                                     \
    VARDLP_RCAT(kCategoryTopologyPTT_NoCut,vars::category_topology,cuts::no_cut);                               \
    VARDLP_RCAT(kCategoryTopologyPTT_All1mu1pCut,vars::category_topology,cuts::all_1mu1p_cut);                  \
    VARDLP_RCAT(kCategoryTopologyPTT_All1muNpCut,vars::category_topology,cuts::all_1muNp_cut);                  \
    VARDLP_RCAT(kCategoryTopologyPTT_All1muXCut,vars::category_topology,cuts::all_1muX_cut);                    \
    VARDLP_RCAT(kCategoryInteractionModePTT_All1mu1pCut,vars::category_interaction_mode,cuts::all_1mu1p_cut);   \
    VARDLP_RCAT(kCategoryInteractionModePTT_All1muNpCut,vars::category_interaction_mode,cuts::all_1muNp_cut);   \
    VARDLP_RCAT(kCategoryInteractionModePTT_All1muXCut,vars::category_interaction_mode,cuts::all_1muX_cut);


/**
 * Preprocessor macro for broadcasting a variable across true interactions
 * with categorical information for all cuts.
 * @param NAME (base) to assign to the variable.
 * @param VAR to broadcast.
*/
#define TCATVAR(NAME,VAR)                                                                                       \
    VARDLP_TCAT(NAME ## _NoCut,vars::VAR,cuts::no_cut);                                                         \
    VARDLP_TCAT(NAME ## _FVCut,vars::VAR,cuts::fiducial_cut);                                                   \
    VARDLP_TCAT(NAME ## _FVConCut,vars::VAR,cuts::fiducial_containment_cut);                                    \
    VARDLP_TCAT(NAME ## _FVConTop1mu1pCut,vars::VAR,cuts::fiducial_containment_topological_1mu1p_cut);          \
    VARDLP_TCAT(NAME ## _FVConTop1muNpCut,vars::VAR,cuts::fiducial_containment_topological_1muNp_cut);          \
    VARDLP_TCAT(NAME ## _FVConTop1muXCut,vars::VAR,cuts::fiducial_containment_topological_1muX_cut);            \
    VARDLP_TCAT(NAME ## _All1mu1pCut,vars::VAR,cuts::all_1mu1p_cut);                                            \
    VARDLP_TCAT(NAME ## _All1muNpCut,vars::VAR,cuts::all_1muNp_cut);                                            \
    VARDLP_TCAT(NAME ## _All1muXCut,vars::VAR,cuts::all_1muX_cut);


/**
 * Preprocessor macro for broadcasting a variable across reco interactions
 * with categorical information for all cuts.
 * @param NAME (base) to assign to the variable.
 * @param VAR to broadcast.
*/
#define RCATVAR(NAME,VAR)                                                                                          \
    VARDLP_PTT(NAME ## _NoCut,vars::VAR,cuts::no_cut,cuts::no_cut);                                                \
    VARDLP_PTT(NAME ## _FVCut,vars::VAR,cuts::no_cut,cuts::fiducial_cut);                                          \
    VARDLP_PTT(NAME ## _FVConCut,vars::VAR,cuts::no_cut,cuts::fiducial_containment_cut);                           \
    VARDLP_PTT(NAME ## _FVConTop1mu1pCut,vars::VAR,cuts::no_cut,cuts::fiducial_containment_topological_1mu1p_cut); \
    VARDLP_PTT(NAME ## _FVConTop1muNpCut,vars::VAR,cuts::no_cut,cuts::fiducial_containment_topological_1muNp_cut); \
    VARDLP_PTT(NAME ## _FVConTop1muXCut,vars::VAR,cuts::no_cut,cuts::fiducial_containment_topological_1muX_cut);   \
    VARDLP_PTT(NAME ## _All1mu1pCut,vars::VAR,cuts::no_cut,cuts::all_1mu1p_cut);                                   \
    VARDLP_PTT(NAME ## _All1muNpCut,vars::VAR,cuts::no_cut,cuts::all_1muNp_cut);                                   \
    VARDLP_PTT(NAME ## _All1muXCut,vars::VAR,cuts::no_cut,cuts::all_1muX_cut);

/**
 * Preprocessor macro for broadcasting a variable across reco interactions
 * passing the three main signal cuts.
 * @param NAME (base) to assign to the variable.
 * @param VAR to broadcast.
*/
#define RECO_SIGNAL_VAR(NAME,VAR)                                         \
    VARDLP_RECO(NAME ## _1mu1p,VAR,cuts::all_1mu1p_data_cut);             \
    VARDLP_RECO(NAME ## _1muNp,VAR,cuts::all_1muNp_data_cut);             \
    VARDLP_RECO(NAME ## _1muX,VAR,cuts::all_1muX_data_cut);               \

#endif