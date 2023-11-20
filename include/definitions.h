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
 * Preprocessor wrapper for looping
 * @param NAME of the resulting SpillMultiVar.
 * @param VAR function to broadcast over the reco particles.
 * @param CAT function that defines the truth category.
 * @param SEL function to select reco particles.
 * @return a vector with the result of VAR called on each reco particle
 * passing SEL that is matched to by the true particle passing category
 * cut CAT.
*/
#define PVAR_TTP(NAME,VAR,CAT,SEL)                                                            \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)                                  \
    {                                                                                         \
        std::vector<double> var;                                                              \
        std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleDLP> *> reco_particles; \
        for(auto const& i : sr->dlp)                                                          \
        {                                                                                     \
            for(auto const& p : i.particles)                                                  \
                reco_particles.insert(std::make_pair(p.id, &p));                              \
        }                                                                                     \
        for(auto const& i : sr->dlp_true)                                                     \
        {                                                                                     \
            for(auto const& p : i.particles)                                                  \
            {                                                                                 \
                if(CAT(p) && p.match.size() > 0 && SEL(*reco_particles[p.match[0]]))          \
                    var.push_back(VAR(*reco_particles[p.match[0]]));                          \
            }                                                                                 \
        }                                                                                     \
        return var;                                                                           \
    })
#endif