/**
 * @file definitions.h
 * @brief Header file for definitions of various helper preprocessor macros.
 * @author justin.mueller@colostate.edu
*/
#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <vector>

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Cut.h"

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

#define VARDLP_TTP(NAME,VAR,CAT,SEL)                                     \
    const SpillMultiVar NAME([](const caf::SRSpillProxy* sr)             \
    {                                                                    \
        std::vector<double> var;                                         \
        for(auto const& i : sr->dlp_true)                                \
        {                                                                \
            if(CAT(i) && i.match.size() > 0 && SEL(sr->dlp[i.match[0]])) \
            var.push_back(VAR(sr->dlp[i.match[0]]));                     \
        }                                                                \
        return var;                                                      \
    })
#endif