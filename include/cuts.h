/**
 * @file cuts.h
 * @brief Header file for definitions of selection cuts.
 * @author justin.mueller@colostate.edu 
*/
#ifndef CUTS_H
#define CUTS_H

#include <functional>
#include <vector>
#include <string>
#include <sstream>

namespace cuts
{
    /**
     * Find the topology of the interaction with cuts applied to each particle.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to find the topology of.
     * @return the topology of the interaction as a string (e.g 0ph0e1mu0pi1p).
     */
    template<class T>
        std::string topology(const T & interaction)
        {
            uint32_t counts[5] = {};
            for(auto &p : interaction.particles)
            {
                if(p.is_primary)
                {
                    double energy(p.csda_ke);
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>) energy = p.energy_deposit;
                    if((p.pid < 4 && energy > 25) || (p.pid == 4 && energy > 50))
                    ++counts[p.pid];
                }
            }
            std::stringstream ss;
            ss  << counts[0] << "ph"
                << counts[1] << "e"
                << counts[2] << "mu"
                << counts[3] << "pi"
                << counts[4] << "p";
            return ss.str();
        }

    /**
     * Apply no cut (all interactions passed).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true (always).
     */
    template<class T>
        bool no_cut(const T & interaction) { return true; }

    /**
     * Apply a fiducial volume cut. Interaction vertex must be within 25 cm of
     * x and y detector faces, 50 cm of downstream (+) z face, and 30 cm of
     * upstream (-) z face.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    template<class T>
        bool fiducial_cut(const T & interaction) { return interaction.is_fiducial; }
    
    /**
     * Apply a containment volume cut. All points within the interaction must be
     * at least 5 cm from the detector boundaries.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the vertex is contained.
     */
    template<class T>
        bool containment_cut(const T & interaction) { return interaction.is_contained; }

    /**
     * Apply a 1mu1p topological cut. The interaction must have a topology
     * matching 1mu1p as defined by the conditions in the topology() function.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has a 1mu1p topology.
     */
    template<class T>
        bool topological_cut(const T & interaction) { return topology<T>(interaction) == "0ph0e1mu0pi1p"; }
    
    /**
     * Apply a flash time cut. The interaction must be matched to an in-time
     * flash. The in-time definition is valid for BNB simulation.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction has been matched to an in-time flash.
     */
    template<class T>
        bool flash_cut(const T & interaction) { return ((interaction.fmatched == 1) && (interaction.flash_time >= 0) && (interaction.flash_time <= 1.6)); }

    /**
     * Apply a fiducial and containment cut (logical "and" of both).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial and containment cut.
     */
    template<class T>
        bool fiducial_containment_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Apply a fiducial, containment, and topological cut (logical "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment, and
     * topological cut.
     */
    template<class T>
        bool fiducial_containment_topological_cut(const T & interaction) { return fiducial_cut<T>(interaction) && containment_cut<T>(interaction) && topological_cut<T>(interaction); }
    
    /**
     * Apply a fiducial, containment, topological, and flash time cut (logical
     * "and" of each).
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction passes the fiducial, containment,
     * topological, and flash time cut.
     */
    template<class T>
        bool all_cut(const T & interaction) { return topological_cut<T>(interaction) && fiducial_cut<T>(interaction) && flash_cut<T>(interaction) && containment_cut<T>(interaction); }

    /**
     * Define the true 1mu1p interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a 1mu1p neutrino interaction.
     */
    template<class T>
        bool signal_1mu1p(const T & interaction) { return topology<T>(interaction) == "0ph0e1mu0pi1p" && interaction.is_neutrino; }

    /**
     * Define the true "other neutrino" interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is an "other neutrino" interaction.
     */
    template<class T>
        bool other_nu(const T & interaction) { return topology<T>(interaction) != "0ph0e1mu0pi1p" && interaction.is_neutrino; }
    
    /**
     * Define the true cosmic interaction classification.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to select on.
     * @return true if the interaction is a cosmic.
     */
    template<class T>
        bool cosmic(const T & interaction) { return !interaction.is_neutrino; }
}
#endif