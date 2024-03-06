/**
 * @file variables.h
 * @brief Header file for definitions of selection variables.
 * @author justin.mueller@colostate.edu
*/
#ifndef VARIABLES_H
#define VARIABLES_H

#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813

#include <algorithm>

namespace vars
{

    /**
     * Variable for counting interactions/particles.
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return 1.0 (always).
     */
    template<class T>
        double count(const T & obj) { return 1.0; }

    /**
     * Variable for image_id (unique identifier for the event).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the image_id of the interaction/particle.
    */
    template<class T>
        double image_id(const T & obj) { return obj.image_id; }

    /**
     * Variable for id (unique identifier for the object).
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the id of the interaction/particle.
    */
    template<class T>
        double id(const T & obj) { return obj.id; }

    /**
     * Variable for the cryostat of the object.
     * @tparam T the type of object (true or reco, interaction or particle).
     * @param interaction/particle to apply the variable on.
     * @return the cryostat of the interaction/particle.
    */
    template<class T>
        double cryostat(const T & obj) { return obj.volume_id; }

    /**
     * Variable for enumerating interaction categories. This is a basic
     * categorization using only signal, neutrino background, and cosmic
     * background as the three categories.
     * 0: 1mu1p (contained and fiducial)
     * 1: 1mu1p (not contained or fiducial)
     * 2: 1muNp (N > 1, contained and fiducial)
     * 3: 1muNp (N > 1, not contained or fiducial)
     * 4: 1muX (not 1muNp, contained and fiducial)
     * 5: 1muX (not 1muNp, not contained or fiducial)
     * 6: Other nu
     * 7: cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    template<class T>
        double category(const T & interaction)
        {
            double cat(7);
            if(cuts::signal_1mu1p(interaction))
            {
                if(cuts::fiducial_containment_cut(interaction)) cat = 0;
                else cat = 1;
            }
            else if(cuts::signal_1muNp(interaction))
            {
                if(cuts::fiducial_containment_cut(interaction)) cat = 2;
                else cat = 3;
            }
            else if(cuts::signal_1muX(interaction))
            {
                if(cuts::fiducial_containment_cut(interaction)) cat = 4;
                else cat = 5;
            }
            else if(cuts::other_nu_1muX(interaction)) cat = 6;
            return cat;
        }

    /**
     * Variable for enumerating interaction categories. This classifies the
     * interactions based on the visible final states.
     * 0: 1mu1p, 1: 1mu0h, 2: 1muNp (N>1), 3: 1mu1p1pi, 4: nu_mu CC Other, 5: NC, 6: Cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    template<class T>
        double category_topology(const T & interaction)
        {
            uint16_t cat(6);
            if(interaction.is_neutrino)
            {
                std::vector<uint32_t> counts(cuts::count_primaries(interaction));
                if(counts[0] == 0 && counts[1] == 0 && counts[2] == 1)
                {
                    if(counts[3] == 0 && counts[4] == 1 && interaction.is_contained && interaction.is_fiducial) cat = 0;
                    else if(counts[3] == 0 && counts[4] == 1) cat = 7;
                    else if(counts[3] == 0 && counts[4] == 0) cat = 1;
                    else if(counts[3] == 0 && counts[4] > 1 && interaction.is_contained && interaction.is_fiducial) cat = 2;
                    else if(counts[3] == 0 && counts[4] > 1) cat = 7;
                    else if(counts[3] == 1 && counts[4] == 1) cat = 3;
                    else if(interaction.nu_current_type == 0) cat = 4;
                }
                else if(interaction.nu_current_type == 0) cat = 4;
                else if(interaction.nu_current_type == 1) cat = 5;
            }
            return cat;
        }

    /**
     * Variable for enumerating interaction categories. This categorization
     * uses the interaction type (generator truth) classify the interactions
     * 0: nu_mu CC QE, 1: nu_mu CC Res, 2: nu_mu CC MEC, 3: nu_mu CC DIS, 4: nu_mu CC Coh, 5: nu_e CC, 6: NC, 7: Cosmic
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the enumerated category of the interaction.
    */
    template<class T>
        double category_interaction_mode(const T & interaction)
        {
            double cat(7);

            if(interaction.is_neutrino)
            {
                if(interaction.nu_current_type == 0)
                {
                    if(abs(interaction.nu_pdg_code) == 14)
                    {
                        if(interaction.nu_interaction_mode == 0) cat = 0;
                        else if(interaction.nu_interaction_mode == 1) cat = 1;
                        else if(interaction.nu_interaction_mode == 10) cat = 2;
                        else if(interaction.nu_interaction_mode == 2) cat = 3;
                        else if(interaction.nu_interaction_mode == 3) cat = 4;
                        else cat = 8;
                    }
                    else cat = 5;
                }
                else cat = 6;
            }

            return cat;
        }

    /**
     * Variable for counting particles in interactions.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the number of particles in the interaction.
     */
    template<class T>
        double count_particles(const T & interaction) { return interaction.num_particles; }

    /**
     * Variable for counting primaries in interactions.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the number of primaries in the interaction.
     */
    template<class T>
        double count_primaries(const T & interaction) { return interaction.num_primaries; }
    
    /**
     * Variable for total visible energy of interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the total visible energy of the interaction.
    */
    template<class T>
        double visible_energy(const T & interaction)
        {
            double energy(0);
            for(const auto & p : interaction.particles)
            {
                if(p.is_primary)
                {
                    if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    {
                        energy += p.energy_deposit;
                    }
                    else
                    {
                        if(p.pid < 2) energy += p.calo_ke;
                        else energy += p.csda_ke;
                    }
                    if(p.pid == 2) energy += MUON_MASS;
                    else if(p.pid == 3) energy += PION_MASS;
                }
            }
            return energy;
        }
    /**
     * Variable for energy of the neutrino primary of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the neutrino energy.
    */
    template<class T>
        double neutrino_energy(const T & interaction) { return 1000*interaction.nu_energy_init; }

    /**
     * Variable for matched interaction flash time.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the matched flash time of the interaction.
    */
    template<class T>
        double flash_time(const T & interaction)
        {
            if(!cuts::valid_flashmatch(interaction))
                return -100000.0;
            else
                return interaction.flash_time;
        }

    /**
     * Variable for particle primary categorizations.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the primary/non-primary designation of the particle.
    */
    template<class T>
        double primary(const T & particle) { return particle.is_primary ? 1 : 0; }

    /**
     * Variable for particle PID.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the PID of the particle.
    */
    template<class T>
        double pid(const T & particle) { return particle.pid; }

    /**
     * Variable for particle PID + primary.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the PID+primary information for the particle.
    */
    template<class T>
        double primary_pid(const T & particle) { return particle.pid + (particle.is_primary ? 5 : 0); }

    /**
     * Variable for particle csda_ke.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the csda_ke of the particle.
    */
    template<class T>
        double csda_ke(const T & particle) { return particle.csda_ke; }

    /**
     * Variable for particle calo_ke.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the calo_ke of the particle.
    */
    template<class T>
        double calo_ke(const T & particle) { return particle.calo_ke; }

    /**
     * Variable for particle csda_ke (muons only).
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the csda_ke of the particle (if a muon).
    */
    template<class T>
        double csda_ke_muon(const T & particle) { return (cuts::muon(particle)) ? csda_ke(particle) : -1; }

    /**
     * Variable for true particle energy deposited.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the energy deposited by the particle (truth only).
    */
    template<class T>
        double energy_deposit(const T & particle) { return particle.energy_deposit; }

    /**
     * Variable for true particle energy starting kinetic energy.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the starting kinetic energy of the particle.
    */
    template<class T>
        double ke_init(const T & particle)
        {
            double energy(particle.energy_init);
            switch (particle.pid)
            {
            case 1:
                energy -= ELECTRON_MASS;
                break;
            case 2:
                energy -= MUON_MASS;
                break;
            case 3:
                energy -= PION_MASS;
                break;
            case 4:
                energy -= PROTON_MASS;
                break;
            default:
                break;
            }
            return energy;
        }
    
    /**
     * Variable for particle overlap (IoU) of best match.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the particle overlap of the best match (IoU).
    */
    template<class T>
        double overlap(const T & particle) { return particle.match.size() > 0 ? (double)particle.match_overlap[0] : 0.0; }

    /**
     * Variable for the particles lowest x-coordinate.
     * @tparam T the type of particle (true or reco)
     * @param particle to apply the variable on.
     * @return the lowest x-coordinate of the particle start/end points.
    */
    template<class T>
        double lowx(const T & particle)
        {
            if(std::isnan(particle.start_point[0]) || std::isnan(particle.end_point[0]) || std::isinf(particle.start_point[0]) || std::isinf(particle.end_point[0]))
                return -100000.;
            else
                return std::min(particle.start_point[0], particle.end_point[0]);
        }

    /**
     * Finds the index corresponding to the leading particle of the specifed
     * particle type.
     * @tparam T the type of intearction (true or reco).
     * @param interaction to operate on.
     * @param pid of the particle type.
     * @return the index of the leading particle (highest KE). 
    */
    template <class T>
        size_t leading_particle_index(const T & interaction, uint16_t pid)
        {
            double leading_ke(0);
            size_t index(0);
            for(size_t i(0); i < interaction.particles.size(); ++i)
            {
                const auto & p = interaction.particles[i];
                double energy(csda_ke(p));
                if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                    energy = ke_init(p);
                if(p.pid == pid && energy > leading_ke)
                {
                    leading_ke = energy;
                    index = i;
                }
            }
            return index;
        }
    
    /**
     * Variable for finding the leading muon kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the kinetic energy of the leading muon.
    */
    template<class T>
        double leading_muon_ke(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            double energy(csda_ke(interaction.particles[i]));
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                energy = ke_init(interaction.particles[i]);
            return energy;
        }
    
    /**
     * Variable for finding the leading proton kinetic energy.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the kinetic energy of the leading muon.
    */
    template<class T>
        double leading_proton_ke(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            double energy(csda_ke(interaction.particles[i]));
            if constexpr (std::is_same_v<T, caf::SRInteractionTruthDLPProxy>)
                energy = ke_init(interaction.particles[i]);
            return energy;
        }

    /**
     * Variable for the transverse momentum of a particle.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the transverse momentum of the particle
    */
    template<class T>
        double transverse_momentum(const T & particle)
        {
            if constexpr (std::is_same_v<T, caf::SRParticleTruthDLPProxy>)
                return std::sqrt(std::pow(particle.truth_momentum[0], 2) + std::pow(particle.truth_momentum[1], 2));
            else
                return std::sqrt(std::pow(particle.momentum[0], 2) + std::pow(particle.momentum[1], 2));
        }
    
    /**
     * Variable for the transverse momentum of the leading muon.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the transverse momentum of the leading muon.
    */
    template<class T>
        double leading_muon_pt(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            return transverse_momentum(interaction.particles[i]);
        }

    /**
     * Variable for the transverse momentum of the leading proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the transverse momentum of the leading proton.
    */
    template<class T>
        double leading_proton_pt(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return transverse_momentum(interaction.particles[i]);
        }
    
    /**
     * Variable for the transverse momentum of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the transverse momentum of the primary particles.
    */
    template<class T>
        double interaction_pt(const T & interaction)
        {
            double px(0), py(0);
            for(const auto & p : interaction.particles)
                if(p.is_primary)
                {
                    px += p.momentum[0];
                    py += p.momentum[1];
                }
            return std::sqrt(std::pow(px, 2) + std::pow(py, 2));
        }

    /**
     * Variable for the cosine of the track angle within the XZ plane.
     * @tparam T the type of particle (true or reco)
     * @param particle to apply the variable on.
     * @return the cosine of the track angle within the XZ plane ()
    */
    template<class T>
        double cosine_theta_xz(const T & particle)
        {
            return particle.start_dir[2] / std::sqrt(std::pow(particle.start_dir[0], 2) + std::pow(particle.start_dir[2], 2));
        }

    /**
     * Variable for cosine theta_xz (transverse) of the leading muon.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the cosine theta_xz of the leading muon.
    */
    template<class T>
        double leading_muon_cosine_theta_xz(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 2));
            return cosine_theta_xz(interaction.particles[i]);
        }

    /**
     * Variable for cosine theta_xz (transverse) of the leading proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the cosine theta_xz of the leading proton.
    */
    template<class T>
        double leading_proton_cosine_theta_xz(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            return cosine_theta_xz(interaction.particles[i]);
        }
    
    /**
     * Variable for the cosine of the opening angle between leading muon and
     * proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the cosine of the opening angle between the leading muon and
     * proton.
    */
    template<class T>
        double cosine_opening_angle(const T & interaction)
        {
            auto & m(interaction.particles[leading_particle_index(interaction, 2)]);
            auto & p(interaction.particles[leading_particle_index(interaction, 4)]);
            double num(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
            return num;
        }

    /**
     * Variable for the cosine of the opening angle between leading muon and
     * proton in the transverse plane.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the cosine of the opening angle between the leading muon and
     * proton in the transverse plane.
    */
    template<class T>
        double cosine_opening_angle_transverse(const T & interaction)
        {
            auto & m(interaction.particles[leading_particle_index(interaction, 2)]);
            auto & p(interaction.particles[leading_particle_index(interaction, 4)]);
            double num(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1]);
            num /= std::sqrt((1-m.start_dir[2]*m.start_dir[2])*(1-p.start_dir[2]*p.start_dir[2]));
            return num;
        }

    /**
     * Variable for the softmax score of the leading muon.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the softmax score of the leading muon.
    */
    template<class T>
        double leading_muon_softmax(const T & interaction)
        {
            auto & m(interaction.particles[leading_particle_index(interaction, 2)]);
            return m.pid_scores[2];
        }
    
    /**
     * Variable for the softmax score of the leading proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the softmax score of the leading proton.
    */
    template<class T>
        double leading_proton_softmax(const T & interaction)
        {
            auto & p(interaction.particles[leading_particle_index(interaction, 4)]);
            return p.pid_scores[4];
        }

    /**
     * Variable for the cosine of the scattering angle between the primary
     * proton and any attached proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the cosine of the scattering angle between the primary proton
     * and any attached proton.
     * @note This is only valid for interactions with a primary proton.
    */
    template<class T>
        double proton_scattering_cosine(const T & interaction)
        {
            double cos(-2);
            size_t i(leading_particle_index(interaction, 4));
            if(interaction.particles.size() <= i)
                return cos;
            auto & p(interaction.particles[i]);
            if(p.pid == 4 && p.is_primary)
            {
                /**
                 * Check if the proton has a daughter proton by comparing the end
                 * of the primary proton to the start of each other particle.
                */
                for(const auto & part : interaction.particles)
                {
                    if(part.pid == 4 && !part.is_primary && part.start_point[0] == p.end_point[0] && part.start_point[1] == p.end_point[1] && part.start_point[2] == p.end_point[2])
                    {
                        double dot(p.truth_momentum[0] * part.truth_momentum[0] + p.truth_momentum[1] * part.truth_momentum[1] + p.truth_momentum[2] * part.truth_momentum[2]);
                        double mag1(std::sqrt(std::pow(p.truth_momentum[0], 2) + std::pow(p.truth_momentum[1], 2) + std::pow(p.truth_momentum[2], 2)));
                        double mag2(std::sqrt(std::pow(part.truth_momentum[0], 2) + std::pow(part.truth_momentum[1], 2) + std::pow(part.truth_momentum[2], 2)));
                        cos = dot / (mag1 * mag2);
                    }
                }

            }
            return cos;
      }
    
    /**
     * Variable for the overlap fraction of the leading proton.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return the overlap fraction of the leading proton.
    */
    template<class T>
        double leading_proton_overlap(const T & interaction)
        {
            size_t i(leading_particle_index(interaction, 4));
            /**
             * Check that the leading particle is actually a proton.
            */
            if(interaction.particles[i].pid == 4)
                return overlap(interaction.particles[i]);
            else
                return -1;
        }
}

#endif