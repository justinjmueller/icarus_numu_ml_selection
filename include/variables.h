/**
 * @file variables.h
 * @brief Header file for definitions of selection variables.
 * @author justin.mueller@colostate.edu
*/
#ifndef VARIABLES_H
#define VARIABLES_H

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
     * Variable for particle PID.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the PID of the particle.
    */
    template<class T>
        double pid(const T & particle) { return particle.pid; }

    /**
     * Variable for particle csda_ke.
     * @tparam T the type of particle (true or reco).
     * @param particle to apply the variable on.
     * @return the csda_ke of the particle.
    */
    template<class T>
        double csda_ke(const T & particle) { return particle.csda_ke; }
}

#endif