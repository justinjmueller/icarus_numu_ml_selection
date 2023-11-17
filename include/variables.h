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
     * Variable for counting interactions.
     * @tparam T the type of interaction (true or reco).
     * @param interaction to apply the variable on.
     * @return 1.0 (always).
     */
    template<class T>
        double count(const T & interaction) { return 1.0; }

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
}

#endif