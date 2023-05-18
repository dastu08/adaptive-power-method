#ifndef _RNG_H_
#define _RNG_H_

#include "general.h"

namespace apm {

/* Random number generator object interface

Initializes the underlying random number generator when the object is created.
Easy interface for getting different kind of random numbers: real, interger...
*/
class Rng {
   private:
    unsigned int seed;

   public:
    /* Initialize the RNG.

    **Parameter**
        - logging: flag for controlling if the constructer logs the seed.

    **Description**
         Use the current time since the epoch as a seed.
    */
    Rng(bool logging = false);

    /* Get the seed used.

    **Return**
        The integer seed that was used during initialization.
    */
    unsigned int getSeed() {
        return this->seed;
    }

    /* Get a random number between 0 and 1.

    **Return**
        A random number between 0 and 1.
    */
    double real();

    /* Get a random integer between 0 and n-1.

    **Parameter**
        - n: number of possible choices (0 to n-1)

    **Return**
        Index of the choice that was choosen. Value is from 0 to n - 1.

    **Description**
        Generate a random integer and return it modulo n.
    */
    int integer(int n);

    /* Choose one possibility with certain probability.

    **Parameter**
        - probabilities: vector of probabilities for each choice.
        Need to sum to 1.

    **Return**
        Index of the choice that was choosen. Value is from 0 to
        `probabilities.size() - 1`.

    **Description**
        Generate a random number between 0 and 1 and search in which interval
        determined by the probabilities this number is. Return the index of the
        interval.
    */
    int choice(vDouble &probabilities);
};

}  // namespace apm

#endif  // _RNG_H_