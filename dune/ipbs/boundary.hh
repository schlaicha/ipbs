/** \file
    \brief A class containing all information about the particles.
    \todo Doc me!
*/

#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#include <vector>

class Boundary
{
  public:
    void set_epsilons(double _epsilonIn, double _epsilonOut);
    void set_charge_density(double _value);
    void set_type (int);
    void set_potential(double);
    double get_charge_density();
    double get_dielectric_factor();
    double get_epsilon();
    double get_potential();
    int get_type();

  private:
    double factor;
    double epsilon;
    double charge_density;
    int type;
    double potential;
};

#endif  // _BOUNDARY_H
