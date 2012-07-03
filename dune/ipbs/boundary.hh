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
    void set_sigma_max(double);
    void set_pK(double);
    void set_Y(double);
    void set_ifShift(bool);

    double get_charge_density();
    double get_dielectric_factor();
    double get_epsilon();
    double get_potential();
    double get_sigma_max();
    double get_pK();
    double get_Y();

    bool doShift();

    int get_type();

  private:
    double factor;
    double epsilon;
    double charge_density;
    int type;
    double potential;
    
    double sigma_max;
    double pK;
    double Y;

    bool shift;
};

#endif  // _BOUNDARY_H
