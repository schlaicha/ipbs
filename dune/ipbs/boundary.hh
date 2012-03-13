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
    double get_charge_density();
    double get_dielectric_factor();
    double get_epsilon();
    void set_res_surface_pot(double _value);
    double get_res_surface_pot();
    void set_isPlane(bool value);
    bool get_isPlane();

  private:
    double factor;
    double epsilon;
    double charge;
    double charge_density;
    double res_surface_pot;
    bool isPlane;
};

#endif  // _BOUNDARY_H
