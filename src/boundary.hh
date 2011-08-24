/** \file
    \brief A class containing all information about the particles.
    \todo Doc me!
*/

#include <vector>

class Boundary
{
  public:
    void set_epsilons(double _epsilonIn, double _epsilonOut);
    void set_charge_density(double _value);
    double get_charge_density();
    double get_dielectric_factor();
    void set_isPlane(bool value);
    bool get_isPlane();

  private:
    double factor;
    double charge;
    double charge_density;
    bool isPlane;
};
