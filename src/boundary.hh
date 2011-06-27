/** \file
    \brief A class containing all information about the particles.
    \todo Doc me!
*/

#include <vector>

class Boundary
{
  public:
    void set_radius(double _radius);
    void set_charge_density(double _value);
    
    double get_charge_density();

  private:
    double espilon;
    double radius;
    double charge;
    double charge_density;
};
