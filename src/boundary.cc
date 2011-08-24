#include "boundary.hh"

void Boundary::set_isPlane(bool value)
{
  isPlane = value;
}

bool Boundary::get_isPlane()
{
  return isPlane;
}
void Boundary::set_epsilons(double epsilonIn, double epsilonOut)
{
    factor = 2. * epsilonIn / (epsilonIn + epsilonOut) ;
}

void Boundary::set_charge_density(double _value)
{
    charge_density = _value;
}

double Boundary::get_charge_density()
{
    return charge_density;
}

double Boundary::get_dielectric_factor()
{
  return factor;
}
