#include "boundary.hh"
#include <iostream>

void Boundary::set_epsilons(double epsilonIn, double epsilonOut)
{
    epsilon = epsilonIn;
    factor = 2. * epsilonIn / (epsilonIn + epsilonOut) ;
}

void Boundary::set_potential(double _value)
{
    potential = _value;
}

double Boundary::get_potential()
{
    return potential;
}

void Boundary::set_type(int _type)
{
    type = _type;
}

int Boundary::get_type()
{
    return type;
}
void Boundary::set_charge_density(double _value)
{
    charge_density = _value;
}

double Boundary::get_epsilon()
{
    return epsilon;
}

double Boundary::get_charge_density()
{
    return charge_density;
}

double Boundary::get_dielectric_factor()
{
  return factor;
}
