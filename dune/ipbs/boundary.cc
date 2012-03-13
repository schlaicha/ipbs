#include "boundary.hh"
#include <iostream>

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
    epsilon = epsilonIn;
    factor = 2. * epsilonIn / (epsilonIn + epsilonOut) ;
    // std::cout << "Epsilon_material = " << epsilonIn << " espilon_solution = " << epsilonOut <<" Dielectric factor is " << factor << std::endl;
}

void Boundary::set_res_surface_pot(double _value)
{
    res_surface_pot = _value;
}

double Boundary::get_res_surface_pot()
{
    return res_surface_pot;
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
