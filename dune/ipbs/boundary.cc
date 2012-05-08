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

double Boundary::get_sigma_max() 
{
    return sigma_max;
}

void Boundary::set_sigma_max(double sigma_max_) 
{
    sigma_max = sigma_max_;
}

double Boundary::get_Y()
{
    return Y;
}

void Boundary::set_Y(double Y_) 
{
    Y=Y_;
}

double Boundary::get_pK()
{
    return pK;
}

void Boundary::set_pK(double pK_)
{
    pK = pK_;
}
