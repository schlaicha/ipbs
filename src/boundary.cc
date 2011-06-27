#include "boundary.hh"

void Boundary::set_radius(double _radius)
{
    radius = _radius;
}

void Boundary::set_charge_density(double _value)
{
    charge_density = _value;
}

double Boundary::get_charge_density()
{
    return charge_density;
}
