#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

// Implementation of SysParams functionality

// Constructor
// Debye and Bjerrum length in [nm], colloid charge in [e], radius in Bjerrum length
// SysParams::SysParams(double _lambda=1.0, double _bjerrum=0.7, 
//                      double _charge=55, double _epsilon=80.0, 
//                      double _radius=1.0)
//                     :lambda(_lambda),bjerrum(_bjerrum),charge(_charge),
//                      epsilon(_epsilon),radius(_radius)
SysParams::SysParams()
{
	totalError = 1E8;
    // Set the systems symmetry
    // 1 is "2D_cylinder"
    // 2 is "2D_sphere"
    // 3 is "3D"  - not verified!
    
    //charge_density = charge / (4 * pi * radius * radius);
}

void SysParams::set_symmetry(int value)
{
  symmetry = value;
}

void SysParams::set_charge_density(double value)
{
  charge_density = value;
}

void SysParams::set_verbose(int value)
{
  verbose = value;
}

int SysParams::get_verbose()
{
  return verbose;
}

void SysParams::set_bjerrum(double value)
{
  bjerrum = value;
}

void SysParams::set_radius(double value)
{
  radius = value;
}

double SysParams::get_radius()
{
	return radius;
}

double SysParams::get_sphere_pos()
{
	return pos;
}

void SysParams::add_error(double error)
{
	if (error > totalError)
	{
	  totalError = error;
	}
}

void SysParams::set_alpha(double alpha_in)
{
	alpha=alpha_in;
}

double SysParams::get_alpha()
{
	return alpha;
}

void SysParams::reset_error()
{
	totalError = 0;
}

double SysParams::get_error()
{
		return totalError;
}

double SysParams::get_epsilon()
{
	return epsilon;
}

double SysParams::get_bjerrum()
{
	return bjerrum;
}

double SysParams::get_lambda2i()
{
	return lambda2i;
}

double SysParams::get_charge_density()
{
	return charge_density;
}

void SysParams::set_lambda(double value)
{
	lambda = value;
	lambda2i = 1 / (lambda * lambda);
}

int SysParams::get_symmetry()
{
    return symmetry;
}

void SysParams::set_meshfile(std::string filename)
{
  meshfile = filename;
}

std::string SysParams::get_meshfile()
{
  return meshfile;
}

void SysParams::set_refinement(int level)
{
  refinement = level;
}

int SysParams::get_refinement()
{
  return refinement;
}

SysParams sysParams;
