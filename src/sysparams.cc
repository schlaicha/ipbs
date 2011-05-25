#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

// Implementation of SysParams functionality

// Constructor
// Debye and Bjerrum length in [nm], colloid charge in [e], radius in Bjerrum length
SysParams::SysParams()
{
	totalError = 1E8;
}

double SysParams::get_r2i()
{
  return r2i;
}

void SysParams::set_symmetry(int value)
{
  // Set the systems symmetry
  // 1 is "2D_cylinder"
  // 2 is "2D_sphere"
  // 3 is "3D"  - not verified!
  symmetry = value;
}

void SysParams::set_boxLength(double value)
{
  boxLength = value;
}

void SysParams::set_tolerance(double value)
{
  tolerance = value;
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
  r2i = 1.0 / (radius * radius);
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

void SysParams::reset_error(double error)
{
  totalError = error;
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

double SysParams::get_tolerance()
{
	return tolerance;
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

double SysParams::get_boxLength()
{
    return boxLength;
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
