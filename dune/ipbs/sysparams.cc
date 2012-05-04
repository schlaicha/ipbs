#include "sysparams.hh"

// Implementation of SysParams functionality

/** Constructor
/  Debye and Bjerrum length in [nm], colloid charge in [e], radius in Bjerrum length */

SysParams::SysParams()
{
	totalError = 1E8;
  epsilon  = 1.;
}

int SysParams::get_npart()
{
  return npart;
}

void SysParams::set_npart(int value)
{
  npart = value;
}

int SysParams::get_salt()
{
  return salt;
}

void SysParams::set_salt(int value)
{
  // Define if we use symmetric salt (case 0) or counterions only (case 1)
  salt = value;
}

void SysParams::set_symmetry(int value)
{
  // Set the systems symmetry
  // 1 is "2D_cylinder"
  // 2 is "2D_sphere"
  // 3 is "3D"  - not verified!
  symmetry = value;
}

void SysParams::set_maxiter(int _maxiter)
{
  maxiter = _maxiter;
}

int SysParams::get_maxiter()
{
  return maxiter;
}

int SysParams::get_refinementSteps()
{
  return refinementSteps;
}

void SysParams::set_refinementSteps(int value)
{
  refinementSteps = value;
}

double SysParams::get_refinementFraction()
{
  return refinementFraction;
}

void SysParams::set_refinementFraction(double value)
{
  refinementFraction = value;
}

void SysParams::set_boxLength(double value)
{
  boxLength = value;
}

void SysParams::set_tolerance(double value)
{
  tolerance = value;
}

void SysParams::set_newton_tolerance(double value)
{
  newton_tolerance = value;
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

void SysParams::add_error(double error)
{
	if (error > totalError)
	{
	  totalError = error;
	}
}

void SysParams::set_alpha_ipbs(double alpha_in)
{
	alpha_ipbs=alpha_in;
}

double SysParams::get_alpha_ipbs()
{
	return alpha_ipbs;
}

void SysParams::set_alpha_ic(double alpha_in)
{
  alpha_ic=alpha_in;
}

double SysParams::get_alpha_ic()
{
	return alpha_ic;
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

double SysParams::get_newton_tolerance()
{
	return newton_tolerance;
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

void SysParams::set_pH(double pH_) {
    pH=pH_;
}

double SysParams::get_pH() 
{
    return pH;
}
