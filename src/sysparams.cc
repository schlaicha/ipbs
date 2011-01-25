#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

// Implementation of SysParams functionality

// Constructor
// Debye and Bjerrum length in [nm], colloid charge in [e], radius in Bjerrum length
SysParams::SysParams(double _lambda=1, double _bjerrum=0.7, int _charge=100, double _epsilon=80.0, double _radius=2.0):lambda(_lambda),bjerrum(_bjerrum),charge(_charge),epsilon(_epsilon),radius(_radius)
{
	lambda2i = 1 / (lambda * lambda);
	sigma_init = bjerrum * charge / (epsilon * radius * radius);
}

double SysParams::get_lambda2i()
{
	return lambda2i;
}

double SysParams::get_sigma_init()
{
	return sigma_init;
}

void SysParams::set_lambda(double value)
{
	lambda = value;
	lambda2i = 1 / (lambda * lambda);
}

SysParams sysParams;
