/** \file
    \brief A class containing all information for the systems setup.
    \todo Doc me!
*/

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H

#include <string>

class SysParams {
  public:
    SysParams();	// We now only need a default constructor

	const double pi;
  
    // Return parameters needed during runtime
	double get_bjerrum();
	double get_epsilon();
	double get_lambda2i();
	double get_lambda();
	double get_charge_density();
	double get_error();
	double get_alpha_ic();
	double get_alpha_ipbs();
	unsigned int get_maxiter();
	double get_tolerance();
 	int get_symmetry();
    int get_verbose();
    int get_salt();
    int get_outStep();
    size_t get_npart();
    std::string get_meshfile();
    std::string get_outname();
    double get_pH();
    double get_integration_d();
    double get_integration_l();

    // Functions setting the private members
    void add_error(double);
    void reset_error();
    void reset_error(double error);
    void set_lambda (double value);
    void set_bjerrum (double value);
    void set_alpha_ic(double);
    void set_alpha_ipbs(double);
    void set_maxiter(unsigned int);
    void set_meshfile(std::string filename);
    void set_outname(std::string filename);
    void set_refinement(int level);
    void set_tolerance(double value);
    void set_verbose(int value);
    void set_charge_density (double value);
    void set_symmetry (int value);
    void set_boxLength (double value);
    void set_refinementFraction (double value);
    void set_refinementSteps (int value);
    void set_salt (int value);
    void set_outStep (int value);
    void set_npart (size_t value);
    void set_pH(double value);
    void set_integration_d(double value);
    void set_integration_l(double value);
	
  private:
	double lambda;
	double lambda2i;
	double bjerrum;
	double charge_density;
	double charge;
	double epsilon;
    double tolerance;
    double totalError;
	int symmetry;
    int npart;
    int salt;
    int outstep;
	double alpha_ic;	// SOR parameter
	double alpha_ipbs;	// SOR parameter
    int maxiter;
    std::string meshfile;
    std::string outname;
    int verbose;
    double integration_d;
    double integration_l;
    double refinementFraction;
    int refinementSteps;
    double pH;
};

#endif  // _SYSPARAMS_H
