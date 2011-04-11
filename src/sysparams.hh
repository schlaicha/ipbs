/** \file
    \brief A class containing all information for the systems setup.
    \todo Doc me!
*/

#include <string>

class SysParams {
 public:
	SysParams();	// We now only need a default constructor

	int counter;
	static const double pi = 3.14159265358979323846;
  
  // Return parameters needed during runtime
	double get_bjerrum();
  double get_radius();
	double get_epsilon();
	double get_lambda2i();
	double get_charge_density();
	double get_sphere_pos();
	double get_error();
	double get_alpha();
 	int get_symmetry();
  int get_refinement();
  std::string get_meshfile();

  // Functions setting the private members
  void add_error(double);
	void reset_error();
	void set_lambda (double value);
	void set_bjerrum (double value);
	void set_alpha(double);
  void set_meshfile(std::string filename);
  void set_refinement(int level);
  void set_radius(double value);
	
  private:
	double lambda;
	double lambda2i;
	double bjerrum;
	double charge_density;
	double charge;
	double epsilon;
	double radius;  // TODO: still needed?
	double totalError;
	int symmetry;
	double alpha;	// SOR parameter
	double pos;
  std::string meshfile;
  int refinement;
};

// Global access to this class via global instance
// TODO is there a better way for parallelization?
extern SysParams sysParams;
