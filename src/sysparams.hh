#include<math.h>

class SysParams {
   public:
	SysParams(double _lambda, double _bjerrum, double _charge_density, double _epsilon, double _radius);	// Constructor
	double get_bjerrum();
	int counter;
	double get_radius();
	double get_charge();
	double get_epsilon();
	double get_lambda2i();
	double get_charge_density();
	double get_E_init();
	void add_error(double);
	void reset_error();
	double get_error();
	void set_lambda (double value);
	void set_E_init (double value);
	static const double pi = 3.14159265358979323846;
	double get_alpha();
	void set_alpha(double);
	static const double error_cut = 4.0E-3;
	
	template <typename Iterator>
	double compute_pbeq(const double &u, const Iterator &it)
	{
	  //if (it->position().two_norm() < 1.01)
	  //  return - 1.0 * (lambda2i * sinh(u) - (charge * bjerrum / (radius * radius)));
	  //else
	    return - 1.0 * (lambda2i * sinh(u));
	}
	
	
   private:
	double lambda;
	double lambda2i;
	double charge_density;
	double E_init;
	double bjerrum;
	int charge;
	double epsilon;
	double radius;
	double totalError;
	double oldValue;
	double alpha;	// SOR parameter
};

extern SysParams sysParams;
