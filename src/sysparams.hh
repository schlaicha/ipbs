#include<math.h>

class SysParams {
   public:
	SysParams(double _lambda, double _bjerrum, double _charge_density, double _epsilon, double _radius);	// Constructor
	double get_bjerrum();
	double get_sphere_pos();
	int counter;
	double get_radius();
	double get_epsilon();
	double get_lambda2i();
	double get_charge_density();
	void add_error(double);
	void reset_error();
	double get_error();
	void set_lambda (double value);
	static const double pi = 3.14159265358979323846;
	double get_alpha();
	void set_alpha(double);
    int get_symmetry();
	
    // TODO: Change this!
	template <typename Iterator>
	double compute_pbeq(const double &u, const Iterator &it)
	{
	    return - 1.0 * (lambda2i * sinh(u));
	}
	
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
};

extern SysParams sysParams;
