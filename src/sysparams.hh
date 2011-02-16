class SysParams {
   public:
	SysParams(double _lambda, double _bjerrum, int _charge, double _epsilon, double _radius);	// Constructor
	double get_bjerrum();
	bool init;
	double get_radius();
	double get_charge();
	double get_epsilon();
	double get_lambda2i();
	double get_sigma_sphere();
	double get_phi_init();
	void add_error(double);
	void reset_error();
	double get_error();
	void set_lambda (double value);
	void set_phi_init (double value);
	static const double pi = 3.14159265358979323846;
	double get_alpha();
	void set_alpha(double);
	static const double error_cut = 4.0E-3;
   private:
	double lambda;
	double lambda2i;
	double sigma_sphere;
	double phi_init;
	double bjerrum;
	int charge;
	double epsilon;
	double radius;
	double totalError;
	double oldValue;
	double alpha;	// SOR parameter
};

extern SysParams sysParams;
