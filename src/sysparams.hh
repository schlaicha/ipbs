class SysParams {
   public:
	SysParams(double _lambda, double _bjerrum, int _charge, double _epsilon, double _radius);	// Constructor
	double get_bjerrum();
	double get_charge();
	double get_epsilon();
	double get_lambda2i();
	double get_sigma_init();
	double get_phi_init();
	void add_error(double);
	void reset_error();
	double get_error();
	void set_lambda (double value);
	void set_phi_init (double value);
	static const double alpha = 0.7;	// SOR parameter
   private:
	double lambda;
	double lambda2i;
	double sigma_init;
	double phi_init;
	double bjerrum;
	int charge;
	double epsilon;
	double radius;
	double totalError;
	double oldValue;
};

extern SysParams sysParams;
