class SysParams {
   public:
	SysParams(double _lambda, double _bjerrum, int _charge, double _epsilon, double _radius);	// Constructor
	double get_lambda2i();
	double get_sigma_init();
	void set_lambda (double value);
   private:
	double lambda;
	double lambda2i;
	double sigma_init;
	double bjerrum;
	int charge;
	double epsilon;
	double radius;
};

extern SysParams sysParams;
