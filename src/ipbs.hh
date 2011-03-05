void solver(NEWTON &newton, SLP &slp);
void save(const DGF &udgf, const U &u, const GV &gv, const std::string filename);

// Evaluate elliptic integrals of the 1st kind

double eval_elliptic(double dr, double dz, double rr_j);

//template <typename PositionVector>
//double compute_pbeq(const double &u, const PositionVector &r);

// container for commandline arguments
typedef struct {
	double alpha_sor; 
	int RefinementLevel;
	std::string GridName;
} Cmdparam;
