// container for commandline arguments

typedef struct {
	int NewtonMaxIteration; 
	int RefinementLevel;
	std::string GridName;
} Cmdparam;

template <typename PositionVector>
double compute_pbeq(const double &u, const PositionVector &r);
