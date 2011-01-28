// std includes
#include<math.h>
#include<iostream>
#include<vector>
#include<string>

// dune includes
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/grid/io/file/gmshreader.hh>

// we use UG
#include<dune/grid/uggrid.hh>

// pdelab includes
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/instationary/onestep.hh>   // Filenamehelper
#include<dune/pdelab/finiteelementmap/p1fem.hh>	// P1 in 1,2,3 dimensions
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

const int dimgrid = 2;
typedef Dune::UGGrid<dimgrid> GridType; 	// 2d mesh
typedef GridType::LeafGridView GV;
typedef GV::Grid::ctype Coord;
typedef double Real;
const int dim = GV::dimension;
typedef Dune::PDELab::P1LocalFiniteElementMap<Coord,Real,dim> FEM;
typedef Dune::PDELab::ConformingDirichletConstraints CON; 
typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
typedef GFS::ConstraintsContainer<Real>::Type CC;
typedef GFS::VectorContainer<Real>::Type U;
typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
typedef Regions<GV,double,std::vector<int>> M;
typedef BCType<GV,std::vector<int>> B;
typedef BCExtension<GV,double,std::vector<int>> G;
typedef BoundaryFlux<GV,double,std::vector<int> > J;
typedef PBLocalOperator<M,B,J> LOP;
typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE,true> GOS;
typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
typedef Dune::PDELab::Newton<GOS,LS,U> NEWTON;
typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
typedef DGF::Traits Traits;

void solver(NEWTON &newton, SLP &slp);
void save(const DGF &udgf, const GV &gv);
void calculate_phi(const GV &gv, const DGF &udgf);

template <typename PositionVector>
double compute_pbeq(const double &u, const PositionVector &r);

// container for commandline arguments
typedef struct {
	int NewtonMaxIteration; 
	int RefinementLevel;
	std::string GridName;
} Cmdparam;
