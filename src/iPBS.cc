// iPBS.cc Read a gmsh file and solve PB eq.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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

// Input/Output
#include <dune/grid/io/file/gnuplot.hh>

// we use UG
#include<dune/grid/uggrid.hh>

// pdelab includes
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
// #include<dune/pdelab/instationary/onestep.hh>   // Filenamehelper
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
typedef Dune::UGGrid<dimgrid> GridType;         // 2d mesh
typedef GridType::LeafGridView GV;
typedef GV::Grid::ctype Coord;
typedef double Real;
const int dim = GV::dimension;
typedef Dune::PDELab::P1LocalFiniteElementMap<Coord,Real,dim> FEM; // DEPRECATED
typedef Dune::PDELab::ConformingDirichletConstraints CON;
typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
typedef GFS::ConstraintsContainer<Real>::Type CC;
typedef GFS::VectorContainer<Real>::Type U;
typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

#include "boundaries.hh"
typedef Regions<GV,double,std::vector<int>> M;
typedef BCType<GV,std::vector<int>> B;
typedef BCExtension_init<GV,double,std::vector<int>> G_init;
typedef BCExtension_iterate<GV,double,std::vector<int>> G;
typedef BoundaryFlux<GV,double,std::vector<int> > J;

#include"PB_operator.hh"
typedef PBLocalOperator<M,B,J> LOP;
typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE,true> GOS;
typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
typedef Dune::PDELab::Newton<GOS,LS,U> NEWTON;
typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
typedef DGF::Traits Traits;

// include application heaeders
#include "ipbs.hh"

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

#include "solve.hh"

//===============================================================
// Main programm
//===============================================================
int main(int argc, char** argv)
{
 try{
  // Initialize Mpi
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  std::cout << "Hello World! This is iPBS." << std::endl;
  if(Dune::MPIHelper::isFake)
    std::cout<< "This is (at the moment) a sequential program!" << std::endl;
  else
  {
    if(helper.rank()==0)
    std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
  }

  // check arguments
  if (argc!=4)
  {
    if (helper.rank()==0)
    {
	std::cout << "usage: ./iPBS <meshfile> <refinement level> <MaxNewtonIterations>" << std::endl;
	return 1;
    }
  }

  // Read in comandline arguments
  Cmdparam cmdparam;
  cmdparam.GridName=argv[1];
  sscanf(argv[2],"%d",&cmdparam.RefinementLevel);
  sscanf(argv[3],"%d",&cmdparam.NewtonMaxIteration);
  std::cout << "Using " << cmdparam.RefinementLevel << " refinement levels." << std::endl;

  
  // <<<1>>> Setup the problem from mesh file

  // instanciate ug grid object
  GridType grid(400);		// heapSize: The size of UG's internal memory in megabytes for this grid.
				// This Constructor is DEPRECATED

  // define vectors to store boundary and element mapping
  std::vector<int> boundaryIndexToEntity;
  std::vector<int> elementIndexToEntity;

  // read a gmsh file
  Dune::GmshReader<GridType> gmshreader;
  gmshreader.read(grid, cmdparam.GridName, boundaryIndexToEntity, elementIndexToEntity, true, false);

  // refine grid
  grid.globalRefine(cmdparam.RefinementLevel);

  // get a grid view
  const GV& gv = grid.leafView();

  //define boundaries
  // inner region
  M m(gv, elementIndexToEntity);
  // boundary
  B b(gv, boundaryIndexToEntity);

  // Create finite element map
  FEM fem;

  // <<<2>>> Make grid function space
  GFS gfs(gv,fem);

  // Create coefficient vector (with zero values)
  U u(gfs,0.0);

  // Get initial solution
  // Dirichlet
  G_init g_init(gv, boundaryIndexToEntity);
  // boundary fluxes
  J j(gv, boundaryIndexToEntity);
  // get initial coefficient vector
  get_solution(u, gv, gfs, m, b, g_init, j);

  // Create Discrete Grid Function Space
  DGF udgf(gfs,u);
  // graphical output
  std::string vtk_filename = "step_0";
  save(udgf, u, gv, vtk_filename);

  // HERE ITERATION STARTS HERE
  int iterationCounter = 1;
  while (sysParams.get_error() > 0.01 || iterationCounter == 1)
  {
    std::cout << std::endl << "IN ITERATION " << iterationCounter <<"\t actual error is: " << sysParams.get_error() << std::endl << std::endl;
    // Create DGF for constructor of G
    DGF udgf_it(gfs,u);
    // Reset error for new iteration
    sysParams.reset_error();
    // use updated Dirichlet
    G g(gv, boundaryIndexToEntity, udgf_it);
    get_solution(u, gv, gfs, m, b, g, j);
    std::stringstream out;
    out << "step_" << iterationCounter;
    vtk_filename = out.str();
    save(udgf_it, u, gv, vtk_filename);
    ++iterationCounter;
  }

  // done
  return 0;
 }
 catch (Dune::Exception &e){
  std::cerr << "Dune reported error: " << e << std::endl;
 }
 catch (...){
  std::cerr << "Unknown exception thrown!" << std::endl;
 }
} 


// ============================================================================

void calculate_phi(const GV &gv, const DGF &udgf)
{
  // Get the potential at the particle's surface
  // std::cout << "Particle Boundaries at:" << std::endl;
  // define iterators
  typedef GV::Codim<0>::Iterator ElementLeafIterator;
  typedef GV::IntersectionIterator IntersectionIterator;
  int iicount=0;
  double phi1=0;
  // double phi2=0;
 for (ElementLeafIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
  {
    //if (it->hasBoundaryIntersections()==true)
    //{
	//for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
	//{
		Traits::DomainType xlocal;
		Traits::RangeType y1;
		Traits::RangeType y2;
		//if (ii->boundary()==true)
		//{
		   Dune::FieldVector<double,dim> x1 = it->geometry().global(xlocal);
		   //Dune::FieldVector<double,dim> x2 = it->geometry().center();
		   udgf.evaluate(*it,x1,y1);
		   //udgf.evaluate(*it,x2,y2);
		   //if (x1.two_norm() < 4.7)				
		   //{    
		     ++iicount;
		     std::cout << "x= " << x1[0] << "\t y = " << x1[1] << "\t Phi = " << y1 << std::endl;
		     //std::cout << "Boundary \t x= " << x2[0] << "\t y = " << x2[1] << "\t Phi = " << y2 << std::endl;
		     phi1 += y1;
		     //phi2 += y2;
		   //}
		//}
	//}
    //}
  }
  std::cout << std::endl << "Phi_init = " << sysParams.get_phi_init() << "\t Phi_S = " << phi1/iicount << std::endl;
//   << "\t Phi_S (Boundary) = " << phi2/iicount << std::endl;
} 

// ============================================================================

template <typename PositionVector>
double compute_pbeq(const double &u, const PositionVector &r)
{
	return (- sysParams.get_lambda2i() * std::sinh(u));
}

void solver (NEWTON &newton, SLP &slp)
{
	newton.apply();
	slp.apply();
}

// ============================================================================

void save(const DGF &udgf, const U &u, const GV &gv, const std::string filename)
{
  //Dune::PDELab::FilenameHelper fn("output");
  // {
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
    vtkwriter.write(filename,Dune::VTK::appendedraw);
  // }
  // Gnuplot output
  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
  gnuplotwriter.addVertexData(u,"solution");
  gnuplotwriter.write(filename + ".dat");
}


