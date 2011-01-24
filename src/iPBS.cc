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
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/gnuplot.hh>

// we use UG
#include<dune/grid/uggrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

// pdelab includes
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>	// P1 in 1,2,3 dimensions
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/instationarygridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

// include application heaeders
#include "extensions.hh"
#include"PB_operator.hh"
#include "solver.hh"
#include "parameters.hh"

//===============================================================
// Main program with grid setup
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

//===============================================================
// Setup the problem from mesh file
//===============================================================

  // instanciate ug grid object
  const int dimgrid = 2;
  typedef Dune::UGGrid<dimgrid> GridType; 	// 2d mesh
  GridType grid(400);		// heapSize: The size of UG's internal memory in megabytes for this grid. 

  // define vectors to store boundary and element mapping
  std::vector<int> boundaryIndexToEntity;
  std::vector<int> elementIndexToEntity;

  // read a gmsh file
  Dune::GmshReader<GridType> gmshreader;
  gmshreader.read(grid, cmdparam.GridName, boundaryIndexToEntity, elementIndexToEntity, true, false);

  // refine grid
  grid.globalRefine(cmdparam.RefinementLevel);

  // get a grid view
  typedef GridType::LeafGridView GV;
  const GV& gv = grid.leafView();

  // Visit boundary elements and say hello
  //boundary_info(grid);

  // inner region, i.e. solve
  typedef Regions<GV,double,std::vector<int>> M;
  M m(gv, elementIndexToEntity);

  // boundary condition
  typedef BCType<GV,std::vector<int>> B;
  B b(gv, boundaryIndexToEntity);
  typedef BCExtension<GV,double,std::vector<int>> G;
  G g(gv, boundaryIndexToEntity);

  // boundary fluxes
  typedef BoundaryFlux<GV,double,std::vector<int> > J;
  J j(gv, boundaryIndexToEntity);

  // call the problem solver
  solver(gv, m, b, g, j, cmdparam);
 
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


