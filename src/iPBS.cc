// meshtest.cc Read a gmsh file, and solve PB eq.
// adapted from dunepdelab-howto/doc/howto/src_examples/
// cadsample.cc 

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
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

// pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>	// We use piecewise linear Lagrange elements
#include<dune/pdelab/finiteelementmap/p1fem.hh>	// P1 in 1,2,3 dimensions
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

// include application heaeders
//#include"cadsample_operator.hh"
//#include"cadsample_parameter.hh"
//#include"cadsample_P1.hh"


//===============================================================
// Setup the problem from mesh file
//===============================================================

void setup_mesh(std::string gridName, int level)
{
  // instanciate ug grid object
  typedef Dune::UGGrid<2> GridType; 	// 2d mesh
  GridType grid(400);			// what's this number good for?

  // vectors for boundary and material conditions
  std::vector<int> boundaryIndexToPhysicalEntity;
  std::vector<int> elementIndexToPhysicalEntity;

  // read a gmsh file
  //std::string gridName = "meshtest.msh";
  Dune::GmshReader<GridType> gmshreader;
  gmshreader.read(grid, gridName, boundaryIndexToPhysicalEntity,
      elementIndexToPhysicalEntity, true, false);

  // refine grid
  grid.globalRefine(level);

  // get a grid view
  typedef GridType::LeafGridView GV;
  const GV& gv = grid.leafView();

  // material conditions
  //typedef CrankDiffusion<GV,double,std::vector<int> > M;
  //M m(gv, elementIndexToPhysicalEntity);
  
  // boundary conditions
  //typedef CrankBCType<GV,std::vector<int> > B;
  //B b(gv, boundaryIndexToPhysicalEntity);
  //typedef CrankBCExtension<GV,double,std::vector<int> > G;
  //G g(gv, boundaryIndexToPhysicalEntity);

  // boundaries fluxes
  //typedef CrankFlux<GV,double,std::vector<int> > J;
  //J j(gv, boundaryIndexToPhysicalEntity);

  // call driver with parameters
  //gridName.erase(0, gridName.rfind("/")+1);
  //gridName.erase(gridName.find(".",0), gridName.length());
  //cadsample_P1(gv, m, b, g, j, gridName);
}

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
 try{
  // Initialize Mpi
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  std::cout << "Hello World! This is meshtest." << std::endl;
  if(Dune::MPIHelper::isFake)
    std::cout<< "This is not a sequential program!" << std::endl;
  else
  {
    if(helper.rank()==0)
    std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
  }

  // check arguments
  if (argc!=3)
  {
    std::cout << "usage: ./cadsample <meshfile> <level>" << std::endl;
    return 1;
  }

  // scan mesh name
  std::string gridName = argv[1];

  // refinement level
  int level = 0;
  sscanf(argv[2],"%d",&level);
  std::cout << "Using " << level << " refinement levels." << std::endl;

  // run simulation
  // setup_mesh(gridName, level);
  
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
