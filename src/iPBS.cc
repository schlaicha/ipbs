// iPBS.cc Read a gmsh file and solve PB eq.
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
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>

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
#include"PB_operator.hh"
//#include"bcextension.hh"
//#include"bctype.hh"
#include "solver.hh"


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
    std::cout<< "This is not a sequential program!" << std::endl;
  else
  {
    if(helper.rank()==0)
    std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
  }

  // check arguments
  if (argc!=3)
  {
    std::cout << "usage: ./iPBS <meshfile> <refinement level>" << std::endl;
    return 1;
  }

  // scan mesh name
  std::string gridName = argv[1];

  // refinement level
  int level = 0;
  sscanf(argv[2],"%d",&level);
  std::cout << "Using " << level << " refinement levels." << std::endl;

//===============================================================
// Setup the problem from mesh file
//===============================================================

  // instanciate ug grid object
  typedef Dune::UGGrid<2> GridType; 	// 2d mesh
  //GridType grid();            // Standard constructor reserver 200MB - but
                                //gmshreader.read can not be found
  GridType grid(400);			// what's this number good for?
                                // heapSize: The size of UG's internal memory in megabytes for this grid. 

  // read a gmsh file
  Dune::GmshReader<GridType> gmshreader;
  gmshreader.read(grid, gridName, true, false);

  // refine grid
  grid.globalRefine(level);

  // get a grid view
  typedef GridType::LeafGridView GV;
  const GV& gv = grid.leafView();

  // call the problem solver
  solver(gv, gridName, level);

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
