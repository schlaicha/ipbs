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

// we use UG
#include<dune/grid/uggrid.hh>

// include application heaeders
#include "ipbs.hh"
#include"PB_operator.hh"
#include "solver.hh"
#include "parameters.hh"

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

template <typename PositionVector>
//double compute_pbeq(const double &u, Dune::FieldVector<double, dim> &r)
double compute_pbeq(const double &u, const PositionVector &r)
{
	return (- sysParams.get_lambda2i() * std::sinh(u));
}


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


  sysParams.set_lambda(2.0);

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

//===============================================================
// Solve
//===============================================================

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


