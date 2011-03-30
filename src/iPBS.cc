/** \file

    \brief IPBS - An Iterative Poisson Boltzmann implementation using DUNE

    Here goes some explanation on what is done :-)
    \todo { Doc Me ! }   
*/

// iPBS.cc Read a gmsh file and solve PB eq.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// std includes
#include<math.h>
#include<iostream>
#include<vector>
#include<string>

// DUNE includes
#include<dune/common/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/timer.hh>
// Input/Output
#include <dune/grid/io/file/gnuplot.hh>
#include<dune/grid/io/file/gmshreader.hh>
// we use UG
#include<dune/grid/uggrid.hh>
// pdelab includes
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
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

// global typedefs
typedef double Real;

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

#include "boundaries.hh"
#include "RefLocalOperator.hh"
#include "ipbs_P1.hh"
#include "ref_P1.hh"

/** \brief container for commandline arguments

    \todo { Implement user interface }
*/
typedef struct {
	double alpha_sor; 
	int RefinementLevel;
	std::string GridName;
} Cmdparam;




//===============================================================
// Main programm
//===============================================================
int main(int argc, char** argv)
{
 try{
  // Initialize Mpi
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
   if(Dune::MPIHelper::isFake)
    std::cout<< "This is (at the moment) a sequential program!" << std::endl;
  else
  {
    if(helper.rank()==0)
    {
       std::cout << "Hello World! This is iPBS." << std::endl;
       std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }
  }

  // check arguments
  if (argc!=4)
  {
    if (helper.rank()==0)
    {
	std::cout << "usage: ./iPBS <meshfile> <refinement level> <SOR Parameter>" << std::endl;
	return 1;
    }
  }

  // Read in comandline arguments
  Cmdparam cmdparam;
  cmdparam.GridName=argv[1];
  sscanf(argv[2],"%d",&cmdparam.RefinementLevel);
  double alpha;
  sscanf(argv[3],"%lf", &alpha);
  if(helper.rank()==0)
    std::cout << "Using " << cmdparam.RefinementLevel << " refinement levels. alpha = " << alpha << std::endl;
  sysParams.set_alpha(alpha);
  
  
//===============================================================
// Setup problem
//===============================================================
  
  
  // <<<1>>> Setup the problem from mesh file
  const int dimgrid = 2;
  typedef Dune::UGGrid<dimgrid> GridType;         // 2d mesh

  // instanciate ug grid object
  GridType grid;

  // define vectors to store boundary and element mapping
  std::vector<int> boundaryIndexToEntity;
  std::vector<int> elementIndexToEntity;

  // read a gmsh file
  Dune::GmshReader<GridType> gmshreader;
  //if (helper.rank()==0)
     gmshreader.read(grid, cmdparam.GridName, boundaryIndexToEntity, elementIndexToEntity, true, false);

  // refine grid
  grid.globalRefine(cmdparam.RefinementLevel);
  
  // get a grid view on the leaf grid
  grid.loadBalance();

  typedef GridType::LeafGridView GV;
  const GV& gv = grid.leafView();
 
  ref_P1(gv, elementIndexToEntity, boundaryIndexToEntity);
  ipbs_P1(gv, elementIndexToEntity, boundaryIndexToEntity);
  
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
