/** \file

    \brief IPBS - An Iterative Poisson Boltzmann implementation using DUNE

    Here goes some explanation on what is done :-)
    \todo { Doc Me ! }   
*/

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
// #include<dune/common/collectivecommunication.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/timer.hh>

// // Adaptivity
// #include <dune/pdelab/adaptivity/adapt.hh>

// Single Geometry Single Codim Mapper
#include <dune/grid/common/scsgmapper.hh>
// quadrature
#include<dune/grid/common/quadraturerules.hh>
// Input/Output
#include <dune/grid/io/file/gnuplot.hh>
#include<dune/grid/io/file/gmshreader.hh>

// we use UG
#ifdef UGGRID
  #include<dune/grid/uggrid.hh>
  #include<dune/grid/uggrid/uggridfactory.hh>
#endif
#ifdef ALUGRID
  #include<dune/grid/alugrid.hh>
  #include<dune/grid/alugrid/2d/alu2dgridfactory.hh>
#endif
#if !(UGGRID || ALUGRID)
  #error It looks like dunecontrol could not detect your UG installation properly.
  #error At the moment, iPBS *STRICTLY* depends on UG!
  #error Compilation will be aborted.
#endif

// pdelab includes
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/finiteelementmap/p1fem.hh>	// P1 in 1,2,3 dimensions
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>	// P1 in 1,2,3 dimensions
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/common/dynmatrix.hh>

// global typedefs
typedef double Real;

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif
#ifndef _PARTICLE_H
#define _PARTICLE_H
#include "boundary.hh"
#endif

// global access to particles
std::vector<Boundary*> boundary;

//#include "ipbsgridview.hh"

#include "parser.hh"
#include "functors.hh"
#include "integrateentity.hh"
//#include "eval_elliptic.hh"
#include "boundaries.hh"
//#include "gradient.hh"
//#include "maxwelltensor.hh"
//#include "force.hh"
//#include "ipbs_boundary.hh"
#include "PBLocalOperator.hh"
//#include "ipbs_ref_P1.hh"
//#include "ipbs_prepare.hh"
//#include "ipbs_P1.hh"
//#include "ipbs_P2.hh"
//#include "ref_P1.hh"
#include "test_driver.hh"
#include "test_P2.hh"

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
  if (argc!=2)
  {
    if (helper.rank()==0)
    {
	std::cout << "usage: ./iPBS <configuration file>" << std::endl;
	return 1;
    }
  }
  
  // Parse configuration file.
  std::string config_file(argv[1]);
  parser(config_file);

  
//===============================================================
// Setup problem
//===============================================================
  
  
  // <<<1>>> Setup the problem from mesh file
  const int dimgrid = 2;         // 2d mesh
  
  // define vectors to store boundary and element mapping
  std::vector<int> boundaryIndexToEntity;
  std::vector<int> elementIndexToEntity;
  
#ifdef UGGRID
  typedef Dune::UGGrid<dimgrid> GridType;
#elif ALUGRID
  typedef Dune::ALUConformGrid< dimgrid, dimgrid > GridType;
#endif
  Dune::GridFactory<GridType> factory;

 
  if(helper.rank() == 0)
  {
    // read a gmsh file
    Dune::GmshReader<GridType> gmshreader;
    gmshreader.read(factory, sysParams.get_meshfile(), boundaryIndexToEntity, elementIndexToEntity, true, true);
  }

 // Setup Dune Collective Communication
 Dune::CollectiveCommunication<MPI_Comm> collCom(helper.getCommunicator());

 // Communicate boundary vector

 int size = boundaryIndexToEntity.size();
 collCom.broadcast (&size, 1, 0);
 if (helper.rank() > 0)
   boundaryIndexToEntity.reserve(size);
 collCom.broadcast(&boundaryIndexToEntity[0],size,0);

 
 // create the grid
 GridType* grid = factory.createGrid();

 // refine grid
  if(helper.rank()==0) {
    std::cout << "Using " << sysParams.get_refinement() << "global refinement steps and" << std::endl;
    std::cout << sysParams.get_refinementSteps() << " adaptive refinement steps with "
      << sysParams.get_refinementFraction() << " percent refinement." << std::endl;
  }

 grid->globalRefine(sysParams.get_refinement());
 
 grid->loadBalance();

 // Call problem drivers
 // ref_P1(grid, elementIndexToEntity, boundaryIndexToEntity, collCom);
 // ipbs_P1(grid, elementIndexToEntity, boundaryIndexToEntity, helper);
 // ipbs_P2(grid, elementIndexToEntity, boundaryIndexToEntity, collCom);
 // test_P1(grid, elementIndexToEntity, boundaryIndexToEntity, helper);
 test_P2(grid, elementIndexToEntity, boundaryIndexToEntity, helper);
  
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
