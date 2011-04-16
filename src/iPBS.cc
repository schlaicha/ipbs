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
#include<dune/common/collectivecommunication.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/timer.hh>
// Global Universal Mapper
// #include <dune/grid/common/universalmapper.hh>
// Single Geometry Single Codim Mapper
#include <dune/grid/common/scsgmapper.hh>
// Multiple Geometry Multiple Codim Mapper
#include <dune/grid/common/mcmgmapper.hh>
// quadrature
#include<dune/grid/common/quadraturerules.hh>
// Input/Output
#include <dune/grid/io/file/gnuplot.hh>
#include<dune/grid/io/file/gmshreader.hh>
// we use UG
#include<dune/grid/uggrid.hh>
#include<dune/grid/uggrid/uggridfactory.hh>
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

#include <dune/grid/common/gridenums.hh>

// global typedefs
typedef double Real;

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

#include "parser.hh"
#include "functors.hh"
#include "integrateentity.hh"
#include "eval_elliptic.hh"
#include "boundaries.hh"
#include "gradient.hh"
#include "ipbs_boundary.hh"
#include "RefLocalOperator.hh"
#include "PBLocalOperator.hh"
#include "ipbs_P1.hh"
#include "ref_P1.hh"


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
  
  typedef Dune::UGGrid<dimgrid> GridType;
  Dune::GridFactory<GridType> factory;

 
  if(helper.rank() == 0)
  {
    // read a gmsh file
    Dune::GmshReader<GridType> gmshreader;
    gmshreader.read(factory, sysParams.get_meshfile(), boundaryIndexToEntity, elementIndexToEntity, true, true);
  }

  // for (int i=0;i<elementIndexToEntity.size();i++)
  //   std::cout << boundaryIndexToEntity[i] << std::endl;

  // Setup Dune Collective Communication
  Dune::CollectiveCommunication<MPI_Comm> collCom(helper.getCommunicator());

  // Communicate boundary vector
  int size = boundaryIndexToEntity.size();
  collCom.broadcast (&size, 1, 0);
  int* boundaryIndexToEntity_carray = (int*) malloc(size*sizeof(int));
  if (sysParams.get_verbose() > 4)
    std::cout << "size is now " << size << "on node " << helper.rank() << "\n";
  if (helper.rank() == 0 ) {
    for (int i =0; i<size; i++)
      boundaryIndexToEntity_carray[i]=boundaryIndexToEntity[i];
  }
  collCom.broadcast(boundaryIndexToEntity_carray,size,0);
  if (sysParams.get_verbose() > 4)
    std::cout << "array bcasted " << "on node " << helper.rank() << "\n";
  if (helper.rank() != 0) {
    for (int i =0; i<size; i++)
      boundaryIndexToEntity.push_back(boundaryIndexToEntity_carray[i]);
    if (sysParams.get_verbose() > 4)
      std::cout << "vector was created " << "on node " << helper.rank() << "\n";
 }
 free(boundaryIndexToEntity_carray);

 // Communicate element vector
 size = elementIndexToEntity.size();
 collCom.broadcast (&size, 1, 0);
 int* elementIndexToEntity_carray = (int*) malloc(size*sizeof(int));
 if (sysParams.get_verbose() > 4)
   std::cout << "size is now " << size << "on node " << helper.rank() << "\n";
 if (helper.rank() == 0 ) {
   for (int i =0; i<size; i++)
     elementIndexToEntity_carray[i]=elementIndexToEntity[i];
 }
 collCom.broadcast(elementIndexToEntity_carray,size,0);
 if (sysParams.get_verbose() > 4)
    std::cout << "array bcasted " << "on node " << helper.rank() << "\n";
 if (helper.rank() != 0) {
   for (int i =0; i<size; i++)
    elementIndexToEntity.push_back(elementIndexToEntity_carray[i]);
    if (sysParams.get_verbose() > 4)
      std::cout << "vector was created " << "on node " << helper.rank() << "\n";
 } 
 free(elementIndexToEntity_carray);
  
 // create the grid
 GridType* grid = factory.createGrid();

 // refine grid
 if(helper.rank()==0)
   std::cout << "Using " << sysParams.get_refinement() << " refinement levels." << std::endl;
 grid->globalRefine(sysParams.get_refinement());
 
 grid->loadBalance();

 // get a grid view on the leaf grid
 typedef GridType::LeafGridView GV;
 const GV& gv = grid->leafView();

 // Call problem drivers
 ref_P1(gv, elementIndexToEntity, boundaryIndexToEntity);
 ipbs_P1(gv, elementIndexToEntity, boundaryIndexToEntity, collCom);
  
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
