/** \file

    \brief IPBS - An Iterative Poisson Boltzmann implementation using DUNE

    Here goes some explanation on what is done :-)
    \todo { Doc Me ! }   
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define BCGS_SSORk    1
#define BCGS_NOPREC   2
#define CG_SSORk      3
#define CG_NOPREC     4
#define CG_Jacobi     5
#define CG_AMG_SSOR   6
#define BCGS_AMG_SSOR 7

// default values
#ifndef PDEGREE 
#define PDEGREE 1
#endif
#ifndef LINEARSOLVER
#define LINEARSOLVER BCGS_SSORk
#endif

// std includes
//#include<math.h>
//#include<iostream>
//#include<vector>
//#include<string>

// global DUNE includes
#include<dune/common/mpihelper.hh>
#include<dune/common/collectivecommunication.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/timer.hh>

//TODO: sort out
//#include <dune/pdelab/adaptivity/adapt.hh>
//#include <dune/grid/common/gridenums.hh>
//#include <dune/common/dynmatrix.hh>
//#include<dune/pdelab/gridfunctionspace/constraints.hh>
//#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
//#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
//#include<dune/pdelab/finiteelementmap/p1fem.hh>	// P1 in 1,2,3 dimensions

/* include grid IO */
#include<dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/gridtype.hh>

// pdelab includes
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>


// global typedefs
typedef double Real;

#include "dune/ipbs/sysparams.hh"
#include "dune/ipbs/boundary.hh"
#include "dune/ipbs/parser.hh"
#include "dune/ipbs/ipbs_Pk.hh"

// global access to particles
std::vector<Boundary*> boundary;
SysParams sysParams;

//===============================================================
// Main programm
//===============================================================
int main(int argc, char** argv)
{
 try{
  // Initialize Mpi
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
   if(Dune::MPIHelper::isFake)
    std::cout<< "This is a sequential program!" << std::endl;
  else
  {
    if(helper.rank()==0)
    {
       std::cout << "Hello World! This is IPBS." << std::endl;
       std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }
  }
  
  // check arguments
  if (argc!=2)
  {
    if (helper.rank()==0)
    {
	std::cout << "usage: ./ipbs <configuration file>" << std::endl;
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
  
  // define vectors to store boundary and element mapping
  typedef std::vector<int> PGMap;
  PGMap boundaryIndexToEntity;
  PGMap elementIndexToEntity;
  
  typedef Dune::GridSelector::GridType GridType;
  Dune::GridFactory<GridType> factory;

  if(helper.rank() == 0)
  {
    // read a gmsh file
    Dune::GmshReader<GridType> gmshreader;
    gmshreader.read(factory, sysParams.get_meshfile(), boundaryIndexToEntity, elementIndexToEntity, true, true);
  }

  // MPIHelper ensures that this works for the sequential case
  Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> colCom(helper.getCommunicator());
 
  // Communicate boundary vector
  int size = boundaryIndexToEntity.size();
  colCom.broadcast (&size, 1, 0);
  if (helper.rank() > 0)
    boundaryIndexToEntity.reserve(size);
  colCom.broadcast(&boundaryIndexToEntity[0],size,0);
  
  // create the grid
  GridType* grid = factory.createGrid();
 
//  // Load balance the parallel grid
//  std::cout << "Grid has been modified by load balancing: " << grid->loadBalance() << std::endl;
  
  typedef GridType::LeafGridView myGV;

  myGV gv = grid->leafView();

  typedef typename myGV::Codim<0>::Iterator GVIT;
  typedef typename myGV::IntersectionIterator IntersectionIterator;
  
  GVIT it = gv.begin<0>();

//  for (; it != gv.end<0>(); ++it) {
//    for (IntersectionIterator iit = gv.ibegin(*it); iit!=  gv.iend(*it); ++iit) 
//      if (iit->boundary())
//        std::cout << "" << iit->geometry().center() << " " << iit->boundarySegmentIndex() 
////            << " " << factory.insertionIndex(*iit) 
////            << " " << boundaryIndexToEntity[factory.insertionIndex(*iit)] 
//            << " " << boundaryIndexToEntity[iit->boundarySegmentIndex()]  << std::endl ;
//    }
  

 // Call problem driver
 ipbs_Pk<GridType, PGMap, PDEGREE>(grid, elementIndexToEntity, boundaryIndexToEntity);
  }
  
 // done
 catch (Dune::Exception &e){
  std::cerr << "Dune reported error: " << e << std::endl;
 }
 catch (...){
  std::cerr << "Unknown exception thrown!" << std::endl;
 }
} 

// ============================================================================
