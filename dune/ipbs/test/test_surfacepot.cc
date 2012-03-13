#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// global DUNE includes
#include<dune/common/mpihelper.hh>
#include<dune/common/collectivecommunication.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/timer.hh>

/* include grid IO */
#include<dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/gridtype.hh>

// global typedefs
typedef double Real;

#include <dune/ipbs/ipbs.hh>
#include <dune/ipbs/sysparams.hh>
#include <dune/ipbs/boundary.hh>
#include <dune/ipbs/ipbs_Pk.hh>

// global access to particles and parameters
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
 
    // set the symmetry of the system
  sysParams.set_symmetry(2);
  
  // Parse other options
  sysParams.set_maxiter(100);
  sysParams.set_alpha_ipbs(0.67);
  sysParams.set_alpha_ic(0.2);
  sysParams.set_newton_tolerance(1e-10);
  sysParams.set_bjerrum(0.71);
  sysParams.set_lambda(1);
  sysParams.set_tolerance(1e-5);
  sysParams.set_verbose(4);
  sysParams.set_salt(0);
  double epsilonOut = 1.;

  // Create particles
  int n_particle = 1;
  sysParams.set_npart(n_particle);
  for (int i = 0; i < n_particle; i++)
    boundary.push_back(new Boundary());

  for(int i = 0; i < n_particle; i++)
  {
    boundary[i]->set_charge_density(1e-3);
    double epsilonIn = 1.; 
    boundary[i]->set_epsilons(epsilonIn, epsilonOut);
  }

  
//===============================================================
// Setup problem
//===============================================================
  
  
  // <<<1>>> Setup the problem from mesh file
  
  // define vectors to store boundary and element mapping
  std::vector<int> boundaryIndexToEntity;
  std::vector<int> elementIndexToEntity;
  
  typedef Dune::GridSelector::GridType GridType;
  Dune::GridFactory<GridType> factory;

  if(helper.rank() == 0)
  {
    // read a gmsh file
    Dune::GmshReader<GridType> gmshreader;
#if GRIDDIM == 2
    gmshreader.read(factory, "sphere2d.msh", boundaryIndexToEntity, elementIndexToEntity, true, true);
#endif
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
 
  // refine grid
  if(helper.rank()==0) {
    std::cout << "Using " << sysParams.get_refinement() << " global refinement steps and" << std::endl;
    std::cout << sysParams.get_refinementSteps() << " adaptive refinement steps with "
      << sysParams.get_refinementFraction() << " percent refinement." << std::endl;
  }

  // Load balance the parallel grid
  std::cout << "Grid has been modified by load balancing: " << grid->loadBalance() << std::endl;

 // Call problem driver
 ipbs_Pk<GridType, 1>(grid, elementIndexToEntity, boundaryIndexToEntity, helper);
 
 double dhsurfacepot = 4.*sysParams.pi*sysParams.get_bjerrum()*(10./(1.+10.*1.));
 for (int i = 0; i < sysParams.get_npart(); i++)
 {
    std::cout << "DH solution:\t" << dhsurfacepot << "\tipbs solution: " <<
        boundary[i]->get_res_surface_pot() << std::endl;
 } 
 // done
 return 0;
 }

 catch (Dune::Exception &e){
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
 }
 catch (...){
  std::cerr << "Unknown exception thrown!" << std::endl;
  return 2;
 }
 
}
// ============================================================================
