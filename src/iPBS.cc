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
typedef Dune::PDELab::P1LocalFiniteElementMap<Coord,Real,dim> FEM;
typedef Dune::PDELab::ConformingDirichletConstraints CON;
typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
typedef GFS::ConstraintsContainer<Real>::Type CC;
typedef GFS::VectorContainer<Real>::Type U;
typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;

#ifndef _P0LAYOUT_H
#define _P0LAYOUT_H
#include "p0layout.hh"
#endif

#include <dune/grid/common/mcmgmapper.hh>
typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType, P0Layout> Mapper;

#include "boundaries.hh"
typedef Regions<GV,double,std::vector<int>> M;
typedef BCType<GV,std::vector<int>> B;
typedef BCExtension_init<GV,double,std::vector<int>> G_init;
typedef BoundaryFlux<GV,double,std::vector<int> > J;
typedef BoundaryFluxRef<GV,double,std::vector<int> > J_ref;

#include"PB_operator.hh"
typedef PBLocalOperator<M,B,J> LOP;
typedef PBLocalOperatorRef<M,B,J_ref> LOP_ref;
typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE,true> GOS;
typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP_ref,CC,CC,MBE,true> GOS_ref;
typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
typedef Dune::PDELab::Newton<GOS,LS,U> NEWTON;
typedef Dune::PDELab::Newton<GOS_ref,LS,U> NEWTON_ref;
typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
typedef Dune::PDELab::StationaryLinearProblemSolver<GOS_ref,LS,U> SLP_ref;
typedef DGF::Traits Traits;

// include application heaeders
#include "ipbs.hh"
#include "solve.hh"

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif



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
	std::cout << "usage: ./iPBS <meshfile> <refinement level> <SOR Parameter>" << std::endl;
	return 1;
    }
  }

  // Read in comandline arguments
  Cmdparam cmdparam;
  cmdparam.GridName=argv[1];
  sscanf(argv[2],"%d",&cmdparam.RefinementLevel);
  //sscanf(argv[3],"%f",&cmdparam.alpha_sor);
  double alpha;
  sscanf(argv[3],"%lf", &alpha);
  std::cout << "Using " << cmdparam.RefinementLevel << " refinement levels. alpha" << alpha << std::endl;
  sysParams.set_alpha(alpha);
  
  
//===============================================================
// Setup problem
//===============================================================
  
  
  // <<<1>>> Setup the problem from mesh file

  // instanciate ug grid object
  GridType grid;

  // define vectors to store boundary and element mapping
  std::vector<int> boundaryIndexToEntity;
  std::vector<int> elementIndexToEntity;

  // read a gmsh file
  Dune::GmshReader<GridType> gmshreader;
  gmshreader.read(grid, cmdparam.GridName, boundaryIndexToEntity, elementIndexToEntity, true, false);

  // refine grid
  grid.globalRefine(cmdparam.RefinementLevel);

  const GV& gv = grid.leafView();

  //define boundaries
  
  // inner region
  M m(gv, elementIndexToEntity);
  // boundary
  B b(gv, boundaryIndexToEntity);
    // Dirichlet
  G_init g(gv, boundaryIndexToEntity);
  // boundary fluxes
  J j(gv, boundaryIndexToEntity);
  J_ref j_ref(gv, boundaryIndexToEntity);
  
  // Create finite element map
  FEM fem;

  // <<<2>>> Make grid function space
  GFS gfs(gv,fem);

  // Create coefficient vector (with zero values)
  U u(gfs,0.0);
  
  // <<<3>>> Make FE function extending Dirichlet boundary conditions
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc); 
  //std::cout << "constrained dofs=" << cc.size() 
  //          << " of " << gfs.globalSize() << std::endl;

  // interpolate coefficient vector
  Dune::PDELab::interpolate(g,gfs,u);
  
  
     // get the Neumann solution as reference
  
  // construct discrete grid function for access to solution
  U u_ref(gfs,0.0);
  const DGF udgf_ref(gfs, u_ref);
  
  LOP_ref lop_ref(m,b,j_ref);
  GOS_ref gos_ref(gfs,cc,gfs,cc,lop_ref);
  
  // <<<5a>>> Select a linear solver backend
  LS ls_ref(5000,true);
    
  // <<<5b>>> Instantiate solver for nonlinear problem
  NEWTON_ref newton_ref(gos_ref,u_ref,ls_ref); 
  newton_ref.setReassembleThreshold(0.0);
  newton_ref.setVerbosityLevel(1);
  newton_ref.setReduction(1e-10);
  newton_ref.setMinLinearReduction(1e-4);
  newton_ref.setMaxIterations(20);
  newton_ref.setLineSearchMaxIterations(10);
    
  // <<<5c>>> Instantiate Solver for linear problem
  SLP_ref slp_ref(gos_ref,u_ref,ls_ref,1e-10); 
  
  // <<<6>>> Solve Problem
  newton_ref.apply();
  slp_ref.apply();
  
  DGF udgf_refSave(gfs,u_ref);
  save(udgf_refSave, u_ref, gv, "reference");
  
  
  // Call function exectuting the iterative procedure
  get_solution(u, gv, gfs, cc, grid, m, b, j);
  
 
  
  U u_relError(gfs), u_absError(gfs);
  std::vector<double> tmp1, tmp2, relError, absError;
  u.std_copy_to(tmp1);
  u_ref.std_copy_to(tmp2);
  relError.resize(tmp1.size());
  absError.resize(tmp1.size());
  for (unsigned int i = 0; i < tmp1.size(); i++) {
    relError[i] = fabs((tmp1[i] - tmp2[i]) / tmp2[i]);
    absError[i] = fabs(tmp1[i] - tmp2[i]);
  }
  u_relError.std_copy_from(relError);
  u_absError.std_copy_from(absError);
  DGF udgf_relErrorSave(gfs,u_relError);
  DGF udgf_absErrorSave(gfs,u_absError);
  save(udgf_relErrorSave, u_relError, gv, "relError");
  save(udgf_absErrorSave, u_absError, gv, "absError");
  
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


