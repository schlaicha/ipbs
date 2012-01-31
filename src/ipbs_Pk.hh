/** \brief Driver for solving the IPBS problem using P1 piecewise linear Lagrange elements

    Dirichlet boundary conditions are used for domain (outer) boundaries, Neumann b.c. on
    symmetry axis and IPBS Neumann b.c. where we want the IPBS iterative procedure.
*/

/*!
   \param gv the view on the leaf grid
   \param elementIndexToEntity mapper defining the index of inner elements
   \param boundaryIndexToEntity mapper defining the index of boundary elements
*/

#include <dune/pdelab/gridoperator/gridoperator.hh>
#if GRIDDIM == 2
#include<dune/pdelab/finiteelementmap/pk2dfem.hh>	// Pk in 2 dimensions
#endif
#include <dune/grid/io/file/gnuplot.hh>

//#include "../dune/iPBS/datawriter.hh"

#include "ipbsolver.hh"
#include "boundaries.hh"
#include "PBLocalOperator.hh"

template<class GridType, int k>
void ipbs_Pk(GridType* grid, const std::vector<int>& elementIndexToEntity,
             const std::vector<int>& boundaryIndexToEntity,
             Dune::MPIHelper& helper)
{
  // We want to know the total calulation time
  Dune::Timer timer;
  timer.start();
  
  std::ofstream status;
  status.open ("status.dat", std::ios::out | std::ios::out); 

  // get a grid view on the leaf grid
  typedef typename GridType::LeafGridView GV;
  const GV& gv = grid->leafView();

  // some typedef
  typedef typename GV::Grid::ctype ctype;

  // <<<1>>> Setup the problem from mesh file

  // NOTE only flux b.c. are iterated, so rest is instantiated only once
  // inner region
  typedef Regions<GV,double,std::vector<int> > M;
  M m(gv, elementIndexToEntity);
  // boundary condition type
  typedef BCType<GV,std::vector<int> > B;
  B b(gv, boundaryIndexToEntity);
  // Class defining Dirichlet B.C.
  typedef BCExtension<GV,double,std::vector<int> > G;
  G g(gv, boundaryIndexToEntity);

  // Create finite element map
#if GRIDDIM == 2
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, ctype, Real, k> FEM;
#endif
  FEM fem(gv);

  // <<<2>>> Make grid function space
  
  // Define ISTL Vector backend - template argument is the blocksize!
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  
  // Setup constraints
#if HAVE_MPI  // UG and AluGrid are nonoverlapping grids
  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;
#else
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
#endif
  CON con;  // initialize constraints

  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem,con);

#if HAVE_MPI  // Compute Ghost Elements
    con.compute_ghosts(gfs);
#endif

  // Create coefficient vector (with zero values)
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;
  U u(gfs,0.0);
  
  // <<<3>>> Make FE function extending Dirichlet boundary conditions
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc); 
  std::cout << "constrained dofs=" << cc.size() 
            << " of " << gfs.globalSize() << std::endl;

  // interpolate coefficient vector
  // this has to be done only during intitialization as we don't change dirichlet b.c.
  Dune::PDELab::interpolate(g,gfs,u);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,u);

  typedef Ipbsolver<GV, GFS> Ipbs;
  Ipbs ipbs(gv, gfs, helper, boundaryIndexToEntity);
  // instanciate boundary fluxes
  typedef BoundaryFlux<GV,double,std::vector<int>, Ipbs > J;
  J j(gv, boundaryIndexToEntity, ipbs);
  ipbs.updateBC(u);

  // <<<4>>> Make Grid Operator Space
  typedef PBLocalOperator<M,B,J> LOP;
  LOP lop(m,b,j,k+1);   // integration order
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
#if HAVE_MPI    // enable overlapping mode
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,
                                    Real,Real,Real,CC,CC,true> GO;
#else
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,
                                    Real,Real,Real,CC,CC> GO;
#endif
  GO go(gfs,cc,gfs,cc,lop);

  // <<<5a>>> Select a linear solver backend
#if HAVE_MPI
  //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GO> LS;
  typedef Dune::PDELab::ISTLBackend_NOVLP_CG_NOPREC<GFS> LS;
  //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS;
  LS ls(gfs);
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000, true);
#endif

  // <<<5b>>> Solve nonlinear problem
  typedef Dune::PDELab::Newton<GO,LS,U> NEWTON;
  NEWTON newton(go,u,ls);
  newton.setLineSearchStrategy(newton.hackbuschReuskenAcceptBest);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(sysParams.get_verbose());
  newton.setReduction(sysParams.get_newton_tolerance());
  newton.setMinLinearReduction(5e-7); // seems to be low in parallel?
  newton.setMaxIterations(20);
  newton.setLineSearchMaxIterations(10);

  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  
  double inittime = timer.elapsed();
  double solvertime = 0.;
  double itertime = 0.;

  // --- Here the iterative loop starts ---

  while (ipbs.next_step())
  {
    timer.reset();
    try{
        newton.apply();
    }
    catch (Dune::Exception &e){
        status << "Dune reported error: " << e << std::endl;
        std::cerr << "Dune reported error: " << e << std::endl;
        break;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
        break;
    }
    solvertime += timer.elapsed();

//   // save snapshots of each iteration step
//   std::stringstream out;
//   out << "ipbs_step_" << sysParams.counter;
//   std::string filename = out.str();
//   DGF udgf_snapshot(gfs,u);
//   Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
//   vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf_snapshot,"solution"));
//   vtkwriter.write(filename,Dune::VTK::appendedraw);
//   // Prepare filename for sequential Gnuplot output
//   if(helper.size()>1)
//   {
//     std::stringstream s;
//     s << 's' << std::setw(4) << std::setfill('0') << helper.size() << ':';
//     s << 'p' << std::setw(4) << std::setfill('0') << helper.rank() << ':';
//     s << filename;
//     filename = s.str();
//   }
//   // Gnuplot output
//   Dune::GnuplotWriter<GV> gnuplotwriter(gv);
//   gnuplotwriter.addVertexData(u,"solution");
//   gnuplotwriter.write(filename + ".dat"); 

    timer.reset();
    ipbs.updateBC(u);
    ipbs.updateIC();
    itertime += timer.elapsed();
 }

  // --- here the iterative loop ends! ---


  double fluxError, icError;
  int iterations;
  status << "reached convergence criterion: " << ipbs.next_step(fluxError, icError, iterations) << std::endl;
  status << "in iteration " << iterations << std::endl
      << "maximum relative change in boundary condition calculation is " << std::endl << fluxError << std::endl
      << "maximum relative change in induced charge density is " << std::endl << icError << std::endl;
  status.close();

  // <<<6>>> graphical output
  DGF udgf(gfs,u);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
  vtkwriter.write("ipbs_solution",Dune::VTK::appendedraw);

  // Prepare filename for sequential Gnuplot output
  std::string filename = "ipbs_solution";
  std::ostringstream s;
  if(helper.size() > 1)
  {
    s << 's' << std::setw(4) << std::setfill('0') << helper.size() << ':';
    s << 'p' << std::setw(4) << std::setfill('0') << helper.rank() << ':';
  }
  s << filename << ".dat";
  filename = s.str();
  
//  MyDataWriter<GV> mydatawriter(gv);
//  mydatawriter.addVertexData(u, "solution");
//  mydatawriter.write(filename);
  // Gnuplot output  - not for higher order elements
  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
  gnuplotwriter.addVertexData(u,"solution");
  gnuplotwriter.write(filename); 

  // Calculate the forces
  ipbs.forces(u);
//  ipbs.forces2(u);
  
  if (helper.rank() == 0) {
    std::cout << "P " << helper.size() << " N: " << elementIndexToEntity.size() << " M: " << ipbs.get_n() << " init: " << inittime << " solver: " << solvertime/sysParams.counter << " boundary update " << itertime/sysParams.counter << std::endl;
    //std::ofstream runtime;
    //runtime.open ("runtime.dat", std::ios::out | std::ios::app); 
    //runtime << "P " << helper.size() << " N: " << elementIndexToEntity.size() << " M: " << ipbs.get_n() << " init: " << inittime << " solver: " << solvertime/sysParams.counter << " boundary update " << itertime/sysParams.counter << std::endl;
    //runtime.close();
 }
}
