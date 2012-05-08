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
#elif GRIDDIM == 3
#include <dune/pdelab/finiteelementmap/pk3dfem.hh>
#endif

// pdelab includes
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>

#if HAVE_MPI
#include <dune/pdelab/backend/novlpistlsolverbackend.hh>
#endif

#include <dune/ipbs/datawriter.hh>
#include <dune/ipbs/ipbsolver.hh>
#include <dune/ipbs/boundaries.hh>
#include <dune/ipbs/PBLocalOperator.hh>

// test some solvers
//#include<dune/pdelab/stationary/linearproblem.hh>
//#include <dune/pdelab/backend/seqistlsolverbackend.hh>

template<class GridType, int k>
void ipbs_Pk(GridType* grid, const std::vector<int>& elementIndexToEntity,
             const std::vector<int>& boundaryIndexToEntity,
             Dune::MPIHelper& helper)
{
  // We want to know the total calulation time
  Dune::Timer timer;
  timer.start();
  
  std::stringstream status;
  //std::ofstream status;
  //status.open ("status.dat", std::ios::out | std::ios::out); 

  // get a grid view on the leaf grid
  typedef typename GridType::LeafGridView GV;
  const GV& gv = grid->leafView();

  // some typedef
  typedef typename GV::Grid::ctype ctype;
  const int dim = GV::dimension;

  // <<<1>>> Setup the problem from mesh file

  // NOTE only flux b.c. are iterated, so rest is instantiated only once
  // inner region
  typedef Regions<GV,double,std::vector<int> > M;
  M m(gv, elementIndexToEntity);
  // boundary condition type
  typedef BCTypeParam<std::vector<int> > B;
  B b(boundaryIndexToEntity);
  // Class defining Dirichlet B.C.
  typedef BCExtension<GV,double,std::vector<int> > G;
  G g(gv, boundaryIndexToEntity);

  // Create finite element map
#if GRIDDIM == 2
  typedef Dune::PDELab::Pk2DLocalFiniteElementMap<GV, ctype, Real, k> FEM;
#elif GRIDDIM == 3
  typedef Dune::PDELab::Pk3DLocalFiniteElementMap<GV, ctype, Real, k> FEM;
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
  //ipbs.updateBC(u);

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
  typedef Dune::PDELab::ISTLBackend_NOVLP_CG_Jacobi< GFS > LS;
  //typedef Dune::PDELab::ISTLBackend_NOVLP_CG_SSORk< GO > LS;
  //typedef Dune::PDELab::ISTLBackend_NOVLP_CG_NOPREC<GFS> LS;
  //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS;
  //LS ls(gfs);
  LS ls(gfs, 20000, 1);
#else
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_CG_ILU0 LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SOR<GOS> LS;
  LS ls(5000, true);
#endif

  //typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
  //SLP slp(gos, u, ls, 1e-10);
  //slp.apply();

  // <<<5b>>> Solve nonlinear problem
  typedef Dune::PDELab::Newton<GO,LS,U> NEWTON;
  NEWTON newton(go,u,ls);
  newton.setLineSearchStrategy(newton.hackbuschReuskenAcceptBest);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(sysParams.get_verbose());
  newton.setReduction(sysParams.get_newton_tolerance());
  newton.setMinLinearReduction(1e-10); // seems to be low in parallel?
  newton.setMaxIterations(100);
  newton.setLineSearchMaxIterations(50);

  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  
  double inittime = timer.elapsed();
  double solvertime = 0.;
  double itertime = 0.;

  double fluxError, icError;
  int iterations = 0;

  DataWriter<GV,dim> mydatawriter(gv, helper);

  // --- Here the iterative loop starts ---

  do
  {
    timer.reset();
    try{
        newton.apply();
    }
    catch (Dune::Exception &e){
        status << "# Dune reported error: " << e << std::endl;
        std::cerr << "Dune reported error: " << e << std::endl;
        break;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
        break;
    }
    solvertime += timer.elapsed();

   // save snapshots of each iteration step
   //std::stringstream out;
   //out << "ipbs_step_" << iterations;
   //std::string filename = out.str();
   //DGF udgf_snapshot(gfs,u);
   //Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
   //vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf_snapshot,"solution"));
   //vtkwriter.write(filename,Dune::VTK::appendedraw);
   //mydatawriter.writeIpbsCellData(gfs, u, "solution", filename, status);

   timer.reset();
   ipbs.updateBC(u);
   ipbs.updateIC();
   itertime += timer.elapsed();
  }
  while (ipbs.next_step(fluxError,icError,iterations));

  // --- here the iterative loop ends! ---


  status << "# reached convergence criterion: " << std::boolalpha <<
    ipbs.next_step(fluxError, icError, iterations) << std::endl;
  status << "# in iteration " << iterations << std::endl
      << "# maximum relative change in boundary condition calculation is " <<  fluxError << std::endl
      << "# maximum relative change in induced charge density is " << icError << std::endl;
  //status.close();

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
  
  mydatawriter.writeIpbsCellData(gfs, u, "solution", "ipbs_solution", status);

  // Calculate the forces
  ipbs.forces(u);
//  ipbs.forces2(u);
  
  if (helper.rank() == 0) {
    std::cout << "P " << helper.size() << " N: " << elementIndexToEntity.size() << " M: " << ipbs.get_n() << " init: " << inittime << " solver: " << solvertime/iterations << " boundary update " << itertime/iterations << std::endl;
    //std::ofstream runtime;
    //runtime.open ("runtime.dat", std::ios::out | std::ios::app); 
    //runtime << "P " << helper.size() << " N: " << elementIndexToEntity.size() << " M: " << ipbs.get_n() << " init: " << inittime << " solver: " << solvertime/sysParams.counter << " boundary update " << itertime/sysParams.counter << std::endl;
    //runtime.close();
    
#ifdef SURFACE_POT
    std::cout << "Now I would calculate the surface potential :-)" << std::endl;
    typedef typename GV::template Codim<0>::template Partition
              <Dune::Interior_Partition>::Iterator LeafIterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    typedef typename DGF::Traits::RangeType RT;
    // Do the loop for all boundaryIDs > 1 (all colloids)
    for (int i = 2; i < sysParams.get_npart()+2; i++)
    {
      int nElems = 0;
      double sum = 0.;
      for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
               	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
      {
        if(it->hasBoundaryIntersections() == true) {
          for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii) {
            if(ii->boundary() == true) {
              if (boundaryIndexToEntity[ii->boundarySegmentIndex()] == i) // check if IPBS boundary
              {
                Dune::FieldVector<Real, dim> evalPos = ii->geometry().center();
                Dune::FieldVector<double,GFS::Traits::GridViewType::dimension> local =
                    it->geometry().local(evalPos);
                RT value;
                // evaluate the potential
                udgf.evaluate(*it, local, value);
                sum += value;
                nElems++;
              }
            }
          }
        }
      }
      sum /= nElems;
      boundary[i]->set_res_surface_pot(sum);
      std::cout << "Averaged surface potential: " << sum << std::endl;
    }
#endif
    
 }
}
