/** \brief Driver for solving the IPBS problem using P1 piecewise linear Lagrange elements

    Dirichlet boundary conditions are used for domain (outer) boundaries, Neumann b.c. on
    symmetry axis and IPBS Neumann b.c. where we want the IPBS iterative procedure.
*/

/*!
   \param gv the view on the leaf grid
   \param elementIndexToEntity mapper defining the index of inner elements
   \param boundaryIndexToEntity mapper defining the index of boundary elements
*/

#ifndef _IPBSOLVER_H
#define _IPBSOLVER_H
#include "ipbsolver.hh"
#endif

template<class GridType>
void test_P1(GridType* grid, const std::vector<int>& elementIndexToEntity,
             const std::vector<int>& boundaryIndexToEntity,
             Dune::MPIHelper& helper)
{
  // We want to know the total calulation time
  Dune::Timer timer;
  timer.start();

  // get a grid view on the leaf grid
  typedef typename GridType::LeafGridView GV;
  const GV& gv = grid->leafView();
 // const IpbsGridView<GridType>& gv = static_cast<IpbsGridView<GridType> > (grid->leafView());

  // some typedef
  const int dim = GV::dimension;
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
  typedef Dune::PDELab::P1LocalFiniteElementMap<ctype,Real,dim> FEM;
  FEM fem;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;
  CON con;
  // Define ISTL Vector backend - template argument is the blocksize!
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem,con);

  // Compute Ghost Elements
  con.compute_ghosts(gfs);

  // Create coefficient vector (with zero values)
  typedef typename GFS::template VectorContainer<Real>::Type U;
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
  LOP lop(m,b,j);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE,true> GOS;
  GOS gos(gfs,cc,gfs,cc,lop);

  // <<<5a>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GOS,double> LS;
  //typedef Dune::PDELab::ISTLBackend_NOVLP_CG_NOPREC<GFS> LS;
  //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS;
  LS ls(gfs);

  // <<<5b>>> Solve nonlinear problem
  typedef Dune::PDELab::Newton<GOS,LS,U> NEWTON;
  NEWTON newton(gos,u,ls);
  newton.setLineSearchStrategy(newton.hackbuschReuskenAcceptBest);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(sysParams.get_verbose());
  newton.setReduction(sysParams.get_newton_tolerance());
  newton.setMinLinearReduction(5e-1); // seems to be low in parallel?
  newton.setMaxIterations(50);
  newton.setLineSearchMaxIterations(25);

  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  sysParams.counter = 0;
  
  double inittime = timer.elapsed();
  double solvertime = 0.;
  double itertime = 0.;

  // --- Here the iterative loop starts ---

  while (ipbs.next_step())
  {
    timer.reset();
    newton.apply();
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

    sysParams.counter ++;

    timer.reset();
    ipbs.updateBC(u);
    ipbs.updateIC();
    itertime += timer.elapsed();
 }

  // --- here the iterative loop ends! ---

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
  
  // Gnuplot output
  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
  gnuplotwriter.addVertexData(u,"solution");
  gnuplotwriter.write(filename); 

  // Calculate the forces
  ipbs.forces(u);
  ipbs.forces2(u);
  // ipbs.forces3(u);
  
  if (helper.rank() == 0) {
    std::ofstream runtime;
    runtime.open ("runtime.dat", std::ios::out | std::ios::app); 
    runtime << "P " << helper.size() << " N: " << elementIndexToEntity.size() << " M: " << ipbs.get_n() << " init: " << inittime << " solver: " << solvertime/sysParams.counter << " boundary update " << itertime/sysParams.counter << std::endl;
    runtime.close();
  }
}
