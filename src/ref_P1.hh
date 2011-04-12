/** \brief Driver for solving the Poisson Boltzmann Eq. with Neumann b.c.
           as a reference solution to IPBS using P1 piecewise linear Lagrange elements

    Dirichlet boundary conditions are used for domain (outer) boundaries, Neumann b.c. on
    symmetry axis and IPBS Neumann b.c. which match the total charge of the particle.
*/

/*!
   \param gv the view on the leaf grid
   \param elementIndexToEntity mapper defining the index of inner elements
   \param boundaryIndexToEntity mapper defining the index of boundary elements
*/

template<class GV, typename Factory>
void ref_P1(const GV& gv, const std::vector<int>& elementIndexToEntity,
             const std::vector<int>& boundaryIndexToEntity,
             const Factory& factory)
{
  // We want to know the total calulation time
  Dune::Timer timer;

  // some typedef
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype ctype;

  // <<<1>>> Setup the problem from mesh file

  // inner region
  typedef Regions<GV,double,std::vector<int>> M;
  M m(gv, elementIndexToEntity);
  // boundary conditioon type
  typedef BCType<GV,std::vector<int>, Factory> B;
  B b(gv, boundaryIndexToEntity, factory);
  // Class defining Dirichlet B.C.
  typedef BCExtension<GV,double,std::vector<int>> G;
  G g(gv, boundaryIndexToEntity);
  // boundary fluxes - this one is for the reference solution!
  typedef RefBoundaryFlux<GV,double,std::vector<int>, Factory> J;
  J j(gv, boundaryIndexToEntity, factory);

  // Create finite element map
  typedef Dune::PDELab::P1LocalFiniteElementMap<ctype,Real,dim> FEM;
  FEM fem;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  CON con;
  // Define ISTL Vector backend - template argument is the blocksize!
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem,con);

  // Create coefficient vector (with zero values)
  typedef typename GFS::template VectorContainer<Real>::Type U;
  U u(gfs,0.0);
  
  // <<<3>>> Make FE function extending Dirichlet boundary conditions
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;

  Dune::PDELab::constraints(b,gfs,cc); 
  //   std::cout << "constrained dofs=" << cc.size() 
  //             << " of " << gfs.globalSize() << std::endl;

  // interpolate coefficient vector
  Dune::PDELab::interpolate(g,gfs,u);

  // <<<4>>> Make Grid Operator Space
  typedef RefLocalOperator<M,B,J> LOP;
  LOP lop(m,b,j);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE,true> GOS;
  GOS gos(gfs,cc,gfs,cc,lop);

  // <<<5a>>> Select a linear solver backend
  // typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
  // LS ls(gfs,cc,5000,5,1);
  // typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GFS> LS;
  // LS ls(gfs);
  typedef Dune::PDELab::ISTLBackend_OVLP_CG_SSORk<GFS, CC> LS;
  LS ls(gfs,cc);

  // <<<5b>>> Solve nonlinear problem
  typedef Dune::PDELab::Newton<GOS,LS,U> NEWTON;
  NEWTON newton(gos,u,ls);
  newton.setLineSearchStrategy(newton.hackbuschReuskenAcceptBest);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(0);
  newton.setReduction(1e-10);
  newton.setMinLinearReduction(1e-4);
  newton.setMaxIterations(25);
  newton.setLineSearchMaxIterations(10);
  newton.apply();

  // <<<6>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
  vtkwriter.write("reference",Dune::VTK::appendedraw);

  // Gnuplot output
  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
  gnuplotwriter.addVertexData(u,"solution");
  gnuplotwriter.write("reference.dat"); 


  std::cout << "Reference total calculation time=" << timer.elapsed() << std::endl;
}
