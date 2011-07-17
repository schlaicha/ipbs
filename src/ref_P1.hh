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

#include "test.hh"

template<class GridType, class ColCom>
void ref_P1(GridType* grid, const std::vector<int>& elementIndexToEntity,
             const std::vector<int>& boundaryIndexToEntity,
             const ColCom& colCom)
{
  // We want to know the total calulation time
  Dune::Timer timer;

  // get a grid view on the leaf grid
  typedef typename GridType::LeafGridView GV;
  const GV& gv = grid->leafView();

  // some typedef
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype ctype;

  // <<<1>>> Setup the problem from mesh file

  // inner region
  typedef Regions<GV,double,std::vector<int>> M;
  M m(gv, elementIndexToEntity);
  // boundary conditioon type
  typedef BCType<GV,std::vector<int> > B;
  B b(gv, boundaryIndexToEntity);
  // Class defining Dirichlet B.C.
  typedef BCExtension<GV,double,std::vector<int>> G;
  G g(gv, boundaryIndexToEntity);
  // boundary fluxes - this one is for the reference solution!
  typedef RefBoundaryFlux<GV,double,std::vector<int> > J;
  J j(gv, boundaryIndexToEntity);

  // Create finite element map
  typedef Dune::PDELab::P1LocalFiniteElementMap<ctype,Real,dim> FEM;
  FEM fem;

// adaptive refinement loop
for (int step = 0; step <= sysParams.get_refinementSteps(); step++)
{
  std::cout << "Refinement step " << step << std::endl;

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
  //   std::cout << "constrained dofs=" << cc.size() 
  //             << " of " << gfs.globalSize() << std::endl;

  // interpolate coefficient vector
  Dune::PDELab::interpolate(g,gfs,u);

  // <<<4>>> Make Grid Operator Space
  typedef PBLocalOperator<M,B,J> LOP;
  LOP lop(m,b,j);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE,true> GOS;
  GOS gos(gfs,cc,gfs,cc,lop);

  // <<<5a>>> Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GOS,double> LS;
  LS ls(gfs);

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

  if (step > 0) {
    // compute the estimated error using GradientSmoothnessOperator and refine
    typedef Dune::PDELab::ResidualErrorEstimation<GFS,U,Dune::PDELab::
      GradientSmoothnessOperator,true> GradientErrorEstimator;
    GradientErrorEstimator gradientErrorEstimator(gfs);
    typedef Dune::PDELab::EstimationAdaptation<GridType,GFS,U,GradientErrorEstimator> 
      EstimationAdaptor;
    EstimationAdaptor estimationAdaptor(*grid, gfs, gradientErrorEstimator, sysParams.get_refinementFraction());
    typedef Dune::PDELab::L2Projection<GFS,U> L2projection;
    L2projection l2projection(2);
    typedef typename Dune::PDELab::GridAdaptor<GridType,GFS,U,EstimationAdaptor,L2projection> GridAdaptor;
    GridAdaptor gridAdaptor(*grid, gfs, estimationAdaptor, l2projection);
    gridAdaptor.adapt(u);
  }

  std::string filename = "reference";
  std::ostringstream n;
  n << filename << "_step_" << step;
  filename = n.str();
  // <<<6>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
  vtkwriter.write(filename,Dune::VTK::appendedraw);

  // Prepare filename for sequential Gnuplot output
  std::ostringstream s;
  if(colCom.size()>1)
  {
    s << 's' << std::setw(4) << std::setfill('0') << colCom.size() << ':';
    s << 'p' << std::setw(4) << std::setfill('0') << colCom.rank() << ':';
  }
  s << filename << ".dat";
  filename = s.str();
  // Gnuplot output
  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
  gnuplotwriter.addVertexData(u,"solution");
  gnuplotwriter.write(filename); 

  // Calculate the forces
  force(gv, boundaryIndexToEntity, gfs, u);

  ipbs_testField(gv,udgf, gfs, u, boundaryIndexToEntity);
}

  std::cout << "Reference total calculation time on rank " << colCom.rank() << " is " << timer.elapsed() << std::endl;
}
