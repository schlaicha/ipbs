/** \brief Driver for solving the IPBS problem using P1 piecewise linear Lagrange elements

    Dirichlet boundary conditions are used for domain (outer) boundaries, Neumann b.c. on
    symmetry axis and IPBS Neumann b.c. where we want the IPBS iterative procedure.
*/

/*!
   \param gv the view on the leaf grid
   \param elementIndexToEntity mapper defining the index of inner elements
   \param boundaryIndexToEntity mapper defining the index of boundary elements
*/

template<class GV, class ColCom>
void ipbs_P1(const GV& gv, const std::vector<int>& elementIndexToEntity,
             const std::vector<int>& boundaryIndexToEntity,
             const ColCom& colCom)
{
  // We want to know the total calulation time
  // Dune::Timer timer;

  // some typedef
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype ctype;

  // <<<1>>> Setup the problem from mesh file

  // NOTE only flux b.c. are iterated, so rest is instantiated only once
  // inner region
  typedef Regions<GV,double,std::vector<int>> M;
  M m(gv, elementIndexToEntity);
  // boundary conditioon type
  typedef BCType<GV,std::vector<int>> B;
  B b(gv, boundaryIndexToEntity);
  // Class defining Dirichlet B.C.
  typedef BCExtension<GV,double,std::vector<int>> G;
  G g(gv, boundaryIndexToEntity);

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
  // if (helper.rank==0)
  //   std::cout << "constrained dofs=" << cc.size() 
  //             << " of " << gfs.globalSize() << std::endl;

  // interpolate coefficient vector
  // this has to be done only during intitialization as we don't change dirichlet b.c.
  Dune::PDELab::interpolate(g,gfs,u);

  std::cout << "Hello I'm rank " << colCom.rank() << " of " << colCom.size() <<
    " and I have " << gv.size(0) << " codim<0> leaf elements of which " << 
    gv.grid().numBoundarySegments() << " are boundaries." << std::endl;

  int counter = 0;
  typedef typename GV::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
  //typename GV::Grid::template Codim<0>::template Partition<Dune::Interior_Partition>::LeafIterator it;
  for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
            	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
  {
    counter++;
  }
  std::cout << "counted " << counter << " Elements on process " << colCom.rank() << std::endl;


/*  
  // Provide a mapper for storing precomputed values
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV, P0Layout> Mapper;
  Mapper mapper(gv);
  // allocate a vector for the data
  std::vector<double> fluxContainer(mapper.size());

  // --- here the iterative loop starts! ---
  
  while (sysParams.get_error() > 1E-3 && sysParams.counter < 11)
  // while (sysParams.counter < 3)
  {	  
    // construct a discrete grid function for access to solution
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
    const DGF udgf(gfs, u);

    // call the function precomputing the boundary flux values
    ipbs_boundary(gv,udgf, mapper, fluxContainer);

    // boundary fluxes - this one is for the reference solution!
    typedef BoundaryFlux<GV,double,std::vector<int> > J;
    J j(gv, boundaryIndexToEntity);

    // <<<4>>> Make Grid Operator Space
    typedef PBLocalOperator<M,B,J> LOP;
    LOP lop(m,b,j, fluxContainer);
    typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
    typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE,true> GOS;
    GOS gos(gfs,cc,gfs,cc,lop);

    // <<<5a>>> Select a linear solver backend
    typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GFS> LS;
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

    // save snapshots of each iteration step
    std::stringstream out;
    out << "step_" << sysParams.counter;
    std::string filename = out.str();
    DGF udgf_snapshot(gfs,u);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf_snapshot,"solution"));
    vtkwriter.write(filename,Dune::VTK::appendedraw);
    // Gnuplot output
    Dune::GnuplotWriter<GV> gnuplotwriter(gv);
    gnuplotwriter.addVertexData(u,"solution");
    gnuplotwriter.write(filename + ".dat"); 
  }

  // --- here the iterative loop ends! ---

  // <<<6>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
  vtkwriter.write("ipbs_solution",Dune::VTK::appendedraw);

  // Gnuplot output
  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
  gnuplotwriter.addVertexData(u,"solution");
  gnuplotwriter.write("ipbs_solution.dat"); 
  
  // std::cout << "Reference total calculation time=" << timer.elapsed() << std::endl;
  

*/
}
