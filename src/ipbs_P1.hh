/** \brief Driver for solving the IPBS problem using P1 piecewise linear Lagrange elements

    Dirichlet boundary conditions are used for domain (outer) boundaries, Neumann b.c. on
    symmetry axis and IPBS Neumann b.c. where we want the IPBS iterative procedure.
*/

/*!
   \param gv the view on the leaf grid
   \param elementIndexToEntity mapper defining the index of inner elements
   \param boundaryIndexToEntity mapper defining the index of boundary elements
*/

#include <map>

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
  // if (helper.rank==0)
  //   std::cout << "constrained dofs=" << cc.size() 
  //             << " of " << gfs.globalSize() << std::endl;

  // interpolate coefficient vector
  // this has to be done only during intitialization as we don't change dirichlet b.c.
  Dune::PDELab::interpolate(g,gfs,u);


  // ================================================================== 
  // Now detect iterative boundary elements on each processor

  if (sysParams.get_verbose() == 5)
    std::cout << "Hello I'm rank " << colCom.rank() << " of " << colCom.size() <<
      " and I have " << gv.size(0) << " codim<0> leaf elements of which " << 
      gv.grid().numBoundarySegments() << " are boundaries." << std::endl;

  int elemCounter = 0;
  int boundCounter = 0;
  int ipbsCount = 0;
  // provide a mapper for storing positions of iterated boundary elements
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> BoundaryElemMapper;
  // typedef Dune::GlobalUniversalMapper<GV::Grid> BoundaryElemMapper;
  BoundaryElemMapper boundaryElemMapper(gv);
  typedef std::map<int, int> StoreMap;
  StoreMap storeMap;

  // provide a vector storing the positions, will be resized automatically
  std::vector<Dune::FieldVector<double,2> > boundaryElemPositions;
  typedef typename GV::template Codim<0>::template Partition
          <Dune::Interior_Partition>::Iterator LeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
            	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
  {
    elemCounter++;
    if(it->hasBoundaryIntersections() == true)
    {
      for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
      {
        if(ii->boundary() == true)
        {
          boundCounter++;
          if (boundaryIndexToEntity[ii->boundarySegmentIndex()] == 2)
          {
            if (sysParams.get_verbose() == 5) // verbose output
              std::cout << "Boundary Center " << ii->geometry().center() << " on rank "
                << colCom.rank() << std::endl;
            ipbsCount ++;
            // insert elements using insert function
            storeMap.insert(std::pair<int, int>(ipbsCount, boundaryElemMapper.map(*it)));
            boundaryElemPositions.push_back(ii->geometry().center());
          }
        }
      }
    }
  }

  // verbose output
  if (sysParams.get_verbose() > 3)
  {
    std::cout << "counted " << elemCounter << " Elements on process " << colCom.rank()
      << " and " << boundCounter << " Boundaries on process " << colCom.rank()
      << " of which " << ipbsCount << " are of iterated type." << std::endl;
    std::cout << "storeMap has size: " << storeMap.size() << std::endl;
    std::cout << "Stored positions are: " << std::endl;

    for(int i = 0; i < boundaryElemPositions.size(); ++i)
    {
      std::cout << boundaryElemPositions[i] << std::endl;
    }

    std::cout << "Store Map looks like:" << std::endl;
    for(int i = 0; i < storeMap.size(); ++i)
    {
      std::cout << storeMap[i] << std::endl;
    }
  }

  // storeMap now contains the local indices of iterative boundary elements
  // and boundaryElemPositions the position vectors on each processor
  // ================================================================== 
  // Now communicate the boundary positions
  
  // master process receives length of vectors on each processor
  MPI_Status status;
  std::vector<int> length_on_processor;
  if (colCom.rank() == 0)
  {
    length_on_processor.push_back(boundaryElemPositions.size());
  }
  for (int i = 1; i < colCom.size(); i++) 
  {
    if (colCom.rank() == i)
    {
      int temp = boundaryElemPositions.size();
      MPI_Send(&temp,1,MPI_INT,0,colCom.rank(),MPI_COMM_WORLD);
    }
    else if (colCom.rank() == 0)
    {
      int temp;
      MPI_Recv(&temp,1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
      length_on_processor.push_back(temp);
    }
  }
  int countBoundElems = 0;
  if (colCom.rank() == 0)
  {
    for (int i=0; i<length_on_processor.size(); i++)
      countBoundElems += length_on_processor[i];
    std::cout << "Master node is now aware that we want to sent " << countBoundElems
      << " boundary elements." << std::endl;
  }
  colCom.broadcast(&countBoundElems,1,0);   // Communicate to all processes how many boundary elements exist

  // now master process receives all positions from other processes
  // first every processor creates an array of its position vectors
  double* my_positions = (double*) malloc(dim*boundaryElemPositions.size()*sizeof(double));
  for (int i = 0; i<boundaryElemPositions.size(); i++){
    for (int j = 0; j<dim; j++)
      my_positions[i*dim+j] = boundaryElemPositions.at(i).vec_access(j);
  }
  
  for (int i = 0; i<storeMap.size(); i++){
    for (int j = 0; j<dim; j+=2)
      std::cout << my_positions[i*dim+j] << my_positions[i*dim+j+1] << std::endl;
  }

  int indexcounter;   // Count how many positions we already got
  double* all_positions = (double*) malloc(dim*countBoundElems*sizeof(double)); // allocate on all processors
  if( colCom.rank() !=0)  // other nodes send their positions to master node
    MPI_Send(my_positions,dim*boundaryElemPositions.size(),MPI_DOUBLE,0,colCom.rank(),MPI_COMM_WORLD);
  if (colCom.rank() == 0) {
    // Write positions of master node
    for (int i = 0; i<boundaryElemPositions.size(); i++){
      for (int j = 0; j<dim; j++)
        all_positions[i*dim+j] = boundaryElemPositions[i][j];
    }
    indexcounter = dim*length_on_processor[0];
    if (sysParams.get_verbose() == 5)
      std::cout << "After performing processor 0 indexcounter is: " << indexcounter << std::endl;
    // get positions from other nodes
    for (int i = 1; i < colCom.size(); i++) {
      std::cout << "Length on processor " << i << " is " << length_on_processor[i] << std::endl;
      MPI_Recv(&all_positions[indexcounter],dim*length_on_processor[i],MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
      indexcounter += dim*length_on_processor[i];
    }
    
    if (sysParams.get_verbose() > 3)
      std::cout << "full iterative boundary elements vectors:" << std::endl;
    for (int i = 0; i<countBoundElems; i++){
      for (int j = 0; j<dim; j+=2)
        std::cout << all_positions[i*dim+j] << all_positions[i*dim+j+1] << std::endl;
    }
  }

  // deploy the iterative boundary positions on all processors
  colCom.broadcast(all_positions,countBoundElems*dim,0); // communicate array

  // ================================================================== 
  // Now every processor knows the iterative boundary elements, i.e. we can start computing

  // allocate a vector for the data whose size is according to the number of iterative boundary elements
  std::vector<double> fluxContainer(countBoundElems);

/*  
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

*/

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
  
}
