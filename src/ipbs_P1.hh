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

template<class GridType, class ColCom>
void ipbs_P1(GridType* grid, const std::vector<int>& elementIndexToEntity,
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
      " and I have " << gv.size(0) << " codim<0> leaf elements." << std::endl;

  int elemCounter = 0;    // debugging variable - count elements on this processor
  int boundCounter = 0;   // debugging variable - count boundary elements on this processor
  int ipbsCount = 0;      // count IPBS boundary elements on this processor
  
  // we store element indices in the order they are found on each processor,
  // so we do not need a mapper but can use an array of appropriate size
  // all other containers use that order and finally a global one with elements
  // of each processor appended is used

  // provide a mapper for getting indices of iterated boundary elements
  typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> BoundaryElemMapper;
  BoundaryElemMapper boundaryElemMapper(gv);
  
  // store ipbs_elements
  std::vector<int> ipbsIndices;
  // provide a vector storing the positions, will be resized automatically
  std::vector<Dune::FieldVector<double,dim> > boundaryElemPositions;
  // provide a vector storing the normals, will be resized automatically
  std::vector<Dune::FieldVector<double,dim> > boundaryElemNormals;
  // provide a vector storing the type of the iterative boundary
  std::vector<int> ipbsType;

  // typedef iterators
  typedef typename GV::template Codim<0>::template Partition
          <Dune::Interior_Partition>::Iterator LeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  
  // loop over elements on this processor
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
          if (boundaryIndexToEntity[ii->boundarySegmentIndex()] > 1) // check if IPBS boundary
          {
            if (sysParams.get_verbose() == 5) // verbose output
              std::cout << "Boundary Center " << ii->geometry().center() << " on rank "
                << colCom.rank() << std::endl;
            ipbsCount ++;
            ipbsIndices.push_back(boundaryElemMapper.map(*it));
            boundaryElemPositions.push_back(ii->geometry().center());
            // remember that the normals point outwards, i.e. when using them we'll have to turn around
            boundaryElemNormals.push_back(ii->centerUnitOuterNormal());
            ipbsType.push_back(boundaryIndexToEntity[ii->boundarySegmentIndex()]);
          }
        }
      }
    }
  }

  // verbose output
  if (sysParams.get_verbose() > 1)
  {
    std::cout << "counted " << elemCounter << " Elements on process " << colCom.rank()
      << " and " << boundCounter << " Boundaries on process " << colCom.rank()
      << " of which " << ipbsCount << " are of iterated type." << std::endl;
    std::cout << "ipbsIndices has size: " << ipbsIndices.size() << std::endl;
    std::cout << "Stored positions are: " << std::endl;
    for(unsigned int i = 0; i < boundaryElemPositions.size(); ++i)
    {
      std::cout << boundaryElemPositions[i] << std::endl;
    }
    std::cout << "Stored normals are: " << std::endl;
    for(unsigned int i = 0; i < boundaryElemNormals.size(); ++i)
    {
      std::cout << boundaryElemNormals[i] << std::endl;
    }
    std::cout << "Store Map looks like:" << std::endl;
    for(unsigned int i = 0; i < ipbsIndices.size(); ++i)
    {
      std::cout << i << "\t" << ipbsIndices[i] << "\t" << std::endl;
    }
  }

  // ipbsIndices now contains the local indices of iterative boundary elements
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
  unsigned int countBoundElems = 0;
  if (colCom.rank() == 0)
  {
    for (unsigned int i=0; i<length_on_processor.size(); i++)
      countBoundElems += length_on_processor[i];
    if (sysParams.get_verbose() > 3)
      std::cout << "Master node is now aware that we want to sent " << countBoundElems
        << " boundary elements." << std::endl;
  }
  colCom.broadcast(&countBoundElems,1,0);   // Communicate to all processes how many boundary elements exist

  // ================================================================== 
  
  // now master process receives all positions from other processes
  // first every processor creates an array of its position vectors
  double* my_positions = (double*) malloc(dim*boundaryElemPositions.size()*sizeof(double));
  double* my_normals = (double*) malloc(dim*boundaryElemNormals.size()*sizeof(double));
  int* my_types = (int*) malloc (ipbsType.size()*sizeof(int));
  for (unsigned int i = 0; i<boundaryElemPositions.size(); i++){
    for (int j = 0; j<dim; j++)
    {
      my_positions[i*dim+j] = boundaryElemPositions.at(i).vec_access(j);
      my_normals[i*dim+j] = boundaryElemNormals.at(i).vec_access(j);
    }
    my_types[i] = ipbsType[i];
  }
  
  if (sysParams.get_verbose() > 4)
  {
    std::cout << "On rank " << colCom.rank() << " positions are:" << std::endl;
    for (unsigned int i = 0; i<boundaryElemPositions.size(); i++){
      for (int j = 0; j<dim; j+=2)
        std::cout << my_positions[i*dim+j] << my_positions[i*dim+j+1] << std::endl;
    }
    std::cout << "On rank " << colCom.rank() << " normals are:" << std::endl;
    for (unsigned int i = 0; i<boundaryElemNormals.size(); i++){
      for (int j = 0; j<dim; j+=2)
        std::cout << my_normals[i*dim+j] << my_normals[i*dim+j+1] << std::endl;
    }
  }

  int indexcounter;   // Count how many positions we already got
  double* all_positions = (double*) malloc(dim*countBoundElems*sizeof(double)); // allocate on all processors
  double* all_normals = (double*) malloc(dim*countBoundElems*sizeof(double)); // allocate on all processors
  int* all_types = (int*) malloc(countBoundElems*sizeof(int));
  if( colCom.rank() !=0) // other nodes send their positions to master node
  {
    MPI_Send(my_positions,dim*boundaryElemPositions.size(),MPI_DOUBLE,0,0,MPI_COMM_WORLD); // pos sent on slot 0
    MPI_Send(my_normals,dim*boundaryElemPositions.size(),MPI_DOUBLE,0,1,MPI_COMM_WORLD); // normals sent on slot 1
    MPI_Send(my_types,boundaryElemPositions.size(),MPI_INT,0,2,MPI_COMM_WORLD); // types sent on slot 2
  }
  if (colCom.rank() == 0)
  {
    // Write positions of master node
    for (unsigned int i = 0; i<boundaryElemPositions.size(); i++){
      for (int j = 0; j<dim; j++)
      {
        all_positions[i*dim+j] = boundaryElemPositions[i][j];
        all_normals[i*dim+j] = boundaryElemNormals[i][j];
      }
      all_types[i] = ipbsType[i];
    }
    // get the right offset
    indexcounter = dim*length_on_processor[0];
    
    if (sysParams.get_verbose() == 5)
      std::cout << "After performing processor 0 indexcounter is: " << indexcounter << std::endl;
    
    // get positions from other nodes
    for (int i = 1; i < colCom.size(); i++) {
      if (sysParams.get_verbose() > 3)
        std::cout << "Length on processor " << i << " is " << length_on_processor[i] << std::endl;
      MPI_Recv(&all_positions[indexcounter],dim*length_on_processor[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
      MPI_Recv(&all_normals[indexcounter],dim*length_on_processor[i],MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);
      MPI_Recv(&all_types[indexcounter/dim],length_on_processor[i],MPI_INT,i,2,MPI_COMM_WORLD,&status);
      indexcounter += dim*length_on_processor[i];
    }
   }

  // deploy the iterative boundary positions on all processors
  colCom.broadcast(all_positions,countBoundElems*dim,0); // communicate array
  colCom.broadcast(all_normals,countBoundElems*dim,0); // communicate array
  colCom.broadcast(all_types,countBoundElems,0); // communicate array
  if (sysParams.get_verbose() > 3)
  {
    std::cout << "on rank " << colCom.rank() << " full iterative boundary elements positions vectors:" << std::endl;
    for (unsigned int i = 0; i<countBoundElems; i++){
      for (int j = 0; j<dim; j+=2)
        if (sysParams.get_verbose() > 3)
          std::cout << all_positions[i*dim+j] << "\t" << all_positions[i*dim+j+1] << std::endl;
    }
    std::cout << "on rank " << colCom.rank() << " full iterative boundary elements normals vectors:" << std::endl;
    for (unsigned int i = 0; i<countBoundElems; i++){
      for (int j = 0; j<dim; j+=2)
        if (sysParams.get_verbose() > 3)
          std::cout << all_normals[i*dim+j] << "\t" << all_normals[i*dim+j+1] << std::endl;
    }

  }

  // ================================================================== 
  // Finally, we need a map so we can reverse lookup the positions of the ipbs_indices
  // note that this map is only needed local on each processor
  typedef std::map<int, int> IndexLookupMap;
  IndexLookupMap indexLookupMap;
  for (unsigned int i=0; i<ipbsIndices.size(); i++)
    indexLookupMap.insert(std::pair<int, int>(ipbsIndices[i],i));


  // ================================================================== 
  // Now every processor knows the iterative boundary elements, i.e. we can start computing

  // allocate array for the data whose size is according to the number of iterative boundary elements
  double* fluxContainer = (double*) malloc(countBoundElems*sizeof(double)); // allocate on all processors
  // for determining the error we need to store the old fluxes
  double* fluxContainerStored = (double*) malloc(countBoundElems*sizeof(double)); // allocate on all processors
  if (colCom.rank() == 0) 
  {
    srand((unsigned)time(0));
    if (sysParams.get_salt() == 0)
      for(unsigned int i=0; i<countBoundElems; i++)
        {
          fluxContainer[i] = 0;   // initialize with zero
          fluxContainerStored[i] = ((float)rand()/RAND_MAX - 0.5) * 0.1;   // random initial b.c.
        }
    else
      for(unsigned int i=0; i<countBoundElems; i++)
        {
          fluxContainer[i] = 0;
          fluxContainerStored[i] = ((float)rand()/RAND_MAX - 1) * 0.5;
        }
  }
  colCom.broadcast(fluxContainerStored,countBoundElems,0); // communicate array
  colCom.broadcast(fluxContainer,countBoundElems,0); // communicate array

  // ================================================================== 
  
  // --- here the iterative loop starts! ---
  
  while (sysParams.get_error() > sysParams.get_tolerance() && sysParams.counter < 15)
  {	  
    // construct a discrete grid function for access to solution
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
    const DGF udgf(gfs, u);

    // Reset error for new iteration
    sysParams.reset_error();

    // call the function precomputing the boundary flux values
    ipbs_boundary(gv,udgf, all_positions, all_normals, fluxContainer, countBoundElems, boundaryIndexToEntity,
        indexLookupMap, boundaryElemMapper);
    // make sure each processor has finished calculations before proceeding
    MPI_Barrier(MPI_COMM_WORLD);
    if (sysParams.get_verbose() > 2)
    {
      std::cout << "Before summing fluxes on rank " << colCom.rank() << ":" << std::endl;
      for (unsigned int i = 0; i<countBoundElems; i++)
        std::cout << i << "\t" << fluxContainer[i] << std::endl;
    }
  
    // now each processor has computed the values, i.e. we must ALL_REDUCE
    // colCom.sum(fluxContainer,countBoundElems);
    MPI_Allreduce(MPI_IN_PLACE, fluxContainer, countBoundElems, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // do SOR step determine the error on each processor
    for (unsigned int i = 0; i < countBoundElems; i++)
    {
      double fluxCoulomb = boundary[all_types[i]-2]->get_charge_density()  * sysParams.get_bjerrum() * (4.0 * sysParams.pi);

      fluxContainer[i] = sysParams.get_alpha() * (fluxContainer[i] + fluxCoulomb)
                          + ( 1 - sysParams.get_alpha()) * fluxContainerStored[i];
      double error = fabs(2.0*(fluxContainer[i]-fluxContainerStored[i])
                      /(fluxContainer[i]+fluxContainerStored[i]));
      sysParams.add_error(error);
    }
    // set the error to the global maximum of errors so we have a unique stop criterion
    // double error = sysParams.get_error();
    // std::cout << "Error on rank " << colCom.rank() << ": " << error << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    // colCom.max(error);
    // std::cout << "Max error: " << error << std::endl;
    // sysParams.reset_error(error);
    // print error on rank 0
    if (colCom.rank()==0)
      std::cout << std::endl << "In iteration " << sysParams.counter <<" the relative error is: " 
        << sysParams.get_error() << std::endl << std::endl;

    int offset;
    // determine the offset for each processor to access its fluxes
    if(colCom.rank()==0)
    {
      offset = 0;
      int sendOffset = 0;
      for (unsigned int i=1;i<length_on_processor.size();i++)
      {
        sendOffset += length_on_processor[i];
        MPI_Send(&sendOffset,1,MPI_INT,i,0,MPI_COMM_WORLD);
      }
    }
    if (colCom.rank()!=0)
      MPI_Recv(&offset,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
    
    if (sysParams.get_verbose() > 2)
    {
      std::cout << "fluxes on rank " << colCom.rank() << ":" << std::endl;
      for (unsigned int i = 0; i<countBoundElems; i++)
        std::cout << i << "\t" << fluxContainer[i] << std::endl;
    }
 
    // instanciate boundary fluxes
    typedef BoundaryFlux<GV,double,std::vector<int>, IndexLookupMap> J;
    J j(gv, boundaryIndexToEntity, fluxContainer,indexLookupMap,offset);

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
    newton.setVerbosityLevel(sysParams.get_verbose());
    newton.setReduction(1e-10);
    newton.setMinLinearReduction(1e-4);
    newton.setMaxIterations(25);
    newton.setLineSearchMaxIterations(10);
    newton.apply();

    // save snapshots of each iteration step
    std::stringstream out;
    out << "ipbs_step_" << sysParams.counter;
    std::string filename = out.str();
    DGF udgf_snapshot(gfs,u);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf_snapshot,"solution"));
    vtkwriter.write(filename,Dune::VTK::appendedraw);
    // Prepare filename for sequential Gnuplot output
    if(colCom.size()>1)
    {
      std::stringstream s;
      s << 's' << std::setw(4) << std::setfill('0') << colCom.size() << ':';
      s << 'p' << std::setw(4) << std::setfill('0') << colCom.rank() << ':';
      s << filename;
      filename = s.str();
    }
    // Gnuplot output
    Dune::GnuplotWriter<GV> gnuplotwriter(gv);
    gnuplotwriter.addVertexData(u,"solution");
    gnuplotwriter.write(filename + ".dat"); 

    sysParams.counter ++;
    // copy flux array to backup one for next iteration step
    for (unsigned int i = 0; i < countBoundElems; i++)
    {
      fluxContainerStored[i] = fluxContainer [i];
    }
  }

  // --- here the iterative loop ends! ---


  // <<<6>>> graphical output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
  vtkwriter.write("ipbs_solution",Dune::VTK::appendedraw);

  // Prepare filename for sequential Gnuplot output
  std::string filename = "ipbs_solution";
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
  
  std::cout << "iPBS calculation time on rank " << colCom.rank() << " is: " << timer.elapsed() << std::endl;
  
}
