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
  BoundaryElemMapper boundaryElemMapper(gv);
  typedef std::map<int, int> StoreMap;
  typedef std::vector<int> StoreMapInverse;
  StoreMap storeMap;
  StoreMapInverse storeMapInverse;

  // provide a vector storing the positions, will be resized automatically
  std::vector<Dune::FieldVector<double,dim> > boundaryElemPositions;
  // provide a vector storing the normals, will be resized automatically
  std::vector<Dune::FieldVector<double,dim> > boundaryElemNormals;

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
            storeMap.insert(std::pair<int, int>(boundaryElemMapper.map(*it),ipbsCount));
            storeMapInverse.push_back(boundaryElemMapper.map(*it));
            boundaryElemPositions.push_back(ii->geometry().center());
            // remember that the normals point outwards, i.e. when using them we'll have to turn around
            boundaryElemNormals.push_back(ii->centerUnitOuterNormal());
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
      std::cout << storeMapInverse[i] << "\t" << storeMap.find(storeMapInverse[i])->second << std::endl;
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
    if (sysParams.get_verbose() > 3)
      std::cout << "Master node is now aware that we want to sent " << countBoundElems
        << " boundary elements." << std::endl;
  }
  colCom.broadcast(&countBoundElems,1,0);   // Communicate to all processes how many boundary elements exist

  // ================================================================== 
  // get a global store map
  StoreMap globalStoreMap;
  if(colCom.rank() == 0)
  {
    globalStoreMap = storeMap; // On rank 0 we don't have an offset
    // Now receive the indices from the other processors
    for (int i=1; i<colCom.rank(); i++)
    {
      int* index = (int*) malloc(length_on_processor[i]*sizeof(int));
      MPI_Recv(index,length_on_processor[i],MPI_INT,i,2,MPI_COMM_WORLD,&status); // send indices on slot 2
      for (int j=0; j<length_on_processor[i]; j++)
        globalStoreMap.insert(std::pair<int, int>(j+length_on_processor[i-1], index[j]));
      free(index);
    }
  }
  if (colCom.rank()!=0)
  {
    // Create an array containing the index and sent to master
    int* index = (int*) malloc(storeMap.size()*sizeof(int));
    for(int i=0; i<storeMap.size(); i++)
      index[i] = storeMap.find(i)->second;
    //  std::cout << "StoreMapInverse: " << storeMapInverse[i] << std::endl;
    MPI_Send(index, storeMap.size(), MPI_INT, 0, 2, MPI_COMM_WORLD);
    free(index);
  }

  // Now send the globalStoreMap
  if (colCom.rank()==0)
  {
    int* index = (int*) malloc(globalStoreMap.size()*sizeof(int));
    for(int i=0; i<globalStoreMap.size(); i++)
      index[i] = globalStoreMap.find(i)->second;
    for (int i=1; i<colCom.size(); i++)
      MPI_Send(index, globalStoreMap.size(), MPI_INT, i, 3, MPI_COMM_WORLD);
    free(index);
  }
  if (colCom.rank()!= 0)
  {
    int* index = (int*) malloc(countBoundElems*sizeof(int));
    MPI_Recv(index,countBoundElems,MPI_INT,0,3,MPI_COMM_WORLD,&status); // recieve global indices on slot 3
    for(int i=0; i<countBoundElems; i++)
      globalStoreMap.insert(std::pair<int, int>(i, index[i]));
    free(index);
  }

  // ================================================================== 
  // now master process receives all positions from other processes
  // first every processor creates an array of its position vectors
  double* my_positions = (double*) malloc(dim*boundaryElemPositions.size()*sizeof(double));
  double* my_normals = (double*) malloc(dim*boundaryElemPositions.size()*sizeof(double));
  for (int i = 0; i<boundaryElemPositions.size(); i++){
    for (int j = 0; j<dim; j++)
    {
      my_positions[i*dim+j] = boundaryElemPositions.at(i).vec_access(j);
      my_normals[i*dim+j] = boundaryElemPositions.at(i).vec_access(j);
    }
  }
  
  if (sysParams.get_verbose() > 4)
  {
    std::cout << "On rank " << colCom.rank() << " positions are:" << std::endl;
    for (int i = 0; i<storeMap.size(); i++){
      for (int j = 0; j<dim; j+=2)
        std::cout << my_positions[i*dim+j] << my_positions[i*dim+j+1] << std::endl;
    }
  }

  int indexcounter;   // Count how many positions we already got
  double* all_positions = (double*) malloc(dim*countBoundElems*sizeof(double)); // allocate on all processors
  double* all_normals = (double*) malloc(dim*countBoundElems*sizeof(double)); // allocate on all processors
  if( colCom.rank() !=0)  // other nodes send their positions to master node
    MPI_Send(my_positions,dim*boundaryElemPositions.size(),MPI_DOUBLE,0,0,MPI_COMM_WORLD); // pos sent on slot 0
    MPI_Send(my_normals,dim*boundaryElemPositions.size(),MPI_DOUBLE,0,1,MPI_COMM_WORLD); // normals sent on slot 1
  if (colCom.rank() == 0) {
    // Write positions of master node
    for (int i = 0; i<boundaryElemPositions.size(); i++){
      for (int j = 0; j<dim; j++)
      {
        all_positions[i*dim+j] = boundaryElemPositions[i][j];
        all_normals[i*dim+j] = boundaryElemNormals[i][j];
      }
    }
    
    indexcounter = dim*length_on_processor[0];
    
    if (sysParams.get_verbose() == 5)
      std::cout << "After performing processor 0 indexcounter is: " << indexcounter << std::endl;
    
    // get positions from other nodes
    for (int i = 1; i < colCom.size(); i++) {
      if (sysParams.get_verbose() > 3)
        std::cout << "Length on processor " << i << " is " << length_on_processor[i] << std::endl;
      MPI_Recv(&all_positions[indexcounter],dim*length_on_processor[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
      MPI_Recv(&all_normals[indexcounter],dim*length_on_processor[i],MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);
      indexcounter += dim*length_on_processor[i];
    }
    
    if (sysParams.get_verbose() > 3)
      std::cout << "full iterative boundary elements vectors:" << std::endl;
    for (int i = 0; i<countBoundElems; i++){
      for (int j = 0; j<dim; j+=2)
        if (sysParams.get_verbose() > 3)
          std::cout << all_positions[i*dim+j] << "\t" << all_positions[i*dim+j+1] << std::endl;
    }
  }

  // deploy the iterative boundary positions on all processors
  colCom.broadcast(all_positions,countBoundElems*dim,0); // communicate array
  colCom.broadcast(all_normals,countBoundElems*dim,0); // communicate array

  // ================================================================== 
  // Now every processor knows the iterative boundary elements, i.e. we can start computing

  // allocate a vector for the data whose size is according to the number of iterative boundary elements
  // std::vector<double> fluxContainer(countBoundElems);
  double* fluxContainer = (double*) malloc(countBoundElems*sizeof(double)); // allocate on all processors

  // ================================================================== 
  
  // --- here the iterative loop starts! ---
  
  while (sysParams.get_error() > 1E-3 && sysParams.counter < 11)
  // while (sysParams.counter < 3)
  {	  
    // construct a discrete grid function for access to solution
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
    const DGF udgf(gfs, u);

    // Reset error for new iteration
    sysParams.reset_error();

    // call the function precomputing the boundary flux values
    ipbs_boundary(gv,udgf, all_positions, all_normals, fluxContainer, countBoundElems, boundaryIndexToEntity);
    if (colCom.rank()==0)
      std::cout << std::endl << "In iteration " << sysParams.counter <<" the relative error is: " 
        << sysParams.get_error() << std::endl << std::endl;

    std::cout << "Before summing fluxes on rank " << colCom.rank() << ":" << std::endl;
    for (int i = 0; i<countBoundElems; i++)
      std::cout << i << "\t" << fluxContainer[i] << std::endl;


    // now each processor has computed the values, i.e. we must ALL_REDUCE
    colCom.sum(fluxContainer,countBoundElems);

    int offset;
    // determine the offset for each processor to access its fluxes
    if(colCom.rank()==0)
    {
      offset = 0;
      int sendOffset = 0;
      for (int i=1;i<length_on_processor.size();i++)
      {
        sendOffset += length_on_processor[i];
        MPI_Send(&sendOffset,1,MPI_INT,i,0,MPI_COMM_WORLD);
      }
    }
    if (colCom.rank()!=0)
      MPI_Recv(&offset,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
    
    // TODO Now the boundary flux for a given element can be obtained by using
    // storeMapInverse[mapper.map(*it)+offset] where offset is given for each procesor

    std::cout << "fluxes on rank " << colCom.rank() << ":" << std::endl;
    for (int i = 0; i<globalStoreMap.size(); i++)
      std::cout << i << "\t" << fluxContainer[i] << std::endl;

    // instanciate boundary fluxes
    typedef BoundaryFlux<GV,double,std::vector<int>, StoreMap> J;
    J j(gv, boundaryIndexToEntity, fluxContainer,globalStoreMap,offset);

    // <<<4>>> Make Grid Operator Space
    typedef PBLocalOperator<M,B,J> LOP;
    LOP lop(m,b,j);
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

    sysParams.counter ++;
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
  
}
