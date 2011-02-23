// solve.hh - solve the PB equation

typedef GFS::LocalFunctionSpace LFSU;
// extract some types
typedef LFSU::Traits::FiniteElementType::
Traits::LocalBasisType::Traits::DomainFieldType DF;
typedef LFSU::Traits::FiniteElementType::
Traits::LocalBasisType::Traits::RangeFieldType RF;
typedef LFSU::Traits::FiniteElementType::
Traits::LocalBasisType::Traits::JacobianType JacobianType;
typedef LFSU::Traits::FiniteElementType::
Traits::LocalBasisType::Traits::RangeType RangeType;
typedef LFSU::Traits::SizeType size_type;


template <typename U>
void get_solution(U &u, const GV &gv, const GFS &gfs, const CC &cc, const GridType& grid, const M &m, const B &b, const J&j)
{
  
  //LFSU lfsu(gfs);
  //lfsu.setup(gfs);
  
  typedef GV::Codim<0>::Iterator ElementLeafIterator;
  typedef GV::IntersectionIterator IntersectionIterator;
  const int dimw = 2;
  
  // Provide a mapper for storing gradients
  Mapper mapper (grid);
  // allocate a vector for the data
  std::vector<double> fluxContainer(mapper.size());
  std::vector<double> fluxBackupContainer(mapper.size());
  //std::vector<Dune::FieldVector<RF,dim>> gradientContainer(mapper.size());
  
  //while (sysParams.get_error() > 2E-3)
  while(sysParams.counter < 10)
  {
    if (sysParams.counter == 1)
    {	
      // Get gradients at boundaries
      for (ElementLeafIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
      {      
	/*if (it->hasBoundaryIntersections()==true && it->geometry().center().two_norm() < 4.7) 
	{
	  for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
	  {
	    if (ii->boundary()==true)
	    {	
	      // Bind Local Function Space to this element  
	      lfsu.bind(*it);
	      
	      // local reference coordinate is the center of the intersection :-)
	      Dune::FieldVector<DF,dim> local = ii->geometry().center();
	      
	      // Vector container for storing the coefficients on the actual element (same as U...)
	      GFS::VectorContainer<Real>::Type x_s(gfs,0.0);
	      // get coefficients of this element
	      lfsu.vread(u,x_s);
	      
	      // evaluate gradient of basis functions on reference element
	      std::vector<JacobianType> js(lfsu.size());
	      lfsu.finiteElement().localBasis().evaluateJacobian(ii->geometry().center(),js);
	      
	      // transform gradients from reference element to real element
	      const Dune::FieldMatrix<DF,dimw,dim> 
	      jac = it->geometry().jacobianInverseTransposed(ii->geometry().center());
	      std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
	      for (size_type i=0; i<lfsu.size(); i++)
		jac.mv(js[i][0],gradphi[i]);
	      
	      // compute gradient of u
	      Dune::FieldVector<RF,dim> gradu(0.0);
	      for (size_type i=0; i<lfsu.size(); i++)
		gradu.axpy(x_s[i],gradphi[i]);
		
	      // Store gradient vector of this element for later access
	      //gradientContainer[mapper.map(*it)] = gradu;
	      
	    }	    
	  }	  
	}*/
	fluxContainer[mapper.map(*it)] = sysParams.get_E_init();
	fluxBackupContainer[mapper.map(*it)] = 0.0;
      }
    }
    
    std::cout << std::endl << "IN ITERATION " << sysParams.counter << std::endl << std::endl;
    // Reset error for new iteration
    sysParams.reset_error();
    
    // construct discrete grid function for access to solution
    const DGF udgf(gfs, u);
    
    // <<<4>>> Make grid operator space
    LOP lop(m,b,j, udgf, gv, mapper, fluxContainer, fluxBackupContainer);
    GOS gos(gfs,cc,gfs,cc,lop);
    
    // <<<5a>>> Select a linear solver backend
    LS ls(5000,true);
    
    // <<<5b>>> Instantiate solver for nonlinear problem
    NEWTON newton(gos,u,ls); 
    newton.setReassembleThreshold(0.0);
    newton.setVerbosityLevel(1);
    newton.setReduction(1e-10);
    newton.setMinLinearReduction(1e-4);
    newton.setMaxIterations(40);
    newton.setLineSearchMaxIterations(20);
    
    // <<<5c>>> Instantiate Solver for linear problem
    SLP slp(gos,u,ls,1e-10); 
    
    // <<<6>>> Solve Problem
    newton.apply();
    slp.apply();
    
    std::stringstream out;
    out << "step_" << sysParams.counter;
    std::string vtk_filename = out.str();
    
    DGF udgf_save(gfs,u);
    save(udgf_save, u, gv, vtk_filename);
    
    std::cout << std::endl << "actual error is: " << sysParams.get_error() << std::endl << std::endl;
    sysParams.counter ++;
    fluxBackupContainer = fluxContainer;
    //gradientBackupContainer = gradientContainer;
  }
}