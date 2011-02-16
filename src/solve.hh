// solve.hh - solve the PB equation

// Get LeafLocalFunctionSpaceNode - needed for gradient calculation
//typedef typename Dune::PDELab::LeafLocalFunctionSpaceNode< GFS >::LeafLocalFunctionSpaceNode LFSU;
typedef typename GFS::LocalFunctionSpace LFSU;
// extract some types
typedef typename LFSU::Traits::FiniteElementType::
  Traits::LocalBasisType::Traits::DomainFieldType DF;
typedef typename LFSU::Traits::FiniteElementType::
  Traits::LocalBasisType::Traits::RangeFieldType RF;
typedef typename LFSU::Traits::FiniteElementType::
  Traits::LocalBasisType::Traits::JacobianType JacobianType;
typedef typename LFSU::Traits::FiniteElementType::
  Traits::LocalBasisType::Traits::RangeType RangeType;
typedef typename LFSU::Traits::SizeType size_type;


template <typename U>
void get_solution(U &u, const GV &gv, const GFS &gfs, const CC &cc, const M &m, const B &b, const J&j)
{
  
  /*LFSU lfsu(gfs);
  lfsu.debug();
  typedef GV::Codim<0>::Iterator ElementLeafIterator;
  typedef GV::IntersectionIterator IntersectionIterator;
  const int dimw = 2;
  */
  
  
  int iterationCounter = 0;
  sysParams.init = true;
  //while (sysParams.get_error() > 1E-2)
  while(iterationCounter < 10)
  {
    /*
    // Try to get gradients here
    for (ElementLeafIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
    {
      if (it->hasBoundaryIntersections()==true && it->geometry().center().two_norm() < 4.7) 
      {
	for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
	{
	  if (ii->boundary()==true)
	{
	std::cout << "Visited element at surface with";
	 
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
          gradu.axpy(u[i],gradphi[i]);
	
	std::cout <<  "\tGradient u is: " << gradu << std::endl;
	}
	}
      }
    }*/
    
    
    std::cout << std::endl << "IN ITERATION " << iterationCounter << std::endl << std::endl;
    // Reset error for new iteration
    sysParams.reset_error();
    
    // construct discrete grid function for access to solution
    //const U storeCoefficientVector = u;
    //const DGF udgf(gfs, storeCoefficientVector);
    
    const DGF udgf(gfs, u);
    
    // <<<4>>> Make grid operator space
    LOP lop(m,b,j, udgf, gv);
    GOS gos(gfs,cc,gfs,cc,lop);
    
    // <<<5a>>> Select a linear solver backend
    LS ls(5000,true);
    
    // <<<5b>>> Instantiate solver for nonlinear problem
    NEWTON newton(gos,u,ls); 
    newton.setReassembleThreshold(0.0);
    newton.setVerbosityLevel(1);
    newton.setReduction(1e-10);
    newton.setMinLinearReduction(1e-4);
    newton.setMaxIterations(20);
    newton.setLineSearchMaxIterations(10);
    
    // <<<5c>>> Instantiate Solver for linear problem
    SLP slp(gos,u,ls,1e-10); 
    
    // <<<6>>> Solve Problem
    //solver(newton, slp);
    newton.apply();
    slp.apply();
    
    std::stringstream out;
    out << "step_" << iterationCounter;
    std::string vtk_filename = out.str();
    
    DGF udgf_save(gfs,u);
    save(udgf_save, u, gv, vtk_filename);
    ++iterationCounter;
    std::cout << std::endl << "actual error is: " << sysParams.get_error() << std::endl << std::endl;
    sysParams.init = false;
  }
  
  
  
  
}
