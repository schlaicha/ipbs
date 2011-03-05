// solve.hh - solve the PB equation

#include "functors.hh"
#include "integrateentity.hh"

template <typename U>
void get_solution(U &u, const GV &gv, const GFS &gfs, const CC &cc, const GridType& grid, const M &m, const B &b, const J&j)
{
  //typedef GV::Traits Traits;
  typedef double ctype;		// TODO: Find correct traits and implement correct!
  
  typedef GV::Codim<0>::Iterator ElementLeafIterator;
  typedef GV::IntersectionIterator IntersectionIterator;
  
  // Provide a mapper for storing precomputed values
  Mapper mapper(grid);
  // allocate a vector for the data
  std::vector<double> fluxContainer(mapper.size());
  
  typename Traits::RangeType value;	// store potential during integration (for calculating sinh-term)
  
  // Initialize functor for integrating coulomb flux
  CoulombFlux<ctype,dim> f;
  
  while (sysParams.get_error() > 1E-3)
  //while (sysParams.counter < 10)    
  {
    std::cout << std::endl << "IN ITERATION " << sysParams.counter << std::endl << std::endl;
    // Reset error for new iteration
    sysParams.reset_error();
  
    // Precompute fluxes
    for (ElementLeafIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
    {                
      for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
      {
	if (ii->boundary()==true && it->geometry().center().two_norm() < 4.7) 
	{	  
	  // Now we have selected the elements for which we want to evaluate the fluxes
	  // So now loop over all elements but the surface ones and integrate sinh term
	  
	  double fluxIntegrated = 0.0;
	  double fluxCoulomb = 0.0;
	  
	  // construct discrete grid function for access to solution
	  const DGF udgf(gfs, u);
	  
	  // Get the vector of the actual intersection
	  Dune::FieldVector<ctype,dim> r = ii->geometry().center();
	  
	  // Get the unit normal vector of the surface element
	  Dune::FieldVector<ctype,dim> unitNormal = ii->centerUnitOuterNormal();
	  unitNormal*=-1.0;	// turn around unit vector as it is outer normal
	  
	  for (ElementLeafIterator integrationIterator = gv.begin<0>(); integrationIterator!=gv.end<0>(); ++integrationIterator)
	  {
	    if (integrationIterator->hasBoundaryIntersections() == false || integrationIterator->geometry().center().two_norm() > 4.7)
	    {
	      // integrating sinh over all elements but all the surface ones, where
	      // the densiy of counterions is zero by definition)
	      
	      Dune::FieldVector<ctype,dim> r_prime = integrationIterator->geometry().center();
	      Dune::FieldVector<ctype,dim> dist = r -r_prime;
	      
	      // Evaluate the potential at the elements center
	      udgf.evaluate(*integrationIterator,integrationIterator->geometry().center(),value);
	      
          switch( sysParams.get_symmetry() ){
            case 1:     // "2D_cylinder"
		        fluxIntegrated += std::sinh(value) / (dist.two_norm() * dist.two_norm())
				    * integrationIterator->geometry().volume() * (dist * unitNormal)
				    * sysParams.get_bjerrum()*sysParams.get_lambda2i()*2.0 ;
            break;
		    
            case 2:
                fluxIntegrated += sysParams.get_lambda2i()/(4.0*sysParams.pi)
                        * eval_elliptic(dist[0],dist[1],(r*r_prime))
                        * std::sinh(value);
		    break;

		    case 3:     // "3D"
		        fluxIntegrated += std::sinh(value) / (dist.two_norm() *dist.two_norm() * dist.two_norm())
				  * integrationIterator->geometry().volume() * (dist * unitNormal)
				  * sysParams.get_bjerrum()*sysParams.get_lambda2i();
		    break;
		    default: // TODO: put some check here and in the other swith(dim) ! 
		    break;
	      }
	    }
	    
	    else	// we are on a surface element and do integration for coulomb flux
	      // add surface charge contribution from all other surface elements but this one
	      // (using standard coulomb field formula)
	      {
		// NOTE: For algorithm validation we use the pillowbox conribution
		if (integrationIterator  == it)
		  fluxCoulomb += 1.0 * sysParams.get_charge_density() *  sysParams.get_bjerrum() * 2.0 *sysParams.pi;
		
		if (integrationIterator != it && 0 ) // NOTE: This one is to be used later
		{
		  // loop over all boundary intersections of this surface element
		  for (IntersectionIterator intersectionIntegrator = gv.ibegin(*integrationIterator); intersectionIntegrator != gv.iend(*integrationIterator); ++intersectionIntegrator)
		  {
		    if (intersectionIntegrator->boundary()==true && intersectionIntegrator->neighbor()==false)
		      // integrate coulomb flux along this intersection
		    fluxCoulomb += 1.0 * sysParams.get_bjerrum()*sysParams.get_charge_density() * integrateentity(ii,f,2,r,unitNormal);
		  }
		}
	      }
	  }
	  

	  double flux = fluxCoulomb + fluxIntegrated;
	  //if (sysParams.counter != 0)
	  //{
	    // Do SOR step and add error
	    flux = sysParams.get_alpha() * flux + (1.0 - sysParams.get_alpha()) * fluxContainer[mapper.map(*it)];
	    double error = fabs(2.0*(flux-fluxContainer[mapper.map(*it)])/(flux+fluxContainer[mapper.map(*it)]));
	    sysParams.add_error(error);
	  //}
	  // Store new flux
	  fluxContainer[mapper.map(*it)] = flux;
	}
      }
      
    }
    
    // <<<4>>> Make grid operator space
    LOP lop(m,b,j,fluxContainer);
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
  }
}
