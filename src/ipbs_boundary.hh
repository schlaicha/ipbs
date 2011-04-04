/** \brief Precompute the IPBS boundary flux values.

    Here goes some explanation on what is done :-)
    \todo { Doc Me ! }   
*/

template <class GV, class DGF, typename Mapper>
void ipbs_boundary(const GV& gv, const DGF& udgf, const Mapper& mapper,
			std::vector<double>& fluxContainer)
{
  // some typedef
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype ctype;
  typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename DGF::Traits::RangeType RT;	// store potential during integration 
  						// (for calculating sinh-term)
  
  // Initialize functor for integrating coulomb flux
  CoulombFlux<ctype,dim> f;
  
  std::cout << std::endl << "IN ITERATION " << sysParams.counter << std::endl << std::endl;
  // Reset error for new iteration
  sysParams.reset_error();
  
 // Precompute fluxes
 for (ElementLeafIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
 {                
   for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
   {
     if (ii->boundary()==true) 
     {	  
       // Now we have selected the elements for which we want to evaluate the fluxes
       // So now loop over all elements but the surface ones and integrate sinh term
	  
       double fluxIntegrated = 0.0;
       double fluxCoulomb = 0.0;
	  
       // Get the vector of the actual intersection
       Dune::FieldVector<ctype,dim> r = ii->geometry().center();
	  
       // Get the unit normal vector of the surface element
       Dune::FieldVector<ctype,dim> unitNormal = ii->centerUnitOuterNormal();
       unitNormal*=-1.0;	// turn around unit vector as it is outer normal

       for (ElementLeafIterator integrationIterator = gv.template begin<0>();
            	integrationIterator!=gv.template end<0>(); ++integrationIterator)
       {
         if (integrationIterator->hasBoundaryIntersections() == false	// inner element
		 || integrationIterator->geometry().center().two_norm() > 7.8) // outer boundary
         {
           // integrate sinh over all elements but all the surface ones, where
           // the densiy of counterions is zero by definition
	      
           Dune::FieldVector<ctype,dim> r_prime = integrationIterator->geometry().center();
           Dune::FieldVector<ctype,dim> dist = r -r_prime;
           double a,b;	// Arguments for elliptic

           // Evaluate the potential at the elements center
	   RT value;
           udgf.evaluate(*integrationIterator,integrationIterator->geometry().center(),value);

           switch( sysParams.get_symmetry() )
	   {
             case 1: // "2D_cylinder"
	     {
               fluxIntegrated += std::sinh(value) / (dist.two_norm() * dist.two_norm())
	                      * integrationIterator->geometry().volume() * (dist * unitNormal)
			      * sysParams.get_bjerrum()*sysParams.get_lambda2i()*1.0;
	       break;
	     }

	     case 2:
	     {
               a = (dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1]);
	       b = 2 * r[1] * r_prime[1];
	       fluxIntegrated += sysParams.get_lambda2i()/(4.0*sysParams.pi)
                              * eval_elliptic(a,b)
                              * std::sinh(value)*integrationIterator->geometry().volume();
	       //std::cout << "fluxIntegrated: " << fluxIntegrated << std::endl; 
	       break;
	     }

	     case 3:	// "3D"
	     {
	       fluxIntegrated += std::sinh(value) / (dist.two_norm() *dist.two_norm() * dist.two_norm())
                              * integrationIterator->geometry().volume() * (dist * unitNormal)
                              * sysParams.get_bjerrum()*sysParams.get_lambda2i();
               break;
	     }
	     
	     //default: // TODO: put some check here and in the other swith(dim) ! 
	   }
	 }

	 else	// we are on a surface element and do integration for coulomb flux
	      	// add surface charge contribution from all other surface elements but this one
	      	// (using standard coulomb field formula)
	 {
	   // NOTE: For algorithm validation we use the pillowbox conribution
	   if (integrationIterator  == it)
              switch ( sysParams.get_symmetry() )
              {
                case 1:	// "2D_cylinder"
		{	
                  fluxCoulomb += 1.0 * sysParams.get_charge_density()
		              * sysParams.get_bjerrum() * 2.0 * sysParams.pi;
                  break;
		}

                case 2: // "2D_sphere"
		{    
                  fluxCoulomb += 1.0 * sysParams.get_charge_density()  * sysParams.get_bjerrum();
                  break;
                }

                default:
		{
                  fluxCoulomb += 0.0;
                  std::cerr << "IPBS WARNING:\tNo Symmetry specified (in flux evaluation). Will use 0." 
                            << std::endl;
                  break;
		}
	      }

           // TODO Adapt to geometry !!!
           else if (integrationIterator != it && 0 ) // NOTE: This one is to be used later
	   {
             // loop over all boundary intersections of this surface element
             for (IntersectionIterator intersectionIntegrator = gv.ibegin(*integrationIterator);
	     		intersectionIntegrator != gv.iend(*integrationIterator); ++intersectionIntegrator)
             {
               if (intersectionIntegrator->boundary()==true && intersectionIntegrator->neighbor()==false)
                  // integrate coulomb flux along this intersection
                  fluxCoulomb += 1.0 * sysParams.get_bjerrum()*sysParams.get_charge_density() 
		              * integrateentity(ii,f,2,r,unitNormal);
             }
           }
	 }
       }	 

       double flux = fluxCoulomb + fluxIntegrated;
       // Do SOR step and add error
       flux = sysParams.get_alpha() * flux + (1.0 - sysParams.get_alpha()) * fluxContainer[mapper.map(*it)];
       double error = fabs(2.0*(flux-fluxContainer[mapper.map(*it)])/(flux+fluxContainer[mapper.map(*it)]));
       sysParams.add_error(error);
       // Store new flux
       fluxContainer[mapper.map(*it)] = flux;
     }
   }
 }
    
    
    /*DGF udgf_save(gfs,u);
    save(udgf_save, u, gv, vtk_filename);
    */
    std::cout << std::endl << "actual error is: " << sysParams.get_error() << std::endl << std::endl;
    sysParams.counter ++;
}
