/** \brief Precompute the IPBS boundary flux values.

    Here goes some explanation on what is done :-)
    \todo { Doc Me ! }   
*/

template <class GV, class DGF>
void ipbs_boundary(const GV& gv, const DGF& udgf,
      const double positions[], const double normals[],
			double fluxContainer[], const int &countBoundElems,
      const std::vector<int>& boundaryIndexToEntity)
{
  // some typedef
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype ctype;
  typedef typename GV::template Codim<0>::template Partition
          <Dune::Interior_Partition>::Iterator LeafIterator;
  typedef typename DGF::Traits::RangeType RT;	// store potential during integration 
                                  						// (for calculating sinh-term)

  
  // // show me the distributed vectors
  // for(int i = 0; i < fluxContainer.size(); i++)
  // {
  //   std::cout << positions[i*dim] << "\t" << positions[i*dim+1] << std::endl;
  // }
  // // show me the distributed vectors
  // for(int i = 0; i < fluxContainer.size(); i++)
  // {
  //   std::cout << normals[i*dim] << "\t" << normals[i*dim+1] << std::endl;
  // }

  // Initialize functor for integrating coulomb flux
  CoulombFlux<ctype,dim> f;
  
  // Precompute fluxes
  for(int i = 0; i < countBoundElems; i++)
  {
    int testcounter = 0;

    // Get the unit normal vector of the surface element
    Dune::FieldVector<ctype,dim> unitNormal;
    Dune::FieldVector<ctype,dim> r_prime;  // vector of iterative surface boundary center
    for(int j=0; j<dim; j++)
    {
        unitNormal[j] = normals[i*dim+j];
        r_prime[j] = positions[i*dim+j];
    }
    unitNormal*=-1.0;	// turn around unit vector as it is outer normal
   
    double fluxIntegrated = 0.0;
    double fluxCoulomb = 0.0;
	   
    for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
            	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
    {
      // Now we go over all elements on this processor
      // and sum up their flux contributions
	  
      // Get the position vector of the this element's center
      Dune::FieldVector<ctype,dim> r = it->geometry().center();
      Dune::FieldVector<ctype,dim> dist = r - r_prime;
	  
      // integrate sinh over all elements but all the surface ones, where
      // the densiy of counterions is zero by definition
      bool isIPBS_Elem = false;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      // check that we are not on an IPBS surface element
      if(it->hasBoundaryIntersections() == true)
      for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
        if(ii->boundary() == true && boundaryIndexToEntity[ii->boundarySegmentIndex()] == 2)
          isIPBS_Elem = true;

      // ================================================================== 
      // Integral over all elements but the surface ones
      // ================================================================== 
     
      if (isIPBS_Elem == false)
      {
        // Evaluate the potential at the elements center
	      RT value;
        udgf.evaluate(*it,it->geometry().center(),value);

        // integration depends on symmetry
        switch( sysParams.get_symmetry() )
        {
          case 1: // "2D_cylinder"
            {
              fluxIntegrated += std::sinh(value) / (dist.two_norm() * dist.two_norm())
                * it->geometry().volume() * (dist * unitNormal)
                / sysParams.get_bjerrum()*sysParams.get_lambda2i()*1.0/4.0/sysParams.pi;
            }
            break;
        }

      }

      // ================================================================== 
      // Integral over surface elements
      // ================================================================== 
      else if (isIPBS_Elem == true)
      {
        // we are on a surface element and do integration for coulomb flux
	      // add surface charge contribution from all other surface elements but this one
	      // (using standard coulomb field formula)
        // NOTE: For algorithm validation we use the pillowbox conribution

        switch ( sysParams.get_symmetry() )
          {
            case 1:	// "2D_cylinder"
              {
                fluxCoulomb += 1.0 * sysParams.get_charge_density()
                  * sysParams.get_bjerrum() * 2.0 * sysParams.pi;
              }
              break;
            case 2: // "2D_sphere"
              {
                fluxCoulomb += 1.0 * sysParams.get_charge_density()  * sysParams.get_bjerrum();
              }
            break;
          }
      }
    }
      
    // ================================================================== 
    // Integral over surface elements
    // ================================================================== 
    double flux = fluxCoulomb + fluxIntegrated;
    // Do SOR step and add error
    flux = sysParams.get_alpha() * flux + (1.0 - sysParams.get_alpha()) * fluxContainer[i];
    double error = fabs(2.0*(flux-fluxContainer[i])/(flux+fluxContainer[i]));
    sysParams.add_error(error);
    // Store new flux
    fluxContainer[i] = flux;
  }


  
  /*
  double a,b;	// Arguments for elliptic

           
           switch( sysParams.get_symmetry() )
	   {
             case 1: // "2D_cylinder"
	     {
               fluxIntegrated += std::sinh(value) / (dist.two_norm() * dist.two_norm())
	                      * integrationIterator->geometry().volume() * (dist * unitNormal)
                              / sysParams.get_bjerrum()*sysParams.get_lambda2i()*1.0/4.0/sysParams.pi;
             }
	     break;
	     

	     case 2:
	     {
               a = (dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1]);
	       b = 2 * r[1] * r_prime[1];
	       fluxIntegrated += sysParams.get_lambda2i()/(4.0*sysParams.pi)
                              * eval_elliptic(a,b)
                              * std::sinh(value)*integrationIterator->geometry().volume();
	       //std::cout << "fluxIntegrated: " << fluxIntegrated << std::endl; 
	     }
	     break;

	     case 3:	// "3D"
	     {
	       fluxIntegrated += std::sinh(value) / (dist.two_norm() *dist.two_norm() * dist.two_norm())
                              * integrationIterator->geometry().volume() * (dist * unitNormal)
                              * sysParams.get_bjerrum()*sysParams.get_lambda2i();
	     }
             break;
	     
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
		}
                break;

                case 2: // "2D_sphere"
		{    
                  fluxCoulomb += 1.0 * sysParams.get_charge_density()  * sysParams.get_bjerrum();
                }
                break;

                default:
		{
                  fluxCoulomb += 0.0;
                  std::cerr << "IPBS WARNING:\tNo Symmetry specified (in flux evaluation). Will use 0." 
                            << std::endl;
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

*/
   
}
