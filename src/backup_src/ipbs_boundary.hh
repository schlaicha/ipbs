/** \brief Precompute the IPBS boundary flux values.

    Here goes some explanation on what is done :-)
    \todo { Doc Me ! }   
*/

#include <gsl/gsl_sf_ellint.h>

#ifndef _P0LAYOUT_H
#define _P0LAYOUT_H
#include "p0layout.hh"
#endif 

template <class GV, class DGF, class ElemPointer,
         typename IndexLookupMap, typename BoundaryElemMapper>

void ipbs_boundary(const GV& gv, const DGF& udgf,
      const std::vector<ElemPointer> ipbsElems,
			double fluxContainer[], double inducedChargeContainer[], const int countBoundElems,
      const std::vector<int>& boundaryIndexToEntity,
      const IndexLookupMap& indexLookupMap, const BoundaryElemMapper& boundaryElemMapper)
{
  // some typedef
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype ctype;
  typedef typename GV::template Codim<0>::template Partition
          <Dune::Interior_Partition>::Iterator LeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename DGF::Traits::RangeType RT;	// store potential during integration 
                                  						// (for calculating sinh-term)
                                              
  // // prepare a mapper - ALL THIS SHOULD GET INTO A CLASS !!!
  // typedef typename Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
  // mapper(gv);
  // const int offset = 0;   //EVIL!!!
  
  // Initialize functor for integrating coulomb flux
  CoulombFlux<ctype,dim> f;
  double dielectric_factor = 0;   // = 2 eps_in / (eps_in + eps_out)
  double epsilon_in= 0;
  double epsilon_out = sysParams.get_epsilon();
  double my_charge = 0;
  
  std::cout << "in ipbs_boundary, countBoundElems = " << countBoundElems << std::endl;
  // Precompute fluxes
  for(int i = 0; i < countBoundElems; i++)
  {
    // Get the unit normal vector of the surface element
    Dune::FieldVector<ctype,dim> unitNormal, r;
    for (IntersectionIterator ii = gv.ibegin(*ipbsElems[i]); ii != gv.iend(*ipbsElems[i]); ++ii)
      if (ii->boundary() == true)
      {
        if (boundaryIndexToEntity[ii->boundarySegmentIndex()] > 1)
        {
          unitNormal = ii->centerUnitOuterNormal();
          r = ii->geometry().center();  // vector of iterative surface boundary center
          dielectric_factor = boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_dielectric_factor();
          epsilon_in =  boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_epsilon();
          my_charge = boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
        }
      }
    unitNormal*=-1.0;	// turn around unit vector as it is outer normal
   
    double summed_flux = 0;
    double summed_volumeElem_flux = 0;
    double summed_surfaceElem_flux = 0;
	   
    for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
            	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
    {
      // Now we go over all elements on this processor
      // and sum up their flux contributions
      Dune::FieldVector<ctype,dim> E_field(0);
	  
      // Get the position vector of the this element's center
      Dune::FieldVector<ctype,dim> r_prime = it->geometry().center();
      Dune::FieldVector<ctype,dim> dist = r - r_prime;
	  
      // integrate sinh over all volume elements 
      double a, b, k, K, E;  // parameters for elliptic
      double itself_flux, volumeElem_flux, surfaceElem_flux;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      
      // ================================================================== 
      // Integral over all volume elements
      // ================================================================== 
      
      // Evaluate the potential at the elements center
	    RT value;
      udgf.evaluate(*it,it->geometry().local(it->geometry().center()),value);
      // Now calculate the flux through surface element caused by this element

      // integration depends on symmetry
      switch( sysParams.get_symmetry() )
      {
        case 1:   // infinite cylinder
          {
            // use minimum image convention
            if (dist[0] > sysParams.get_boxLength()/2.0) dist[0]-=sysParams.get_boxLength();
            else if (dist[0] < -sysParams.get_boxLength()/2.0) dist[0] +=sysParams.get_boxLength();
            a = dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1] + 2.0 * r[1] * r_prime[1];
            b = 4.0 * r[1] * r_prime[1];
            k = sqrt(b/a);
            E = gsl_sf_ellint_Ecomp (k, GSL_PREC_DOUBLE);
            K = gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE);
            E_field[0] = -2.0 * dist[0] / ((a-b)*sqrt(a)) * E;
            E_field[1] = 2.0 * ( 2.0 * r_prime[1] * a - b * r_prime[1] - b * r[1])
                           / ((a-b)*sqrt(a)*b) * E
                         + 2.0 * (-2.0 * r_prime[1] * a + 2.0 * b * r_prime[1]) /  ((a-b)*sqrt(a)*b) * K;
            volumeElem_flux = E_field * unitNormal;
            volumeElem_flux *= -2.0 / (4.0*sysParams.pi) * sysParams.get_lambda2i() * it->geometry().volume() * r_prime[1];
          }
          break;
        
        case 2:   // spherical reduced symmetry
	        {
            a = dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1] + 2.0 * r[1] * r_prime[1];
            b = 4.0 * r[1] * r_prime[1];
            k = sqrt(b/a);
            E = gsl_sf_ellint_Ecomp (k, GSL_PREC_DOUBLE);
            K = gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE);
            E_field[0] = -2.0 * dist[0] / ((a-b)*sqrt(a)) * E;
            E_field[1] = 2.0 * ( 2.0 * r_prime[1] * a - b * r_prime[1] - b * r[1])
                           / ((a-b)*sqrt(a)*b) * E
                         + 2.0 * (-2.0 * r_prime[1] * a + 2.0 * b * r_prime[1]) /  ((a-b)*sqrt(a)*b) * K;
            volumeElem_flux = E_field * unitNormal;
            volumeElem_flux *= -2.0 / (4.0*sysParams.pi) * sysParams.get_lambda2i() * it->geometry().volume() * r_prime[1];
          }
	        break;

        case 3: // An ininite plane in 2d cartesion coordinates 
          {
            // use minimum image convention
            // if (dist[0] > sysParams.get_boxLength()/2.0) dist[0]-=sysParams.get_boxLength();
            // else if (dist[0] < -sysParams.get_boxLength()/2.0) dist[0] +=sysParams.get_boxLength();
           
            E_field = dist;
            E_field /= dist.two_norm() * dist.two_norm();
                      
            volumeElem_flux = E_field * unitNormal;
            volumeElem_flux *= 2. * sysParams.get_lambda2i() / (4.0 * sysParams.pi) * it->geometry().volume();
          }
          break;
     }

        // Now get the counterion distribution and calculate flux
        switch (sysParams.get_salt())
        {
          case 0:
            volumeElem_flux *= std::sinh(value);
            break;
          case 1:
            volumeElem_flux *= std::exp(value); // Counterions have opposite sign!
            break;
        }
        summed_volumeElem_flux += volumeElem_flux;

      // ================================================================== 
      // Integral over surface elements
      // ================================================================== 
 
      if (it->hasBoundaryIntersections() == true)
      {
        for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
          {
             if(ii->boundary() == true)
             {
                if (boundaryIndexToEntity[ii->boundarySegmentIndex()] > 1)
                {
                  if (it != ipbsElems[i])
                  {
                    r_prime = ii->geometry().center();
                    dist = r-r_prime;
                  
                    // Calculate the position of this particular surface element in the inducedChargeContainer
                    int mappedIndex = indexLookupMap.find(boundaryElemMapper.map(*ii->inside()))->second;
                    
                    // we are on a surface element and do integration for coulomb flux
                    // add surface charge contribution from all other surface elements but this one
                    // (using standard coulomb field formula)

                    // Integration depends on symmetry!
                    switch ( sysParams.get_symmetry() )
                    {
                      case 1: // an infinite cylinder
                        {
                          // use minimum image convention
                          if (dist[0] > sysParams.get_boxLength()/2.0) dist[0]-=sysParams.get_boxLength();
                          else if (dist[0] < -sysParams.get_boxLength()/2.0) dist[0] +=sysParams.get_boxLength();
                          a = dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1] + 2.0 * r[1] * r_prime[1];
                          b = 4.0 * r[1] * r_prime[1];
                          k = sqrt(b/a);
                          k = sqrt(b/a);
                          E = gsl_sf_ellint_Ecomp (k, GSL_PREC_DOUBLE);
                          K = gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE);
                          E_field[0] = -2.0 * dist[0] / ((a-b)*sqrt(a)) * E;
                          E_field[1] = 2.0 * ( 2.0 * r_prime[1] * a - b * r_prime[1] - b * r[1])
                                                  / ((a-b)*sqrt(a)*b) * E
                                  + 2.0 * (-2.0 * r_prime[1] * a + 2.0 * b * r_prime[1]) /  ((a-b)*sqrt(a)*b) * K;
                          surfaceElem_flux = E_field * unitNormal;
                          surfaceElem_flux *= 2.0 * sysParams.get_bjerrum() * ii->geometry().volume() * r_prime[1] * 
                                              (boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density() + inducedChargeContainer[mappedIndex]);
                        }
                      break;
                      case 2: // "spherical symmetry"
                      {
                        a = dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1] + 2.0 * r[1] * r_prime[1];
                        b = 4.0 * r[1] * r_prime[1];
                        k = sqrt(b/a);
                        k = sqrt(b/a);
                        E = gsl_sf_ellint_Ecomp (k, GSL_PREC_DOUBLE);
                        K = gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE);
                        E_field[0] = -2.0 * dist[0] / ((a-b)*sqrt(a)) * E;
                        E_field[1] = 2.0 * ( 2.0 * r_prime[1] * a - b * r_prime[1] - b * r[1])
                                                / ((a-b)*sqrt(a)*b) * E
                                + 2.0 * (-2.0 * r_prime[1] * a + 2.0 * b * r_prime[1]) /  ((a-b)*sqrt(a)*b) * K;
                        surfaceElem_flux = E_field * unitNormal;
                        surfaceElem_flux *= 2.0 * sysParams.get_bjerrum() * ii->geometry().volume() * r_prime[1] * 
                                            (boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density() + inducedChargeContainer[mappedIndex]);
                      }
                      break;
                      case 3:
                      {
                        // use minimum image convention
                        // if (dist[0] > sysParams.get_boxLength()/2.0) dist[0]-=sysParams.get_boxLength();
                        // else if (dist[0] < -sysParams.get_boxLength()/2.0) dist[0] +=sysParams.get_boxLength();
                       
                        E_field = dist;
                        E_field /= dist.two_norm() * dist.two_norm();
                                  
                        surfaceElem_flux = E_field * unitNormal;
                        surfaceElem_flux *= -2.0 * sysParams.get_bjerrum() * ii->geometry().volume() *
                                            (boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density() + inducedChargeContainer[mappedIndex]) ;
                      }
                      break;
                    }
                  }
                  else  // this is the flux through the surface element that we impose
                  {
                    if (boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_isPlane() == true)
                      itself_flux = 4 * sysParams.get_bjerrum() * sysParams.pi *boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
                    else
                      itself_flux = -2. * sysParams.get_bjerrum() * sysParams.pi *boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
                    summed_flux += itself_flux;
                  }
                  summed_surfaceElem_flux += surfaceElem_flux;

                  
                }
            }
          }
      }
      // std::cout << std::endl <<" summed_flux is " << summed_flux << " surfaceElem_flux = " << surfaceElem_flux << " volumeElem_flux = " << volumeElem_flux << " itself_flux = " << itself_flux<< std::endl;
    }

    
    // Include dielectrics! Multiply with dielectric factor...
     double diel_flux = 2/(epsilon_in+epsilon_out) * summed_flux + 
                         2*epsilon_in/(epsilon_in+epsilon_out) * (summed_surfaceElem_flux + summed_volumeElem_flux);
     fluxContainer[i] = diel_flux;
    
    // calculate the induced charge in this surface element
    double induced_charge = - (epsilon_in-epsilon_out)/(epsilon_in+epsilon_out) *
                                        ( my_charge + epsilon_out/(2.*sysParams.pi*sysParams.get_bjerrum())*(summed_volumeElem_flux+summed_surfaceElem_flux));
    // double induced_charge = -(epsilon_out-epsilon_in)/(epsilon_out+2.*epsilon_in) * my_charge
    //                              + (epsilon_out-epsilon_in)/(epsilon_out+epsilon_in)*epsilon_out/(2.*sysParams.pi*sysParams.get_bjerrum())*(summed_volumeElem_flux+summed_surfaceElem_flux);
    // fluxContainer[i] = (summed_surfaceElem_flux+summed_volumeElem_flux) + 2.*sysParams.pi*sysParams.get_bjerrum()*(my_charge+induced_charge);
    // Update induced charge on surface elements
    inducedChargeContainer[i] = induced_charge;
  }
}
