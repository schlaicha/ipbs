/** \brief Calculate the field using the Coulomb integral 
    and compare it to the field obtain via gradient of potential ;-) 

 Here goes some explanation on what is done :-)
    \todo { Doc Me ! }   
*/

#include <gsl/gsl_sf_ellint.h>

template <class GV, class DGF, class GFS, class U>
void ipbs_testField(const GV& gv, const DGF& udgf, const GFS& gfs, const U& u,
    std::vector<int> boundaryIndexToEntity)
{
  // some typedef
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype ctype;
  typedef typename GV::template Codim<0>::template Partition
          <Dune::Interior_Partition>::Iterator LeafIterator;
  typedef typename DGF::Traits::RangeType RT;	// store potential during integration 
                                  						// (for calculating sinh-term)
   typedef typename GV::IntersectionIterator IntersectionIterator;
  
  // Initialize functor for integrating coulomb flux
  CoulombFlux<ctype,dim> f;
  
  std::ofstream coulomb_file, gradient_file;
  coulomb_file.open ("efield_coulomb.dat", std::ios::out );
  gradient_file.open ("efield_gradient.dat", std::ios::out );
  for (LeafIterator elemIt = gv.template begin<0,Dune::Interior_Partition>();
            	elemIt!=gv.template end<0,Dune::Interior_Partition>(); ++elemIt)
  {
    if (elemIt->hasBoundaryIntersections() == true)
     {
         for (IntersectionIterator elemii = gv.ibegin(*elemIt); elemii != gv.iend(*elemIt); ++elemii)
         {
          if(elemii->boundary() == true)
          {
             if (boundaryIndexToEntity[elemii->boundarySegmentIndex()] > 1)
             {
    // Loop over all elements and calculate the field
  
    // Get the unit normal vector of the surface element
    Dune::FieldVector<ctype,dim> r = elemii->geometry().center();
    Dune::FieldVector<ctype,dim> summed_field(0);
	   
    for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
            	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
    {
      Dune::FieldVector<ctype,dim> surfaceElem_field(0);
      Dune::FieldVector<ctype,dim> volumeElem_field(0);
      // Now we go over all elements on this processor
      // and sum up their flux contributions
	  
      // Get the position vector of the this element's center
      Dune::FieldVector<ctype,dim> r_prime = it->geometry().center();
      Dune::FieldVector<ctype,dim> dist = r - r_prime;
	  
      double a, b, k, K, E;  // parameters for elliptic

      // ================================================================== 
      // Integral over all elements 
      // ================================================================== 
      
      //if (it != elemIt)
      //{
        // Evaluate the potential at the elements center
	      RT value;
        udgf.evaluate(*it,it->geometry().center(),value);
        // Calculate the flux through surface element caused by this element

        // integration depends on symmetry
        // spherical reduced 3d symmetry
	   
        // Calculate the integrated sinh term for finite geometry
        //  a = (dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1]);
	      //  b = 2.0 * r[1] * r_prime[1];
        //  k = sqrt(2.0*b/(a+b));
         
        a = dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1] + 2.0 * r[1] * r_prime[1];
        b = 4.0 * r[1] * r_prime[1];
        k = sqrt(b/a);

        E = gsl_sf_ellint_Ecomp (k, GSL_PREC_DOUBLE);
        K = gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE);

        volumeElem_field[0] = -1.0 * dist[0] / ((a-b)*sqrt(a)) * E;
        volumeElem_field[1] = ( 2.0 * r_prime[1] * a - b * r_prime[1] - b * r[1]) / ((a-b)*sqrt(a)*b) * E
                    + 2.0 * (-2.0 * r_prime[1] * a + 2.0 * b * r_prime[1]) /  ((a-b)*sqrt(a)*b) * K;

        // volumeElem_field[0] = dist[0]/(dist.two_norm()*dist.two_norm()*dist.two_norm());
        // volumeElem_field[1] = dist[1]/(dist.two_norm()*dist.two_norm()*dist.two_norm());
        // volumeElem_field *= 4.0 / sysParams.pi * sysParams.get_lambda2i() * it->geometry().volume();
         
        //  E_field[0] = sqrt(a+b) * E / (a-b) * dist[0];
         //  E_field[1] = sqrt(a+b) * ((-a * r_prime[1] + b * r[1]) * E
         //                  + (a-b) * r_prime[1] * K) / (a-b) / b;
         // tmpFlux = E_field * unitNormal;

         //  E_field[0] = 4.0*sqrt((a-b)/(a+b)) * E / sqrt((a-b)*(a-b)*(a-b))* dist[0];
         //  E_field[1] = 4.0*sqrt((a-b)/(a+b)) * ((-a * r_prime[1] + b * r[1]) * E
         //                  + (a-b) * r_prime[1] * K / sqrt((a-b)*(a-b)*(a-b))) / b;
         //          
         //  tmpFlux = E_field * unitNormal;
         //  tmpFlux *= -1.0 / (4.0 * sysParams.pi)
         //        * sysParams.get_lambda2i() * it->geometry().volume() * r_prime[1];

        // std::cout << "field: " << volumeElem_field << " value: " << value << " sinh: " << std::sinh(value) << " product: ";
        // Now get the counterion distribution and calculate flux
         switch (sysParams.get_salt())
         {
           case 0:
             volumeElem_field *= std::sinh(value);
             break;
           case 1:
             volumeElem_field *= std::exp(value); // Counterions have opposite sign!
             break;
         }
        volumeElem_field *= 1.0 / (2.0*sysParams.pi) * sysParams.get_lambda2i() * it->geometry().volume() * r_prime[1];
      //}

      // Now also get the surface contribution
      if (it->hasBoundaryIntersections() == true)
      {
          for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
          {
           if(ii->boundary() == true)
           {
              if (boundaryIndexToEntity[ii->boundarySegmentIndex()] > 1)
              {
                if (it != elemIt)
                {
                  r_prime = ii->geometry().center();
                  dist = r-r_prime;

                  // Spherical symmetry - working!
                  a = dist[0]*dist[0] + r[1]*r[1] + r_prime[1]*r_prime[1] + 2.0 * r[1] * r_prime[1];
                  b = 4.0 * r[1] * r_prime[1];
                  if (a > 1e-3)
                  {
                    k = sqrt(b/a);
                    E = gsl_sf_ellint_Ecomp (k, GSL_PREC_DOUBLE);
                    K = gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE);
                  }
                  else
                  {
                    E = sysParams.pi / 2.0;
                    K = sysParams.pi / 2.0;
                  }
                  surfaceElem_field[0] = -2.0 * dist[0] / ((a-b)*sqrt(a)) * E;
                  surfaceElem_field[1] = 2.0 * ( 2.0 * r_prime[1] * a - b * r_prime[1] - b * r[1]) / ((a-b)*sqrt(a)*b) * E
                           + 2.0 * (-2.0 * r_prime[1] * a + 2.0 * b * r_prime[1]) /  ((a-b)*sqrt(a)*b) * K;
                  surfaceElem_field *= 4 * sysParams.get_bjerrum() * ii->geometry().volume() * r_prime[1] * boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
                  // // Cylindrical symmetry - working!
                  // surfaceElem_field[0] = dist[0]/(dist.two_norm()*dist.two_norm()*dist.two_norm());
                  // surfaceElem_field[1] = dist[1]/(dist.two_norm()*dist.two_norm()*dist.two_norm());
                  // surfaceElem_field *=  4.0 / sysParams.pi * sysParams.get_bjerrum() * ii->geometry().volume() * r_prime[1] * boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
                }
                else
                {
                  surfaceElem_field = ii->centerUnitOuterNormal();
                  surfaceElem_field *= 0.5 * sysParams.get_bjerrum() * boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
                }
                // {
                //   // The contribution of the surface charge density to this element
                //   // j_s = 4 * pi * l_B * sigma
                //   surfaceElem_field = ii->centerUnitOuterNormal();
                //   surfaceElem_field *= 2.0 * sysParams.pi * sysParams.get_bjerrum() * boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
                //   r_prime = ii->geometry().center();
                //   dist = r-r_prime;
                //   // We want the field at the center of the element for testing purposes (there it equals the gradient)
                //   surfaceElem_field /= dist.two_norm() * dist.two_norm();
                // }
              }
           }
          }
      }
      summed_field += surfaceElem_field;
      // summed_field += volumeElem_field;
    } // End of loop over all elements contributing to this one

   coulomb_file << elemii->geometry().center() << " " << summed_field << std::endl;
   
   // Get the gradient for comparing the field
   Dune::FieldVector<Real,dim> E = gradient(gfs, elemIt, u);
   E *= -1.0;
   gradient_file << elemIt->geometry().center() << " " << E << std::endl;
  }
 }}}}
   gradient_file.close();   
   coulomb_file.close();                
}
