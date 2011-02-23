#ifndef _IPBS_COULOMB_FLUX_INTEGRAL
#define _IPBS_COULOMB_FLUX_INTEGRAL
#include <dune/common/fvector.hh>


// Function calculating Coulomb part of the surface flux
template <typename ct, int dim>
class CoulombFlux {
public:
  CoulombFlux() {}
  double operator() (const Dune::FieldVector<ct,dim>& dist,
		     const Dune::FieldVector<ct,dim>& normal) const
  {
    //std::cout << "dist: " << dist << "\tnormal: " << normal << "\tdist*normal" << (dist*normal) << std::endl;
    switch(dim){
      case 2:
          return 2*((dist * normal) / (dist.two_norm()*dist.two_norm()));
          break;
      case 3:
	  return ((dist * normal) / (dist.two_norm()*dist.two_norm()*dist.two_norm()));
	  break;
      default:
	  break;
    }
  }
};

#endif // _IPBS_COULOMB_FLUX_INTEGRAL