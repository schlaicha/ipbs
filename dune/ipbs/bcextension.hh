#ifndef _BCEXTENSION_HH
#define _BCEXTENSION_HH

#include <dune/common/exceptions.hh>

#include "p0layout.hh"
#include "ipbsolver.hh"
#include "boundary.hh"


/** \brief Set Dirichlet B.C.
*
*   This sets the potential (Dirichlet B.C.) at the boundaries */

template<typename GV, typename RF, typename PGMap>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
           BCExtension<GV,RF,PGMap> >
{
public :

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename GV::ctype ctype;

  //! construct from grid view
  BCExtension(const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  const inline bool global_on_intersection(Dune::FieldVector<ctype, GV::dimensionworld> 
          integrationPointGlobal, IntersectionIterator& ii ) const 
  {

    if (GV::dimensionworld == 2) {
      Dune::FieldVector<ctype,GV::dimensionworld> p_vec = 
          ii->geometry().corner(1) - ii->geometry().corner(0);
      p_vec/=p_vec.two_norm();
      Dune::FieldVector<ctype,GV::dimensionworld> p2 = p_vec;
      p2 *= ((integrationPointGlobal - ii->geometry().corner(0))*p_vec);
      Dune::FieldVector<ctype,GV::dimensionworld> dist_vec = 
          p2 - (integrationPointGlobal - ii->geometry().corner(0));

      if (dist_vec.two_norm() < 1e-9)
        return true;
      return false;
    }
  //  else DUNE_THROW(Dune::NotImplemented,"Dirichlet interpolation for 3d still has to be done");
    else 
        return false;
  }


  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
     //! Set value for potential at outer domain boundaries

    Dune::FieldVector<ctype,GV::dimensionworld> integrationPointGlobal =  e.geometry().global(xlocal);
    int physgroup_index = -1;
    for (IntersectionIterator ii = gv.ibegin(e); ii != gv.iend(e) ; ++ii) {
      if (ii->boundary()) {
        if (global_on_intersection(integrationPointGlobal, ii)) {
          physgroup_index = pg[ii->boundarySegmentIndex()];
        }
      } 
      else if (ii->neighbor())     // !ii->boundary(), so we know this is a valid intersection
                                    // with inner, overlap or ghost entity
      {
        typename GV::Traits::Grid::template Codim<0>::EntityPointer o( ii->outside() );
        for  (IntersectionIterator ii2 = gv.ibegin(*o); ii2 != gv.iend(*o) ; ++ii2) {
          if (ii2->boundary()) {
            if (global_on_intersection(integrationPointGlobal, ii2)) {
              physgroup_index = pg[ii2->boundarySegmentIndex()];
            }
          }
        }
      }
    }
    if (physgroup_index > 0 )
      if (boundary[physgroup_index]->get_type() == 0)
        y = boundary[physgroup_index]->get_potential();
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView() { return gv; }

private :

  const GV& gv;
  const PGMap& pg;
};

#endif // _BCEXTENSION_HH
