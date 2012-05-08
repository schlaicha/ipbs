#ifndef _BCTYPE_HH
#define _BCTYPE_HH

#include "p0layout.hh"
#include "ipbsolver.hh"
#include "boundary.hh"


/** \brief boundary grid function selecting boundary conditions
 * 
 * \arg 0 means Neumann
 * \arg 1 means Dirichlet
 */

template<typename GV, typename PGMap>
class BCType : public Dune::PDELab::BoundaryGridFunctionBase<
        Dune::PDELab::BoundaryGridFunctionTraits<GV,int,1,
        Dune::FieldVector<int,1> >,BCType<GV,PGMap> >
{
public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,int,1,Dune::FieldVector<int,1> > Traits;

  //! construct from grid view
  BCType (const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! return bc type at point on intersection
  template<typename I>
  inline void evaluate (I& i, const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {

    /** use physical index to determine B.C.
    *   values are specified in .geo file and compared to boundary type in the configuration
    *   \arg 0 is for Dirichlet surfaces
    *   \arg 1 for Neumann
    *   \arg 2 for IPBS iterated boundaries  */

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    if (boundary[physgroup_index]->get_type() == 0)
      y = 1;
    else
      y = 0;
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

private:

  const GV&    gv;
  const PGMap& pg;
};

#endif // _BCTYPE_HH
