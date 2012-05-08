#ifndef _BCTYPE_HH
#define _BCTYPE_HH

#include <dune/common/fvector.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/constraintsparameters.hh>

#include "boundary.hh"


//! \brief Parameter class selecting boundary conditions
template<typename PGMap>
class BCTypeParam
  : public Dune::PDELab::DirichletConstraintsParameters
{
public:
  //! Test whether boundary is Dirichlet-constrained
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
     /** use physical index to determine B.C.
    *   values are specified in .geo file and compared to boundary type in the configuration
    *   \arg 0 is for Dirichlet surfaces
    *   \arg 1 for Neumann
    *   \arg 2 for IPBS iterated boundaries  */

    int physgroup_index = pg[intersection.intersection().boundarySegmentIndex()];
    if (boundary[physgroup_index]->get_type() == 0)
      return true;
    return false;
  }
 
  //! constructor for obtaining BoundaryIndexToEntity
  BCTypeParam (const PGMap& pg_) : pg(pg_) {}

private:
  const PGMap& pg;

};

#endif // _BCTYPE_HH
