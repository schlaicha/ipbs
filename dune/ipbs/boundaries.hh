/** \file boundaries.hh - Adjust the boundary conditions for the iPBS scheme 

    \brief    determine the boundary conditions
    
          define what are inner, boundary (Neumann)
          and extended boundary (Dirichlet) elements
          and set to IPBS / reference cases
    */

#ifndef _BOUNDARIES_HH
#define _BOUNDARIES_HH

#include "p0layout.hh"
#include "ipbsolver.hh"
#include "boundary.hh"

// BC Type
#include "bctype.hh"

// (Dirichlet) Boundary Extension
#include "bcextension.hh"

#include <dune/common/exceptions.hh>
class BC_Error : public Dune::IOError {};

// ============================================================================
/** \brief class defining the inner elements */
template<typename GV, typename RF, typename PGMap>
class Regions
  : public Dune::PDELab::GridFunctionBase<
        Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
        Regions<GV,RF,PGMap> >
{
public:

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,1,
      Dune::FieldVector<RF,1> >, Regions<GV,RF,PGMap> > BaseT;

  // constructor
  Regions(const GV& gv_, const PGMap& pg_) : gv(gv_), mapper(gv), pg(pg_) {}

  // evaluate the inner elements
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // retrieve element index and corresponding material index on level 0
    typename GV::template Codim<0>::EntityPointer ep(e);
    while (ep->level() != 0) ep = ep->father();
    const int ei              = mapper.map(e);
    const int physgroup_index = pg[ei];

    // evaluate physical group map and set values accordingly
    switch ( physgroup_index )
    {
      default : y = 1.0;  break; // only one material here
                                 // i.e. only one domain of computation
    }
  }

  inline const typename Traits::GridViewType& getGridView() { return gv; }

private:

  const GV& gv;
  const Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
  const PGMap& pg;
};

// ============================================================================

//! \brief class defining radiation and Neumann boundary conditions

template<typename GV, typename RF, typename PGMap, class Ipbssolver>
class BoundaryFlux
  : public Dune::PDELab::BoundaryGridFunctionBase<
           Dune::PDELab::BoundaryGridFunctionTraits<GV,RF,1,
           Dune::FieldVector<RF,1> >, BoundaryFlux<GV,RF,PGMap, Ipbssolver> >
{
public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename Traits::GridViewType::Grid::ctype ctype;

  //! constructor
  BoundaryFlux(const GV& gv_, const PGMap& pg_, const Ipbssolver& ipbsolver_)
    : gv(gv_), pg(pg_), ipbsolver(ipbsolver_) {}

  //! evaluate flux boundary condition
  template<typename I>
  inline void evaluate(I& i, typename Traits::RangeType& y) const
  {
    /** use physical index to determine B.C.
    *   values are specified in .geo file
    *   \arg 0 is for Dirichlet surfaces
    *   \arg 1 for Neumann
    *   \arg 2 for iPBS iterated boundaries  */

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    if (boundary[physgroup_index]->get_type() == 2)
      y = -ipbsolver.get_flux(i);
    else if (boundary[physgroup_index]->get_type() == 1)
      y = -4. * sysParams.get_bjerrum() * sysParams.pi
          * boundary[physgroup_index]->get_charge_density();
    else {
      DUNE_THROW(BC_Error, "Something went wrong in assigning the BC!");
      y = 0;}
    return;
  }

private:

  const GV&    gv;
  const PGMap& pg;
  const Ipbssolver& ipbsolver;
};

#endif  // _BOUNDARIES_HH
