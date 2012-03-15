/** \file boundaries.hh - Adjust the boundary conditions for the iPBS scheme 

    \brief    determine the boundary conditions
    
          define what are inner, boundary (Neumann)
          and extended boundary (Dirichlet) elements
          and set to IPBS / reference cases
    */

#ifndef _BOUNDARIES_HH
#define _BOUNDARIES_HH

#include<dune/common/fvector.hh>
#include "p0layout.hh"
#include "ipbsolver.hh"

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
    }
  }

  inline const typename Traits::GridViewType& getGridView() { return gv; }

private:

  const GV& gv;
  const Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
  const PGMap& pg;
};

// ============================================================================

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
    *   values are specified in .geo file
    *   \arg 0 is for Dirichlet surfaces
    *   \arg 1 for Neumann
    *   \arg 2 for iPBS iterated boundaries  */

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    if (physgroup_index == 0)// || physgroup_index > 1)
      y = 1;
    else
      y = 0;
    // if (physgroup_index > 0)
    //   y = 0;  // All others are Neumann boundaries
    // else
    //   y = 1;  // only zero is Dirichlet
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

private:

  const GV&    gv;
  const PGMap& pg;
};

// ============================================================================

/** \brief Set Dirichlet B.C.
*
*   This sets the potential to zero (Dirichlet B.C.) at the system outer boundaries */

template<typename GV, typename RF, typename PGMap>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
           BCExtension<GV,RF,PGMap> >
{
public :

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename GV::IntersectionIterator IntersectionIterator;

  //! construct from grid view
  BCExtension(const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    //! Set value for potential at outer domain boundaries
    // const int dim = Traits::GridViewType::Grid::dimension;
    // typedef typename Traits::GridViewType::Grid::ctype ctype;
    // Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal) ;
    // if (x.two_norm() < 21.)
    //   y = 1.0;
    // else
      y = 0.0;
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView() { return gv; }

private :

  const GV& gv;
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
    if (physgroup_index > 1)
      y = ipbsolver.get_flux(i);
    else
      y = 0;
    // std::cout << i.geometry().center() << " booundary " << physgroup_index << " flux " << y << std::endl;
    return;
  }

private:

  const GV&    gv;
  const PGMap& pg;
  const Ipbssolver& ipbsolver;
};


// ============================================================================
//! \brief class defining radiation and Neumann boundary conditions for reference solution

template<typename GV, typename RF, typename PGMap>
class RefBoundaryFlux
  : public Dune::PDELab::BoundaryGridFunctionBase<
           Dune::PDELab::BoundaryGridFunctionTraits<GV,RF,1,
           Dune::FieldVector<RF,1> >, RefBoundaryFlux<GV,RF,PGMap> >
{
public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename Traits::GridViewType::Grid::ctype ctype;

  //! constructor
  RefBoundaryFlux(const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! evaluate flux boundary condition
  template<typename I>
  inline void evaluate(I& i,
                       typename Traits::RangeType& y) const
  {
     /** use physical index to determine B.C.
    *   values are specified in .geo file
    *   \arg 0 is for Dirichlet surfaces
    *   \arg 1 for Neumann
    *   \arg 2 for surface flux given by charge density */

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    if (physgroup_index > 1)
           y = boundary[physgroup_index-2]->get_charge_density()  * 4.0 * sysParams.pi * sysParams.get_bjerrum() / boundary[physgroup_index-2]->get_epsilon();
    else y = 0;
    return;
  }
  
  private:

  const GV&    gv;
  const PGMap& pg;
};

// ============================================================================

//! \brief class defining radiation and Neumann boundary conditions

template<typename GV, typename RF, typename PGMap, typename IndexLookupMap>
class OldBoundaryFlux
  : public Dune::PDELab::BoundaryGridFunctionBase<
           Dune::PDELab::BoundaryGridFunctionTraits<GV,RF,1,
           Dune::FieldVector<RF,1> >, OldBoundaryFlux<GV,RF,PGMap, IndexLookupMap> >
{
public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename Traits::GridViewType::Grid::ctype ctype;

  //! constructor
  OldBoundaryFlux(const GV& gv_, const PGMap& pg_, const double fluxValues[],
      const IndexLookupMap& indexLookupMap_, const int offset_) : gv(gv_), pg(pg_),
      indexLookupMap(indexLookupMap_), offset(offset_), mapper(gv)
  {fluxContainer = fluxValues;}

  //! evaluate flux boundary condition
  template<typename I>
  inline void evaluate(I& i, typename Traits::RangeType& y) const
  {
    /** use physical index to determine B.C.
* values are specified in .geo file
* \arg 0 is for Dirichlet surfaces
* \arg 1 for Neumann
* \arg 2 for iPBS iterated boundaries */

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    if (physgroup_index > 1)
    {
      int mappedIndex = indexLookupMap.find(mapper.map(*i.inside()))->second + offset;
y = fluxContainer[mappedIndex];
    }
    else y = 0;
    return;
  }

private:

  const GV& gv;
  const PGMap& pg;
  const double* fluxContainer;
  const IndexLookupMap& indexLookupMap;
  const int offset;
  const Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
};


#endif  // _BOUNDARIES_HH
