/** \file boundaries.hh - Adjust the boundary conditions for the iPBS scheme 

    \brief    determine the boundary conditions
    
          define what are inner, boundary (Neumann)
          and extended boundary (Dirichlet) elements
          and set to IPBS / reference cases
    */

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

#ifndef _P0LAYOUT_H
#define _P0LAYOUT_H
#include "p0layout.hh"
#endif

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
    switch ( physgroup_index )
    {
      case 0:   y = 1.0;  break; // Set Dirichlet
      case 1:   y = 0.0;  break; // Set Neumann
      case 2:   y = 0.0;  break; // Set Neumann for iPBS
      default : y = 1.0;  break; // Default is Dirichlet
    }
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

template<typename GV, typename RF, typename PGMap, typename StoreMap>
class BoundaryFlux
  : public Dune::PDELab::BoundaryGridFunctionBase<
           Dune::PDELab::BoundaryGridFunctionTraits<GV,RF,1,
           Dune::FieldVector<RF,1> >, BoundaryFlux<GV,RF,PGMap,StoreMap> >
{
public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename Traits::GridViewType::Grid::ctype ctype;

  //! constructor
  BoundaryFlux(const GV& gv_, const PGMap& pg_, const double fluxValues[],
      StoreMap storeMap_, const int offset_) : gv(gv_), pg(pg_), storeMap(storeMap_),
      offset(offset_), mapper(gv)
  {fluxContainer = fluxValues;}

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
    switch ( physgroup_index )
    {
      case 1:   y = 0.0;  break; // Set Neumann for symmetry
      case 2:   // Set Neumann for iPBS
      {
        int mappedIndex = storeMap.find(mapper.map(*i.inside()))->second + offset;
		    y = fluxContainer[mappedIndex]; 
 //       std::cout << "At element " << mapper.map(*i.inside()) << " stored at pos " << storeMap.find(mapper.map(*i.inside()))->second+offset << " y= " << y << std::endl;
        break;
		  }
      //! /todo missing check
      default : y = 0.0; std::cerr << "Boundary flux detection failed!" << std::endl; break;
    }
    return;
  }

private:

  const GV&    gv;
  const PGMap& pg;
  const double* fluxContainer;
  StoreMap storeMap;
  const int offset;
  const Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
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
  template<typename I, typename E>
  inline void evaluate(I& i, E& e,
                       typename Traits::RangeType& y) const
  {
     /** use physical index to determine B.C.
    *   values are specified in .geo file
    *   \arg 0 is for Dirichlet surfaces
    *   \arg 1 for Neumann
    *   \arg 2 for surface flux given by charge density */

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    switch ( physgroup_index )
    {
      case 1:   y = 0.0;  break; // Set Neumann
      case 2:   {// Set Neumann for Reference
		  switch (sysParams.get_symmetry() )
		  {
		    case 1:  y = 1.0 * sysParams.get_charge_density()
		               * sysParams.get_bjerrum() * 2.0 * sysParams.pi;
			     break;
		    case 2:  y = 1.0 * sysParams.get_charge_density()  * sysParams.get_bjerrum();
		             break;
        //! \todo missing check
        default : y = 0.0; std::cerr << "Boundary flux detection failed!" << std::endl; break;
		  }
		break;}
      //! \todo missing check
      default : y = 0.0; std::cerr << "Boundary flux detection failed!" << std::endl; break;
    }
    return;
  }
  
  private:

  const GV&    gv;
  const PGMap& pg;
};
