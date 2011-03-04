// boundaries.hh - Adjust the boundary conditions for the iPBS scheme

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

#ifndef _P0LAYOUT_H
#define _P0LAYOUT_H
#include "p0layout.hh"
#endif

//#include <limits>
//#include <cstdlib>


//===============================================================
// boundary classes - define what are inner, boundary (Neumann)
// and extended boundary (Dirichlet)elements
//===============================================================

// function defining the inner elements
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
 * 0 means Neumann
 * 1 means Dirichlet
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

    // use physical index to determine B.C.
    // values are specified in .geo file
    // 0 is for Dirichlet surfaces
    // 1 for Neumann
    // 2 for iPBS

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    //std::cout << "Set boundary type" << physgroup_index << " at " << i.geometry().center() << std::endl;
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

// This sets the potential to zero (Dirichlet B.C.) at the system outer boundaries

template<typename GV, typename RF, typename PGMap>
class BCExtension_init
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
           BCExtension_init<GV,RF,PGMap> >
{
public :

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;

  //! construct from grid view
  BCExtension_init(const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    // Set value for potential at outer domain boundaries
      y = 0.0;				// set domain boundary
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView() { return gv; }

private :

  const GV& gv;
  const PGMap& pg;
};


// ============================================================================

// function for defining radiation and Neumann boundary conditions

template<typename GV, typename RF, typename PGMap>
class BoundaryFlux
  : public Dune::PDELab::BoundaryGridFunctionBase<
           Dune::PDELab::BoundaryGridFunctionTraits<GV,RF,1,
           Dune::FieldVector<RF,1> >, BoundaryFlux<GV,RF,PGMap> >
{
public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename Traits::GridViewType::Grid::ctype ctype;

  // constructor
  BoundaryFlux(const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  // evaluate flux boundary condition
  template<typename I>
  inline void evaluate(I& i, typename Traits::RangeType& y,
		       std::vector<double>& fluxContainer,
		       const Mapper& mapper) const
  {
    // use physical index to determine B.C.
    // values are specified in .geo file
    // 0 is for Dirichlet surfaces
    // 1 for Neumann
    // 2 for iPBS

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    switch ( physgroup_index )
    {
      case 1:   y = 0.0;  break; // Set Neumann
      case 2:   // Set Neumann for iPBS
                y = fluxContainer[mapper.map(*i.inside())];
                break;
      default : y = 0.0;  break;
    }
    return;
  }

private:

  const GV&    gv;
  const PGMap& pg;
  //const Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
  //const Mapper mapper;
};


// ============================================================================

// function for defining radiation and Neumann boundary conditions

template<typename GV, typename RF, typename PGMap>
class BoundaryFluxRef
  : public Dune::PDELab::BoundaryGridFunctionBase<
           Dune::PDELab::BoundaryGridFunctionTraits<GV,RF,1,
           Dune::FieldVector<RF,1> >, BoundaryFluxRef<GV,RF,PGMap> >
{
public:

  typedef Dune::PDELab::BoundaryGridFunctionTraits<
          GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename Traits::GridViewType::Grid::ctype ctype;

  // constructor
  BoundaryFluxRef(const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  // evaluate flux boundary condition
  template<typename I, typename E>
  inline void evaluate(I& i, E& e,
                       typename Traits::RangeType& y) const
  {
    // use physical index to determine B.C.
    // values are specified in .geo file
    // 0 is for Dirichlet surfaces
    // 1 for Neumann
    // 2 for iPBS

    int physgroup_index = pg[i.intersection().boundarySegmentIndex()];
    switch ( physgroup_index )
    {
      case 1:   y = 0.0;  break; // Set Neumann
      case 2:   // Set Neumann for Reference
        {
         std::cout << "Set Neumann boundary at " << i.geometry().center() << std::endl;
         if ( sysParams.get_symmetry()  == 1) // "2D_cylinder"
                y = 1.0 * sysParams.get_charge_density()  * sysParams.get_bjerrum() * 2 * sysParams.pi ;
         else if (sysParams.get_symmetry() == 2) // "2D_sphere"
                y = 1.0 * sysParams.get_charge_density()  * sysParams.get_bjerrum();
         else 
         {
                y = 0;
                std::cerr << "IPBS WARNING:\tNo Symmetry specified (in flux evaluation). Will use 0." 
                        << std::endl;
         }
         break;
        }
      default : 
         {
                y = 0.0;  
                std::cerr << "IPBS WARNING:\tCould not evaluate flux B.C. Will use 0." 
                        << std::endl;
         }
    }
    return;
  }
  
  private:

  const GV&    gv;
  const PGMap& pg;
};
