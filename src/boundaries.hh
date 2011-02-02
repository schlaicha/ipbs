// boundaries.hh - Adjust the boundary conditions for the iPBS scheme

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

#include <limits>

// ============================================================================

// layout for codim0 data
template <int dim>
struct P0Layout
{
  bool contains(Dune::GeometryType gt)
  {
    if ( gt.dim() == dim ) return true;
    return false;
  }
};

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

  // evaluate the inner elementa
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
      case 1  : y = 1.0;  break;
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
    // use with global ccordinates
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = i.geometry().global(xlocal);
    
    // outer borders should be set to Zero Dirichlet BC (Open Boundary)
    //if (x.two_norm() < 4.7)  y = 0; // Neumann in inner region
    //else
	    y=1;
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

private:

  const GV&    gv;
  const PGMap& pg;
};

// ============================================================================

// This sets the initial potential at the particle's surface and zero Dirichlet at the domain boundaries

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
  BCExtension_init(const GV& gv_, const PGMap& pg_) : gv(gv_), mapper(gv), pg(pg_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    // evaluate with global cordinates
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

    //if (x.two_norm() < 4.7 && e.hasBoundaryIntersections() == true)
    //if (e.hasBoundaryIntersections() == true)
    if (x.two_norm() < 1.0+2E-1)
	y = sysParams.get_phi_init();	// set particles potential
    else
      y = 0.0;				// set domain boundary
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView() { return gv; }

private :

  const GV& gv;
  const Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
  const PGMap& pg;
};

// ============================================================================

// This sets the potential at the particle's surface and zero Dirichlet at the domain boundaries during iteration

template<typename GV, typename RF, typename PGMap>
class BCExtension_iterate
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
           BCExtension_iterate<GV,RF,PGMap> >
{
public :

  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;

  //! construct from grid view
  BCExtension_iterate(const GV& gv_, const PGMap& pg_, const DGF &udgf_) : gv(gv_), mapper(gv), pg(pg_), udgf(udgf_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    // evaluate with global ccordinates
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    typename Traits::RangeType y_old;	// store potential at actual element for SOR step
    udgf.evaluate(e,xlocal,y_old);



    //if (x.two_norm() < 3.0 && e.hasBoundaryIntersections() == true)	// i.e. is particle
    if (x.two_norm() < 1.0+2E-1)
    {
	// Take a SOR step towards the correct solution (see paper)
	// calculate B.C. according to eq. 2 (paper)
	y = 0;
	//calculate_phi(gv, udgf);

    typedef typename GV::template Codim<0>::Iterator LeafIterator;
	for (LeafIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
	{
	  typename Traits::RangeType value;
	  Dune::FieldVector<ctype,dim> x_prime = it->geometry().global(xlocal);
	  if (x != x_prime)	// make sure we have no zero division
	  {
	    udgf.evaluate(*it,it->geometry().local(it->geometry().center()),value);
	    Dune::FieldVector<ctype,dim> dist = x - x_prime;
	    //std::cout << "value = " << value << std::endl;

  	    y += std::sinh(value) / dist.two_norm() * it->geometry().volume();
	  }
	}
	y *= sysParams.get_lambda2i() / sysParams.get_bjerrum() / (4*3.14);
	// Currently we use a fixed point charge at origin (sphere case)
	y += sysParams.get_bjerrum() * sysParams.get_charge() / (sysParams.get_epsilon() * x.two_norm());
	// SOR step
	y = sysParams.alpha * y + (1.0 - sysParams.alpha) * y_old;
	// Calculate error
	sysParams.add_error(y);
	//if (y!=y || y == std::numeric_limits<double>::infinity())
	//    y = -1.0;
    //    std::cout << "x[1] = " << x.vec_access(0) << "\tx[2] = " << x.vec_access(1) << "\ty = " << y << std::endl;
    }
    else
      y = 0.0;			// Dirichlet B.C. for domain
    // std::cout << "x[1] = " << x.vec_access(0) << "\tx[2] = " << x.vec_access(1) << "\ty = " << y << std::endl;
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView() { return gv; }

private :

  const GV& gv;
  const Dune::MultipleCodimMultipleGeomTypeMapper<GV,P0Layout> mapper;
  const PGMap& pg;
  const DGF& udgf;
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

  // constructor
  BoundaryFlux(const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  // evaluate flux boundary condition
  template<typename I>
  inline void evaluate(I& i, const typename Traits::DomainType& xlocal,
                       typename Traits::RangeType& y) const
  {
    y = sysParams.get_sigma_init();
    return;
  }

private:

  const GV&    gv;
  const PGMap& pg;
};
