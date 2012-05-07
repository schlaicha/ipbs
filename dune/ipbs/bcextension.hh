#ifndef _BCEXTENSION_HH
#define _BCEXTENSION_HH

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

  //! construct from grid view
  BCExtension(const GV& gv_, const PGMap& pg_) : gv(gv_), pg(pg_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    for (IntersectionIterator ii = gv.ibegin(e); ii != gv.iend(e); ++ii) {
      if(ii->boundary() == true) {
        int physgroup_index = pg[ii->boundarySegmentIndex()];
        if (boundary[physgroup_index]->get_type() == 0)
          y = boundary[physgroup_index]->get_potential();
        else 
          y=0;
      }
    }
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView() { return gv; }

private :

  const GV& gv;
  const PGMap& pg;
};

#endif // _BCEXTENSION_HH
