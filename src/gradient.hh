/** \brief Compute the gradient of the potential at the actual element

**/

template <typename GFS, class Iterator, typename U>
Dune::FieldVector<double,GFS::Traits::GridViewType::dimension>
gradient(const GFS& gfs, const Iterator& it, const U& u,
    const Dune::FieldVector<double,GFS::Traits::GridViewType::dimension>& global)
{
  const int dim = GFS::Traits::GridViewType::dimension;
  const int dimw = GFS::Traits::GridViewType::dimensionworld;
  
  typedef typename GFS::LocalFunctionSpace LFSU;
  typedef typename LFSU::Traits::FiniteElementType::
	  Traits::LocalBasisType::Traits::RangeFieldType RF;
  typedef typename LFSU::Traits::FiniteElementType::
  	Traits::LocalBasisType::Traits::JacobianType JacobianType;
  typedef typename LFSU::Traits::FiniteElementType::
	  Traits::LocalBasisType::Traits::DomainFieldType DF;
  typedef typename LFSU::Traits::SizeType size_type;
  
  // get local function space
  LFSU lfsu(gfs);
  lfsu.setup(gfs);

  // Bind Local Function Space to this element
  lfsu.bind(*it);

  // Vector container for storing the coefficients on the actual element (same as U...)
  typename GFS::template VectorContainer<Real>::Type x_s(gfs,0.0);
  // get coefficients of this element
  lfsu.vread(u,x_s);

  // evaluate gradient of basis functions on reference element
  std::vector<JacobianType> js(lfsu.size());
  //lfsu.finiteElement().localBasis().evaluateJacobian(it->geometry().center(),js);
  lfsu.finiteElement().localBasis().evaluateJacobian(global,js);

  // transform gradients from reference element to real element
  const Dune::FieldMatrix<DF,dimw,dim>
  //        jac = it->geometry().jacobianInverseTransposed(it->geometry().center());
          jac = it->geometry().jacobianInverseTransposed(global);
  std::vector<Dune::FieldVector<RF,dim> > gradphi(lfsu.size());
  for (size_type i=0; i<lfsu.size(); i++)
    jac.mv(js[i][0],gradphi[i]);

  // compute gradient of u
  Dune::FieldVector<RF,dim> gradu(0.0);
  for (size_type i=0; i<lfsu.size(); i++)
    gradu.axpy(x_s[i],gradphi[i]);

  return gradu;
}
