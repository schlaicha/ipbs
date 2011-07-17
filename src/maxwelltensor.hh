/** \brief Get the Maxwell stress tensor
 
    Some more description :-)
*/

/*!
   \param gfs the grid function space
   \param it iterator refering to the element where the evaluation should be done
   \param u solution vector
*/


template <typename GFS, typename Iterator, typename Intersection, typename U>
Dune::FieldMatrix<Real, GFS::Traits::GridViewType::dimension,
  GFS::Traits::GridViewType::dimension>
  maxwelltensor(const GFS& gfs, const Iterator& it, const Intersection& ii, const U& u)
{
  const int dim = GFS::Traits::GridViewType::dimension;
  
  // construct a discrete grid function for access to solution                                                                                                                                                              
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;                                                                                                                                                                    
  const DGF udgf(gfs, u);  
  typedef typename DGF::Traits::RangeType RT;
  RT value;
  // evaluate the potential
  udgf.evaluate(*it, ii->geometry().center(), value);

  // Get the E-Field (grad Phi)
  Dune::FieldVector<Real,dim> E = gradient(gfs, it, u);
  E *= -1.0;

  // Build the stress tensor
  Dune::FieldMatrix<Real, dim, dim> res;
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      res[i][j] = E[i] * E[j];

  // for the correct force on the particles not only electrostatic but also osmotic pressure has to be included!
  // this is in principle proportional to n_+ + n_- - 2c_s = 2 * ( cosh (phi) - 1 )
  for (int i = 0; i < dim; i++)
    res[i][i] -= -2 * sysParams.get_lambda2i() * std::cosh(value) + .5 * E.two_norm2();
  return res;
}
