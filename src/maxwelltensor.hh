/** \brief Get the Maxwell stress tensor
 
    Some more description :-)
*/

/*!
   \param gfs the grid function space
   \param it iterator refering to the element where the evaluation should be done
   \param u solution vector
*/


template <typename GFS, typename Iterator, typename U>
Dune::FieldMatrix<Real, GFS::Traits::GridViewType::dimension,
  GFS::Traits::GridViewType::dimension>
  maxwelltensor(const GFS& gfs, const Iterator& it, const U& u)
{
  const int dim = GFS::Traits::GridViewType::dimension;
  // Get the E-Field (grad Phi)
  Dune::FieldVector<Real,dim> E = gradient(gfs, it, u);
  E *= -1.0;

  // Build the stress tensor
  Dune::FieldMatrix<Real, dim, dim> res;
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      res[i][j] = E[i] * E[j];
  for (int i = 0; i < dim; i++)
    res[i][i] -= .5 * E.two_norm2();
  return res;
}

