#ifndef GRADIENT_H
#define GRADIENT_H
#include "gradient.hh"
#endif

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
  maxwelltensor(const GFS& gfs, const Iterator& it, 
      const Dune::FieldVector<double,GFS::Traits::GridViewType::dimension>& global, const U& u)
{
  const int dim = GFS::Traits::GridViewType::dimension;
  
  Dune::FieldVector<double,GFS::Traits::GridViewType::dimension> local = it->geometry().local(global);

  // construct a discrete grid function for access to solution
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  const DGF udgf(gfs, u);
  typedef typename DGF::Traits::RangeType RT;
  RT value;
  // evaluate the potential
  udgf.evaluate(*it, local, value);


  // Get the E-Field (-grad Phi)
  //Dune::FieldVector<Real,dim> gradphi = gradient(gfs, it, u, global);
  Dune::FieldVector<Real,dim> gradphi;
  typename Dune::PDELab::DiscreteGridFunctionGradient< GFS, U > grads(gfs,u);
  grads.evaluate(*it, local, gradphi);

  // Build the stress tensor
  Dune::FieldMatrix<Real, dim, dim> res;
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
    {
      res[i][j] = 1.0/(4.0 * sysParams.pi * sysParams.get_bjerrum()) * gradphi[i] * gradphi[j];
    }

  // for the correct force on the particles not only electrostatic but also osmotic pressure has to be included!
  // this is in principle proportional to n_+ + n_- - 2c_s = 2 * ( cosh (phi) - 1 )
  for (int i = 0; i < dim; i++)
  {
    //res[i][i] -= 1.0/(4.0 * sysParams.pi * sysParams.get_bjerrum()) * .5 * gradphi.two_norm2();
    res[i][i] -= 1.0/(4.0 * sysParams.pi * sysParams.get_bjerrum()) * ( sysParams.get_lambda2i() * (std::cosh(value) - 1.0) + .5 * gradphi.two_norm2());
  }
  return res;
}
