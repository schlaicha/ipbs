/** \brief Calculate the force between colloids

    Some more description :-)
*/

/*!
   \param
*/

#include <fstream>

template <typename GV, typename GFS, typename U>
void force(const GV& gv, const std::vector<int>& boundaryIndexToEntity,
    const GFS& gfs, const U& u)
{
  const int dim = GFS::Traits::GridViewType::dimension;

  // Here we once more loop trough all elements (as for the measurement of force performance is
  // not important anyway) and integrate Maxwell stress tensor over the particles surface
  // (see Hsu06a, eq. 61)

  typedef typename GV::template Codim<0>::template Partition
          <Dune::Interior_Partition>::Iterator LeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  
  // Open output file for force on particles
  std::ofstream force_file, force2_file, vector_force_file, vector_force2_file;
  force_file.open ("forces.dat", std::ios::out);
  force2_file.open ("forces2.dat", std::ios::out);
  vector_force_file.open ("vector_forces.dat", std::ios::out);
  vector_force2_file.open ("vector_forces2.dat", std::ios::out);

  
  // Do the loop for all boundaryIDs > 1 (all colloids)
  for (int i = 2; i < sysParams.get_npart()+2; i++)
  {
    std::cout << "Calculating the force actin on particle " << i << std::endl;
    
    //double summed_force_x = 0;
    //double summed_force2 = 0;
    Dune::FieldVector<Real, dim> F(0);
    Dune::FieldVector<Real, dim> F2(0);
    // loop over elements on this processor
    for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
              	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
    {
     if(it->hasBoundaryIntersections() == true)
      {
        for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
        {
          if(ii->boundary() == true)
          {
            if (boundaryIndexToEntity[ii->boundarySegmentIndex()] == i) // check if IPBS boundary
            {
              Dune::FieldVector<Real, dim> normal = ii->centerUnitOuterNormal();
              Dune::FieldVector<Real, dim> forcevec;
              normal *= -1.0 * ii->geometry().volume() * ii->geometry().center()[1]
                       * 2.0 * sysParams.pi; // Surface normal
              Dune::FieldMatrix<Real, GFS::Traits::GridViewType::dimension, 
                  GFS::Traits::GridViewType::dimension>
                      sigma = maxwelltensor(gfs, it, ii, u);
              sigma.umv(normal, F);
              sigma.mv(normal, forcevec);
              forcevec *= 2.0 * sysParams.pi * ii->geometry().center()[1];
              //double force_x = forcevec[0];
              //summed_force_x += force_x;
              // Get the E-Field (-grad Phi / (4*pi))
              Dune::FieldVector<Real,dim> tmp = gradient(gfs, it, u, ii->geometry().center());
              tmp *= -1. / (4. * sysParams.pi);
              // Calculate F = q * E
              tmp *= ii->geometry().volume() * ii->geometry().center()[1] 
                      * boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
              // Integrate in theta
              tmp *= 2.0 * sysParams.pi;              
              vector_force_file << ii->geometry().center() << " " << forcevec << std::endl;
              vector_force2_file << ii->geometry().center() << " " << tmp << std::endl;
              F2 += tmp;
            }
          }
        }
      }
    }

  // To get SI units divide by 4*PI (E = - 1/(4*pi) grad phi)
  // F /= 4*sysParams.pi;
  // F2 /= 4*sysParams.pi;

  // MPI_Allreduce(MPI_IN_PLACE, &F, dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

    force_file << i << " " << F << std::endl;
    force2_file << i << " " << F2 << std::endl;
  }
  force_file.close();
  force2_file.close();
  vector_force_file.close();
  vector_force2_file.close();
}
