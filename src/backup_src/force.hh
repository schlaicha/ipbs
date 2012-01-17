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
  std::ofstream force_file, vector_force_file;
  std::ofstream force2_file, vector_force2_file;
  std::ofstream surface_potential_file;;
  std::ofstream surface_field_file;;
  force_file.open ("forces.dat", std::ios::out);
  force2_file.open ("forces2.dat", std::ios::out);
  surface_potential_file.open ("surface_potential.dat", std::ios::out);
  surface_field_file.open ("surface_field.dat", std::ios::out);
  vector_force_file.open ("vector_forces.dat", std::ios::out);
  vector_force2_file.open ("vector_forces2.dat", std::ios::out);
  surface_potential_file << "# Averaged surface potential over boundary of type" << std::endl;
  surface_field_file << "# Averaged surface field over boundary of type" << std::endl;
  
  // Do the loop for all boundaryIDs > 1 (all colloids)
  for (int i = 2; i < sysParams.get_npart()+2; i++)
  {
    std::cout << "Calculating the force actin on particle " << i << std::endl;
    
    double summed_potential = 0;
    double elemCounter = 0;
    Dune::FieldVector<Real, dim> F(0);
    Dune::FieldVector<Real, dim> F2(0);
    //Dune::FieldVector<Real, dim> osmotic_force(0);

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
              normal *= -1.0 * ii->geometry().volume(); // Surface normal
              Dune::FieldMatrix<Real, GFS::Traits::GridViewType::dimension, 
                  GFS::Traits::GridViewType::dimension>
                      sigma = maxwelltensor(gfs, it, ii, u);
              //sigma.umv(normal, F);
              sigma.mv(normal, forcevec);
              forcevec *= 2.*sysParams.pi*ii->geometry().center()[1]; // integration in theta
              F += forcevec;

              // Testing the osmotic part
              // -----------------------
              // typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
              // const DGF testdgf(gfs, u);
              // typedef typename DGF::Traits::RangeType RT;
              // RT phi_local;
              // // evaluate the potential
              // testdgf.evaluate(*it, it->geometry().local(ii->geometry().center()), phi_local);

              // osmotic_force=ii->centerUnitOuterNormal();
              // osmotic_force*=-1;
              // osmotic_force*=cosh(phi_local)-1.;
              // osmotic_force*=ii->geometry().center()[1];
              // F3 += osmotic_force;

              // double force_x = forcevec[0];
              // summed_force_x += force_x;
              // Get the E-Field (-grad Phi))

              // Testing the gradient
              // --------------------
              Dune::FieldVector<Real,dim> tmp = gradient(gfs, it, u, ii->geometry().center());
              Dune::FieldVector<Real,dim> tmp2 = normal;
              // Forget about tmp2 here!!! We don't know!!!!
              tmp *= -1.; // E = -grad phi
              tmp2 *= -2*sysParams.pi*boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
              tmp -= tmp2;
              // surface_field_file << ii->geometry().center() << tmp << std::endl;
              // Calculate F = q * E = sigma * A * E
              tmp *= ii->geometry().volume() * boundary[boundaryIndexToEntity[ii->boundarySegmentIndex()]-2]->get_charge_density();
              // Integrate in theta
              tmp *= 2.0 * sysParams.pi * ii->geometry().center()[1];               

              vector_force_file << ii->geometry().center() << " " << forcevec << std::endl;
              vector_force2_file << ii->geometry().center() << " " << tmp << std::endl;
              F2 += tmp;
              
              // Evaluate the averaged surface potential
              // ---------------------------------------
              typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
              DGF udgf(gfs,u);
              typedef typename DGF::Traits::RangeType RT;
              RT value;
              udgf.evaluate(*it, it->geometry().local(ii->geometry().center()), value);
              elemCounter++;
              summed_potential += value;
            }
          }
        }
      }
    }

  // Build arithmetic average of surface potential
  double surface_potential = summed_potential / elemCounter;
  surface_potential_file << i << " " << surface_potential << std::endl;

  // MPI_Allreduce(MPI_IN_PLACE, &F, dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

    force_file << i << " " << F << std::endl;
    force2_file << i << " " << F2 << std::endl;
  }
  force_file.close();
  force2_file.close();
  vector_force_file.close();
  vector_force2_file.close();
  surface_potential_file.close();
  surface_field_file.close();
}
