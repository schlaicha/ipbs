#ifndef _IPBSANALYSIS_HH
#define _IPBSANALYSIS_HH

/** \file 
 * Analysis stuff for IPBS
 * \todo Doc Me!
*/


#include <dune/common/fvector.hh>

#include "maxwelltensor.hh"
#include "sysparams.hh"
#include "boundary.hh"

extern SysParams sysParams;
extern std::vector<Boundary*> boundary;

template <class GV, class GFS, typename PGMap, class IPBSolver>
class IpbsAnalysis
{

  // Some typedef
  typedef typename GV::Grid::ctype ctype;
  static const int dim = GV::dimension;
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;


  public:

    IpbsAnalysis(const GV& _gv, const GFS& _gfs, const PGMap& _pgmap, const IPBSolver& ipbsolver_)
      : gv(_gv), gfs(_gfs), pgmap(_pgmap), communicator( gv.comm() ), ipbsolver(ipbsolver_) {} 

     typedef typename GV::template Codim<0>::template Partition
              <Dune::Interior_Partition>::Iterator LeafIterator;
     typedef typename GV::IntersectionIterator IntersectionIterator;


    // ------------------------------------------------------------------------
    /// Force computation
    // ------------------------------------------------------------------------
    void forces(const U& u, std::string fname) const
    {
      // Here we once more loop over all elements on this node (need of the element information
      // for the gradient calculation) and integrate Maxwell stress tensor over the particles surface
      // (see Hsu06a, eq. 61)
  
      // Open output file for force on particles
      std::ofstream force_file, vector_force_file;
  
      vector_force_file.open (filename_helper(sysParams.get_outname() + "_forceVec").c_str(), std::ios::out);
      if (communicator.rank() == 0) {
        force_file.open (fname.c_str(), std::ios::out);
      }
       
      // Do the loop for boundary type 2 (iterated b.c.)
      // TODO: implement that!
      for (size_t i = 0; i < sysParams.get_npart(); i++)
      {
        if (communicator.rank() == 0)
          std::cout << "Calculating the force acting on particle id " << i << std::endl;
        Dune::FieldVector<Real, dim> F(0);
  
        for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
                 	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
        {
          if(it->hasBoundaryIntersections() == true) {
            for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii) {
              if(ii->boundary() == true) {
                if ( pgmap[ii->boundarySegmentIndex()] == (int)i ) // check if IPBS boundary
                {
                  Dune::FieldVector<Real, dim> normal = ii->centerUnitOuterNormal();
                  Dune::FieldVector<Real, dim> forcevec;
                  normal *= -1.0 * ii->geometry().volume(); // Surface normal
                  Dune::FieldVector<Real, dim> evalPos = ii->geometry().center();
                  Dune::FieldMatrix<Real, dim, dim> sigma = maxwelltensor(gfs, it, evalPos, u);
                  //sigma.umv(normal, F);
                  sigma.mv(normal, forcevec);
                  if (sysParams.get_symmetry() > 0) {
                    // integration in theta
                    forcevec *= 2.*sysParams.pi*evalPos[1]; 
                  }
                  F += forcevec;
                  vector_force_file << evalPos << " " << forcevec << std::endl;
                }
              }
            }
          }
        }
        // Sum up force of all nodes
        communicator.barrier();
        communicator.sum(&F[0], F.dim());
        if (communicator.rank() == 0)
          force_file << i << " " << F << std::endl;        
      }
    
      // close output
      if (communicator.rank() == 0) {
        force_file.close();
        vector_force_file.close();
      }
    }

    // ------------------------------------------------------------------------
    /// Average surface potential (useful for many testcases)
    // ------------------------------------------------------------------------
    
    void surfacepot(const U& u, std::string filename) const
    {
      std::ofstream pot_file;

      if (gv.comm().rank() ==0) {
        pot_file.open( filename.c_str() );
        
        typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
        DGF udgf(gfs,u);
        
        typedef typename DGF::Traits::RangeType RT;
  
        for (unsigned int i = 0; i < sysParams.get_npart(); i++)
        {
          int nElems = 0;
          double sum = 0.;
          for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
                 	it!=gv.template end<0,Dune::Interior_Partition>(); ++it) {
            if(it->hasBoundaryIntersections() == true) {
              for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii) {
                if(ii->boundary() == true) {
                  if (pgmap[ii->boundarySegmentIndex()] == int(i)) {
                    Dune::FieldVector<Real, dim> evalPos = ii->geometry().center();
                    Dune::FieldVector<double,GFS::Traits::GridViewType::dimension> local =
                      it->geometry().local(evalPos);
                    RT value;
                    // evaluate the potential
                    udgf.evaluate(*it, local, value);
                    sum += value;
                    nElems++;
                  }
                }
              }
            }
          }
          sum /= nElems;
          //boundary[i]->set_res_surface_pot(sum);
          pot_file << i << " " << sum << std::endl;
        }
      }
    }

    // ------------------------------------------------------------------------
    /// Determine the system's electrostatic energy
    // ------------------------------------------------------------------------
    double energy(const U&u, std::string filename) const {

      std::ofstream en_file;

      if (gv.comm().rank() ==0) {
        en_file.open( filename.c_str() );
      }
      
      typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
        DGF udgf(gfs,u);
      typedef typename DGF::Traits::RangeType RT;

      double energy = 0;
      double senergy = 0;
      
      for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
        it!=gv.template end<0,Dune::Interior_Partition>(); ++it) {
        
        //Dune::FieldVector<Real, dim> evalPos = it->geometry().center();
        //Dune::FieldVector<Real, dim> local = it->geometry().local(evalPos);
        //RT value;
        //    
        //// evaluate the potential
        //udgf.evaluate(*it, local, value);

        //double charge_density = - ( sysParams.get_lambda2i() / 
        //      (4.*sysParams.pi * sysParams.get_bjerrum()) 
        //      * std::sinh( value) );
        //double local_energy = charge_density * value * it->geometry().volume();
        //if (sysParams.get_symmetry() > 0) {
        //  local_energy *= 2. * sysParams.pi * evalPos[1];
        //}
        //energy += local_energy;

        if(it->hasBoundaryIntersections() == true) {
          for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii) {
            if(ii->boundary() == true) {
                Dune::FieldVector<Real, dim> sevalPos = ii->geometry().center();
                Dune::FieldVector<double,GFS::Traits::GridViewType::dimension> slocal =
                  it->geometry().local(sevalPos);
                RT svalue;
                // evaluate the potential
                udgf.evaluate(*it, slocal, svalue);

                double local_senergy = 0;

                if (boundary[ pgmap[ii->boundarySegmentIndex()] ]->get_type() == 0) { 
                  // Dirichlet
                  DUNE_THROW(Dune::NotImplemented,"Dirichlet boundaries are not yet supported for energy calculations"); 
                  // 1/4/pi grad(phi) * n * phi
                  // just implement :P
                }
                else if (boundary[ pgmap[ii->boundarySegmentIndex()] ]->get_type() == 1) { 
                  // Neumann
                  local_senergy =  boundary[ pgmap[ii->boundarySegmentIndex()] ]
                    ->get_charge_density() * svalue * ii->geometry().volume();
                }
                else if (boundary[ pgmap[ii->boundarySegmentIndex()] ]->get_type() == 2) { 
                  // IPBS
                  local_senergy = ipbsolver.get_lcd(*ii) * svalue * ii->geometry().volume();
                }
                
                if (sysParams.get_symmetry() > 0) {
                  local_senergy *= 2. * sysParams.pi * sevalPos[1];
                }
                senergy += local_senergy;
            }
          }
        }

      }
      communicator.sum(&senergy,1);
      if (communicator.rank() == 0) {
        en_file << senergy << std::endl;
        en_file.close();
      }

      //std::cout << "Volume: " << energy << " Surface: " << senergy << " sum: " << senergy+energy << std::endl;
    }
    
    // ------------------------------------------------------------------------
    /// Print the external field at IPBS boundaries to file
    // ------------------------------------------------------------------------
    void E_ext(const U& u, std::string filename) const {
        std::ofstream E_ext_file;
        E_ext_file.open( filename.c_str() );
        for (unsigned int i=0; i<ipbsolver.ipbsPositions.size(); i++) {
            E_ext_file << ipbsolver.ipbsPositions[i] << " " << ipbsolver.E_ext[i] << std::endl;
        }
    }
  
  private:
     
    // ------------------------------------------------------------------------
    /// Filenamehelper for parallel output
    // ------------------------------------------------------------------------

    std::string filename_helper(const std::string &name) const
    {
      // generate filename for process data
      std::ostringstream pieceName;
      if( communicator.size() > 1 )
      {
        pieceName << "s" << std::setfill( '0' ) << std::setw( 4 ) << communicator.size() << ":";
        pieceName << "p" << std::setfill( '0' ) << std::setw( 4 ) << communicator.rank() << ":";
      }
      pieceName << name << ".dat";
      return pieceName.str();
    }

    const GV& gv;
    /// The grid function space
    const GFS& gfs;
    /// BoundaryIndexToEntity Mapper
    const PGMap &pgmap;

    /// The communicator decides weither to use MPI or fake
    typedef typename GV::Traits::CollectiveCommunication CollectiveCommunication;
    const CollectiveCommunication & communicator;

    const IPBSolver ipbsolver;
};


#endif  // _IPBSANALYSIS_HH
