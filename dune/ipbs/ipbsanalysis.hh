#ifndef _IPBSANALYSIS_HH
#define _IPBSANALYSIS_HH

/** \file 
 * Analysis stuff for IPBS
 * \todo Doc Me!
*/


#include <dune/common/fvector.hh>

#include "maxwelltensor.hh"
#include "sysparams.hh"

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
    void forces(const U& u) const
    {
      // Here we once more loop over all elements on this node (need of the element information
      // for the gradient calculation) and integrate Maxwell stress tensor over the particles surface
      // (see Hsu06a, eq. 61)
  
      // Open output file for force on particles
      std::ofstream force_file, vector_force_file;
  
      vector_force_file.open (filename_helper("vector_forces").c_str(), std::ios::out);
      if (communicator.rank() == 0) {
        force_file.open ("forces.dat", std::ios::out);
      }
       
      // Do the loop for boundary type 2 (iterated b.c.)
      // TODO: implement that!
      for (size_t i = 0; i < sysParams.get_npart(); i++)
      {
        std::cout << std::endl << "Calculating the force acting on particle id " << i << std::endl;
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

    void surfacepot(const U& u, std::string filename) const
    {
      std::ofstream pot_file;

      if (gv.comm().rank() ==0) {
        pot_file.open(filename);
        
        typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
        DGF udgf(gfs,u);
        
        typedef typename DGF::Traits::RangeType RT;
  
        for (int i = 0; i < sysParams.get_npart(); i++)
        {
          int nElems = 0;
          double sum = 0.;
          for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
                 	it!=gv.template end<0,Dune::Interior_Partition>(); ++it) {
            if(it->hasBoundaryIntersections() == true) {
              for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii) {
                if(ii->boundary() == true) {
                  if (pgmap[ii->boundarySegmentIndex()] == i) {
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

    void E_ext(const U& u, std::string filename) const {
        std::ofstream E_ext_file;
        E_ext_file.open(filename);
        for (int i=0; i<ipbsolver.ipbsPositions.size(); i++) {
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
