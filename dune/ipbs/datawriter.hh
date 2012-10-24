#ifndef IPBS_DATAWRITER_HH
#define IPBS_DATAWRITER_HH

/** @file
    @author Alexander Schlaich
    @brief Provides gnuplot readable output
    based on the Dune::Gnuplotwriter,
    but iterating over elements for element-data
*/

#include <string>
#include <iostream>
#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>

#include "sysparams.hh"
#include "boundary.hh"
extern SysParams sysParams;
extern std::vector<Boundary*> boundary;

/** \brief Writer for grid data in gnuplot format
    \tparam GridType the grid
    \tparam GridView Level- or LeafGridView
*/

template<class GridView, int codim=0>
class DataWriter {
    
  typedef typename GridView::Grid::ctype ctype;
  enum {dimworld = GridView::dimensionworld};
  enum {dim = GridView::dimension};

  public:

    DataWriter (const GridView & _gv) : gv(_gv), 
      communicator(gv.comm()) {}
    
    /** \brief Add cell data
        \param data An ISTL compliant vector type
        \param name associated with the data
        \param filename Filename for the output (*.dat)
    */

    template <class GFS, class DataContainer>
    void writeIpbsCellData(const GFS& gfs, const DataContainer& data, 
            const std::string& name, const std::string &filename,
            const std::stringstream& statusMessage = "")
      {
        std::string myFilename = filename_helper(filename);
        std::ofstream out;
        out.open (myFilename.c_str(), std::ios::out);
        out.precision(5);

        out << "# IPBS solution for a system given by " 
          << sysParams.get_meshfile() << std::endl << std::endl
          << "# kappa = " << sqrt(1./sysParams.get_lambda2i()) << std::endl
          << "# bjerrum = " << sysParams.get_bjerrum() << std::endl
          << "# epsilon = " << sysParams.get_epsilon() << std::endl;
        for (size_t i = 0; i < sysParams.get_npart(); i++)
          out << "# boundary " << i << ":\tsigma = " 
            << boundary[i]->get_charge_density() 
            << "\tepsilon = " << boundary[i]->get_epsilon() << std::endl;
        out << std::endl << statusMessage.str();
        //out << "final residual update:\t" 
        //  << "b.c. update " << sysParams.get_bcError()
        //  << " IC update " << sysParams.get_fluxError()
        //  << " in " << sysParams.get_n << " iterations" <<std::endl;
        out << std::endl;
        out << std::left << std::setw(dim*12+1) << "#coordinates\t";
        out << std::left << std::setw(12) << "potential\t";
        out << std::left << std::setw(dim*12+1) << "electric field\t";
        out << std::left << std::setw(12) << "|E|\t" << "total ion density";
        out << std::left << std::setw(12) << "element volume";
        out << std::endl << std::endl;

        typedef typename GridView::template Codim<0>::template Partition
            <Dune::Interior_Partition>::Iterator LeafIterator;
    
        // do output
        for (LeafIterator it = 
            gv.template begin<0,Dune::Interior_Partition>();
            it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
        {
            Dune::FieldVector<Real, dim> evalPos = 
              it->geometry().center();
            Dune::FieldVector<Real, dim> local = 
              it->geometry().local(evalPos);
            // construct a discrete grid function for access to solution
            typedef Dune::PDELab::DiscreteGridFunction
              <GFS,DataContainer> DGF;
            const DGF udgf(gfs, data);
            typedef typename DGF::Traits::RangeType RT;
            RT value;
            // evaluate the potential
            udgf.evaluate(*it, local, value);
            Dune::FieldVector<Real,dim> gradphi;
            typename Dune::PDELab::DiscreteGridFunctionGradient
              < GFS, DataContainer > grads(gfs,data);
            grads.evaluate(*it, local, gradphi);

            double density = ( sysParams.get_lambda2i() / 
                  (4.*sysParams.pi * sysParams.get_bjerrum()) 
                  * std::sinh(- value) );

            out << std::left << std::scientific << evalPos << "\t";
            out << std::left << value << "\t";
            out << std::left << gradphi << "\t";
            out << std::left << gradphi.two_norm() << "\t";
            out << std::left << density << "\t";
            out << std::left << it->geometry().volume() << "\n";
        }
        
        out.close();
      }
  
  private:
    const GridView &gv;
    
    /// The communicator decides weither to use MPI or fake
    typedef typename GridView::Traits::CollectiveCommunication CollectiveCommunication;
    const CollectiveCommunication & communicator;

    std::string filename_helper(const std::string &name)
    {
      // generate filename for process data
      std::ostringstream pieceName;
      if( communicator.size() > 1 )
      {
        pieceName << "s" << std::setfill( '0' ) << std::setw( 4 ) << communicator.size() << "-";
        pieceName << "p" << std::setfill( '0' ) << std::setw( 4 ) << communicator.rank() << "-";
      }
      pieceName << name << ".dat";
      return pieceName.str();
    }

};

#endif // IPBS_DATAWRITER_HH
