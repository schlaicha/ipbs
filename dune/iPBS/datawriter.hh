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

/** \brief Writer for grid data in gnuplot format
    \tparam GridType the grid
    \tparam GridView Level- or LeafGridView
*/

template<class GridView>
class DataWriter {
    
  typedef typename GridView::Grid::ctype ctype;
  enum {dimworld = GridView::dimensionworld};

  public:
    DataWriter (const GridView & _gv, Dune::MPIHelper & _helper)
      : gv(gv), helper(_helper)
      {}
    
    /** \brief Add cell data
        \param data An ISTL compliant vector type
        \param name associated with the data
        \param filename Filename for the output (*.dat)
    */
    template <class DataContainer>
    void writeIpbsCellData(const DataContainer& data, const std::string & name, const std::string &filename)
      {
        std::string myFilename = filename_helper(filename);
      }
  
  private:
    typedef typename GridView::template Codim<0>::template Partition
      <Dune::Interior_Partition>::Iterator LeafIterator;
    const GridView &gv;
    const Dune::MPIHelper &helper;

    std::string filename_helper(const std::string &name)
    {
      // generate filename for process data
      std::ostringstream pieceName;
      if( helper.size() > 1 )
      {
        pieceName << "s" << std::setfill( '0' ) << std::setw( 4 ) << helper.size() << ":";
        pieceName << "p" << std::setfill( '0' ) << std::setw( 4 ) << helper.rank() << ":";
      }
      pieceName << name << ".dat";
      return pieceName.str();
    }

};

#endif // IPBS_DATAWRITER_HH
