// some container for extensions where there is no clear place to put


typedef struct {
	int NewtonMaxIteration; 
	int RefinementLevel; 
} Cmdparam;



// EXPLAIN WHAT YOU ARE DOING!!!
template<class G>
void boundary_info(const G& grid)
{
  typedef typename G::LeafGridView GV;
  const GV& gv = grid.leafView();

  const int dim = G::dimension;
  typedef typename G::ctype ctype;
  //Dune::FieldVector<ctype,dim> x = geo.global(xlocal);
 

  // make a mapper for codim 0 entities so we can decide if we deal with solvent, particles, or end-of-world
//  Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout> mapper(grid);

  // allocate a vector for storing this information
  std::vector<int> innerBoundary;

  // Iterate through all elements of codim 0
  typedef typename GV::template Codim<0>::Iterator ElementLeafIterator;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  int itcount=0;
  int iicount=0;
  //Dune::FieldVector<ctype,dim> globalpos;

  // loop over all entities (elements) of the grid
  for (ElementLeafIterator it = gv.template begin<0>(); it != gv.template end<0>(); ++it)
  {
	  // We are only interested in access to boundary segments
	  if (it->hasBoundaryIntersections()==true)
	  {
		//std::cout << "Boundary detected at Element " << itcount << std::endl;
	  	for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
	  	{
			if (ii->boundary()==true)
			{
				std::cout << "Boundary Intersection detected, ID " << ii->boundarySegmentIndex() << std::endl;
			
	  			//globalpos = ii->geometry();
				++iicount;
			}
	  	}
	  ++itcount;
	  }
  }
  std::cout << "visited " << itcount << " elements with boundaries" << std::endl;
  std::cout << "visited " << iicount << " intersections with boundaries" << std::endl;
}

