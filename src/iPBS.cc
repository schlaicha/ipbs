// iPBS.cc Read a gmsh file and solve PB eq.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// include application heaeders
#include"PB_operator.hh"
#include "parameters.hh"
#include "ipbs.hh"

#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

template <typename PositionVector>
double compute_pbeq(const double &u, const PositionVector &r)
{
	return (- sysParams.get_lambda2i() * std::sinh(u));
}

void solver (NEWTON &newton, SLP &slp)
{
	newton.apply();
	slp.apply();
}

void save(const DGF &udgf, const GV &gv)
{
  Dune::PDELab::FilenameHelper fn("output");
   {
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
    vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
  }
}

void calculate_phi(const GV &gv, const DGF &udgf)
{
  // Get the potential at the particle's surface
  std::cout << "Particle Boundaries at:" << std::endl;
  // define iterators
  typedef GV::Codim<0>::Iterator ElementLeafIterator;
  typedef GV::IntersectionIterator IntersectionIterator;
  int iicount=0;
  double phi=0;
 for (ElementLeafIterator it = gv.begin<0>(); it != gv.end<0>(); ++it)
  {
	if (it->hasBoundaryIntersections()==true)
	{
		for (IntersectionIterator ii = gv.ibegin(*it); ii != gv.iend(*it); ++ii)
		{
  			Traits::DomainType xlocal;
  			Traits::RangeType y;
			if (ii->boundary()==true)
			{
			   udgf.evaluate(*it,xlocal,y);
			   Dune::FieldVector<double,dim> x = it->geometry().global(xlocal);
			   if (x.two_norm() < 4.7)				
			   {    
				//std::cout << "x= " << x[0] << "\t y = " << x[1] << "\t Phi = " << y << std::endl;
			        ++iicount;
				phi += y;
			   }
			}
		}
	}
  }
  std::cout << std::endl << "Phi_init = " << sysParams.get_phi_init() << "\t Phi_S = " << phi/iicount << std::endl;
}

//===============================================================
// Main programm
//===============================================================
int main(int argc, char** argv)
{
 try{
  // Initialize Mpi
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  std::cout << "Hello World! This is iPBS." << std::endl;
  if(Dune::MPIHelper::isFake)
    std::cout<< "This is (at the moment) a sequential program!" << std::endl;
  else
  {
    if(helper.rank()==0)
    std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
  }

  // check arguments
  if (argc!=4)
  {
    if (helper.rank()==0)
    {
	std::cout << "usage: ./iPBS <meshfile> <refinement level> <MaxNewtonIterations>" << std::endl;
	return 1;
    }
  }

  // Read in comandline arguments
  Cmdparam cmdparam;
  cmdparam.GridName=argv[1];
  sscanf(argv[2],"%d",&cmdparam.RefinementLevel);
  sscanf(argv[3],"%d",&cmdparam.NewtonMaxIteration);
  std::cout << "Using " << cmdparam.RefinementLevel << " refinement levels." << std::endl;

  
  // <<<1>>> Setup the problem from mesh file

  // instanciate ug grid object
  GridType grid(400);		// heapSize: The size of UG's internal memory in megabytes for this grid. 

  // define vectors to store boundary and element mapping
  std::vector<int> boundaryIndexToEntity;
  std::vector<int> elementIndexToEntity;

  // read a gmsh file
  Dune::GmshReader<GridType> gmshreader;
  gmshreader.read(grid, cmdparam.GridName, boundaryIndexToEntity, elementIndexToEntity, true, false);

  // refine grid
  grid.globalRefine(cmdparam.RefinementLevel);

  // get a grid view
  const GV& gv = grid.leafView();

  // inner region, i.e. solve
  M m(gv, elementIndexToEntity);

  // boundary condition
  B b(gv, boundaryIndexToEntity);
  G g(gv, boundaryIndexToEntity);

  // boundary fluxes
  J j(gv, boundaryIndexToEntity);

  // <<<2>>> Make grid function space
  FEM fem;
  CON con;
  GFS gfs(gv,fem);
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc);                          // assemble constraints
  std::cout << "constrained dofs=" << cc.size() 
            << " of " << gfs.globalSize() << std::endl;

  // <<<3>>> Make FE function extending Dirichlet boundary conditions
  U u(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,u);                           // interpolate coefficient vector
  
  // <<<4>>> Make grid operator space
  LOP lop(m,b,j);
  GOS gos(gfs,cc,gfs,cc,lop);

  // <<<5a>>> Select a linear solver backend
  LS ls(5000,true);
  
  // <<<5b>>> Instantiate solver for nonlinear problem
  NEWTON newton(gos,u,ls); 
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setReduction(1e-10);
  newton.setMinLinearReduction(1e-4);
  newton.setMaxIterations(cmdparam.NewtonMaxIteration);
  newton.setLineSearchMaxIterations(10);

  // <<<5c>>> Instantiate Solver for linear problem
  SLP slp(gos,u,ls,1e-10); 

  // <<<6>>> Solve Problem
  solver(newton, slp);

  // Create Grid Function Space
  DGF udgf(gfs,u);

  // graphical output
  save(udgf, gv);

  calculate_phi(gv, udgf);
  
  // Try what we need in the iteration circle
  sysParams.set_phi_init(5.0);
  U u_new(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,u_new);
  LOP lop_new(m,b,j);
  GOS gos_new(gfs,cc,gfs,cc,lop_new);
  NEWTON newton_new(gos_new,u_new,ls); 
  newton_new.setReassembleThreshold(0.0);
  newton_new.setVerbosityLevel(2);
  newton_new.setReduction(1e-10);
  newton_new.setMinLinearReduction(1e-4);
  newton_new.setMaxIterations(cmdparam.NewtonMaxIteration);
  newton_new.setLineSearchMaxIterations(10);
  SLP slp_new(gos_new,u_new,ls,1e-10);
  solver(newton_new, slp_new);
  DGF udgf_new(gfs,u_new);
  calculate_phi(gv, udgf_new);

  // done
  return 0;
 }
 catch (Dune::Exception &e){
  std::cerr << "Dune reported error: " << e << std::endl;
 }
 catch (...){
  std::cerr << "Unknown exception thrown!" << std::endl;
 }
}  
