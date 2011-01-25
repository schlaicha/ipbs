// Input/Output
#include <dune/grid/io/file/gnuplot.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

// pdelab includes
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include<dune/pdelab/instationary/onestep.hh>   // Filenamehelper
#include<dune/pdelab/finiteelementmap/p1fem.hh>	// P1 in 1,2,3 dimensions
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/constraints.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/backend/istlvectorbackend.hh>
#include<dune/pdelab/backend/istlmatrixbackend.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

template<class GV, typename M, typename B, typename G, typename J>
void solver (const GV& gv, const M& m, const B& b, const G& g, const J& j, const Cmdparam &cmdparam)
{
  // <<<1>>> Choose domain and range field type
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
  const int dim = GV::dimension;

  // <<<2>>> Make grid function space
  typedef Dune::PDELab::P1LocalFiniteElementMap<Coord,Real,dim> FEM;
  FEM fem;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;     // constraints class
  //typedef Dune::PDELab::NonoverlappingConformingDirichletConstraints CON;     // constraints class
  CON con;
  typedef Dune::PDELab::ISTLVectorBackend<1> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);

  //con.compute_ghosts(gfs);	// con stores indices of ghost nodes

  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc);                          // assemble constraints
  std::cout << "constrained dofs=" << cc.size() 
            << " of " << gfs.globalSize() << std::endl;

  // <<<3>>> Make FE function extending Dirichlet boundary conditions
  typedef typename GFS::template VectorContainer<Real>::Type U;
  U u(gfs,0.0);
  Dune::PDELab::interpolate(g,gfs,u);                           // interpolate coefficient vector
  //Dune::PDELab::set_nonconstrained_dofs(cc,0.0,u);

  
  // <<<4>>> Make grid operator space
  typedef PBLocalOperator<M,B,J> LOP;
  LOP lop(m,b,j);
  typedef Dune::PDELab::ISTLBCRSMatrixBackend<1,1> MBE;
  typedef Dune::PDELab::GridOperatorSpace<GFS,GFS,LOP,CC,CC,MBE,true> GOS;
  GOS gos(gfs,cc,gfs,cc,lop);

   // <<<5a>>> Select a linear solver backend
   typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
   LS ls(5000,true);
   //typedef Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GFS> LS; // parallel backend
   //LS ls(gfs, 5000, 1);

  // <<<5b>>> Solve nonlinear problem
  Dune::PDELab::Newton<GOS,LS,U> newton(gos,u,ls);                        
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setReduction(1e-10);
  newton.setMinLinearReduction(1e-4);
  newton.setMaxIterations(cmdparam.NewtonMaxIteration);
  newton.setLineSearchMaxIterations(10);
  newton.apply();
  
  // <<<6>>> assemble and solve linear problem
  typedef Dune::PDELab::StationaryLinearProblemSolver<GOS,LS,U> SLP;
  SLP slp(gos,u,ls,1e-10);
  slp.apply();

  // <<<7>>> graphical output
  Dune::PDELab::FilenameHelper fn("output");
   {
    typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
    DGF udgf(gfs,u);
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTKOptions::conforming);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf,"solution"));
    vtkwriter.write(fn.getName(),Dune::VTKOptions::binaryappended);
  }

  // Gnuplot output
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
  DGF udgf(gfs,u);
  Dune::GnuplotWriter<GV> gnuplotwriter(gv);
  gnuplotwriter.addVertexData(u,"solution");
  gnuplotwriter.write("gnuplot.dat");
}
