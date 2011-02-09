// solve.hh - solve the PB equation

template <typename U, typename G>
void get_solution(U &u, const GV &gv, const GFS &gfs, const M &m, const B &b, const G &g, const J &j)
{
  // <<<3>>> Make FE function extending Dirichlet boundary conditions
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc); 
  //std::cout << "constrained dofs=" << cc.size() 
  //          << " of " << gfs.globalSize() << std::endl;

  // interpolate coefficient vector
  Dune::PDELab::interpolate(g,gfs,u);
  
      
  // construct discrete grid function for access to solution
  const DGF udgf(gfs, u);

  // <<<4>>> Make grid operator space
  LOP lop(m,b,j, udgf);
  GOS gos(gfs,cc,gfs,cc,lop);

  // <<<5a>>> Select a linear solver backend
  LS ls(5000,true);

  // <<<5b>>> Instantiate solver for nonlinear problem
  NEWTON newton(gos,u,ls); 
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(1);
  newton.setReduction(1e-10);
  newton.setMinLinearReduction(1e-4);
  newton.setMaxIterations(20);
  newton.setLineSearchMaxIterations(10);

  // <<<5c>>> Instantiate Solver for linear problem
  SLP slp(gos,u,ls,1e-10); 

  // <<<6>>> Solve Problem
  //solver(newton, slp);
  newton.apply();
  slp.apply();
}
