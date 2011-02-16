// solve.hh - solve the PB equation

template <typename U>
void get_solution(U &u, const GV &gv, const GFS &gfs, const CC &cc, const M &m, const B &b, const J&j)
{
  int iterationCounter = 0;
  sysParams.init = true;
  while (sysParams.get_error() > 1E-2)
  //while(iterationCounter < 10)
  {
    std::cout << std::endl << "IN ITERATION " << iterationCounter << std::endl << std::endl;
    // Reset error for new iteration
    sysParams.reset_error();
    
    // construct discrete grid function for access to solution
    const U storeCoefficientVector = u;
    const DGF udgf(gfs, storeCoefficientVector);
  
    //const DGF udgf(gfs, u);
    
    // <<<4>>> Make grid operator space
    LOP lop(m,b,j, udgf, gv);
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
    
    
    std::stringstream out;
    out << "step_" << iterationCounter;
    std::string vtk_filename = out.str();
    DGF udgf_save(gfs,u);
    save(udgf_save, u, gv, vtk_filename);
    ++iterationCounter;
    std::cout << std::endl << "actual error is: " << sysParams.get_error() << std::endl << std::endl;
    sysParams.init = false;
  }
  
  
  
  
}
