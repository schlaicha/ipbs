
#ifndef _PARSER_H
#define _PARSER_H

/** \file
    
    \brief Read in IPBS configuration parameters from a config file 
    
    \todo Doc me!
*/

#include<string>
#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>
#include "sysparams.hh"
#include "boundary.hh"

extern SysParams sysParams;
extern std::vector<Boundary*> boundary;

void parser(std::string config_file)
{
  Dune::ParameterTree configuration;
  Dune::ParameterTreeParser parser;
  
  try{
      parser.readINITree( config_file, configuration );
  }
  catch(...){
      std::cerr << "Could not read config file \"" 
                << config_file << "\"!" << std::endl;
      exit(1);
  }

  // Mesh file must be specified
  sysParams.set_meshfile(configuration.get<std::string>("mesh.filename"));

  // default values
  const double alpha_sor = 0.7;
  const int level = 0;
  const int verbose = 4;
  const double lambda = 1.0;
  const double bjerrum = 0.7;
  const double newton_tolerance = 1e-10;

  // set the symmetry of the system
  sysParams.set_symmetry(configuration.get<double>("mesh.symmetry"));
  sysParams.set_boxLength(configuration.get<double>("mesh.boxLength",0.0));
  
  // Parse other options
  sysParams.set_maxiter(configuration.get<int>("solver.maxiter",100));
  sysParams.set_alpha_ipbs(configuration.get<double>("solver.alpha_ipbs",alpha_sor));
  sysParams.set_alpha_ic(configuration.get<double>("solver.alpha_ic",alpha_sor));
  sysParams.set_newton_tolerance(configuration.get<double>("solver.newton_tolerance",newton_tolerance));
  sysParams.set_refinement(configuration.get<int>("mesh.global_refinement_level",level));
  sysParams.set_refinementFraction(configuration.get<double>("mesh.adaptive_refinement_fraction",0));
  sysParams.set_refinementSteps(configuration.get<int>("mesh.adaptive_refinement_steps",1));
  sysParams.set_bjerrum(configuration.get<double>("system.bjerrum",bjerrum));
  sysParams.set_lambda(configuration.get<double>("system.lambda",lambda));
  sysParams.set_tolerance(configuration.get<double>("solver.tolerance"));
  sysParams.set_verbose(configuration.get<int>("system.verbose",verbose));
  sysParams.set_salt(configuration.get<int>("system.salt"));
  double epsilonOut = configuration.get<double>("system.epsilon");

  // Create particles
  int n_particle = configuration.get<int>("system.NPart");
  sysParams.set_npart(n_particle);
  for (int i = 0; i < n_particle; i++)
    boundary.push_back(new Boundary());

  for(int i = 0; i < n_particle; i++)
  {
    std::string p_name = "boundary_";
    std::ostringstream s;
    s << p_name << i+2; // we don't need to set 0 and 1
                        // remember that vector starts with 0, so access the boundaries with -2
                        // TODO find a clever way!
    p_name = s.str();
    boundary[i]->set_charge_density(configuration.get<double>(p_name+".charge_density"));
    double epsilonIn = configuration.get<double>(p_name+".epsilon"); 
    boundary[i]->set_epsilons(epsilonIn, epsilonOut);
    boundary[i]->set_isPlane(configuration.get<bool>(p_name+".isPlane", false));
  }

}

#endif  // _PARSER_H
