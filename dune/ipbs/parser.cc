
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
 
  /// default values
  const double alpha_sor = 0.7;
  const int verbose = 1;
  const double lambda = 1.0;
  const double bjerrum = 0.7;


  /** \brief A parser to read in the system setup
   *
   *  \param filename   Configuration to read in
   */
  
  // Mesh file must be specified
  sysParams.set_meshfile(configuration.get<std::string>("mesh.filename"));
  /** set the symmetry of the system
   * PARAMS */
  sysParams.set_symmetry(configuration.get<double>("mesh.symmetry"));
  

  // Parse other options
  sysParams.set_maxiter(configuration.get<int>("solver.maxiter",100));
  sysParams.set_alpha_ipbs(configuration.get<double>("solver.alpha_ipbs",alpha_sor));
  sysParams.set_alpha_ic(configuration.get<double>("solver.alpha_ic",alpha_sor));
  sysParams.set_bjerrum(configuration.get<double>("system.bjerrum",bjerrum));
  sysParams.set_lambda(configuration.get<double>("system.lambda",lambda));
  sysParams.set_tolerance(configuration.get<double>("solver.tolerance"));
  sysParams.set_verbose(configuration.get<int>("system.verbose",verbose));
  sysParams.set_salt(configuration.get<int>("system.salt"));
  sysParams.set_pH(configuration.get<double>("system.pH", 7.));
  double epsilonOut = configuration.get<double>("system.epsilon");
  
  sysParams.set_integration_l(configuration.get<double>("solver.l", 0.15*sysParams.get_lambda()));
  sysParams.set_integration_d(configuration.get<double>("solver.d", 0.075*sysParams.get_lambda()));
  sysParams.set_integration_maxintorder(configuration.get<double>("solver.maxintorder", 10));

  // Output
  sysParams.set_outStep(configuration.get<int>("output.steps",0));

  // Create particles
  size_t n_particle = configuration.get<size_t>("system.NPart");
  sysParams.set_npart(n_particle);
  for (size_t i = 0; i < n_particle; i++)
    boundary.push_back(new Boundary());

  for(size_t i = 0; i < n_particle; i++)
  {
    std::string p_name = "boundary_";
    std::ostringstream s;
    s << p_name << i; 
                        
    p_name = s.str();
    boundary[i]->set_charge_density(configuration.get<double>(p_name+".charge_density",0));
    double epsilonIn = configuration.get<double>(p_name+".epsilon",1); 
    boundary[i]->set_epsilons(epsilonIn, epsilonOut);
    boundary[i]->set_type(configuration.get<int>(p_name+".type",0));
    boundary[i]->set_potential(configuration.get<double>(p_name+".potential",0));
    boundary[i]->set_sigma_max(configuration.get<double>(p_name+".sigma_max",0));
    boundary[i]->set_Y(configuration.get<double>(p_name+".Y",0));
    boundary[i]->set_pK(configuration.get<double>(p_name+".pK",0));
    boundary[i]->set_ifShift(configuration.get<bool>(p_name+".shifted",true));
  }

}

#endif  // _PARSER_H
