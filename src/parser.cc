/** \file
    
    \brief Read in IPBS configuration parameters from a config file 
    
    \todo Doc me!
*/

#include<string>
#include<dune/common/parametertree.hh>
#include<dune/common/parametertreeparser.hh>
#ifndef _SYSPARAMS_H
#define _SYSPARAMS_H
#include "sysparams.hh"
#endif

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
  static const double alpha_sor = 0.7;
  static const int level = 0;
  static const double lambda = 1.0;
  static const double bjerrum = 0.7;
  static const double radius = 1.0;

  // Parse other options
  sysParams.set_alpha(configuration.get<double>("solver.alpha_sor",alpha_sor));
  sysParams.set_refinement(configuration.get<int>("mesh.global_refinement_level",level));
  sysParams.set_bjerrum(configuration.get<double>("system.bjerrum",bjerrum));
  sysParams.set_lambda(configuration.get<double>("system.lambda",lambda));
  sysParams.set_radius(configuration.get<double>("system.radius",radius));
}