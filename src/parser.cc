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
  const double alpha_sor = 0.7;
  const int level = 0;
  const int verbose = 4;
  const double lambda = 1.0;
  const double bjerrum = 0.7;
  const double radius = 1.0;

  // set the symmetry of the system
  sysParams.set_symmetry(configuration.get<double>("mesh.symmetry"));
  sysParams.set_boxLength(configuration.get<double>("mesh.boxLength",0.0));
  
  // Parse other options
  sysParams.set_alpha(configuration.get<double>("solver.alpha_sor",alpha_sor));
  sysParams.set_refinement(configuration.get<int>("mesh.global_refinement_level",level));
  sysParams.set_refinementFraction(configuration.get<double>("mesh.adaptive_refinement_fraction",0));
  sysParams.set_refinementSteps(configuration.get<int>("mesh.adaptive_refinement_steps",1));
  sysParams.set_bjerrum(configuration.get<double>("system.bjerrum",bjerrum));
  sysParams.set_lambda(configuration.get<double>("system.lambda",lambda));
  sysParams.set_radius(configuration.get<double>("system.radius",radius));
  sysParams.set_tolerance(configuration.get<double>("solver.tolerance"));
  sysParams.set_verbose(configuration.get<int>("system.verbose",verbose));
  sysParams.set_salt(configuration.get<int>("system.salt"));;

  if (sysParams.get_symmetry() == 2)
    // TODO check the charge for 2d case!!!
  {
    double charge_density = configuration.get<double>("system.charge") / (4.0 * sysParams.pi
        * sysParams.get_radius() * sysParams.get_radius());
    sysParams.set_charge_density(charge_density);
  }
  else 
    sysParams.set_charge_density(configuration.get<double>("system.charge_density"));
}
