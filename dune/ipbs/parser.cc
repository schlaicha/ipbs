
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
  const double bjerrum = 0.7;
  const double lambda = 1.0;
  const double pH = 7.;
  const int verbose = 1;

  /** \brief A parser to read in the system setup
   *
   *  \param filename   Configuration to read in
   */
  
  // ===================================================================
  // MESH section
  // ===================================================================

  // Mesh file must be specified
  sysParams.set_meshfile(configuration.get<std::string>("mesh.filename"));
  // set the symmetry of the system
  sysParams.set_symmetry(configuration.get<double>("mesh.symmetry"));
  
  // ===================================================================
  // SYSTEM section
  // ===================================================================

  int n_particle = configuration.get<int>("system.NPart");
  double epsilonOut = configuration.get<double>("system.epsilon");
  sysParams.set_npart(n_particle);
  sysParams.set_bjerrum(configuration.get<double>("system.bjerrum",bjerrum));
  sysParams.set_lambda(configuration.get<double>("system.lambda",lambda));
  sysParams.set_pH(configuration.get<double>("system.pH", pH));
  sysParams.set_salt(configuration.get<int>("system.salt"));
  sysParams.set_verbose(configuration.get<int>("system.verbose",verbose));
  
  // TODO Remove or use!
  //sysParams.set_maxiter(configuration.get<int>("solver.maxiter",100));

  // ===================================================================
  // SOLVER section
  // ===================================================================
  
  sysParams.set_alpha_ipbs(configuration.get<double>("solver.alpha_ipbs",alpha_sor));
  sysParams.set_alpha_ic(configuration.get<double>("solver.alpha_ic",alpha_sor));
  sysParams.set_tolerance(configuration.get<double>("solver.tolerance"));

  // TODO Remove or use!
  // sysParams.set_newton_tolerance(configuration.get<double>("solver.newton_tolerance",newton_tolerance));

  // ===================================================================
  // BOUNDARY[...] section
  // ===================================================================

  for(int i = 0; i < n_particle; i++)
  {
    boundary.push_back(new Boundary());
    
    std::string p_name = "boundary_";
    std::ostringstream s;
    s << p_name << i; 
    p_name = s.str();

    double epsilonIn = configuration.get<double>(p_name+".epsilon",1); 
    boundary[i]->set_charge_density(configuration.get<double>(p_name+".charge_density",0));
    boundary[i]->set_epsilons(epsilonIn, epsilonOut);
    boundary[i]->set_type(configuration.get<int>(p_name+".type",0));
    boundary[i]->set_pK(configuration.get<double>(p_name+".pK",0));
    boundary[i]->set_potential(configuration.get<double>(p_name+".potential",0));
    boundary[i]->set_sigma_max(configuration.get<double>(p_name+".sigma_max",0));
    boundary[i]->set_Y(configuration.get<double>(p_name+".Y",0));
  }
}

#endif  // _PARSER_H
