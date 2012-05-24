#ifndef _IPBSOLVER_HH
#define _IPBSOLVER_HH

#include <dune/grid/common/scsgmapper.hh> // Single Geometry Single Codim Mapper
#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>

#include <gsl/gsl_sf_ellint.h>
#include "sysparams.hh"
#include "boundary.hh"
#include "e_field.hh"

#include <time.h>

extern SysParams sysParams;
extern std::vector<Boundary*> boundary;



template <class GV, class GFS>
class Ipbsolver

/** \brief Encapsulation of the algorithm determining the iterated boundary values

    Here goes some explanation on what is done :-)
    \todo { Doc Me ! }   
*/

{

  // Some typedef
  typedef typename GV::Grid::ctype ctype;
  static const int dim = GV::dimension;
  typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;

  public:
    Ipbsolver(const GV& gv_, const GFS& gfs_,
        const std::vector<int>& boundaryIndexToEntity_, const int intorder_=1,
        const bool use_guess=true) :
      gv(gv_), gfs(gfs_), boundaryIndexToEntity(boundaryIndexToEntity_),
      intorder(intorder_),  boundaryElemMapper(gv), communicator(gv.comm()), 
      my_offset(0), my_len(0), iterationCounter(0), fluxError(0)
     
    /*!
       \param gv the view on the leaf grid
       \param boundaryIndexToEntity physical property of boundary elements
    */
    {
      // Resize the containers specific to each physical surface
      physArea.resize( sysParams.get_npart() );
      physQTot.resize( sysParams.get_npart() );
      physQIpbs.resize( sysParams.get_npart() );
      physShift.resize( sysParams.get_npart() );

      init(); // Detect iterative elements
      communicateIpbsData(); 

      bContainer.resize(ipbsPositions.size(),0);
      inducedChargeDensity.resize(ipbsPositions.size(),0);
      E_ext.resize(ipbsPositions.size(),0);
      /// Prepare container for storing the potential at each intersection
      regulatedChargeDensity.resize(ipbsPositions.size(), 0.);

      if (use_guess) initial_guess();
    }

    int get_n()
    {
      return ipbsPositions.size();
    }

    bool converged()
    {
      iterationCounter++;
      std::cout << "in iteration " << iterationCounter << " the relative fluxError is " << fluxError
        << " relative error in induced charge density is " << icError << std::endl;
      if ( std::max(fluxError, icError) < sysParams.get_tolerance() ) {
        return true;
      }
      else
        return false;
    }
    
    /*!
     * \param _fluxError return the current maximum relative change in boundary condition calulation
     * \param _icError return the current maximum relative change in induced charge computation
     */

    bool converged(double& _fluxError, double& _icError, int& _iterations)
    {
      iterationCounter++;
      _fluxError = fluxError;
      _icError = icError;
      _iterations = iterationCounter;
      std::cout << "in iteration " << iterationCounter << " the relative fluxError is " << fluxError
        << " relative error in induced charge density is " << icError << std::endl;
      if ( std::max(fluxError, icError) < sysParams.get_tolerance() ) {
        return true;
      }
      else
        return false;
    }
    
    template <typename I>
    double get_flux(const I& i) const
    {
      int mappedIndex = indexLookupMap.find(boundaryElemMapper.map(*i.inside()))->second + my_offset;
      double y = bContainer[mappedIndex];
      return y;
    }

    // ------------------------------------------------------------------------
    /// Induced charge computation
    // ------------------------------------------------------------------------
    void updateIC()
    {
      icError = 0;  // reset the fluxError for next iteration step
      //std::cout << "in updateIC() my_offset = " << my_offset << " my_len = " << my_len << std::endl;
      ContainerType ic(inducedChargeDensity.size(), 0.);
      double eps_out = sysParams.get_epsilon();  
      unsigned int target = my_offset + my_len;
      for (unsigned int i = my_offset; i < target; i++) {
        double eps_in = boundary[ipbsType[i]]->get_epsilon();

        // include charge regulation
        double my_charge = boundary[ipbsType[i]]->get_charge_density() + regulatedChargeDensity[i];

        // calculate the induced charge in this surface element
        ic[i] = -(eps_in - eps_out) / (eps_in+eps_out)
                                * ( my_charge + eps_out/(2.*sysParams.pi*sysParams.get_bjerrum())
                                    * E_ext[i] );
      }
      communicator.barrier();
      communicator.sum(&ic[0], ic.size());

      // Do the SOR 
      for (unsigned int i = 0; i < inducedChargeDensity.size(); i++) {
        inducedChargeDensity[i] = sysParams.get_alpha_ic() * ic[i]
                          + ( 1 - sysParams.get_alpha_ic()) * inducedChargeDensity[i];
        double local_icError = fabs(2.0*(ic[i]-inducedChargeDensity[i])
                      /(ic[i]+inducedChargeDensity[i]));
        icError = std::max(icError, local_icError);
      }
      communicator.barrier();
      communicator.max(&icError, 1);
    }

    // ------------------------------------------------------------------------
    /// This method will do the update on the boundary conditions
    // ------------------------------------------------------------------------

    void updateBC(const U& u)
    {
      fluxError = 0;  // reset the fluxError for next iteration step
      
      /// Construct a discrete grid function space for access to solution
      typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
      typedef typename DGF::Traits::RangeType RT;
      typedef typename DGF::Traits::DomainFieldType DF;

      DGF udgf(gfs,u);

      typedef typename GV::template Codim<0>::template Partition
              <Dune::Interior_Partition>::Iterator LeafIterator;

      /// Store the new calculated values
      ContainerType fluxes(ipbsPositions.size(),0.);
      ContainerType physFluxShift( sysParams.get_npart(), 0 );

      for (size_t i =0; i<E_ext.size(); i++) {
          E_ext[i] = 0;
      }
      for (size_t i = 0; i < sysParams.get_npart(); i++){
          physQIpbs[i] = 0;
      }

      // Loop over all elements and calculate the volume integral contribution
      for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
               	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
      {

        // For each element on this processor calculate the contribution to volume integral part of the flux
        for (size_t i = 0; i < ipbsPositions.size(); i++)
        {

          // select quadrature rule
          Dune::GeometryType gt = it->geometry().type();
          const Dune::QuadratureRule<DF,dim>& 
            rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

          // loop over quadrature points
          for (typename Dune::QuadratureRule<DF,dim>::const_iterator 
                     q_it=rule.begin(); q_it!=rule.end(); ++q_it)
          {
            double E_ext_ions = 0.;
            // Get the position vector of the this element's center
            Dune::FieldVector<ctype,dim> r_prime = it->geometry().global(q_it->position());
            //Dune::FieldVector<ctype,dim> r_prime = it->geometry().center();
                      
            RT value;
            udgf.evaluate(*it,q_it->position(),value);
            //udgf.evaluate(*it,it->geometry().local(r_prime),value);
            
            double weight = q_it->weight() * it->geometry().integrationElement(q_it->position());
            //double weight = it->geometry().volume();

            Dune::FieldVector<ctype,dim> r (ipbsPositions[i]);
            Dune::FieldVector<ctype,dim> dist = r - r_prime;
            Dune::FieldVector<ctype, dim> unitNormal(ipbsNormals[i]);
            unitNormal *= -1.;

            Dune::FieldVector<ctype,dim> e_field(0.);

            e_field = E_field<Dune::FieldVector<ctype,dim> ,Dune::FieldVector<ctype,dim> > 
                        (r, r_prime, sysParams.get_symmetry());

            e_field *= weight;
            E_ext_ions = e_field * unitNormal;
            E_ext_ions *= 1./ (4.0*sysParams.pi) * sysParams.get_lambda2i();

            if ( sysParams.get_symmetry() > 0) {
                E_ext_ions *= r_prime[1];
            }

            // Now get the counterion distribution and calculate flux
            switch (sysParams.get_salt())
            {
                case 0:
                    E_ext_ions *= -std::sinh(value);
                    break;
                case 1:
                    E_ext_ions *= std::exp(value); // Counterions have opposite sign!
                    break;
            }
            E_ext[i] += E_ext_ions;
            if (sysParams.get_symmetry() == 0)
              physFluxShift[ ipbsType[i] ] -= E_ext_ions*ipbsVolumes[i];
            else
              physFluxShift[ ipbsType[i] ] -= E_ext_ions*ipbsVolumes[i]*2*sysParams.pi*r[1];
          }
        }
      }


      unsigned int target = my_offset + my_len;
      // For each element on this processor calculate the contribution to surface integral part of the flux
      for (unsigned int i = my_offset; i < target; i++)
      {
        Dune::FieldVector<ctype,dim> r (ipbsPositions[i]);
        Dune::FieldVector<ctype, dim> unitNormal(ipbsNormals[i]);
        unitNormal *= -1.;
        for (size_t j = 0; j < ipbsPositions.size(); j++)
        {
          double surfaceElem_flux = 0.;
          double lcd = boundary[ipbsType[j]]->get_charge_density() 
                          + inducedChargeDensity[j]
                          + regulatedChargeDensity[j]; /**< The local charge density 
                                                       on this particular surface element */
          if (i!=j)
          { 
            Dune::FieldVector<ctype,dim> r_prime (ipbsPositions[j]);
            Dune::FieldVector<ctype,dim> e_field(1.);

            e_field = E_field<Dune::FieldVector<ctype,dim> ,Dune::FieldVector<ctype,dim> > 
                        (r, r_prime, sysParams.get_symmetry());

            e_field *= ipbsVolumes[j];
      
            surfaceElem_flux = e_field * unitNormal;
            if ( sysParams.get_symmetry() > 0 ) {  // TODO This is better using a mirror switch and a cylinder switch
              surfaceElem_flux *= r_prime[1];
            }
            surfaceElem_flux *= lcd;
          } 
          else // if i==j - only contributes in case of mirroring
            if (sysParams.get_symmetry() == 2) {
              Dune::FieldVector<ctype,dim> r_prime2(ipbsPositions[i]);
              r_prime2[0] *= -1;
              
              Dune::FieldVector<ctype,dim> e_field(0.);
              e_field = E_field<Dune::FieldVector<ctype,dim> ,Dune::FieldVector<ctype,dim> > (r, r_prime2, 1);
              e_field *= ipbsVolumes[j];
              surfaceElem_flux = e_field * unitNormal;
              // cylindrical coordinates
              surfaceElem_flux *= r_prime2[1];
              surfaceElem_flux *= lcd;
          }
          E_ext[i] += surfaceElem_flux;
          if (sysParams.get_symmetry() == 0)
            physFluxShift[ ipbsType[i] ] -= surfaceElem_flux*ipbsVolumes[i];
          else
            physFluxShift[ ipbsType[i] ] -= surfaceElem_flux*ipbsVolumes[i]*2*sysParams.pi*r[1];
        } // end of j loop
      } // end of i loop

      
      // Collect results from all nodes
      communicator.barrier();
      communicator.sum(&E_ext[0], fluxes.size());
      
      for (unsigned int i = my_offset; i < target; i++) {
        
        double thisChargeDensity = boundary[ipbsType[i]]->get_charge_density() 
                                  + inducedChargeDensity[i] + regulatedChargeDensity[i];
      
        /* Equation 3.4.4 DA Schlaich */
        fluxes[i] = E_ext[i] + 2. * sysParams.pi*sysParams.get_bjerrum() * thisChargeDensity;

        double thisCharge = 0;
        if (sysParams.get_symmetry() > 0)
          thisCharge = 2 * sysParams.pi * ipbsPositions[i][1] * ipbsVolumes[i] * thisChargeDensity;
        else
          thisCharge = ipbsVolumes[i] * thisChargeDensity;

        physQIpbs[ ipbsType[i] ] += thisCharge;
        physFluxShift[ ipbsType[i] ] += 2. * sysParams.get_bjerrum()* sysParams.pi * thisCharge; 
      }  
        
      // Collect results from all nodes
      communicator.barrier();
      communicator.sum( &physQIpbs[0], physQIpbs.size() );
      communicator.sum( &physFluxShift[0], physFluxShift.size() );

      if (communicator.rank() == 0 && sysParams.get_verbose() > 0) {
        for (size_t i = 0; i < sysParams.get_npart(); i++) {
          if (boundary[i]->get_type() == 2)
            std::cout << "Charge on surface " << i << ": " << physQIpbs[i] << ", flux shifted by " 
              << physFluxShift[i] / physArea[ ipbsType[i] ] << std::endl;
        }
      }
      
      // Do the SOR 
      for (unsigned int i = my_offset; i < target; i++) {

        // do the shift
        //fluxes[i] += physFluxShift[ ipbsType[i] ] / physArea[ ipbsType[i] ];

        bContainer[i] = sysParams.get_alpha_ipbs() * fluxes[i]
                        + ( 1 - sysParams.get_alpha_ipbs()) * bContainer[i];
        double local_fluxError = fabs(2.0*(fluxes[i]-bContainer[i])
                     /(fluxes[i]+bContainer[i]));
        fluxError = std::max(fluxError, local_fluxError);
      }
      communicator.barrier();
      communicator.max(&fluxError, 1);
 
    }


    // ------------------------------------------------------------------------
    /// Charge regulation calculation
    // ------------------------------------------------------------------------
    void updateChargeRegulation(const U& u)
    { 
      /// Construct a discrete grid function space for access to solution
      typedef Dune::PDELab::DiscreteGridFunction<GFS,U> DGF;
      typedef typename DGF::Traits::RangeType RT;
      typedef typename DGF::Traits::DomainFieldType DF;

      DGF udgf(gfs,u);

      for ( size_t i = 0; i < ipbsElemPointers.size(); i++)
      {
        // Calculate regulated charge density
        RT value;
        ElemPointer it = ipbsElemPointers[i];
        Dune::FieldVector<ctype, dim> local = it->geometry().local( ipbsPositions[i] );
        udgf.evaluate(*it, local,value);
        double pK = boundary[ ipbsType[i] ]->get_pK();
        double pH = sysParams.get_pH();
        double exponential = std::exp ( std::log10(10)*( pH - pK ) - value );
        regulatedChargeDensity[i] = boundary[ ipbsType[i] ]->get_sigma_max() * exponential / ( 1. + exponential);
        //std::cout << "i: " << i << " pH " << pH << " pK " << pK << " Value: " << value << " exponential: " << exponential << " regulatedChargeDensity: " << regulatedChargeDensity[i] << std::endl;
      }
    }



  // ------------------------------------------------------------------------


  private:
    void initial_guess()
    {
      /** \brief Get an inital guess for the iterative boundaries */

      srand ( time(NULL) );
      unsigned int target = my_offset + my_len;
      // For each element on this processor calculate the contribution to surface integral part of the flux
      for (unsigned int i = my_offset; i < target; i++)
      {
        // initialize with constant surface charge density
        bContainer[i] = 4. * sysParams.get_bjerrum() * sysParams.pi * boundary[ipbsType[i]]->get_charge_density();
        //bContainer[i] = 4. * sysParams.get_bjerrum() * sysParams.pi * boundary[ipbsType[i]]->get_charge_density() * ( rand()/RAND_MAX*.4 + .8);
      }
    }
 
    // ------------------------------------------------------------------------
    /// This method executes everything related to initialization
    // ------------------------------------------------------------------------

    void init()
    {
       /** \brief Traverse the grid and determine which boundaries shall be iterated
        *
        *   Assign a local list of iterative boundary elements and communicate this to all processors.
        */
      typedef typename GV::template Codim<0>::template Partition
              <Dune::Interior_Partition>::Iterator LeafIterator;
      typedef typename GV::IntersectionIterator IntersectionIterator;
    
      int counter = 0;
      int elemcounter = 0;
      /// loop over elements on this processor and get information on the iterated boundaries
      for (LeafIterator it = gv.template begin<0,Dune::Interior_Partition>();
               	it!=gv.template end<0,Dune::Interior_Partition>(); ++it)
      {
        elemcounter++;
        if (it->hasBoundaryIntersections() == true)
          for (IntersectionIterator ii = it->ileafbegin(); ii!=it->ileafend(); ++ii) {
            if(ii->boundary()==true && 
                boundary[ boundaryIndexToEntity[ii->boundarySegmentIndex()] ]->get_type() == 2)
              {
                int physType = boundaryIndexToEntity[ii->boundarySegmentIndex()];
                double volume = ii->geometry().volume();
                Dune::FieldVector<ctype, dim> center = ii->geometry().center();

                ipbsPositions.push_back(center);
                ipbsNormals.push_back(ii->centerUnitOuterNormal());
                indexLookupMap.insert(std::pair<int, int>(boundaryElemMapper.map(*it),counter));
                ipbsType.push_back( physType );
                ipbsVolumes.push_back(volume);
                ipbsElemPointers.push_back(*it);   // Point to the IPBS elements local on each node
                counter++;
                // charge neutralisation stuff
                if (sysParams.get_symmetry() > 0)
                  volume *= 2. * sysParams.pi * center[1];
                physArea[physType] += volume;
                physQTot[physType] += boundary[physType]->get_charge_density() * volume;
              }
          }
      }

      communicator.sum( &physArea[0], physArea.size() );
      communicator.sum( &physQTot[0], physQTot.size() );
      if (communicator.rank() == 0 && sysParams.get_verbose() > 0) {
        std::cout << "Init detected the following physical iterative surfaces:" << std::endl;
        for (size_t i = 0; i < sysParams.get_npart(); i++) {
          if (boundary[i]->get_type() == 2)
            std::cout << i << "\t area: " << physArea[i] << "\tQ: " << physQTot[i] << std::endl;
        }
      }
    }

    int communicateIpbsData() /** The containers now contain information on the iterative boundary elements on each processor.
                                 Collect the information into one container that is distributed on all processors.*/
    {
#if HAVE_MPI
      std::vector<unsigned int> length_on_processor;
      length_on_processor.resize(communicator.size(),0);
      // Collect on rank 0  - we seem to have to use the standard MPI calls, this is why we NEED to check isFake() before
      MPI_Status status;
      if (communicator.rank() == 0)
      {
        length_on_processor[0] = ipbsVolumes.size();
        for (int i = 1; i < communicator.size(); i++)
        {
          int temp;
          MPI_Recv(&temp,1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
          length_on_processor[i]=temp;
        }
      }
      else
      {
        int temp = ipbsVolumes.size();
        MPI_Send(&temp,1,MPI_INT,0,communicator.rank(),MPI_COMM_WORLD); // Ensures we recover the correct order
      }
      int countBoundElems = 0;
      if (communicator.rank() == 0)
      {
         for (std::vector<unsigned int>::iterator i=length_on_processor.begin(); i != length_on_processor.end(); ++i)
         {
           // std::cout << std::endl << *i;
           countBoundElems += *i;
         }
         std::cout << "Detected " << countBoundElems << " boundary elements." << std::endl;
      }
      communicator.broadcast(&countBoundElems, 1, 0);
      communicator.broadcast(&length_on_processor[0], communicator.size(), 0);
      // Calculate the offset and length of data for each property-vector on each node
      for (int i = 0; i < communicator.rank(); i++) {
        my_offset += length_on_processor[i];
      }
      my_len = length_on_processor[communicator.rank()]; 
      // Now we send/receive the boundary data
      ipbsVolumes.resize(countBoundElems,0.);
      ipbsType.resize(countBoundElems,0.);
      ipbsNormals.resize(countBoundElems,Dune::FieldVector<ctype,dim>(0.));
      ipbsPositions.resize(countBoundElems,Dune::FieldVector<ctype,dim>(0.));
  
      // for communicating position-vectors we temporally store into arrays
      int indexcounter; // Count how many position-vectors we already got
      double* all_positions = (double*) malloc(dim*countBoundElems*sizeof(double)); // allocate on all processors
      double* all_normals = (double*) malloc(dim*countBoundElems*sizeof(double)); // allocate on all processors
      if( communicator.rank() !=0) // other nodes send their positions to master node
      {
        double* my_positions = (double*) malloc(dim*length_on_processor[communicator.rank()]*sizeof(double));
        double* my_normals = (double*) malloc(dim*length_on_processor[communicator.rank()]*sizeof(double));
        for (unsigned int i = 0; i<length_on_processor[communicator.rank()]; i++){
          for (int j = 0; j<dim; j++)
          {
            my_positions[i*dim+j] = ipbsPositions.at(i).vec_access(j);
            my_normals[i*dim+j] = ipbsNormals.at(i).vec_access(j);
          }
        }
        MPI_Send(my_positions,dim*length_on_processor[communicator.rank()],MPI_DOUBLE,0,0,MPI_COMM_WORLD); // pos sent on slot 0
        MPI_Send(my_normals,dim*length_on_processor[communicator.rank()],MPI_DOUBLE,0,1,MPI_COMM_WORLD); // normals sent on slot 1
      }
      else
      {   // Write positions of master node
        for (unsigned int i = 0; i<length_on_processor[communicator.rank()]; i++){
          for (int j = 0; j<dim; j++)
          {
            all_positions[i*dim+j] = ipbsPositions[i][j];
            all_normals[i*dim+j] = ipbsNormals[i][j];
          }
        }
        // get the right offset
        indexcounter = dim*length_on_processor[communicator.rank()];
        // get positions from other nodes
        for (int i = 1; i < communicator.size(); i++) {
          MPI_Recv(&all_positions[indexcounter],dim*length_on_processor[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
          MPI_Recv(&all_normals[indexcounter],dim*length_on_processor[i],MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);
          indexcounter += dim*length_on_processor[i];
        }
      }
      communicator.broadcast(all_positions,countBoundElems*dim,0); // communicate array
      communicator.broadcast(all_normals,countBoundElems*dim,0); // communicate array
      for(int i = 0; i < countBoundElems; i++) {
        for(int j=0; j<dim; j++)
        {
          ipbsNormals[i][j] = all_normals[i*dim+j];
          ipbsPositions[i][j] = all_positions[i*dim+j];
        }
      }

      // Commicate everything else
      /// determine the offset in the boundary proerty vectors on each node
      int offset = 0;
      for (int i = 1; i < communicator.rank(); i++) {
        offset += length_on_processor[i];
      }
      
      if( communicator.rank() !=0) // other nodes send their positions to master node
      {
        MPI_Send(&ipbsType[0],length_on_processor[communicator.rank()],MPI_INT,0,1,MPI_COMM_WORLD);
        MPI_Send(&ipbsVolumes[0],length_on_processor[communicator.rank()],MPI_DOUBLE,0,2,MPI_COMM_WORLD);
      }
      else
      {
        int tmpcounter = length_on_processor[0];
        // get positions from other nodes
        for (int i = 1; i < communicator.size(); i++) {
          MPI_Recv(&ipbsType[tmpcounter],length_on_processor[i],MPI_INT,i,1,MPI_COMM_WORLD,&status);
          MPI_Recv(&ipbsVolumes[tmpcounter],length_on_processor[i],MPI_DOUBLE,i,2,MPI_COMM_WORLD,&status);
          tmpcounter += length_on_processor[i];
        }
      }

      //if (communicator.rank() == 0)
      //{  for(unsigned int i=0;i<ipbsType.size();i++)
      //    std::cout << ipbsType[i] << " " << ipbsVolumes[i] << " " << ipbsPositions[i] << std::endl;
      //}
      communicator.broadcast(&ipbsType[0], countBoundElems, 0);
      communicator.broadcast(&ipbsVolumes[0], countBoundElems, 0);
      communicator.barrier();
      std::cout << "Myrank: " << communicator.rank() << " Mylen: " << my_len << std::endl;
      return 0;
#else
      my_len = ipbsType.size();
      my_offset = 0;
      return -1;
#endif
    }
 
    const GV& gv;
    /// The grid function space
    const GFS& gfs;
    const std::vector<int>& boundaryIndexToEntity;
    // provide a mapper for getting indices of iterated boundary elements
    typedef Dune::SingleCodimSingleGeomTypeMapper<GV, 0> BoundaryElemMapper;
    BoundaryElemMapper boundaryElemMapper;
    /// The communicator decides weither to use MPI or fake
    typedef typename GV::Traits::CollectiveCommunication CollectiveCommunication;
    const CollectiveCommunication & communicator;

    /// Store the center of boundary intersections of iterative type
    std::vector<Dune::FieldVector<ctype,dim> > ipbsPositions;
    /// Store the normal vector at the center of boundary intersections of iterative type
    std::vector<Dune::FieldVector<ctype,dim> > ipbsNormals;
    /// Store the volume of boundary intersections of iterative type
    std::vector<Real> ipbsVolumes;
    /// Provide a vector storing the type of the iterative boundary @todo Maybe it's more reasonable to store its surface charge density
    std::vector<int> ipbsType;
    /// Element pointers to local IPBS elements
    typedef typename GV::template Codim<0>::EntityPointer ElemPointer;
    std::vector<ElemPointer> ipbsElemPointers;

    /// For fast access to the precomputed fluxes we use binary tree  - only local index is needed :-)
    typedef std::map<int, int> IndexLookupMap;
    IndexLookupMap indexLookupMap;

    typedef std::vector<double> ContainerType;
    ContainerType bContainer;  /// Store the flux on the surface elements
    ContainerType inducedChargeDensity;  // Store the induced charge density
    ContainerType E_ext;   /**< Store the electric field on boundary elements, 
                                 which is calculated during updateBC(),
                                 \f[ \vec{E}(\vec{r}) \propto \int_V \sinh(\Phi(\vec{r}))
                                 \textnormal{d}\vec{r} \f] */
    ContainerType regulatedChargeDensity;  /**< Store the regulatedChargeDensity
                                                CAVE: Only local on each node! */
    
    /// Offset and length of data stream on each node
    unsigned int my_offset, my_len;
    unsigned int iterationCounter;
    double fluxError,  icError;

    const int intorder; 

    /// Charge neutrallity constraint
    std::vector<double> physArea, physQTot, physQIpbs, physShift;
};


#endif  /* _IPBSOLVER_HH */
