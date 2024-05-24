module setup_flow
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines generate initial conditions, either from restart files or
  !! analytic forms.
  use kind_parameters
  use common_parameter
  use common_vars
#ifdef mp
  use mpi_transfers
#endif    
  use turbulence
  implicit none
  
  integer(ikind), parameter :: n_modes=8
  real(rkind),dimension(n_modes) :: mode_amp,mode_phase
  real(rkind), parameter :: ptbn_size = zero!half   !! perturbation size in multiples of fl_thck
  
contains
!! ------------------------------------------------------------------------------------------------
  subroutine initial_solution
     !! Allocates arrays for primary and secondary properties, and calls routines to populate arrays
     !! for initial conditions.
     use mirror_boundaries
     use derivatives
     use thermodynamics
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro     
     
     !! Allocate arrays for properties - primary
     allocate(rou(np),rov(np),row(np),ro(np),roE(np),divvel(np))
     allocate(Yspec(np,nspec))
     rou=zero;rov=zero;row=zero;ro=one;roE=one;Yspec=one;divvel=zero

     !! Secondary properties
     allocate(T(np));T=T_ref
     allocate(p(np));p=zero
     allocate(u(np),v(np),w(np));u=zero;v=zero;w=zero
     allocate(hrr(npfb));hrr = zero !! Array for heat release rate  
     
     !! Transport properties
     allocate(visc(np));visc = visc_ref
#ifndef isoT
     allocate(lambda_th(np))
#endif
     allocate(roMdiff(np,nspec))
     allocate(cp(np),Rgas_mix(np))
     
     !! Allocate the boundary temperatures
     if(nb.ne.0) then
        allocate(T_bound(nb));T_bound = T_ref  
        allocate(u_inflow_local(nb));u_inflow_local = u_char   
        allocate(dudt_inflow_local(nb));dudt_inflow_local=zero            
     end if

     !! Evaluate mass fractions of reactants and products assuming complete combustion
     call initialise_composition     

     !! =======================================================================
     !! Choose initial conditions
#ifndef restart     

     !! Initialise the perturbation
     call initialise_perturbation            

     !! Make the base flow
     call make_baseflow   
     
     !! FLOW TYPE CHOICE =====================================       
     if(flag_flow_type.eq.1) then        
        call make_1d_flame
     else if(flag_flow_type.eq.2) then    
        call make_gaussian_hotspot   
     else if(flag_flow_type.eq.3) then
        call load_flame_file
     else if(flag_flow_type.eq.4) then
        !! A messy routine to play with for other initial conditions
        call hardcode_initial_conditions     
     else
        !! NOT SUPPORTED
        write(6,*) "Requested flow type not supported, aborting."
        stop
     end if

     !! Density stratification (for, e.g. Rayleigh-Taylor flow)
     if(flag_strat.eq.1) then
        call add_density_stratification
     end if     

     !! END FLOW TYPE CHOICE =================================

#else    
     !! RESTART OPTION. 
     call load_restart_file

     !! Un-comment this to (re-)ignite a restarting simulation
!     call make_gaussian_hotspot  

#endif

     !! Add some turbulence to the velocity field
     if(flag_turbulent.eq.1) then
        call make_turbulent_velocity_field
     end if        
     !! =======================================================================
     
     !! Convert from velocity to momentum and Y to roY
     !$omp parallel do private(ispec)
     do i=1,np
        rou(i) = ro(i)*u(i)
        rov(i) = ro(i)*v(i)
        row(i) = ro(i)*w(i)                
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec(i,ispec)*ro(i)
        end do
     end do
     !$omp end parallel do   
              
     !! Set energy from ro,u,Y,T
     call initialise_energy     
     
     !! Mirrors and halos                        
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchanges_all
#endif          

     !! Obtain velocity from momentum for mirrors
     !$omp parallel do
     do i=npfb+1,np
        u(i)=rou(i)/ro(i)
        v(i)=rov(i)/ro(i)
        w(i)=row(i)/ro(i)                
     end do
     !$omp end parallel do
     
     !! Set the initial velocity divergence
#ifdef dim3        
     call calc_divergence(u,v,w,divvel(1:npfb))
#else
     call calc_divergence(u,v,divvel(1:npfb))     
#endif              
     
     !! Mirrors and halos for divvel                   
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchange_divvel
#endif         
         
     !! Initialise the variable which holds inflow velocity locally
     if(nb.ne.0)then
        !$omp parallel do
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.1) then !! Inflow node
              u_inflow_local(j) = u(i)
           endif
        end do
        !$omp end parallel do        
     end if
     
     !! Initialise the variable holding the inflow species     
     if(nb.ne.0) then
        allocate(Yspec_inflow(nspec))
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.1) then !! Inflow node
              do ispec=1,nspec
                 Yspec_inflow(ispec) = Yspec(i,ispec)/ro(i)
              end do
           end if
        end do
     endif
           
     !! SOME ADDITIONAL INITIALISATION STUFF ------------------------------------------------------
     !! Profiling - re-zero time accumulators
     segment_time_local = zero
     cputimecheck = zero         
     
#ifndef restart
     !! Pre-set the time-stepping (necessary for PID controlled stepping)
     dt = 1.0d-10            
     !! Initialise PID controller variables
     emax_np1=pid_tol;emax_n=pid_tol;emax_nm1=pid_tol
#endif     
     
     !! Initialise time-stepping error normalisation based on expected magnitudes
     ero_norm = one/one !! Unity density
     erou_norm = one/(one*u_char)  !! Characteristic velocity
     eroE_norm = one/p_ref         !! Energy of the order of pressure
     eroY_norm = one/(one*one)     !! Unity mass fraction
     

     !! Initialise PID controller variables for <|u|>
#ifndef restart
     eflow_nm1 = one
     sum_eflow = zero   
     driving_force(:) = zero !! Will only be changed if using PID controller
#endif     
                   
     return
  end subroutine initial_solution   
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
!! N.B. In the routines below here, we are loading or generating initial conditions on the
!! primitive variables (ro,u,v,w,T,Y), and hence Yspec holds Y. Everywhere else in the code, 
!! Yspec holds roY.
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_composition
     !! This subroutine evaluates the reactant and product mass fractions for the given
     !! stoichiometry.
     real(rkind) :: Yin_H2,Yin_O2,Yin_N2,Yout_H2O,Yout_O2
     real(rkind) :: o2n2_ratio,h2o2_stoichiometric,h2o2_ratio
     real(rkind) :: Yin_CH4,Yout_CO2
     real(rkind) :: ch4o2_ratio,ch4o2_stoichiometric,co2h2o_ratio     
  
     !! Allocate two arrays
     allocate(Yspec_reactants(nspec),Yspec_products(nspec))
  
     !! Single-phase
     if(nspec.eq.1) then
        Yspec_reactants(1)=one
        Yspec_products(1)=one
     end if
  
     !! Single-step chemistry
     if(nspec.eq.2) then 
        Yspec_reactants(1) = one
        Yspec_reactants(2) = zero
        Yspec_products(1) = zero
        Yspec_products(2) = one
     end if
    
     !! 9-species H2 chemistry
     !! Assuming lean conditions, so some O2 remains in products
     if(nspec.eq.9) then
        o2n2_ratio = one*molar_mass(2)/(3.76d0*molar_mass(9))
        h2o2_stoichiometric = two*molar_mass(1)/(one*molar_mass(2))
        h2o2_ratio = h2o2_stoichiometric*phi_in

        Yin_O2 = one/(one + h2o2_ratio + one/o2n2_ratio)
        Yin_H2 = Yin_O2*h2o2_ratio
        Yin_N2 = one - Yin_H2 - Yin_O2            
        Yout_O2 = Yin_O2*(one-phi_in)
        Yout_H2O = one - Yin_N2 - Yout_O2   

        Yspec_reactants(1) = Yin_H2
        Yspec_reactants(2) = Yin_O2
        Yspec_reactants(3:8) = zero
        Yspec_reactants(9) = Yin_N2
        
        Yspec_products(1) = zero
        Yspec_products(2) = Yout_O2
        Yspec_products(3) = Yout_H2O
        Yspec_products(4:8) = zero
        Yspec_products(9) = Yin_N2      
     end if
     
     !! 3 species H2 inert
     !! Assuming lean conditions, so some O2 remains in products
     if(nspec.eq.3) then
        o2n2_ratio = one*molar_mass(2)/(3.76d0*molar_mass(3))
        h2o2_stoichiometric = two*molar_mass(1)/(one*molar_mass(2))
        h2o2_ratio = h2o2_stoichiometric*phi_in

        Yin_O2 = one/(one + h2o2_ratio + one/o2n2_ratio)
        Yin_H2 = Yin_O2*h2o2_ratio
        Yin_N2 = one - Yin_H2 - Yin_O2            

        Yspec_reactants(1) = Yin_H2
        Yspec_reactants(2) = Yin_O2
        Yspec_reactants(3) = Yin_N2
        
        Yspec_products(1) = Yin_H2
        Yspec_products(2) = Yin_O2
        Yspec_products(3) = Yin_N2      
     end if     
     
     !! 15 species CH4 chemistry
     if(nspec.eq.16) then
        o2n2_ratio = one*molar_mass(2)/(3.76d0*molar_mass(16))
        ch4o2_stoichiometric = one*molar_mass(1)/(two*molar_mass(2))
        ch4o2_ratio = ch4o2_stoichiometric*phi_in
          
        Yin_O2 = one/(one + ch4o2_ratio + one/o2n2_ratio)
        Yin_CH4 = Yin_O2*ch4o2_ratio
        Yin_N2 = one - Yin_CH4 - Yin_O2
   
        co2h2o_ratio = one*molar_mass(3)/(two*molar_mass(4))  !! Assume complete combustion??
        Yout_CO2 = (one - Yin_N2)/(one + one/co2h2o_ratio)
        Yout_H2O = Yout_CO2/co2h2o_ratio        
     
        Yspec_reactants(1) = Yin_CH4
        Yspec_reactants(2) = Yin_O2
        Yspec_reactants(3:15) = zero
        Yspec_reactants(16) = Yin_N2
        
        Yspec_products(1:2) = zero
        Yspec_products(3) = Yout_CO2
        Yspec_products(4) = Yout_H2O
        Yspec_products(5:15) = zero
        Yspec_products(16) = Yin_N2   
     end if
  
     return
  end subroutine initialise_composition
!! ------------------------------------------------------------------------------------------------ 
  subroutine make_1d_flame
     integer(ikind) :: i,ispec,j
     real(rkind) :: fl_thck_nd,ptbn
     real(rkind) :: c,Rmix_local,x,y,z,ro_inflow

     !! Scale thickness because position vectors are scaled...
     fl_thck_nd = fl_thck/L_char 

     !! Inflow mixture gas constant
     Rmix_local = zero
     do ispec=1,nspec
        Rmix_local = Rmix_local + Yspec_reactants(ispec)*one_over_molar_mass(ispec)
     end do
     Rmix_local = Rmix_local*Rgas_universal     

     !! Inflow density
     ro_inflow = p_ref/(Rmix_local*T_ref)
    
     !! Loop over all nodes and impose values
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec,ptbn)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Make a perturbation
        call make_perturbation(y/(ymax-ymin),ptbn)
        x = x+ptbn        
                
        !! Error function based progress variable
        c = half*(one + erf((x-fl_pos_x)/fl_thck_nd))
        
        !! Temperature profile
        T(i) = T(i) + (T_hot - T(i))*c
        
        !! Composition
        do ispec = 1,nspec
           Yspec(i,ispec) = (one-c)*Yspec_reactants(ispec) + c*Yspec_products(ispec)
        end do
        
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal
        
        !! Density
        ro(i) = p_ref/(Rmix_local*T(i))
        
        !! Adjust the velocity
        u(i) = u(i)*ro_inflow/ro(i)
        v(i) = zero
        w(i) = zero
                        
     end do
     !$omp end parallel do
     
     !! Special apply values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              if(flag_wall_type.eq.1) then
                 ro(i) = ro(i)*T(i)/T_wall
                 T(i) = T_wall
              end if
           end if                 
           T_bound(j) = T(i)
        end do
     end if   
 
  
  
     return
  end subroutine make_1d_flame
!! ------------------------------------------------------------------------------------------------  
  subroutine add_density_stratification
     integer(ikind) :: i,ispec,j
     real(rkind) :: c,Rmix_local,x,y,z
         
     !! Loop over all nodes and impose values
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
                
        !! Add a perturbation
        call make_perturbation(y/(ymax-ymin),c)
        x=x+c
        
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal              
        
        !! Density
        ro(i) = ro(i)*exp(grav(1)*(x-fl_pos_x)/(Rmix_local*T(i)))
        
                        
     end do
     !$omp end parallel do
         
     return
  end subroutine add_density_stratification   
!! ------------------------------------------------------------------------------------------------ 
  subroutine make_gaussian_hotspot
     !! Add a Gaussian hotspot to the flow. The composition and velocity are unchanged.
     !! Temperature follows a prescribed Gaussian, whilst density is modified to ensure pressure
     !! is unchanged.
     integer(ikind) :: i,ispec,j
     real(rkind) :: pR_local
     real(rkind) :: c,u_reactants,Rmix_local,x,y,z
     real(rkind) :: fl_thck_nd

     !! Scale thickness because position vectors are scaled...
     fl_thck_nd = fl_thck/L_char 
   
     !! Loop over all nodes and impose values
     !$omp parallel do private(x,y,z,c,pR_local)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Get local pressure
        pR_local = ro(i)*T(i) 
               
        !! Gaussian progress variable
        c = exp(-((x-fl_pos_x)/fl_thck_nd)**two - ((y-fl_pos_y)/fl_thck_nd)**two) !! 2D
!        c = exp(-((x-fl_pos_x)/fl_thck_nd)**two)                               !! 1D 
                
        !! Adjust temperature profile
        T(i) = T(i)*(one-c) + T_hot*c
        
        !! Composition - unchanged
       
        !! Density
        ro(i) = pR_local/(T(i))
        
        !! Velocity
        !! ...is unchanged.
                        
     end do
     !$omp end parallel do
     
     !! Special apply values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              if(flag_wall_type.eq.1) then
                 ro(i) = ro(i)*T(i)/T_wall
                 T(i) = T_wall
              end if
           end if                 
           T_bound(j) = T(i)
        end do
     end if   
    
     return
  end subroutine make_gaussian_hotspot
!! ------------------------------------------------------------------------------------------------
  subroutine make_baseflow
     !! Initialise the base flow: uniform ro, T, p, Y=Y_reactants, and spatially varying u,v,w
     integer(ikind) :: i,ispec,j
     real(rkind) :: Rmix_local,x,y,z,ro_inflow


#ifndef isoT
     !! Inflow mixture gas constant
     Rmix_local = zero
     do ispec=1,nspec
        Rmix_local = Rmix_local + Yspec_reactants(ispec)*one_over_molar_mass(ispec)
     end do
     Rmix_local = Rmix_local*Rgas_universal     

     !! Inflow density based on reference pressure and temperature
     ro_inflow = p_ref/(Rmix_local*T_ref)
#else
     ro_inflow = rho_char
#endif     
    
     !! Loop over all nodes and impose values
     !$omp parallel do private(x,y,z,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
              
        !! Temperature profile
        T(i) = T_ref
        
        !! Composition - reactants everywhere
        do ispec = 1,nspec
           Yspec(i,ispec) = Yspec_reactants(ispec)
        end do
              
        !! Density
        ro(i) = ro_inflow
        
        !! Velocity
        y=y/(ymax-ymin)
        u(i) = u_inflow_start*(uprof_a0 + uprof_a1*y + uprof_a2*y*y)
        v(i) = zero
        w(i) = zero
                        
     end do
     !$omp end parallel do
     
     !! Special apply values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              if(flag_wall_type.eq.1) then
                 ro(i) = ro(i)*T(i)/T_wall
                 T(i) = T_wall
              end if
           end if                 
           T_bound(j) = T(i)
        end do
     end if   
    
    
     return
  end subroutine make_baseflow  
!! ------------------------------------------------------------------------------------------------  
  subroutine load_flame_file
     !! Load data for a one-dimensional flame generated with the oned branch on this code.
     use thermodynamics
     integer(ikind) :: i,ispec,j,nflamein
     real(rkind) :: cell_pos,flamein_x0,flamein_xN
     real(rkind) :: c,x,y,z,ptbn
     real(rkind) :: roR,roP,uR,uP,SL
     real(rkind),dimension(:,:),allocatable :: flamein_data
     real(rkind) :: dx_flamein

     !! Open a file containing flame profile
     open(unit=19,file='init_flame.in')

     !! Allocate space for the input data
     read(19,*) nflamein
     allocate(flamein_data(nflamein,9+nspec));
     do i=1,nflamein
        read(19,*) flamein_data(i,:)
     end do
     close(19)
     
     !! Flamein data spacing (assumed uniform)
     dx_flamein = flamein_data(2,1)-flamein_data(1,1)
     flamein_x0 = flamein_data(1,1)
     flamein_xN = flamein_data(nflamein,1)
     
     !! Evaluate flame speed:
     roR = flamein_data(1,6)   !! Reactant density
     roP = minval(flamein_data(:,6))  !! Product density
     uR = flamein_data(1,3)          !! Reactant velocity
     uP = maxval(flamein_data(:,3))  !! Product velocity
     SL = (uP-uR)*roP/(roR-roP)   !! (eqn 1)     
       
     !! Loop through all particles. Find the "cell" the particle resides in. Copy data.     
     !$omp parallel do private(j,x,y,z,c,ispec,cell_pos,ptbn)
     do i=1,npfb
        x=rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Make perturbation
        call make_perturbation(y/(ymax-ymin),ptbn)
        x=x+ptbn     
     
        !! Dimensional coordinate of particle relative to desired flame position
        x = (x-fl_pos_x)*L_char 
        
        !! Nearest index in flame-in data (to left of x) and position within "cell"
        if(x.le.flamein_x0) then
           j=1
           ro(i) = flamein_data(j,6)
!           u(i) = flamein_data(j,3);v(i) = flamein_data(j,4);w(i) = flamein_data(j,5)
           roE(i) = flamein_data(j,7)
           do ispec=1,nspec
              Yspec(i,ispec) = flamein_data(j,9+ispec)
           end do
           T(i) = flamein_data(j,8)
           p(i) = flamein_data(j,9)
        else if(x.ge.flamein_xN) then
           j=nflamein
           ro(i) = flamein_data(j,6)
!           u(i) = flamein_data(j,3);v(i) = flamein_data(j,4);w(i) = flamein_data(j,5)
           roE(i) = flamein_data(j,7)
           do ispec=1,nspec
              Yspec(i,ispec) = flamein_data(j,9+ispec)
           end do           
           T(i) = flamein_data(j,8)
           p(i) = flamein_data(j,9)
        else
           j = floor((x-flamein_x0)/dx_flamein) + 1
           cell_pos = (x-flamein_data(j,1))/dx_flamein
           ro(i) = flamein_data(j,6)*(one-cell_pos) + flamein_data(j+1,6)*cell_pos
!           u(i) = flamein_data(j,3)*(one-cell_pos) + flamein_data(j+1,3)*cell_pos
!           v(i) = flamein_data(j,4)*(one-cell_pos) + flamein_data(j+1,4)*cell_pos
!           w(i) = flamein_data(j,5)*(one-cell_pos) + flamein_data(j+1,5)*cell_pos                                 
           roE(i) = flamein_data(j,7)*(one-cell_pos) + flamein_data(j+1,7)*cell_pos           
           do ispec=1,nspec
              Yspec(i,ispec) = flamein_data(j,9+ispec)*(one-cell_pos) + flamein_data(j+1,9+ispec)*cell_pos           
           end do
           T(i) = flamein_data(j,8)*(one-cell_pos) + flamein_data(j+1,8)*cell_pos
           p(i) = flamein_data(j,9)*(one-cell_pos) + flamein_data(j+1,9)*cell_pos
        end if
        

        !! Set velocities in case flame input file has different reactant velocity
        u(i) = SL*(roR-ro(i))/ro(i) + u(i)
        v(i) = zero
        w(i) = zero
        
        !! Convert to conservative variables
        rou(i) = u(i)*ro(i)
        rov(i) = v(i)*ro(i)
        row(i) = w(i)*ro(i)
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec(i,ispec)*ro(i)
        end do

     end do
     !$omp end parallel do
               
     
     !! Free up space
     deallocate(flamein_data)
     
     !! Temporarily copy some energy data to halos and mirrors (it will be later overwritten, but
     !! just prevents the NR solver from crashing at set-up)
     !$omp parallel do
     do i=npfb+1,np
        roE(i) = roE(1)
        ro(i) = ro(1)
        u(i) = u(1);v(i) = v(1);w(i) = w(1)
        rou(i) = rou(1);rov(i)=rov(1);row(i) = row(1)
        Yspec(i,:) = Yspec(1,:)
     end do
     !$omp end parallel do
     
     !! Re-evaluate temperature from energy.  
!     call evaluate_temperature_and_pressure
     
     !! Special apply values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              if(flag_wall_type.eq.1) then
                 ro(i) = ro(i)*T(i)/T_wall
                 T(i) = T_wall
              end if
           end if                 
           T_bound(j) = T(i)
        end do
     end if    
  
     !! Return Yspec to primitive form
     !$omp parallel do private(ispec)
     do i=1,npfb
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec(i,ispec)/ro(i)
        end do
     end do
     !$omp end parallel do
     
  
     return
  end subroutine load_flame_file
!! ------------------------------------------------------------------------------------------------
  subroutine hardcode_initial_conditions
     !! Temporary routine to generate initial conditions from some hard-coded functions.
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro,Rmix_local
     
     !! Values within domain
     !$OMP PARALLEL DO PRIVATE(x,y,z,tmp,ispec,Rmix_local)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)

        !! TG 3D Re1600 as in Cant 2022, Sandam 2017 etc (ish)
        u(i) = half*sqrt(three)*sin(x)*cos(y)*cos(z)!*oosqrt2
        v(i) =-half*sqrt(three)*cos(x)*sin(y)*cos(z)    !!c c
        w(i) = zero!u(i);u(i)=zero
                       
        
#ifndef isoT
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal               
        p(i) = p_ref + (ro(i)*U_char*U_char/16.0d0)*(cos(two*x)+cos(two*y))*(two+cos(two*z))      
        T(i) = T_ref
        ro(i) = p(i)/(Rmix_local*T(i))
#else
        p(i) = ro(i)*csq
#endif        
        
     !1.4028035124466285              
     end do
     !$OMP END PARALLEL DO
     
    
     
     !! Special apply values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              if(flag_wall_type.eq.1) then
                 ro(i) = ro(i)*T(i)/T_wall
                 T(i) = T_wall
              end if
           end if                 
           T_bound(j) = T(i)
        end do
     end if   
  
     return
  end subroutine hardcode_initial_conditions  
!! ------------------------------------------------------------------------------------------------
  subroutine load_restart_file
     !! Load initial conditions from a dump file
     integer(ikind) :: k,i,j
     real(rkind) :: tmp,tmpro
     character(70) :: fname,fname2  

#ifdef mp
     k=10000+iproc
#else
     k=10000
#endif     

     !! Construct the file name:
     write(fname,'(A17,I5)') './restart/fields_',k
     write(fname2,'(A16,I5)') './restart/nodes_',k

     !! Load the "smoothing length" from nodes file
     open(15,file=fname2)
     read(15,*) k
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart. NODES FILE."
     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(15,*) tmp,tmp,tmp,tmp,h(i),k
#else
        read(15,*) tmp,tmp,tmp,h(i),k
#endif        
        if(k.ne.node_type(i)) then
           write(6,*) "ERROR: Problem in restart file. STOPPING."
#ifdef mp
           call MPI_Abort(MPI_COMM_WORLD, k, ierror)
#else
           stop
#endif
        end if
     end do
     close(15)

     !! Open the field files
     open(14,file=fname)
     read(14,*) !! Skip line
     read(14,*) k
     read(14,*) emax_np1,emax_n,emax_nm1,dt
     read(14,*) eflow_nm1,sum_eflow,driving_force
     read(14,*) !! Skip line
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart. FIELDS FILE."

     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(14,*) tmpro,u(i),v(i),w(i),tmp,T(i),p(i),hrr(i),Yspec(i,1:nspec)
#else
        read(14,*) tmpro,u(i),v(i),tmp,T(i),p(i),hrr(i),Yspec(i,1:nspec)
#endif        
        ro(i) = tmpro
        p(i) = p(i) + p_ref !! The output files hold p-p_ref
     end do      
     close(14)
     
     !! Special apply values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              if(flag_wall_type.eq.1) then
                 ro(i) = ro(i)*T(i)/T_wall
                 T(i) = T_wall
              end if
           end if                 
           T_bound(j) = T(i)
        end do
     end if   

  
     return
  end subroutine load_restart_file
!! ------------------------------------------------------------------------------------------------  
  subroutine initialise_perturbation
     !! Create the mode amplitudes and phases for a peturbed interface
     integer(ikind) :: i

     mode_amp = zero
     mode_phase = zero
     do i=1,n_modes
        mode_amp(i) = rand()*fl_thck*ptbn_size/L_char
        mode_phase(i) = two*pi*rand()
     end do
     return
  end subroutine initialise_perturbation
!! ------------------------------------------------------------------------------------------------
  subroutine make_perturbation(y,c)
     !! Return a scalar which adjusts a flame position according to the perturbation
     real(rkind),intent(in) :: y !! y is in [0,1] or [-1/2,1/2]
     real(rkind),intent(out) :: c
     integer(ikind) :: i
     
     c=zero
     do i=1,n_modes
        c = c + mode_amp(i)*sin(two*pi*dble(i)*y + mode_phase(i))
     end do
  
     return
  end subroutine make_perturbation
!! ------------------------------------------------------------------------------------------------  
!! ------------------------------------------------------------------------------------------------
end module setup_flow
