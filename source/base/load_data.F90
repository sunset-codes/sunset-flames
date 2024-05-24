module load_data
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines read in control data.
  use kind_parameters
  use common_parameter
  use common_vars
#ifdef mp
  use mpi_transfers
#endif    
  implicit none
  
contains
!! ------------------------------------------------------------------------------------------------
  subroutine load_control_data_LUonly
     integer(ikind) :: dummy_int
     real(rkind) :: dummy_real
     
     !! Load data from the control file
     open(unit=12,file='control.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)

     !! Length-scale
     read(12,*)
     read(12,*) L_char
     read(12,*)
     
     !! Velocity-scale
     read(12,*) 
     read(12,*) U_char
     read(12,*)
     
     !! Set Z-length-scale and characteristic time-scale
     Time_char = L_char/u_char
     
     
     !! Default set nspec in case of non-multispecies
     nspec =1
     close(12)
          
     return
  end subroutine load_control_data_LUonly
!! ------------------------------------------------------------------------------------------------  
  subroutine load_control_data_all
     integer(ikind) :: dummy_int
     real(rkind) :: dummy_real
     
     !! Load data from the control file
     open(unit=12,file='control.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)

     !! Length-scale
     read(12,*)
     read(12,*) dummy_real
     read(12,*)
     
     !! Velocity-scale
     read(12,*) 
     read(12,*) dummy_real
     read(12,*)
     
     !! Start and end time (in multiples of characteristic time)
     read(12,*)
     read(12,*) time,time_end
     read(12,*)
     time = time*Time_char;time_end = time_end*Time_char
     itime = 0
     
     !! Output frequency (in multiples of characteristic time)
     read(12,*)
     read(12,*) dt_out,dt_out_stats
     read(12,*)
     dt_out = dt_out*Time_char
     dt_out_stats = dt_out_stats*Time_char
        
     !! Gravity
     read(12,*)
     read(12,*) grav(:)
     read(12,*)
          
     !! Inflow equivalence ratio
     read(12,*) 
     read(12,*) phi_in
     read(12,*) 
     
     !! Reference temp
     read(12,*)
     read(12,*) T_ref
     read(12,*) 
     
     !! Reference viscosity
     read(12,*)
     read(12,*) visc_ref
     read(12,*) 
     
     !! Reference pressure
     read(12,*)
     read(12,*) p_ref
     read(12,*) 
         
     !! Prandtl number
     read(12,*)
     read(12,*) Pr
     read(12,*) 
     
     !! Diffusion flag (1 gives mixture averaging)
     read(12,*)
     read(12,*) flag_mix_av
     read(12,*)
     
     !! Mach number (only used for isothermal flows)
     read(12,*)
     read(12,*) Ma
     read(12,*) 
     
     !! set the sound speed squared
#ifdef isoT
     csq = (u_char/Ma)**two     
#endif      

     !! Read in T-exponent for TDTP if required
     read(12,*)
     read(12,*) r_temp_dependence
     read(12,*)     
     
     !! Read in inflow boundary type
     read(12,*)
     read(12,*) flag_inflow_type
     read(12,*)

     !! Read inflow velocity control parameters
     read(12,*)
     read(12,*) flag_uinflow_control,u_inflow_start,u_inflow_end,u_inflow_ramptime
     read(12,*)     
     u_inflow_start=u_inflow_start*u_char
     u_inflow_end = u_inflow_end*u_char
     u_inflow_ramptime = u_inflow_ramptime*Time_char    
     
     !! Read in wall boundary type
     read(12,*)
     read(12,*) flag_wall_type, T_wall
     read(12,*)     

     !! Read in base flow type
     read(12,*)
     read(12,*) flag_base_flow_profile
     read(12,*)  
     
     !! Specify coefficients for base flow profile u=U*(a0 + a1*y + a2*y*y)
     if(flag_base_flow_profile.eq.1) then  !! Uniform flow
        uprof_a0=one;uprof_a1=zero;uprof_a2=zero
     else if(flag_base_flow_profile.eq.2) then  !! Poiseuille
        uprof_a0=1.5d0;uprof_a1=zero;uprof_a2=-six
     else if(flag_base_flow_profile.eq.3) then  !! Half-Poiseuille
        uprof_a0=zero;uprof_a1=zero;uprof_a2=1.5d0
     else                                           !! Zero-flow
        uprof_a0=zero;uprof_a1=zero;uprof_a2=zero        
     end if
     


     !! Read in initial condition type
     read(12,*)
     read(12,*) flag_flow_type
     read(12,*)
     
     !! Read in initial flame position, size and temperature
     read(12,*)
     read(12,*) fl_pos_x,fl_pos_y,fl_thck,T_hot
     read(12,*)
     
     !! Read in turbulence flags/params
     read(12,*)
     read(12,*) flag_turbulent,turb_lengthscale,turb_intensity
     read(12,*)
     
     !! Read in density stratification flag
     read(12,*)
     read(12,*) flag_strat
!     read(12,*)
     
     close(12)
     
     !! Set characteristic density to unity. It is only used for isothermal flows
     rho_char = one
#ifdef isoT
     p_ref = rho_char*csq
#endif     
          
     return
  end subroutine load_control_data_all
!! ------------------------------------------------------------------------------------------------  
  subroutine load_chemistry_data
     integer(ikind) :: ispec,iorder,istep,dummy_int,jspec
     real(rkind) :: dummy_real
     
     !! Load data from the thermochemistry control file
     open(unit=12,file='thermochem.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)
     
     !! Number of chemical species
     read(12,*)
     read(12,*) nspec
     read(12,*)

     !! This section reads in transport data for all species.
     
     !! Number of coefs and order of polynomial for cp(T)
     read(12,*)
     read(12,*) ncoefs_cp
     read(12,*)
     !! ncoefs = (polyorder + 1) + 1 + 1: terms in polynomial,coef for h, coef for s
     polyorder_cp = ncoefs_cp - 3
     
     !! Temperature limits
     read(12,*)
     read(12,*) T_low,T_high
     read(12,*)
     
     !! Reference pressure for Gibbs functions
     read(12,*)
     read(12,*) p_ref_gibbs
     read(12,*)

     !! Allocate space for molar mass, Lewis number, and polynomial fitting for cp(T)   
     allocate(molar_mass(nspec),one_over_Lewis_number(nspec),one_over_molar_mass(nspec))
     allocate(coef_cp(nspec,ncoefs_cp),coef_h(nspec,ncoefs_cp),coef_dcpdT(nspec,ncoefs_cp))
     allocate(coef_gibbs(nspec,ncoefs_cp))
          
     !! Load molar mass, Lewis, and polynomial fits   
     read(12,*) !! Read comment line  
     do ispec = 1,nspec
        read(12,*) !! Read species identifier comment line
        read(12,*) molar_mass(ispec)
        read(12,*) one_over_Lewis_number(ispec)
        one_over_Lewis_number(ispec) = one/one_over_Lewis_number(ispec)
        one_over_molar_mass(ispec) = one/molar_mass(ispec)

        read(12,*)  !! Comment line
        do iorder=1,ncoefs_cp
           read(12,*) coef_cp(ispec,iorder)
        end do
        read(12,*) !! Blank line                            
                           
        !! Pre-divide gibbs coefficients (and leave in molar form). Also include log(p_ref_gibbs/R0)
        !! and a log(T) term, so result of creating gibbs includes all terms required for
        !! backwards rate calculation
        coef_gibbs(ispec,:) = coef_cp(ispec,:)
        do iorder = 2,polyorder_cp+1
           coef_gibbs(ispec,iorder) = coef_gibbs(ispec,iorder)/dble(iorder*(iorder-1))
        end do                           
        coef_gibbs(ispec,1) = coef_cp(ispec,polyorder_cp+3) &
                            - coef_cp(ispec,1) &
                            + log(p_ref_gibbs/Rgas_universal)
        coef_gibbs(ispec,polyorder_cp+3) = coef_cp(ispec,1) - one
                     
                           
        !! Convert cp coefficients from molar to mass based
        coef_cp(ispec,:) = coef_cp(ispec,:)*Rgas_universal*one_over_molar_mass(ispec)

        !! Pre-divide coefs by iorder for h.
        do iorder = 1,polyorder_cp + 1
           coef_h(ispec,iorder) = coef_cp(ispec,iorder)/dble(iorder)
        end do
        coef_h(ispec,polyorder_cp+2) = coef_cp(ispec,polyorder_cp+2)
        
        !! Pre-multiply coefs by iorder-1 for dcp/dT
        do iorder = 1,polyorder_cp + 1
           coef_dcpdT(ispec,iorder) = coef_cp(ispec,iorder)*dble(iorder-1)
        end do
                
     end do     

     !! Next section reads in reaction mechanism. 
     read(12,*) !! Read comment line    
     
     !! Number of steps
     read(12,*)
     read(12,*) nsteps
     read(12,*)
     
     if(nsteps.ne.0)then
        !! Number of different third body efficiencies     
        read(12,*)
        read(12,*) nthirdbodies
        read(12,*)
     end if     
     
     !! Space for rate constants and coefficients etc
     allocate(Arrhenius_coefs(nsteps,3))
     
     !! Numbers of r and p, and lists
     allocate(num_reactants(nsteps),num_products(nsteps))
     allocate(reactant_list(nsteps,3),product_list(nsteps,3))  !! Limit of 3 reactants and 3 products in any given step
     allocate(stepspecies_list(nsteps,6)) !! List of all species in reaction (reactants then products)
     
     !! Stoichiometric coefficients
     allocate(nu_dash(nsteps,nspec),nu_ddash(nsteps,nspec),delta_nu(nsteps,nspec)) 
     nu_dash = zero;nu_ddash = zero;delta_nu = zero
     
     !! Flags, efficiencies and coefficients
     allocate(gibbs_rate_flag(nsteps),lindemann_form_flag(nsteps)) !! Flags for backwards and lindemann
     allocate(third_body_flag(nsteps))
     if(nthirdbodies.ne.0) allocate(third_body_efficiencies(nthirdbodies,nspec))
     allocate(lindemann_coefs(nsteps,4));lindemann_coefs = zero
    
     
     read(12,*) !! Comment line
     if(nsteps.ne.0)then
        do istep = 1,nsteps
           read(12,*) dummy_int  !! Reaction number
           if(dummy_int.ne.istep) write(6,*) "Warning, error reading reaction mech. Expect seg fault."

           read(12,*) num_reactants(istep)            !! Number of reactants
           do ispec = 1,num_reactants(istep)  !! Loop over all reactants
              read(12,*) dummy_int,dummy_real
              reactant_list(istep,ispec) = dummy_int !! Identity of reactant
              stepspecies_list(istep,ispec) = dummy_int
              nu_dash(istep,dummy_int) = dummy_real  !! Stoichiometric coefficient
           end do

           read(12,*) num_products(istep)            !! Number of products
           do ispec = 1,num_products(istep)  !! Loop over all products
              read(12,*) dummy_int,dummy_real
              product_list(istep,ispec) = dummy_int !! Identity of product
              stepspecies_list(istep,num_reactants(istep) + ispec) = dummy_int
              nu_ddash(istep,dummy_int) = dummy_real  !! Stoichiometric coefficient
           end do
        
           !! Calculate deltas
           delta_nu(istep,:) = nu_ddash(istep,:) - nu_dash(istep,:)
        
           !! Coefficients for arrhenius rate constant
           read(12,*) Arrhenius_coefs(istep,1:3)
   
           !! Take logarithm of pre-exponential factor
           Arrhenius_coefs(istep,1) = log(Arrhenius_coefs(istep,1))
           Arrhenius_coefs(istep,3) = Arrhenius_coefs(istep,3)/(Rgas_universal)                     
        
           !! Gibbs based backwards rate?
           read(12,*) gibbs_rate_flag(istep)
           
           !! Lindemann form?
           read(12,*) lindemann_form_flag(istep)
           
           !! Third bodies?
           read(12,*) third_body_flag(istep)        
                 
           read(12,*) !! Blank line      
        end do

        !! Lists of third body efficiencies
        if(nthirdbodies.ne.0) then
        
           do istep = 1,nthirdbodies
              read(12,*) !! Comment line
              do ispec = 1,nspec
                 read(12,*) dummy_int,third_body_efficiencies(istep,ispec)     !! List of efficiencies
              end do
           end do
        end if
        read(12,*) !! Blank line
     
        !! Lists of Lindemann coefficients
        nlindemann = maxval(lindemann_form_flag(1:nsteps))
        if(nlindemann.ne.0) then
           read(12,*) !! Comment line
           do istep = 1,nlindemann
              read(12,*) dummy_int,lindemann_coefs(istep,1:4)
              !! Take logarithm of pre-exponential factor etc
              lindemann_coefs(istep,1) = log(lindemann_coefs(istep,1))
              lindemann_coefs(istep,3) = lindemann_coefs(istep,3)/(Rgas_universal)
              lindemann_coefs(istep,4) = log(lindemann_coefs(istep,4))
           end do         
        end if
     
        !! Build list of species which require gibbs evaluation
        allocate(gibbs_flag_species(nspec));gibbs_flag_species=0     
        do istep = 1,nsteps
           if(gibbs_rate_flag(istep).ne.0) then          !! For gibbs steps
              do jspec = 1, num_reactants(istep) + num_products(istep)  !! Loop over products and reactants
                 ispec = stepspecies_list(istep,jspec)              
                 gibbs_flag_species(ispec) = 1              !! Flag this species as requiring gibbs evaluation
              end do
           end if
        end do
        !! Count the number of species which require gibbs evaluation
        num_gibbs_species=0
        do ispec=1,nspec
           if(gibbs_flag_species(ispec).ne.0)then
              num_gibbs_species = num_gibbs_species  + 1           
              gibbs_flag_species(ispec) = num_gibbs_species  !! Modify flag to become a counter
           end if
        end do

     end if  !! End of "only read if not zero step mechanism" if-statement     
     close(12)               
     
     return
  end subroutine load_chemistry_data    
!! ------------------------------------------------------------------------------------------------
  subroutine load_transport_file
     !! Subroutine loads the transport file transport.in and prepares the coefficients for mixture
     !! average transport properties
     integer(ikind) :: i,j,ispec,jspec
     real(rkind) :: store1
     
     open(unit=13,file='transport.in')     
     
     !! Read header
     read(13,*)
     read(13,*)
     
     !! Reference temperature
     read(13,*)
     read(13,*) T_ref_mxav
     read(13,*)
     
     !! Reference pressure
     read(13,*)
     read(13,*) p_ref_mxav
     read(13,*)
     
     !! Polynomial coefficients for viscosity
     read(13,*)
     allocate(mxav_coef_visc(nspec,4))
     do ispec=1,nspec
        read(13,*) i,mxav_coef_visc(ispec,1),mxav_coef_visc(ispec,2),mxav_coef_visc(ispec,3),mxav_coef_visc(ispec,4)
     end do
     read(13,*)
              
     !! Polynomial coefficients for thermal conductivity
     read(13,*)
     allocate(mxav_coef_lambda(nspec,4))
     do ispec=1,nspec
        read(13,*) i,mxav_coef_lambda(ispec,1),mxav_coef_lambda(ispec,2), &
                     mxav_coef_lambda(ispec,3),mxav_coef_lambda(ispec,4)
     end do
     read(13,*)
     
     !! Polynomial coefficients for molecular diffusivity
     read(13,*)
     allocate(mxav_coef_Diff(nspec,nspec,4));mxav_coef_Diff=zero
     do ispec=1,nspec
        do jspec=1,ispec
           read(13,*) i,j,mxav_coef_Diff(ispec,jspec,1),mxav_coef_Diff(ispec,jspec,2), &
                          mxav_coef_Diff(ispec,jspec,3),mxav_coef_Diff(ispec,jspec,4)
        end do
     end do
     read(13,*)
     
     !! Copy lower triangular matrix to symmetric matrix (though not actually used)
     do ispec=1,nspec
        do jspec=1,ispec
           mxav_coef_Diff(jspec,ispec,1) = mxav_coef_Diff(ispec,jspec,1)
           mxav_coef_Diff(jspec,ispec,2) = mxav_coef_Diff(ispec,jspec,2)
           mxav_coef_Diff(jspec,ispec,3) = mxav_coef_Diff(ispec,jspec,3)
           mxav_coef_Diff(jspec,ispec,4) = mxav_coef_Diff(ispec,jspec,4)                                 
        end do
     end do                   
  
     return
  end subroutine load_transport_file
!! ------------------------------------------------------------------------------------------------
end module load_data
