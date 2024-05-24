module chemistry
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2022 onwards     |Main developer                     
  !! JRCK               |Feb 2023         |Debug and optimise OMP parts
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to evaluate chemical production rates for reacting flows.
  !! If the flow is inert, it does nothing. 
  !! Reaction rates are added to the rhs_Yspec arrays.
  !! N.B. Yspec contains ro*Y
  !!
  !! For some loops in here, the overheads of OMP parallelisation are not worth it, and serial is
  !! always faster.
  
  use kind_parameters
  use common_parameter
  use common_vars
  use thermodynamics
  use omp_lib
  implicit none


contains
!! ------------------------------------------------------------------------------------------------  
  subroutine calculate_chemical_production_rate
#ifdef react
     integer(ikind) :: j,i,ispec,istep,jspec,jstep,ithirdbody,kspec
     real(rkind) :: arrhenius_rate_back,gibbs_tmp,tmpro,nu_dash_local
     real(rkind) :: mass_production_rate,arrhenius_rate0,arrhenius_rate
     real(rkind) :: p_reduced,net_rate,logroYovW,enthalpy,heat_release
     real(rkind) :: delta_nu_local,nu_ddash_local
     real(rkind),dimension(:,:),allocatable :: rateYspec,gibbs
     real(rkind),dimension(:),allocatable :: rate,third_body_conc,backward_rate
         
     segment_tstart = omp_get_wtime()                  
     
     allocate(rateYspec(npfb,nspec));rateYspec = zero
     allocate(rate(npfb));rate=zero
     allocate(third_body_conc(npfb));third_body_conc=one
     allocate(backward_rate(npfb));backward_rate = zero        
   
     
     !! Pre-evaluate gibbs functions as required ==============================
     if(num_gibbs_species.ne.0) then
        !! Allocate space
        allocate(gibbs(npfb,num_gibbs_species));gibbs=zero
        
        !! Loop over all species 
        do ispec = 1,nspec
           if(gibbs_flag_species(ispec).ne.0) then !! If this species needs gibbs evaluation
              jspec = gibbs_flag_species(ispec)  !! jspec is the species index within the gibbs-list
              
              !! Loop over all nodes
              !$omp parallel do private(gibbs_tmp)
              do i=1,npfb
                 !! Evaluate gibbs
                 call evaluate_gibbs_at_node(T(i),log(T(i)),ispec,gibbs_tmp)                 
                 
                 !! Store in gibbs array
                 gibbs(i,jspec) = gibbs_tmp
              end do
              !$omp end parallel do
           
           end if
        end do            
     end if

     !! =======================================================================        
     !! Loop over all steps ===================================================
     do istep = 1,nsteps       

        !! Forward rate =======================================================

        !! Evaluate Arrhenius rate for this step
        !$omp parallel do private(arrhenius_rate)
        do i=1,npfb        
           arrhenius_rate = arrhenius_coefs(istep,1) &
                          + arrhenius_coefs(istep,2)*log(T(i)) &
                          - arrhenius_coefs(istep,3)/T(i)
 
           !! Loop over all reactants and build log of forward rate
           rate(i) = arrhenius_rate
        end do
        !$omp end parallel do
        !! rate contains ln(k_{f,m})

                    
        !! Third bodies =======================================================
        if(third_body_flag(istep).ne.0) then
           ithirdbody = third_body_flag(istep) !! This is the type of third body for step istep
           
           !! Build third-body concentration
           !$omp parallel do
           do i=1,npfb
              third_body_conc(i) = zero
              do ispec = 1,nspec
                 third_body_conc(i) = third_body_conc(i) &
                                    + third_body_efficiencies(ithirdbody,ispec)* &
                                      Yspec(i,ispec)* &
                                      one_over_molar_mass(ispec)
              end do
           end do
           !$omp end parallel do
        else
           third_body_conc(:) = one
        end if   

           
        !! Lindemann steps ====================================================
        if(lindemann_form_flag(istep).ne.0) then
           jstep = lindemann_form_flag(istep)  !! This is the jstep-th Lindemann step

           !$omp parallel do private(arrhenius_rate0,p_reduced)
           do i=1,npfb              
              !! Evaluate k0
              arrhenius_rate0 = lindemann_coefs(jstep,1) + &
                                lindemann_coefs(jstep,2)*log(T(i)) - &
                                lindemann_coefs(jstep,3)/T(i)

              !! Reduced pressure
              p_reduced = exp(arrhenius_rate0 - rate(i))*third_body_conc(i)
              
              !! rate constant for this step
              rate(i) = rate(i) + log(p_reduced/(one+p_reduced)) &
                                + lindemann_coefs(jstep,4)       
              
              !! Reset third body concentrations here
              third_body_conc(i) = one               
           end do
           !$omp end parallel do
        end if                     
           
        !! Store log(k_f) for use in any backward steps
        backward_rate(:) = rate(:)

        !! At this stage, rate contains ln(k_f) for both regular and Lindemann steps
                   
        !! Finalise forward rate ==============================================
        !! Multiply up each reactant contrib (in log space)
        do jspec = 1,num_reactants(istep)
           ispec = reactant_list(istep,jspec)  !! ispec is the jspec-th reactant of step istep

           !! Store nu_dash for this step and species
           nu_dash_local = nu_dash(istep,ispec)

           !$omp parallel do private(logroYovW)
           do i=1,npfb
              
              logroYovW = log(max(Yspec(i,ispec),verysmall)*one_over_molar_mass(ispec))
              
              rate(i) = rate(i) + nu_dash_local*logroYovW

           end do
           !$omp end parallel do
        end do        

        !! Backward rate ======================================================
        if(gibbs_rate_flag(istep).eq.1) then

           !! Evaluate backwards rate                                    
           !! Loop over all species in step
           do jspec = 1,num_reactants(istep) + num_products(istep)
              ispec = stepspecies_list(istep,jspec)
                                 
              !! index of this species in gibbs list
              kspec = gibbs_flag_species(ispec)

              !! Local delta_nu
              delta_nu_local = delta_nu(istep,ispec)

              backward_rate(1:npfb) = backward_rate(1:npfb) + delta_nu_local*gibbs(1:npfb,kspec)              
           end do          
           !! backward_rate contains k_{b,m}
              
           !! Loop over all products and build log of backward_rate
           do jspec = 1,num_products(istep)
              ispec = product_list(istep,jspec)
              
              !! Store nu_ddash for this step and species                 
              nu_ddash_local = nu_ddash(istep,ispec)

              !$omp parallel do private(logroYovW)
              do i=1,npfb                 

                 logroYovW = log(max(Yspec(i,ispec),verysmall)*one_over_molar_mass(ispec))

                 backward_rate(i) = backward_rate(i) + nu_ddash_local*logroYovW           
              end do
              !$omp end parallel do
           end do
        else                    
           !! No backward production otherwise       
           backward_rate(:) = -verylarge
        end if

                
        !! Net production rate ================================================
        !$omp parallel do private(net_rate)
        do i=1,npfb         
           !! net
           net_rate = exp(rate(i)) - exp(backward_rate(i))
           
           !! Add third body concs
           rate(i) = net_rate*third_body_conc(i)
        end do
        !$omp end parallel do
       
        !! Loop over all species in this step, and add to the total rate for that species
        do jspec = 1,num_reactants(istep) + num_products(istep)
           ispec = stepspecies_list(istep,jspec)        

           !! Add net molar production rate for this species 
           rateYspec(1:npfb,ispec) = rateYspec(1:npfb,ispec) + delta_nu(istep,ispec)*rate(1:npfb)           
        end do
              
           
     end do !! End of steps loop ==============================================
     !! =======================================================================                 
         
     !! Add rate onto RHS for each node and each species ====================== 
     hrr = zero     
     do ispec = 1,nspec
        !$omp parallel do
        do i=1,npfb                   
        
           !! Convert from molar to mass production rate
           rateYspec(i,ispec) = rateYspec(i,ispec)*molar_mass(ispec)
        
           !! Augment the RHS of Yspec for species ispec
           rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) + rateYspec(i,ispec)           
                              
           !! Augment the heat release (production rate x enthalpy of formation)
           hrr(i) = hrr(i) - rateYspec(i,ispec)*coef_h(ispec,polyorder_cp+2)

        end do
        !$omp end parallel do
        
     end do


     !! Build contribution to source terms for boundary conditions ============
     if(nb.ne.0) then
        allocate(sumoverspecies_homega(nb));sumoverspecies_homega=zero
        allocate(reaction_rate_bound(nb,nspec));reaction_rate_bound=zero
        !$omp parallel do private(i,ispec,enthalpy)
        do j=1,nb
           i=boundary_list(j)
              
           !! Loop over species
           do ispec=1,nspec                                    
              !! Store the reaction rate on the boundary
              reaction_rate_bound(j,ispec) = rateYspec(i,ispec)                                   
              
              !! Evaluate enthalpy of species ispec
              call evaluate_enthalpy_only_at_node(T(i),ispec,enthalpy)           
           
              !! Evaluate "reduced enthalpy"
              enthalpy = enthalpy - cp(i)*T(i)*Rgas_universal/(Rgas_mix(i)*molar_mass(ispec))
           
              !! Augment sum of reduced_h*omega
              sumoverspecies_homega(j) = sumoverspecies_homega(j) + &
                                      enthalpy*rateYspec(i,ispec)        
              
           end do
        end do
        !$omp end parallel do
     end if                                        
     
     !! De-allocation of arrays
     if(num_gibbs_species.ne.0) deallocate(gibbs)
     deallocate(rateYspec,rate,backward_rate,third_body_conc)
  
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(6) = segment_time_local(6) + segment_tend - segment_tstart  

#else
     !! Dummy routine if we want to compile with react, but suppress chemical reactions (i.e. for
     !! pre-ignition inert simulations.       
     segment_tstart = omp_get_wtime()                  

     !! Build contribution to source terms for boundary conditions ============
     if(nb.ne.0) then
        allocate(sumoverspecies_homega(nb));sumoverspecies_homega=zero
        allocate(reaction_rate_bound(nb,nspec));reaction_rate_bound=zero
     end if           
       
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(6) = segment_time_local(6) + segment_tend - segment_tstart    
#endif  
     return
  end subroutine calculate_chemical_production_rate  
!! ------------------------------------------------------------------------------------------------ 
end module chemistry

