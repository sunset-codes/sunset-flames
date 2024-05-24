module thermodynamics
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to calculate thermodynamic properties (e.g. p,T, from ro,u,roE,Y)
  
  !! Various thermodynamic options ::
  !! 1) isoT      - ISOTHERMAL FLOW. Don't solve an energy equation, p=ro*c*c with c a constant. Many
  !!                arrays are not used, and so not allocated.
  !! 2) not(isoT) - THERMAL FLOW: cp is a polynomial function of T (it can be a polynomial of 
  !!                order 0), and T is obtained from ro,u,roE,Y via solution of a non-linear 
  !!                equation with a Newton-Raphson method. 
  
  !! Evaluation of and direct reference to CHEMKIN polynomials (i.e. use of coef_cp,coef_h) should 
  !! only happen within the routines in this module.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  implicit none


contains
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_temperature_and_pressure
     !! Evaluate temperature numerically with cp = polynomial(T), and the pressure from P=ro*R*T
     !! Also evaluate the mixture gas constant, and the mixture specific heat capacity
     integer(ikind) :: i,ispec,iorder,NRiters,maxiters,sumiters
     real(rkind) :: deltaT,cp_tmp,tmpro
     real(rkind),dimension(:),allocatable :: fT_coef_C,dfT_coef_C
     real(rkind) :: fT_coef_C0
     real(rkind) :: T_tmp,fT,dfT
     real(rkind),parameter :: T_tolerance=1.0d-10
     integer(ikind),parameter :: NRiters_max=100
     logical :: keepgoing
   
     segment_tstart = omp_get_wtime()   

#ifndef isoT
     allocate(fT_coef_C(polyorder_cp+1),dfT_coef_C(polyorder_cp+1))
     
     !! Loop over ALL nodes
     maxiters = 0;sumiters=0
     !$omp parallel do private(fT_coef_C0,fT_coef_C,ispec,iorder,T_tmp,NRiters,fT,dfT,deltaT,dfT_coef_C, &
     !$omp cp_tmp,tmpro,keepgoing) &
     !$omp reduction(max:maxiters) reduction(+:sumiters)
     do i=1,np
                  
        !! Store inverse of density
        tmpro = one/ro(i) 

        !! Evaluate the gas constant for the mixture
        Rgas_mix(i) = zero
        do ispec = 1,nspec
           Rgas_mix(i) = Rgas_mix(i) + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rgas_mix(i) = tmpro*Rgas_mix(i)*Rgas_universal     
                  
        !!Initialise coefficients:
        fT_coef_C0 = half*ro(i)*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i)) - roE(i)  !! K.E. - roE
        fT_coef_C(1) = -Rgas_mix(i)*ro(i)      !! ro*R
        fT_coef_C(2:polyorder_cp+1) = zero
        
        !! Build coefficients
        do iorder = 1,polyorder_cp+1
           do ispec=1,nspec
              fT_coef_C(iorder) = fT_coef_C(iorder) + Yspec(i,ispec)*coef_h(ispec,iorder)
           end do       
        end do
        do ispec = 1,nspec
           fT_coef_C0 = fT_coef_C0 + Yspec(i,ispec)*coef_cp(ispec,polyorder_cp+2)
        end do
        
        !! Make coefficients for dfT
        do iorder = 1,polyorder_cp+1
           dfT_coef_C(iorder) = dble(iorder)*fT_coef_C(iorder)
        end do
                
        !! Initial guess for T is current temperature..
        T_tmp = T(i)
        
        !! Newton-Raphson iterations
        keepgoing = .true.
        NRiters = 0
        do while(keepgoing)
           NRiters = NRiters + 1
     
           !! Evaluate f(T) and f'(T)
           fT = fT_coef_C(polyorder_cp+1)*T_tmp
           dfT = dfT_coef_C(polyorder_cp+1)
           do iorder = polyorder_cp,1,-1
              fT = (fT + fT_coef_C(iorder))*T_tmp
              dfT = dfT*T_tmp + dfT_coef_C(iorder)
           end do
           fT = fT + fT_coef_C0           
        
           !! Calculate new T
           deltaT = - fT/(dfT)
           T_tmp = T_tmp + deltaT

           !! Check for convergence
           if(abs(deltaT).le.T_tolerance) then
              keepgoing = .false.
           end if
           if(NRiters.ge.NRiters_max) then
              keepgoing = .false.
           end if
           
        end do
                
        !! Pass new T back to temperature array
        T(i) = T_tmp
                
        !! Find the maximum number of iterations over this processor and sum of iterations
        maxiters = max(NRiters,maxiters)
        sumiters = sumiters + NRiters

        !! Evaluate the specific heat capacity of the mixture
        cp(i) = zero
        do ispec=1,nspec
           cp_tmp = coef_cp(ispec,polyorder_cp+1)
           do iorder=polyorder_cp,1,-1
              cp_tmp = cp_tmp*T(i) + coef_cp(ispec,iorder)
           end do  
           cp(i) = cp(i) + Yspec(i,ispec)*cp_tmp          
        end do      
        cp(i) = cp(i)*tmpro !! Divide by ro to get cp

        !! Evaluate the pressure        
        p(i) = ro(i)*Rgas_mix(i)*T(i)
     end do
     !$omp end parallel do
     
     deallocate(fT_coef_C,dfT_coef_C)
#else
     !! Isothermal, set constant T, and p proportional to density
     T(:) = zero
     !$omp parallel do private(tmpro)
     do i=1,np    !! N.B. this is over ALL nodes.
        p(i) = csq*ro(i)
     end do
     !$omp end parallel do      
     
#endif     
         
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(9) = segment_time_local(9) + segment_tend - segment_tstart  
       
     return
  end subroutine evaluate_temperature_and_pressure
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_enthalpy_at_node(Temp,ispec,enthalpy,cpispec,dcpdT)  
     !! Evaluate the enthalpy and cp of species ispec based on temperature at one node. Also
     !! evaluate the rate of change of cp with T.
     integer(ikind),intent(in) :: ispec
     real(rkind),intent(in) :: Temp
     real(rkind),intent(out) :: enthalpy,cpispec,dcpdT
     integer(ikind) :: iorder
     
!     segment_tstart = omp_get_wtime()        
     
     !! Enthalpy and cp
     enthalpy = Temp*coef_h(ispec,polyorder_cp+1)    
     cpispec = coef_cp(ispec,polyorder_cp+1)
    
     do iorder=polyorder_cp,1,-1

        enthalpy = Temp*(enthalpy + coef_h(ispec,iorder))   
        cpispec = cpispec*Temp + coef_cp(ispec,iorder)
    
     end do    
     enthalpy = enthalpy + coef_h(ispec,polyorder_cp+2)
    
     !! Rate of change of cp with T
     dcpdT = coef_dcpdT(ispec,polyorder_cp+1)
     
     do iorder=polyorder_cp,2,-1
     
        dcpdT = dcpdT*Temp + coef_dcpdT(ispec,iorder)
        
     end do

     !! Profiling
!     segment_tend = omp_get_wtime()
!     segment_time_local(9) = segment_time_local(9) + segment_tend - segment_tstart                        
                      
     return
  end subroutine evaluate_enthalpy_at_node
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_enthalpy_only_at_node(Temp,ispec,enthalpy)  
     !! Evaluate the enthalpy and cp of species ispec based on temperature at one node. As above, 
     !! but doesn't calculate anything else
     integer(ikind),intent(in) :: ispec
     real(rkind),intent(in) :: Temp
     real(rkind),intent(out) :: enthalpy
     integer(ikind) :: iorder
     
!     segment_tstart = omp_get_wtime()             
     
     !! Enthalpy only
     enthalpy = Temp*coef_h(ispec,polyorder_cp+1)        
     do iorder=polyorder_cp,1,-1
        enthalpy = Temp*(enthalpy + coef_h(ispec,iorder))     
     end do    
     enthalpy = enthalpy + coef_h(ispec,polyorder_cp+2)
    
     !! Profiling
!     segment_tend = omp_get_wtime()
!     segment_time_local(9) = segment_time_local(9) + segment_tend - segment_tstart                
                 
     return
  end subroutine evaluate_enthalpy_only_at_node  
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_gibbs_at_node(Temp,logT,ispec,gibbs)
     !! Evaluate the gibbs function of species ispec at a node, given T,logT and ispec
     !! N.B. actually returns molar_gibbs/(R0*T)
     integer(ikind),intent(in) :: ispec
     real(rkind),intent(in) :: Temp,logT
     real(rkind),intent(out) :: gibbs
     integer(ikind) :: iorder
     
     !! Polynomial (in T) terms            
     gibbs = coef_gibbs(ispec,polyorder_cp+1)
     do iorder=polyorder_cp,1,-1
        gibbs = coef_gibbs(ispec,iorder) + gibbs*Temp
     end do
     gibbs = coef_gibbs(ispec,polyorder_cp+2)/Temp &
           - coef_gibbs(ispec,polyorder_cp+3)*logT &
           - gibbs
     
     return
  end subroutine evaluate_gibbs_at_node
!! ------------------------------------------------------------------------------------------------
  function evaluate_sound_speed_at_node(cp_local,Rgm_local,T_local) result(c)
     !! Sound speed from cp,Rmix and T (or prescribed by csq if isoT)
     real(rkind), intent(in) :: cp_local,Rgm_local,T_local
     real(rkind) ::  c
#ifdef isoT
     c = sqrt(csq)
#else  
     c = dsqrt(cp_local*Rgm_local*T_local/(cp_local-Rgm_local))

#endif      
  end function evaluate_sound_speed_at_node
!! ------------------------------------------------------------------------------------------------  
  subroutine set_energy_on_bound(i,tmpT)
     !! Evaluate roE based on u,ro,T at a specific node
     integer(ikind),intent(in) :: i
     real(rkind),intent(in) :: tmpT
     integer(ikind) :: ispec
     real(rkind) :: enthalpy,Rgas_mix_local,cpispec,dummy_real
     
#ifndef isoT              
        !! Initialise roE with K.E. term
        roE(i) = half*ro(i)*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))

        !! Loop over species
        Rgas_mix_local = zero
        do ispec=1,nspec       
           !! Evaluate local species enthalpy
           call evaluate_enthalpy_at_node(tmpT,ispec,enthalpy,cpispec,dummy_real)
           
           !! Add species enthalpy contribution
           roE(i) = roE(i) + Yspec(i,ispec)*enthalpy
           
           !! Build the local mixture gas constant - actually holds ro*R
           Rgas_mix_local = Rgas_mix_local + Yspec(i,ispec)*Rgas_universal*one_over_molar_mass(ispec)
        end do           

        !! Subtract roRgasT
        roE(i) = roE(i) - Rgas_mix_local*tmpT     
                    
#endif

       
     return
  end subroutine set_energy_on_bound 
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_energy
     !! Evaluate roE based on u,ro,T, over the whole domain. 
     !! This routine is only called at start-up, and it is because loading temperature is a more
     !! intuitive variable to use for input than roE...
     
     !! Additionally calculate the pressure on outflow boundary nodes, and set P_outflow to the
     !! average of this (it should be uniform along bound)
     integer(ikind) :: i,ispec,j,nsum_in,nsum_out
     real(rkind) :: enthalpy,Rgas_mix_local,psum_in,tmpro,cpispec,dummy_real,psum_out
     
#ifndef isoT     
     !! Evaluate the energy (roE) and pressure.
     !$omp parallel do private(ispec,enthalpy,Rgas_mix_local,tmpro,cpispec)
     do i=1,npfb
          
        !! Initialise roE with K.E. term
        roE(i) = half*ro(i)*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        
        !! Store local density
        tmpro = ro(i) 

        !! Loop over species
        Rgas_mix_local = zero
        do ispec=1,nspec       
           !! Evaluate local species enthalpy
           call evaluate_enthalpy_at_node(T(i),ispec,enthalpy,cpispec,dummy_real)
           
           !! Add species enthalpy contribution
           roE(i) = roE(i) + Yspec(i,ispec)*enthalpy
           
           !! Build the local mixture gas constant - actually holds ro*R
           Rgas_mix_local = Rgas_mix_local + Yspec(i,ispec)*Rgas_universal*one_over_molar_mass(ispec)
        end do           

        !! Evaluate the pressure           
        p(i) = Rgas_mix_local*T(i)

        !! Subtract p=ro*R*T
        roE(i) = roE(i) - p(i)
        
     end do
     !$omp end parallel do
          

     !! Pressure on inflow and outflow nodes 
     psum_in = zero;nsum_in = 0
     psum_out = zero;nsum_out = 0     
     if(nb.ne.0)then
#ifdef restart
        !! For restarts, we don't want to update the target pressure
        P_outflow = p_ref
        P_inflow = p_ref
#else             
        !$omp parallel do private(i) reduction(+:psum_in,psum_out,nsum_in,nsum_out)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.1) then
                          
              !! Augment the accumulators for sum of pressure and # inflow nodes
              psum_in = psum_in + p(i)
              nsum_in = nsum_in + 1
           end if
           if(node_type(i).eq.2) then
                          
              !! Augment the accumulators for sum of pressure and # outflow nodes
              psum_out = psum_out + p(i)
              nsum_out = nsum_out + 1
           end if
        end do
        !$omp end parallel do
        
        !! Find the average for processors with in/outflows. Set to zero otherwise
        if(nsum_in.ne.0) then
           P_inflow = psum_in/dble(nsum_in)
        else
           P_inflow = p_ref
        end if
        if(nsum_out.ne.0) then
           P_outflow = psum_out/dble(nsum_out)
        else
           P_outflow = p_ref
        end if
        
#endif
     end if
#else
     P_outflow = csq*rho_char
     P_inflow = csq*rho_char
#endif

       
     return
  end subroutine initialise_energy  
!! ------------------------------------------------------------------------------------------------  
end module thermodynamics
