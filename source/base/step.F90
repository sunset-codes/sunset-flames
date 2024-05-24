module step
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains time-stepping routines, and a routine to set the time-step.
  !! Time-stepping routines start with "step_" and perform one time step. They are
  !! called from the main loop, and they themselves call routines in the rhs module.

  use kind_parameters
  use common_parameter
  use common_vars
  use mirror_boundaries
  use characteristic_boundaries
  use rhs
  use mpi_transfers
#ifdef mp  
  use mpi
#endif  
  implicit none

  private
  public set_tstep,step_rk3_4S_2R,step_rk3_4S_2R_EE, &
         set_tstep_PID
  
  !! Error norms for RK3(2)4S[2R+]C scheme
  real(rkind) :: enrm_ro,enrm_rou,enrm_rov,enrm_E,enrm_row
  real(rkind),dimension(nspec_max) :: enrm_Yspec
  real(rkind), parameter :: outflow_error_scaling=1.0d0 !! Scaling factor for errors on in/out flows


contains
!! ------------------------------------------------------------------------------------------------  
  subroutine step_rk3_4S_2R
     !! 3rd order 4step 2 register Runge Kutta
     !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
     !! 'U' and 'S' in comments relate to p31 of Senga2 user guide
     !! Implemented over three registers, because speed is more
     !! important than memory at present.    
     
     !! Register 1 is ro_reg1,rou_reg1,rov_reg1 etc
     !! Register 2 is ro,u,v etc
     !! Register 3 is rhs_ro, rhs_rou, rhs_rov etc
     use derivatives
     integer(ikind) :: i,k,ispec
     real(rkind) :: time0
     real(rkind),dimension(:),allocatable :: rou_reg1,rov_reg1,row_reg1,ro_reg1,roE_reg1
     real(rkind),dimension(:,:),allocatable :: Yspec_reg1
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb,RKc
         
     !! Set RKa,RKb with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)
     RKc(:) = dt*rk3_4s_2r_c(:)     

     allocate(rou_reg1(npfb),rov_reg1(npfb),ro_reg1(npfb),roE_reg1(npfb),row_reg1(npfb))
     allocate(rhs_rou(npfb),rhs_rov(npfb),rhs_ro(npfb),rhs_roE(npfb),rhs_row(npfb))
     allocate(Yspec_reg1(npfb,nspec),rhs_Yspec(npfb,nspec))

     !! Store primary variables in register 1 (w-register)
     !$omp parallel do private(ispec)
     do i=1,npfb
        ro_reg1(i)=ro(i);rou_reg1(i)=rou(i);rov_reg1(i)=rov(i);row_reg1(i)=row(i);roE_reg1(i)=roE(i)
        do ispec=1,nspec
           Yspec_reg1(i,ispec) = Yspec(i,ispec)
        end do
     end do
     !$omp end parallel do             

     !! Temporary storage of time
     time0=time
     iRKstep=0

     do k=1,3
        iRKstep = iRKstep + 1
        
        !! Set the intermediate time
        time = time0 + RKc(iRKstep)         
        
        !! Calculate the RHS
        call calc_all_rhs
     
        !! Set w_i and new u,v
        !$omp parallel do private(ispec)
        do i=1,npfb
           !! Store next U in register 2
           ro(i) = ro_reg1(i) + RKa(k)*rhs_ro(i)
           rou(i) = rou_reg1(i) + RKa(k)*rhs_rou(i)
           rov(i) = rov_reg1(i) + RKa(k)*rhs_rov(i)
           row(i) = row_reg1(i) + RKa(k)*rhs_row(i)           
#ifndef isoT           
           roE(i) = roE_reg1(i) + RKa(k)*rhs_roE(i)
#endif
           do ispec=1,nspec
              Yspec(i,ispec) = Yspec_reg1(i,ispec) + RKa(k)*rhs_Yspec(i,ispec)
           end do

           !! Store next S in register 1
           ro_reg1(i) = ro_reg1(i) + RKb(k)*rhs_ro(i)
           rou_reg1(i) = rou_reg1(i) + RKb(k)*rhs_rou(i)
           rov_reg1(i) = rov_reg1(i) + RKb(k)*rhs_rov(i) 
           row_reg1(i) = row_reg1(i) + RKb(k)*rhs_row(i)            
#ifndef isoT
           roE_reg1(i) = roE_reg1(i) + RKb(k)*rhs_roE(i)
#endif
           do ispec=1,nspec
              Yspec_reg1(i,ispec) = Yspec_reg1(i,ispec) + RKb(k)*rhs_Yspec(i,ispec)
           end do
        end do
        !$omp end parallel do
              
        !! Apply BCs and update halos
        if(nb.ne.0) call apply_time_dependent_bounds        
        call reapply_mirror_bcs
        call halo_exchanges_all
        
        !! Get velocity from momentum
        call get_velocity_from_momentum
        
        !! Velocity divergence
#ifdef dim3        
        call calc_divergence(u,v,w,divvel(1:npfb))
#else
        call calc_divergence(u,v,divvel(1:npfb))
#endif        
        call reapply_mirror_bcs_divvel_only
        call halo_exchange_divvel
        
     end do
     
     !! Final substep: returns solution straight to ro,u,v,E (register 2)
     !! and doesn't update S
     iRKstep = iRKstep + 1
     
     !! Set the intermediate time
     time = time0 + RKc(iRKstep)      
     
     !! Build the right hand side
     call calc_all_rhs  
     
     !$omp parallel do private(ispec)
     do i=1,npfb
        !! Final values of prim vars
        ro(i) = ro_reg1(i) + RKb(4)*rhs_ro(i)
        rou(i) = rou_reg1(i) + RKb(4)*rhs_rou(i)
        rov(i) = rov_reg1(i) + RKb(4)*rhs_rov(i)
        row(i) = row_reg1(i) + RKb(4)*rhs_row(i)        
#ifndef isoT
        roE(i) = roE_reg1(i) + RKb(4)*rhs_roE(i)
#endif
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec_reg1(i,ispec) + RKb(4)*rhs_Yspec(i,ispec)
        end do
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(rou_reg1,rov_reg1,row_reg1,ro_reg1,roE_reg1,Yspec_reg1)
     deallocate(rhs_rou,rhs_rov,rhs_row,rhs_ro,rhs_roE,rhs_Yspec)     
     
     !! Apply BCs and update halos
     if(nb.ne.0) call apply_time_dependent_bounds     
     call reapply_mirror_bcs
     call halo_exchanges_all
              
     !! Filter the solution 
     call filter_variables

     !! Apply BCs and update halos
     if(nb.ne.0) call apply_time_dependent_bounds
     call reapply_mirror_bcs
     call halo_exchanges_all
     
     !! Get velocity from momentum
     call get_velocity_from_momentum     

     !! Velocity divergence
#ifdef dim3        
        call calc_divergence(u,v,w,divvel(1:npfb))
#else
        call calc_divergence(u,v,divvel(1:npfb))
#endif        
     call reapply_mirror_bcs_divvel_only
     call halo_exchange_divvel     


     return
  end subroutine step_rk3_4S_2R
!! ------------------------------------------------------------------------------------------------
  subroutine get_velocity_from_momentum
     !! Divide momentum by density to get velocities
     integer(ikind) :: i
     real(rkind) :: ooro
     
     !$omp parallel do private(ooro)
     do i=1,np
        ooro = one/ro(i)
        u(i) = rou(i)*ooro
        v(i) = rov(i)*ooro
        w(i) = row(i)*ooro
     end do
     !$omp end parallel do
     
     return
  end subroutine get_velocity_from_momentum  
!! ------------------------------------------------------------------------------------------------
  subroutine step_rk3_4S_2R_EE
     use derivatives  
     !! 3rd order 4step 2 register Runge Kutta, with embedded 2nd order
     !! scheme for error estimation.     
     !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
     !! 'U' and 'S' in comments relate to p31 of Senga2 user guide
     !! Implemented over three registers, because speed is more
     !! important than memory at present.    
     
     !! Register 1 is ro_reg1,rou_reg1,rov_reg1 etc
     !! Register 2 is ro,u,v etc...
     !! Register 3 is rhs_ro,rhs_rou,rhs_rov (only used for RHS)
     !! Register 4 is e_acc_ro,e_acc_rou,e_acc_rov - error accumulator
     integer(ikind) :: i,k,ispec
     real(rkind) :: time0,emax_Y
     real(rkind),dimension(:),allocatable :: rou_reg1,rov_reg1,row_reg1,ro_reg1,roE_reg1
     real(rkind),dimension(:),allocatable :: e_acc_ro,e_acc_rou,e_acc_rov,e_acc_E,e_acc_row
     real(rkind),dimension(:,:),allocatable :: Yspec_reg1,e_acc_Yspec
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb,RKbmbh,RKc
     
     !! Push the max error storage back one
     emax_nm1 = emax_n;emax_n=emax_np1
     
    
     !! Set RKa,RKb,RKbmbh with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)
     RKbmbh(:) = dt*rk3_4s_2r_bmbh(:)
     RKc(:) = dt*rk3_4s_2r_c(:)

     allocate(rou_reg1(npfb),rov_reg1(npfb),ro_reg1(npfb),roE_reg1(npfb),row_reg1(npfb))
     allocate(rhs_rou(npfb),rhs_rov(npfb),rhs_ro(npfb),rhs_roE(npfb),rhs_row(npfb))
     allocate(e_acc_ro(npfb),e_acc_rou(npfb),e_acc_rov(npfb),e_acc_E(npfb),e_acc_row(npfb))
     allocate(Yspec_reg1(npfb,nspec),rhs_Yspec(npfb,nspec),e_acc_Yspec(npfb,nspec))
     e_acc_ro=zero;e_acc_rou=zero;e_acc_rov=zero;e_acc_E=zero;e_acc_Yspec=zero;e_acc_row=zero
     
     !! Store prim vars in register 1 (w-register)
     !$omp parallel do private(ispec)
     do i=1,npfb
        ro_reg1(i)=ro(i);rou_reg1(i)=rou(i);rov_reg1(i)=rov(i);row_reg1(i)=row(i);roE_reg1(i)=roE(i)
        do ispec=1,nspec
           Yspec_reg1(i,ispec)=Yspec(i,ispec) 
        end do
     end do
     !$omp end parallel do             

     !! Temporary storage of time
     time0=time
     iRKstep = 0
     
     do k=1,3
        iRKstep = iRKstep + 1

        !! Set the intermediate time
        time = time0 + RKc(iRKstep)

        !! Calculate the RHS        
        call calc_all_rhs      

        !! Set w_i and new u,v
        !$omp parallel do private(ispec)
        do i=1,npfb
        
           !! Store next U in register 2
           ro(i) = ro_reg1(i) + RKa(k)*rhs_ro(i)
           rou(i) = rou_reg1(i) + RKa(k)*rhs_rou(i)
           rov(i) = rov_reg1(i) + RKa(k)*rhs_rov(i)
           row(i) = row_reg1(i) + RKa(k)*rhs_row(i)           
#ifndef isoT           
           roE(i) = roE_reg1(i) + RKa(k)*rhs_roE(i)
#endif
           do ispec=1,nspec
              Yspec(i,ispec) = Yspec_reg1(i,ispec) + RKa(k)*rhs_Yspec(i,ispec)
           end do

           !! Store next S in register 1
           ro_reg1(i) = ro_reg1(i) + RKb(k)*rhs_ro(i)
           rou_reg1(i) = rou_reg1(i) + RKb(k)*rhs_rou(i)
           rov_reg1(i) = rov_reg1(i) + RKb(k)*rhs_rov(i) 
           row_reg1(i) = row_reg1(i) + RKb(k)*rhs_row(i)            
#ifndef isoT
           roE_reg1(i) = roE_reg1(i) + RKb(k)*rhs_roE(i)
#endif
           do ispec=1,nspec
              Yspec_reg1(i,ispec) = Yspec_reg1(i,ispec) + RKb(k)*rhs_Yspec(i,ispec)
           end do
           
           !! Error accumulation
           e_acc_ro(i) = e_acc_ro(i) + RKbmbh(k)*rhs_ro(i)       
           e_acc_rou(i) = e_acc_rou(i) + RKbmbh(k)*rhs_rou(i)
           e_acc_rov(i) = e_acc_rov(i) + RKbmbh(k)*rhs_rov(i)  
           e_acc_row(i) = e_acc_row(i) + RKbmbh(k)*rhs_row(i)             
#ifndef isoT
           e_acc_E(i) = e_acc_E(i) + RKbmbh(k)*rhs_roE(i)                    
#endif
           do ispec=1,nspec         
              e_acc_Yspec(i,ispec) = e_acc_Yspec(i,ispec) + RKbmbh(k)*rhs_Yspec(i,ispec)
           end do
        end do
        !$omp end parallel do
       
        !! Apply BCs and update halos
        if(nb.ne.0) call apply_time_dependent_bounds        
        call reapply_mirror_bcs
        call halo_exchanges_all
        
        !! Get velocity from momentum
        call get_velocity_from_momentum           
        
        !! Velocity divergence
#ifdef dim3        
        call calc_divergence(u,v,w,divvel(1:npfb))
#else
        call calc_divergence(u,v,divvel(1:npfb))
#endif        
        call reapply_mirror_bcs_divvel_only
        call halo_exchange_divvel        
        
     end do
     
     !! Final substep: returns solution straight to ro,u,v,E (register 2)
     !! and doesn't update S
     iRKstep = iRKstep + 1

     !! Set the intermediate time
     time = time0 + RKc(iRKstep)     
     
     !! Build the right hand sides
     call calc_all_rhs    
          
     enrm_ro=zero;enrm_rou=zero;enrm_rov=zero;enrm_E=zero;enrm_Yspec=zero;enrm_row=zero
     !$omp parallel do private(ispec) reduction(max:enrm_ro,enrm_rou,enrm_rov,enrm_E,enrm_Yspec,enrm_row)
     do i=1,npfb
     
        !! Final values of conservative variables
        ro(i) = ro_reg1(i) + RKb(iRKstep)*rhs_ro(i)
        rou(i) = rou_reg1(i) + RKb(iRKstep)*rhs_rou(i)
        rov(i) = rov_reg1(i) + RKb(iRKstep)*rhs_rov(i)
        row(i) = row_reg1(i) + RKb(iRKstep)*rhs_row(i)        
#ifndef isoT
        roE(i) = roE_reg1(i) + RKb(iRKstep)*rhs_roE(i)
#endif
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec_reg1(i,ispec) + RKb(iRKstep)*rhs_Yspec(i,ispec)
        end do
        
        !! Final error accumulators
        e_acc_ro(i) = e_acc_ro(i) + RKbmbh(iRKstep)*rhs_ro(i)       
        e_acc_rou(i) = e_acc_rou(i) + RKbmbh(iRKstep)*rhs_rou(i)
        e_acc_rov(i) = e_acc_rov(i) + RKbmbh(iRKstep)*rhs_rov(i) 
        e_acc_row(i) = e_acc_row(i) + RKbmbh(iRKstep)*rhs_row(i)         
#ifndef isoT
        e_acc_E(i) = e_acc_E(i) + RKbmbh(iRKstep)*rhs_roE(i)   
#endif
        do ispec=1,nspec
           e_acc_Yspec(i,ispec) = e_acc_Yspec(i,ispec) + RKbmbh(iRKstep)*rhs_Yspec(i,ispec)
        end do
        
        !! Trick for outflow stability - upscale the errors at the outflow. This is because the 
        !! low-order mixed discretisation at the boundary is less accurate than internally, and
        !! requires a more restrictive time step. For non-uniform resolutions, we usually have
        !! the resolution in the outflow boundary larger than around obstacles in the domain, in which
        !! circumstances this isn't necessary (outflow boundary isn't most restrictive). This is controlled by 
        !! "scale_outflow_errors", which is set in setup_domain.
        if(scale_outflow_errors.eq.1) then
           if(node_type(i).eq.2.or.node_type(i).eq.1) then
              e_acc_ro(i) = e_acc_ro(i)*outflow_error_scaling
              e_acc_rou(i) = e_acc_rou(i)*outflow_error_scaling
              e_acc_rov(i) = e_acc_rov(i)*outflow_error_scaling
              e_acc_row(i) = e_acc_row(i)*outflow_error_scaling
#ifndef isoT      
              e_acc_E(i) = e_acc_E(i)*outflow_error_scaling
#endif           
              e_acc_Yspec(i,:) = e_acc_Yspec(i,:)*outflow_error_scaling
           end if
        end if
        
        !! Calculating L_infinity norm of errors, normalised by eX_norm. 
        enrm_ro = max(enrm_ro,abs(e_acc_ro(i))*ero_norm)      
        enrm_rou = max(enrm_rou,abs(e_acc_rou(i))*erou_norm)
        enrm_rov = max(enrm_rov,abs(e_acc_rov(i))*erou_norm)
        enrm_row = max(enrm_row,abs(e_acc_row(i))*erou_norm)
#ifndef isoT
        enrm_E = max(enrm_E,abs(e_acc_E(i))*eroE_norm)
#endif
        do ispec=1,nspec
           enrm_Yspec(ispec) = max(enrm_Yspec(ispec),abs(e_acc_Yspec(i,ispec))*eroY_norm)
        end do

        !! Uncomment this if we want to see the distribution of time-stepping errors (useful for finding
        !! the least-stable nodes)
!        alpha_out(i) = max(abs(e_acc_ro(i))*ero_norm,max(abs(e_acc_rou(i))*erou_norm, &
!                       max(abs(e_acc_rov(i))*erou_norm,max(abs(e_acc_E(i))*eroE_norm,max(&
!                       abs(e_acc_Yspec(i,1))*eroY_norm,abs(e_acc_Yspec(i,2))*eroY_norm)))))


     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(rou_reg1,rov_reg1,row_reg1,ro_reg1,roE_reg1,Yspec_reg1)
     deallocate(rhs_rou,rhs_rov,rhs_row,rhs_ro,rhs_roE,rhs_Yspec) 
     deallocate(e_acc_ro,e_acc_rou,e_acc_rov,e_acc_E,e_acc_Yspec,e_acc_row)
     
     !! Finalise L_infinity error norms: find max and ensure it's >0     
     emax_Y = maxval(enrm_Yspec(1:nspec))     
       
     emax_np1 = max( &
                    max( &
                        max(enrm_ro,enrm_E), & 
                        max( &
                            max(enrm_rou,enrm_rov),&
                            max(enrm_row,emax_Y) &
                            ) &
                        ), &
                    doublesmall)   
                    
                        
     !! Apply BCs and update halos
     if(nb.ne.0) call apply_time_dependent_bounds     
     call reapply_mirror_bcs
     call halo_exchanges_all
          
     !! Filter the solution 
     call filter_variables

     !! Apply BCs and update halos
     if(nb.ne.0) call apply_time_dependent_bounds
     call reapply_mirror_bcs
     call halo_exchanges_all
     
     !! Get velocity from momentum
     call get_velocity_from_momentum        
     
     !! Velocity divergence
#ifdef dim3        
     call calc_divergence(u,v,w,divvel(1:npfb))
#else
     call calc_divergence(u,v,divvel(1:npfb))
#endif        
     call reapply_mirror_bcs_divvel_only
     call halo_exchange_divvel        

     return
  end subroutine step_rk3_4S_2R_EE  
!! ------------------------------------------------------------------------------------------------
  subroutine set_tstep
     use thermodynamics
     use transport
     integer(ikind) :: i
     real(rkind) :: dt_visc,dt_therm,dt_spec  
     real(rkind) :: c,uplusc
         
     call evaluate_temperature_and_pressure
     call evaluate_transport_properties   
     
     !! Find minimum values for cfl, visc, thermal diff terms
     dt_cfl = 1.0d10;dt_visc = 1.0d10;dt_therm=1.0d10;dt_spec=1.0d10
     cmax = zero;umax = zero
!     !$omp parallel do private(c,uplusc) reduction(min:dt_cfl,dt_visc,dt_therm,dt_spec) &
!     !$omp reduction(max:cmax,umax)
     do i=1,npfb
        !! Sound speed 
#ifndef isoT        
        c = evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i))    
#else
        c = sqrt(csq)
#endif                
 
        !! Max velocity and sound speed
        umax = sqrt(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        cmax = max(c,cmax)
        
        !! Max speed of information propagation
        uplusc = umax + c
        
        !! Acoustic:: s/(u+c)
        !! Slightly reduce on outflows for stability
        if(node_type(i).eq.2.or.node_type(i).eq.1) then
           dt_cfl = min(dt_cfl,0.5d0*s(i)/uplusc)
        else
           dt_cfl = min(dt_cfl,s(i)/uplusc)
        endif

        !! Viscous:: s*s*ro/visc
        dt_visc = min(dt_visc,s(i)*s(i)*ro(i)/visc(i))
        
#ifndef isoT        
        !! Thermal:: s*s*ro*cp/lambda_th
        dt_therm = min(dt_therm,s(i)*s(i)*ro(i)*cp(i)/lambda_th(i))
#endif      

        !! Molecular diffusivity::  s*s*ro/roMdiff
        dt_spec = min(dt_spec,s(i)*s(i)*ro(i)/maxval(roMdiff(i,1:nspec)))
         
     end do
!     !$omp end parallel do

     !! Scale by characteristic lengths and coefficients
#ifdef isoT     
     dt_cfl = 0.4d0*dt_cfl*L_char  !! Isothermal seems to need more restrictive cfl at BCs. needs investigating.
#else
     dt_cfl = one*dt_cfl*L_char
#endif     
     dt_visc = 0.3d0*dt_visc*L_char*L_char
     dt_therm = 0.3d0*dt_therm*L_char*L_char
     dt_spec = two*dt_spec*L_char*L_char
                           
     !! Find most restrictive parabolic constraint
     dt_parabolic = min(dt_visc,min(dt_therm,dt_spec)) 
     
     !! Set dt if not reacting. If reacting, it will be set later by PID.
#ifndef react
     dt = min(dt_parabolic,dt_cfl)
#endif     
     
#ifdef mp     
     !! Global cfl-based time-step and parabolic parts based time-step
     call global_reduce_min(dt_cfl)
     call global_reduce_min(dt_parabolic)     
                 
     !! Output time-step (only if not reacting)
#ifndef react
     !! Global time step
     call global_reduce_min(dt)
     if(iproc.eq.0) then
        write(192,*) time,dt,one
        flush(192)        
     end if
#endif     
#else
#ifndef react
     write(192,*) time,dt,one
     flush(192)    
#endif     
#endif     

     return
  end subroutine set_tstep
!! ------------------------------------------------------------------------------------------------  
  subroutine set_tstep_PID
     !! Adapt the time-step using the PID controller based on intergation errors.
     !! This routine is generally only called for reacting flows, and *presumes* that 
     !! the CFL-type time-step constraints have been calculated already.
     real(rkind) :: dtfactor
     real(rkind) :: facA,facB,facC
     real(rkind) :: dt_max
        
                  
#ifdef mp     
     !! Parallel transfer to obtain the global maximum error          
     call global_reduce_max(emax_np1)
#endif     
   
     !! P, I and D factors..  Calculation done in log space...
     facA = pid_a*log(pid_tol/emax_np1)
     facB = pid_b*log(emax_n/pid_tol)
     facC = pid_c*log(pid_tol/emax_nm1)
      
     !! Combined factor
     dtfactor = pid_k*exp(facA+facB+facC)
          
     !! Suppress big changes in time step (especially increases). N.B. this significantly reduces
     !! the stiffness of the PID system, and seems faster so far.
!     dtfactor = one + one*atan((dtfactor-one)/one)

     !! Alternative approach is impose hard limits on dtfactor
     dtfactor = max(half,dtfactor)
     dtfactor = min(1.02d0,dtfactor)         
       
     !! Set time new step
     dt = dt*dtfactor
     
     !! Impose upper limit based on CFL-type constraints
     dt_max = min(dt_cfl,dt_parabolic)
     if(dt.gt.dt_max) dt = dt_max


#ifdef mp     
     !! Find global time-step
     call global_reduce_min(dt)
     if(iproc.eq.0) then
        write(192,*) time,dt,dt/dt_cfl,emax_np1
        flush(192)
     end if
#else
     write(192,*) time,dt,dt/dt_cfl,emax_np1
     flush(192)         
#endif 
 

     return
  end subroutine set_tstep_PID
!! ------------------------------------------------------------------------------------------------  
end module step
