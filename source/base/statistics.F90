module statistics
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2022 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to perform reduction calculations to obtain and output global 
  !! statistics, and write data to standard-out.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  use neighbours
#ifdef mp
  use mpi
  use mpi_transfers
#endif    
  implicit none

contains
!! ------------------------------------------------------------------------------------------------
  subroutine open_stats_files
     !! Opens units for statistic outputting files

     !! Main time out-file, contains time, dt, etc.
     open(unit=21,file='./data_out/time.out')
     
     !! Total cpu time per step
     open(unit=191,file='data_out/statistics/cputime.out')
     
     !! Time, value of time-step, value of time-step divided by CFL-based timestep (if using PID),
     !! and maximum sound speed in domain
     open(unit=192,file='data_out/statistics/dt.out')
     
     !! Time, total mass within the domain, total volume
     open(unit=193,file='data_out/statistics/masscheck.out')
     
     !! Time, net forces on wall boundaries. Needs modifying on case-by-case basis
     open(unit=194,file='data_out/statistics/liftdrag.out')
     
     !! Time, mean (or L2) velocity, and components thereof
     open(unit=195,file='data_out/statistics/velcheck.out')
     
     !! Time, L2 error of velocity field relative to analytic solution for Taylor-Green flow  
     open(unit=196,file='data_out/statistics/l2.out')
     
     !! Time, total energy within domain.
     open(unit=197,file='data_out/statistics/energy_sum.out')  
     
     !! Time, enstrophy
     open(unit=198,file='data_out/statistics/enstrophy.out')

     !! Time, L2 error in mass fraction sum, total mass of each species in domain
     open(unit=199,file='data_out/statistics/species.out')
     
     !! Time, heat release rate
     open(unit=200,file='data_out/statistics/heat_release_rate.out')
     
     !! Time, volume integrated kinetic energy
     open(unit=201,file='data_out/statistics/ke.out')
     
     !! Time, mean pressure, mean outflow pressure, mean inflow pressure
     open(unit=202,file='data_out/statistics/p.out')
     
     ! Time, mean temperature, mean boundary temperature
     open(unit=203,file='data_out/statistics/T.out')
     
     ! Time, furthest upstream flame location
     open(unit=204,file='data_out/statistics/flpos.out')
     
     return     
  end subroutine open_stats_files
!! ------------------------------------------------------------------------------------------------
  subroutine statistics_control(m_out_stats)
     !! This routine controls the statistics calculation and output routines.
     integer(ikind),intent(inout) :: m_out_stats
     integer(ikind) :: i
     integer(ikind),parameter :: istats_freq = 10
     
     
     !! Evaluate the mean velocity and adjust pressure gradient if required
     !! This should be done every step if using PID for pressure gradient
#ifdef pgrad     
     call velocity_control
#endif  
     
     !! All
     if(itime.eq.0.or.time.gt.m_out_stats*dt_out_stats) then
        m_out_stats = m_out_stats+1

#ifndef pgrad
        !! Velocity control
        call velocity_control     
#endif

        !! Check conservation of mass and energy
        call mass_and_energy_check
        
        !! Check pressure
        call pressure_check
        
        !! Check temperature
        call temperature_check
        
        !! Flame positions
        call flame_position
     
        !! Error evaluation for Taylor Green vortices?
!        call error_TG
     
        !! Error evaluation for Poiseuille flow
        call poiseuille_l2norm

        !! Calculate the lift and drag on all solid obstacles
!        call liftdrag

        !! Calculate how well balanced the MPI decomposition is
!        call check_load_balance

        !! Check the conservation of the species equations
        call species_check
        
        !! Check heat release
        call heat_release_rate
     endif

     return
  end subroutine statistics_control  
!! ------------------------------------------------------------------------------------------------   
  subroutine check_load_balance  
     !! N.B. This subroutine is out-of-date with the rest of the code, and needs updating if you
     !! wish to use it.
     real(rkind),dimension(:),allocatable :: prof_tmp,prof_tmp_local,load_n_n
     integer(ikind),dimension(:),allocatable :: sum_n_n,sum_n_n_local
     integer(ikind) :: sum_sum_n_n
#ifdef mp     
     !! Calculate relative amounts of work being done by each processor
     allocate(prof_tmp_local(nprocs),prof_tmp(nprocs));prof_tmp_local=zero
     prof_tmp_local(iproc+1) = sum(segment_time_local(1:8))
     call MPI_ALLREDUCE(prof_tmp_local,prof_tmp,nprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     
     !! Sum of the total neighbours of all nodes on each processor
     allocate(sum_n_n(nprocs),sum_n_n_local(nprocs));sum_n_n_local=0
     sum_n_n_local(iproc+1) = sum(ij_count(1:npfb))
     call MPI_ALLREDUCE(sum_n_n_local,sum_n_n,nprocs,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)
     sum_sum_n_n = sum(sum_n_n(:))
     allocate(load_n_n(nprocs))
     load_n_n = sum_n_n/dble(npfb)!dble(nprocs)*dble(sum_n_n(:))/dble(sum_sum_n_n)
     
     
     !! Output the amounts of work being done
     if(iproc.eq.0) write(6,*) dble(nprocs)*prof_tmp(:)/sum(prof_tmp(:))

     !! Output the expected load (nodes*neighbours)
     if(iproc.eq.0) write(6,*) load_n_n       
     
     deallocate(prof_tmp_local,prof_tmp)
     deallocate(sum_n_n,sum_n_n_local,load_n_n)
#endif
  
     return
  end subroutine check_load_balance   
!! ------------------------------------------------------------------------------------------------
  subroutine liftdrag
     !! TO DO: Update for 3 dimensional simulations  
     !! TO DO: Update for non-isothermal simulations
     use derivatives
#ifdef isoT     
     integer(ikind) :: i,j
     real(rkind),dimension(dims) :: gradu0,gradv0,Fn,force,force_tmp
     real(rkind),dimension(dims,dims) :: Jinv,sigma
     real(rkind) :: xn,yn
   
  
     !! Calculate the velocity gradient 
     allocate(gradu(npfb,dims),gradv(npfb,dims),gradw(npfb,dims))
     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)
     call calc_gradient(w,gradw)     

     force = zero
     !$omp parallel do private(i,Jinv,gradu0,gradv0,xn,yn,sigma,Fn) reduction(+:force)
     do j=1,nb
        i=boundary_list(j)
        if(node_type(i).eq.0)then
           xn = rnorm(i,1);yn=rnorm(i,2)
           Jinv(1,1)=xn;Jinv(1,2)=-yn;Jinv(2,1)=yn;Jinv(2,2)=xn   !! Jacobian for normal-tangent to x-y  
           Jinv(3,:)=zero;Jinv(:,3)=zero;Jinv(3,3)=one         
           gradu0(:) = matmul(Jinv,gradu(i,:))  !! Velocity gradients in x-y FoR
           gradv0(:) = matmul(Jinv,gradv(i,:))           
           
           !! Total stress on surface
           sigma(1,1) = visc(i)*(fourthirds*gradu0(1)-twothirds*gradv0(2)) - csq*(ro(i)-one)
           sigma(1,2) = visc(i)*(gradu0(2)+gradv0(1))
           sigma(2,1) = visc(i)*(gradu0(2)+gradv0(1))           
           sigma(2,2) = visc(i)*(fourthirds*gradv0(2)-twothirds*gradu0(1)) - csq*(ro(i)-one)
          
           Fn(:) = matmul(sigma,rnorm(i,:))   !! Force on surface (sigma.n)
                     
           force(:) = force(:) + Fn(:)*s(i)*L_char  !! Integrate over surface... s(i) is dimensionless node spacing
        end if
     end do
     !$omp end parallel do
     deallocate(gradu,gradv,gradw)

#ifdef mp
     force_tmp = force
     call MPI_ALLREDUCE(force_tmp,force,dims,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)     
     if(iproc.eq.0)then
        write(194,*) time/Time_char,force
        flush(194)
     end if        
#else
     write(194,*) time/Time_char,force
     flush(194)
#endif     
     
#endif     
     return
  end subroutine liftdrag
!! ------------------------------------------------------------------------------------------------  
  subroutine velocity_control
     !! Output the L2 of velocity over the domain
     integer(ikind) :: i
     real(rkind) :: tot_vel,tot_vol,tmpro,dVi,tmpvel
     real(rkind),dimension(dims) :: tot_u
     real(rkind) :: facA,facB,facC,facT,deflowdt
       
     tot_vel = zero
     tot_vol = zero
     tot_u = zero
     !$omp parallel do private(tmpro,tmpvel,dVi) reduction(+:tot_vel,tot_vol,tot_u)
     do i=1,npfb
        tmpro = ro(i)
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif               
        tmpvel = (u(i)*u(i)+v(i)*v(i)+w(i)*w(i))
        tot_vel = tot_vel + tmpvel*dVi
        tot_vol = tot_vol + dVi
        tot_u = tot_u + (/u(i),v(i),w(i)/)*dVi
     end do
     !$omp end parallel do

#ifdef mp
     !! Reduce
     call global_reduce_sum(tot_vel)
     call global_reduce_sum(tot_vol)
     call global_reduce_sum(tot_u(1))
     call global_reduce_sum(tot_u(2))
     call global_reduce_sum(tot_u(3))                    
#endif         
    
     !! Normalise over volume
     tot_vel = sqrt(tot_vel/tot_vol)
     tot_u(:) = tot_u(:)/tot_vol

     !! If we want to P.I.D. control over the mean velocity
#ifdef pgrad     
     !! New error     
     eflow_n = u_char - tot_u(1)!tot_vel
          
     !! Integral term
     sum_eflow = sum_eflow + eflow_n*dt
     
     !! Derivative term
     deflowdt = (eflow_n-eflow_nm1)/dt
    
     !! P, I and D factors..  
     facA = two*one
     facB = facA/1.0d-1
     facC = facA*0.02d0
         
     driving_force(1) = facA*eflow_n + facB*sum_eflow + facC*deflowdt
     !! Impose some upper and lower limits
     driving_force(1) = min(5.0d0,driving_force(1))
     driving_force(1) = max(-5.0d0,driving_force(1))
                       
     !! Pass new eflow to old eflow
     eflow_nm1=eflow_n
#else
     driving_force = zero
#endif
      
#ifdef mp
     if(iproc.eq.0)then
        write(195,*) time/Time_char,tot_vel,tot_u,driving_force(1)
        flush(195)
     end if
#else
     write(195,*) time/Time_char,tot_vel,tot_u,driving_force(1)    
     flush(195)
#endif       
     
     
     return
  end subroutine velocity_control   
!! ------------------------------------------------------------------------------------------------  
  subroutine heat_release_rate
     !! This subroutine calculates the total heat release rate integrated over the domain
     integer(ikind) :: i
     real(rkind) :: tot_hrr,dVi
#ifdef react    
   
     tot_hrr = zero
     !$omp parallel do private(dVi) reduction(+:tot_hrr)
     do i=1,npfb
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif        
        tot_hrr = tot_hrr + hrr(i)*dVi
     end do
     !$omp end parallel do
     
#ifdef mp
     call global_reduce_sum(tot_hrr)
#endif     

     !! Scale by length-scale squared or cubed
#ifdef dim3
     tot_hrr = tot_hrr*L_char**3.0d0
#else
     tot_hrr = tot_hrr*L_char**two
#endif     
        
#ifdef mp
     if(iproc.eq.0)then
        write(200,*) time/Time_char,tot_hrr
        flush(200)
     end if
#else
     write(200,*) time/Time_char,tot_hrr
     flush(200)
#endif

#endif

     return
  end subroutine heat_release_rate 
!! ------------------------------------------------------------------------------------------------
  subroutine mass_and_energy_check
     !! This subroutine calculates the total mass and total energy in the domain.
     integer(ikind) :: i
     real(rkind) :: tot_mass,tot_vol,tmpro,dVi,tot_roE
     real(rkind) :: tot_ke
    
   
     tot_mass = zero;tot_vol = zero;tot_roE = zero;tot_ke=zero
     !$omp parallel do private(tmpro,dVi) reduction(+:tot_mass,tot_vol,tot_roE,tot_ke)
     do i=1,npfb
        tmpro = ro(i)
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif        
        tot_mass = tot_mass + tmpro*dVi
        tot_vol = tot_vol + dVi
#ifndef isoT
        tot_roE = tot_roE + dVi*roE(i)
#endif        
        tot_ke = tot_ke + half*tmpro*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))*dVi
     end do
     !$omp end parallel do
     
     !! Scale to make dimensional
#ifdef dim3
     tot_mass = tot_mass*L_char**three
     tot_vol = tot_vol*L_char**three
     tot_roE = tot_roE*L_char**three         
     tot_ke = tot_ke*L_char**three
#else
     tot_mass = tot_mass*L_char**two
     tot_vol = tot_vol*L_char**two
     tot_roE = tot_roE*L_char**two
     tot_ke = tot_ke*L_char**two          
#endif     
     
#ifdef mp
     call global_reduce_sum(tot_mass)
     call global_reduce_sum(tot_vol)
     call global_reduce_sum(tot_roE)  
     call global_reduce_sum(tot_ke)        
     if(iproc.eq.0)then
        write(193,*) time/Time_char,tot_mass,tot_vol
        flush(193)
        write(197,*) time/Time_char,tot_roE
        flush(197)
        write(201,*) time/Time_char,tot_ke
        flush(201)
        
     end if
#else
     write(193,*) time/Time_char,tot_mass,tot_vol
     flush(193)
     write(197,*) time/Time_char,tot_roE
     flush(197)
     write(201,*) time/Time_char,tot_ke
     flush(201)
#endif

     return
  end subroutine mass_and_energy_check 
!! ------------------------------------------------------------------------------------------------  
  subroutine pressure_check
     !! This subroutine calculates the mean pressure in the domain, and on the inflow/outflows
     integer(ikind) :: i,j
     real(rkind) :: tot_vol,tot_p,tot_pout,tot_pin,dVi,tot_volout,tot_volin
   
   
     tot_vol = zero;tot_p = zero
     !$omp parallel do private(dVi) reduction(+:tot_vol,tot_p)
     do i=1,npfb
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif        
        tot_vol = tot_vol + dVi
        tot_p = tot_p + p(i)*dVi
     end do
     !$omp end parallel do
     tot_p = tot_p - p_ref*tot_vol !! Remove p_ref

     !! Reduce across processors     
#ifdef mp
     call global_reduce_sum(tot_p)
     call global_reduce_sum(tot_vol)     
#endif     
     !! Mean pressure
     tot_p = tot_p/tot_vol

     tot_volout = zero;tot_volin=zero
     tot_pout=zero;tot_pin=zero
     if(nb.ne.0)then
        !$omp parallel do private(i,dVi) reduction(+:tot_volin,tot_volout,tot_pin,tot_pout)
        do j=1,nb
           i=boundary_list(j)        
           dVi = s(i)
        
           if(node_type(i).eq.1) then  !! Inflow
              tot_volin = tot_volin + dVi
              tot_pin = tot_pin + p(i)*dVi
           else if(node_type(i).eq.2) then !! Outflow
              tot_volout = tot_volout + dVi
              tot_pout = tot_pout + p(i)*dVi
           end if
        
        end do
        !$omp end parallel do
        tot_pin = tot_pin - p_ref*tot_volin !! Remove p_ref
        tot_pout = tot_pout - p_ref*tot_volout !! Remove p_ref
     end if

     !! Reduce across processors     
#ifdef mp
     call global_reduce_sum(tot_pin)
     call global_reduce_sum(tot_volin)     
     call global_reduce_sum(tot_pout)
     call global_reduce_sum(tot_volout)          
#endif     
     !! Mean pressures
     tot_pout = tot_pout/tot_volout
     tot_pin = tot_pin/tot_volin    
     
        
#ifdef mp
     if(iproc.eq.0)then
        write(202,*) time/Time_char,tot_p,tot_pin,tot_pout
        flush(202)
     end if
#else
     write(202,*) time/Time_char,tot_p,tot_pin,tot_pout
     flush(202)
#endif


     return
  end subroutine pressure_check  
!! ------------------------------------------------------------------------------------------------  
  subroutine temperature_check
     !! This subroutine calculates the mean temperature in the domain, and the mean temperature on 
     !! solid boundaries (i.e. on a bluff body)
     integer(ikind) :: i,j
     real(rkind) :: tot_vol,tot_t,tot_tb,dVi
   
     !! Global mean temperature
     tot_vol = zero;tot_t = zero
     !$omp parallel do private(dVi) reduction(+:tot_vol,tot_t)
     do i=1,npfb
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif        
        tot_vol = tot_vol + dVi
        tot_t = tot_t + T(i)*dVi
     end do
     !$omp end parallel do

     !! Reduce across processors     
#ifdef mp
     call global_reduce_sum(tot_t)
     call global_reduce_sum(tot_vol)     
#endif     
     !! Mean temperature
     tot_t = tot_t/tot_vol     
   
     !! Mean on solid boundaries
     tot_vol = zero
     tot_tb=zero
     if(nb.ne.0)then
        !$omp parallel do private(i,dVi) reduction(+:tot_tb,tot_vol)
        do j=1,nb
           i=boundary_list(j)        
           dVi = s(i)
#ifdef dim3
           dVi = dVi*dz
#endif        
           if(node_type(i).eq.0) then  !! Wall
              tot_vol = tot_vol + dVi
              tot_tb = tot_tb + T(i)*dVi
           end if
        
        end do
        !$omp end parallel do
     end if

     !! Reduce across processors     
#ifdef mp
     call global_reduce_sum(tot_tb)
     call global_reduce_sum(tot_vol)     
#endif     
     !! Mean pressures
     tot_tb = tot_tb/tot_vol
     
        
#ifdef mp
     if(iproc.eq.0)then
        write(203,*) time/Time_char,tot_t,tot_tb
        flush(203)
     end if
#else
     write(203,*) time/Time_char,tot_t,tot_tb
     flush(203)
#endif


     return
  end subroutine temperature_check    
!! ------------------------------------------------------------------------------------------------  
  subroutine flame_position
     !! This subroutine calculates left-most position at which Y(1)<0.5Y(1,inflow)
     integer(ikind) :: i,j
     real(rkind) :: leftmost_flame,c
   
     !! Global mean temperature
     leftmost_flame=1.0d10
     do i=1,npfb
        c = one - Yspec(i,1)/(ro(i)*Yspec_reactants(1))
        if(c.ge.half.and.rp(i,1).le.leftmost_flame)then
           leftmost_flame = rp(i,1)
        endif                     
     end do

#ifdef mp
     call global_reduce_min(leftmost_flame)
     if(iproc.eq.0)then
        write(204,*) time/Time_char,leftmost_flame
        flush(204)
     end if
#else
     write(204,*) time/Time_char,leftmost_flame
     flush(204)
#endif

     return
  end subroutine flame_position 
!! ------------------------------------------------------------------------------------------------   
  subroutine species_check
     !! This subroutine calculates total mass of each species in the domain
     integer(ikind) :: i,ispec
     real(rkind) :: dVi,tmpY,sumY,tot_error,tmpro,tot_vol
     real(rkind),dimension(:),allocatable :: tot_Yspec
     
     !! Allocate and zero accumulators
     allocate(tot_Yspec(nspec))   
     tot_Yspec = zero     
     tot_error = zero
     tot_vol=zero
     !$omp parallel do private(tmpro,dVi,ispec,tmpY,sumY) reduction(+:tot_Yspec,tot_error,tot_vol)
     do i=1,npfb
        !! Extract one/ro
        tmpro = one/ro(i)
     
        !! size of volume element
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif        

        !! Sum over all species
        sumY = -one
        do ispec=1,nspec
           
           tmpY = Yspec(i,ispec)
           
           !! Total mass of species in domain - integral of roY dV
           tot_Yspec(ispec) = tot_Yspec(ispec) + dVi*tmpY
           
           !! 1- sum_over_species(Y)
           sumY = sumY + tmpY*tmpro
        end do
        
        !! Volume integral
        tot_vol = tot_vol + dVi
        
        !! L2 error in sumY
        tot_error = tot_error + sumY*sumY*dVi
     end do
     !$omp end parallel do
     
     !! Rescale integral quantities to be dimensional
#ifdef dim3
     tot_error = tot_error*L_char**3.0d0
     tot_Yspec = tot_Yspec*L_char**3.0d0
#else
     tot_error = tot_error*L_char**two
     tot_Yspec = tot_Yspec*L_char**two
#endif     
     
#ifdef mp
     do ispec=1,nspec
        call global_reduce_sum(tot_Yspec(ispec))
     end do
     call global_reduce_sum(tot_error)
     call global_reduce_sum(tot_vol)     
     
     !! L2 norm for error           
     tot_error = sqrt(tot_error/tot_vol)
          
     if(iproc.eq.0)then
        write(199,*) time/Time_char,tot_error,tot_Yspec(:)
        flush(199)
     end if
#else
     !! Normalise error by volume. Species masses are not normalised
     write(199,*) time/Time_char,sqrt(tot_error/tot_vol),tot_Yspec(:)
     flush(199)
#endif
     
     deallocate(tot_Yspec)

     return
  end subroutine species_check    
!! ------------------------------------------------------------------------------------------------    
  subroutine poiseuille_l2norm
     !! Output the L2norm of the velocity compared with Poiseuille flow
     !! Umax=1, unit width domain...
     !! N.B. This routine needs updating for multi-processor simulations
     integer(ikind) :: i,j,jj
     real(rkind) :: y,uexact,local_error,sum_e,sum_exact,dVi
     real(rkind) :: N,X1,X2,y1
       
     sum_e = zero
     sum_exact = zero
     !$omp parallel do private(y,uexact,local_error,dVi,j,jj,N,X1,X2,y1) &
     !$omp reduction(+:sum_e,sum_exact)
     do i=1,npfb
        y = rp(i,2)
        uexact = (half-y)*(half+y)*4.0d0
        y1 = y+half


        jj = 50!floor(10.0d0/sqrt(time+0.1d0)); % vary the degree of expansion to avoid NaNs...
        do j=1,jj
           N=(2*j-1)*pi
           X1 = sin(N*y1)/N**3.0d0
           X2 = exp(-N*N*visc(i)*time/Time_char)
           uexact = uexact - 32.0d0*X1*X2
        end do       
        
        dVi = vol(i)
        local_error = u(i)-uexact
        sum_e = sum_e + dVi*local_error**two
        sum_exact = sum_exact + dVi*uexact**two
     end do
     !$omp end parallel do
     
#ifdef mp
     call global_reduce_sum(sum_e)
     call global_reduce_sum(sum_exact)
              
     if(iproc.eq.0) then
        write(196,*) time/Time_char,sqrt(sum_e/sum_exact)         
        flush(196)
     endif
#else     
     write(196,*) time/Time_char,sqrt(sum_e/sum_exact)
     flush(196)  
#endif
     return     
  
     
  end subroutine poiseuille_l2norm  
!! ------------------------------------------------------------------------------------------------
  subroutine error_TG
    !! Error relative to analytic solution for 2D Taylor-Green vortex decay
    implicit none
    integer(ikind) :: i
    real(rkind) :: u_exact,v_exact,x,y,expo
    real(rkind) :: U_ex,U_num,error_sum,L2error,Ex_sum
  
    error_sum = zero;Ex_sum =zero
    !$omp parallel do private(x,y,u_exact,v_exact,U_ex,U_num) reduction(+:Ex_sum,error_sum)
    do i=1,npfb
       x=rp(i,1)
       y=rp(i,2)
       expo = exp(-8.0d0*pi*pi*time/time_char/200.0d0)   !! Hard-coded for Re=200. Need to modify.
       u_exact = -expo*cos(2.0*pi*x)*sin(2.0*pi*y)
       v_exact = expo*sin(2.0*pi*x)*cos(2.0*pi*y)
       U_ex = dsqrt(u_exact**2. + v_exact**2.)
       U_num = dsqrt(u(i)**2. + v(i)**2. + w(i)**2.)
     
       error_sum = error_sum + (U_num-U_ex)**2.
       Ex_sum = Ex_sum + (U_ex)**2.
    end do
    !$omp end parallel do
#ifdef mp
    call global_reduce_sum(error_sum)
    call global_reduce_sum(ex_sum)    
    
    L2error = dsqrt(error_sum/Ex_sum)
    if(iproc.eq.0)then
!       write(6,*) time,L2error,expo,maxval(u(1:npfb))
       write(196,*) time,L2error
       flush(196) 
    end if
#else   
    L2error = dsqrt(error_sum/Ex_sum)
!    write(6,*) time,L2error,expo,maxval(u(1:npfb))
    write(196,*) time/Time_char,L2error
    flush(196) 
#endif    
  
  
    return
  end subroutine error_TG   
!! ------------------------------------------------------------------------------------------------    
  subroutine check_enstrophy
     !! Evaluate volume averaged enstrophy and also the volume averaged kinetic energy
     !! This routine is designed specifically for 3D Taylor-Green vortex tests.
     integer(ikind) :: i
     real(rkind) :: srtnorm,sum_enstrophy,dVi,sum_vol,sum_ke
     
     sum_enstrophy = zero
     sum_vol = zero
     sum_ke = zero
     !$omp parallel do private(srtnorm,dVi) reduction(+:sum_enstrophy,sum_vol,sum_ke)
     do i=1,npfb
     
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif        
        srtnorm = dot_product(gradu(i,:),gradu(i,:)) + &
                  dot_product(gradv(i,:),gradv(i,:)) + &
                  dot_product(gradw(i,:),gradw(i,:))
        sum_enstrophy = sum_enstrophy + visc(i)*srtnorm*dVi
        sum_vol = sum_vol + dVi
        
        !!
        sum_ke = sum_ke + ro(i)*(u(i)**two + v(i)**two + w(i)**two)*dVi
     end do
     !$omp end parallel do
     
#ifdef mp
     call global_reduce_sum(sum_enstrophy)
     call global_reduce_sum(sum_ke)
     call global_reduce_sum(sum_vol)          
     if(iproc.eq.0)then
        write(198,*) time/Time_char,sum_enstrophy/sum_vol,sum_ke/sum_vol
        flush(198)
     end if        
#else
     write(198,*) time/Time_char,sum_enstrophy/sum_vol,sum_ke/sum_vol
     flush(198)
#endif     
           
     return
  end subroutine check_enstrophy  
!! ------------------------------------------------------------------------------------------------ 
end module statistics
