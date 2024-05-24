module turbulence
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |May 2023         |New module to generate initial turbulent velocity field
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to generate initial turbulent velocity fields and turbulent
  !! inflows
  !!
  !! N.B. at present this is only coded up for ISOTROPIC turbulence. (i.e. only the trace of
  !! the Reynolds stress tensor.

  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  implicit none

  real(rkind),dimension(:),allocatable :: u_noturb,v_noturb,w_noturb
  real(rkind),dimension(:),allocatable :: psix,psiy,psiz
  real(rkind) :: L_turb,u_turb
  integer(ikind) :: randseed


contains
!! ------------------------------------------------------------------------------------------------
  subroutine make_turbulent_velocity_field
     !! This is the top-level routine for generating turbulent velocity field, which follows the
     !! diffusion method of Kempf, Klein & Janicka (2005) Flow, Turbulence & Combustion 74:67-84    
     integer(ikind) :: i
     
     !! Copy input parameters to main module variables
     L_turb = turb_lengthscale
     u_turb = turb_intensity
  
  
     !! Store the existing velocity field and zero u,v,w
     allocate(u_noturb(npfb),v_noturb(npfb),w_noturb(npfb))
     !$omp parallel do
     do i=1,npfb
        u_noturb(i) = u(i);u(i) = zero
        v_noturb(i) = v(i);v(i) = zero
        w_noturb(i) = w(i);w(i) = zero    
     end do
     !$omp end parallel do
     
     !! Initialise random number gen
     randseed = 0
#ifdef mp
     randseed = randseed + iproc
#endif     
     call srand(randseed)

  
     !! Apply some diffusion until we get to the required lengthscales
     call noise_and_diffuse_psi
!     call noise_and_diffuse_velocity_field  !! This does just a tiny amount 
     
     !! Scale the turbulent velocity field to give the desired turbulence statistics
     call scale_velocity_field
  
     !! Finalise by adding back the original velocity field
     !$omp parallel do
     do i=1,npfb
        u(i) = u(i) + u_noturb(i)
        v(i) = v(i) + v_noturb(i)
#ifdef dim3
        w(i) = w(i) + w_noturb(i)                
#else
        w(i) = w_noturb(i)
#endif        
     end do
     !$omp end parallel do
     deallocate(u_noturb)
  
     return
  end subroutine make_turbulent_velocity_field
!! ------------------------------------------------------------------------------------------------  
  subroutine scale_velocity_field
     use mpi_transfers
     !! Scale the velocity field to give the desired turbulence characteristics
     integer(ikind) :: i
     real(rkind) :: u_mean,v_mean,w_mean,tot_vol,dVi
     real(rkind) :: u_var
     
     !! Evaluate the mean velocity
     u_mean = zero;v_mean = zero;w_mean = zero;tot_vol = zero
     !$omp parallel do private(dVi) reduction(+:u_mean,v_mean,w_mean,tot_vol)
     do i=1,npfb
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif            
        u_mean = u_mean + u(i)*dVi
        v_mean = v_mean + v(i)*dVi
        w_mean = w_mean + w(i)*dVi
        tot_vol = tot_vol + dVi
     end do
     !$omp end parallel do
#ifdef mp     
     call global_reduce_sum(u_mean)
     call global_reduce_sum(v_mean)
     call global_reduce_sum(w_mean)
     call global_reduce_sum(tot_vol)          
#endif      
     u_mean = u_mean/tot_vol
     v_mean = v_mean/tot_vol
     w_mean = w_mean/tot_vol
     if(iproc.eq.0) write(6,*) "Mean velocities:",u_mean,v_mean,w_mean
     
     !! Remove the mean
     !$omp parallel do 
     do i=1,npfb
        u(i) = u(i) - u_mean
        v(i) = v(i) - v_mean
        w(i) = w(i) - w_mean                
     end do
     !$omp end parallel do
     
     !! Evaluate the variance
     u_var = zero
     !$omp parallel do private(dVi) reduction(+:u_var)
     do i=1,npfb
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz        
#endif        
        u_var = u_var + (u(i)*u(i) + v(i)*v(i) + w(i)*w(i))*dVi
     end do
     !$omp end parallel do
#ifdef mp     
     call global_reduce_sum(u_var)
#endif      
     u_var = u_var/tot_vol
  
     !! Re-scale to give unity variance
     !$omp parallel do
     do i=1,npfb
        u(i) = u_turb*u(i)/sqrt(u_var)
        v(i) = u_turb*v(i)/sqrt(u_var)
        w(i) = u_turb*w(i)/sqrt(u_var)
     end do
     !$omp end parallel do     

     if(iproc.eq.0) write(6,*) "Variance:",u_var     
  
  
     return
  end subroutine scale_velocity_field
!! ------------------------------------------------------------------------------------------------        
  subroutine noise_and_diffuse_psi
     use derivatives
     use mpi_transfers
     !! Applies diffusion to a velocity field
     integer(ikind) :: it_diff,nt_diff,i,j
     real(rkind) :: dt_diff,tot_vol,x_mean,y_mean,z_mean,dVi
     real(rkind),dimension(:),allocatable :: lapx,lapy,lapz
     real(rkind),dimension(:,:),allocatable :: gradx,grady,gradz
     real(rkind) :: randscaling,rnum,oosqrtv
       
     !! Evaluate time-step and number of iterations    
     dt_diff = 0.1d0*(min(smin_global,dz)*L_char)**two
     nt_diff = 1 + floor(L_turb*L_turb/(two*pi*dt_diff))
      
     !! Space for RHS
     allocate(lapx(npfb),lapy(npfb),lapz(npfb))
     allocate(psix(np),psiy(np),psiz(np))
     psix=zero;psiy=zero;psiz=zero
     
     !! Loop over time-steps
     do it_diff = 1,nt_diff+10
                                 
        !! Evaluate Laplacian
#ifdef dim3
        call calc_laplacian(psix,lapx)
        call calc_laplacian(psiy,lapy)
#endif        
        call calc_laplacian(psiz,lapz)                

        !! Set the amplitude of the random noise                 
        if(it_diff.le.nt_diff-1) then
           randscaling = (dble(nt_diff-it_diff)/dble(nt_diff))**1.0d0
        else
           randscaling = zero
        endif
                 
                             
        !! Update - 1st order time-stepping is fine
        x_mean=zero;y_mean=zero;z_mean=zero;tot_vol=zero
        !$omp parallel do private(dVi,rnum,oosqrtv) reduction(+:x_mean,y_mean,z_mean,tot_vol)
        do i=1,npfb
           !! Node volume           
           dVi = vol(i)
#ifdef dim3
           dVi = dVi*dz
#endif            
           oosqrtv = one/sqrt(dVi)
                   
           !! Update
#ifdef dim3
           rnum = rand();psix(i) = psix(i) + dt_diff*lapx(i) + randscaling*(two*rnum-one)*oosqrtv
           rnum = rand();psiy(i) = psiy(i) + dt_diff*lapy(i) + randscaling*(two*rnum-one)*oosqrtv
#endif
           rnum = rand();psiz(i) = psiz(i) + dt_diff*lapz(i) + randscaling*(two*rnum-one)*oosqrtv           

           !! Calculate means
#ifdef dim3           
           x_mean = x_mean + psix(i)*dVi
           y_mean = y_mean + psiy(i)*dVi
#endif           
           z_mean = z_mean + psiz(i)*dVi                                      
           tot_vol = tot_vol + dVi                    

        end do
        !$omp end parallel do            
        
        !! Reduce to find means        
#ifdef mp     
#ifdef dim3
        call global_reduce_sum(x_mean)
        call global_reduce_sum(y_mean)
#endif
        call global_reduce_sum(z_mean)
        call global_reduce_sum(tot_vol)          
#endif      
        x_mean = x_mean/tot_vol
        y_mean = y_mean/tot_vol
        z_mean = z_mean/tot_vol
        
        !! Remove mean from psi
        !$omp parallel do
        do i=1,npfb
#ifdef dim3
           psix(i) = psix(i) - x_mean
           psiy(i) = psiy(i) - y_mean
#endif           
           psiz(i) = psiz(i) - z_mean
        end do
        !$omp end parallel do        
        
        !! Impose boundary conditions
        !$omp parallel do private(i)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! Walls
              psix(i) = zero;psiy(i) = zero;psiz(i) = zero           
           else if(node_type(i).eq.1) then !! Inflow
              psix(i) = zero;psiy(i) = zero;psiz(i) = zero           
           else if(node_type(i).eq.2) then !! Outflow
              psix(i) = zero;psiy(i) = zero;psiz(i) = zero           
           end if

        end do
        !$omp end parallel do            
        
       
        !! Update mirrors and halos
        !$omp parallel do private(i)
        do j=npfb+1,np_nohalo
           i=irelation(j)
#ifdef dim3
           psix(j) = psix(i)
           psiy(j) = psiy(i)
#endif           
           psiz(j) = psiz(i)
        end do        
        !$omp end parallel do
#ifdef mp
#ifdef dim3
        call halo_exchange(psix)
        call halo_exchange(psiy)        
#endif
        call halo_exchange(psiz)
#endif                
        
        if(iproc.eq.0) write(6,*) "Diffusion progress:",100.0d0*dble(it_diff)/dble(nt_diff),"%."
        
      end do
      
      deallocate(lapx,lapy,lapz)
   
     !! Evaluate gradients
     allocate(gradx(npfb,dims),grady(npfb,dims),gradz(npfb,dims))
#ifdef dim3
     call calc_gradient(psix,gradx)
     call calc_gradient(psiy,grady)
#endif     
     call calc_gradient(psiz,gradz)
     
     !! Evaluate the velocity
     !$omp parallel do
     do i=1,npfb
#ifdef dim3
        u(i) = gradz(i,2) - grady(i,3)
        v(i) = gradx(i,3) - gradz(i,1)
        w(i) = grady(i,1) - gradx(i,2)
#else
        u(i) = gradz(i,2)
        v(i) = -gradz(i,1)
        w(i) = zero
#endif
     
     end do
     !$omp end parallel do
     
     !! Set u,v,w on boundaries
     !$omp parallel do private(i)
     do j=1,nb
        i=boundary_list(j)
        if(node_type(i).eq.0) then
           u(i) = zero;v(i) = zero;w(i) = zero
        end if
     end do
     !$omp end parallel do
      
     deallocate(psix,psiy,psiz,gradx,grady,gradz)
  
     return
  end subroutine noise_and_diffuse_psi
!! ------------------------------------------------------------------------------------------------ 
!! ------------------------------------------------------------------------------------------------  
  subroutine noise_and_diffuse_velocity_field
     use derivatives
     use mpi_transfers
     use mirror_boundaries    
     !! Applies diffusion to a velocity field
     integer(ikind) :: it_diff,nt_diff,i,j
     real(rkind) :: dt_diff,tot_vol,u_mean,v_mean,w_mean,dVi
     real(rkind),dimension(:),allocatable :: lapu,lapv,lapw
     real(rkind) :: randscaling
       
     !! Evaluate time-step and number of iterations    
     dt_diff = 0.1d0*(min(smin_global,dz)*L_char)**two
     nt_diff = 1 + floor(L_turb*L_turb/(two*pi*dt_diff))
     nt_diff = 1 + floor((three**two)*10.0d0/(two*pi)) !! three**two sets the cut-off lengthscale = 3*smin
      
     !! Space for RHS
     allocate(lapu(npfb),lapv(npfb),lapw(npfb))
     
     !! Loop over time-steps
     do it_diff = 1,nt_diff
     
        !! Add some random noise to the velocity field
!        randscaling = zero!(one/dble(it_diff))**10.0d0
!        call add_random_to_velocity_field(randscaling)
    
        !! Update mirrors and halos
        call mirror_bcs_vel_only
#ifdef mp        
        call halo_exchange(u)
        call halo_exchange(v)        
        call halo_exchange(w)
#endif        
                          
        !! Evaluate RHS
        call calc_laplacian(u,lapu)
        call calc_laplacian(v,lapv)
        call calc_laplacian(w,lapw)                
               
        !! Update - 1st order time-stepping is fine
        !! Also evaluate mean velocity components
        u_mean = zero;v_mean = zero;w_mean = zero;tot_vol = zero
        !$omp parallel do private(dVi) reduction(+:u_mean,v_mean,w_mean,tot_vol)
        do i=1,npfb
        
           !! Update
           u(i) = u(i) + dt_diff*lapu(i)
           v(i) = v(i) + dt_diff*lapv(i)
           w(i) = w(i) + dt_diff*lapw(i)                      
 
           !! Node volume           
           dVi = vol(i)
#ifdef dim3
           dVi = dVi*dz
#endif            

           !! Augment sums
           u_mean = u_mean + u(i)*dVi
           v_mean = v_mean + v(i)*dVi
           w_mean = w_mean + w(i)*dVi
           tot_vol = tot_vol + dVi


        end do
        !$omp end parallel do
 
        !! Reduce to find means        
#ifdef mp     
        call global_reduce_sum(u_mean)
        call global_reduce_sum(v_mean)
        call global_reduce_sum(w_mean)
        call global_reduce_sum(tot_vol)          
#endif      
        u_mean = u_mean/tot_vol
        v_mean = v_mean/tot_vol
        w_mean = w_mean/tot_vol

        !! Subtract mean
        !$omp parallel do 
        do i=1,npfb
           u(i) = u(i) - u_mean
           v(i) = v(i) - v_mean
           w(i) = w(i) - w_mean                
        end do
        !$omp end parallel do
              
        
        !! Impose boundary conditions
        !$omp parallel do private(i)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(j).eq.0) then !! Walls
           
           else if(node_type(j).eq.1) then !! Inflow
           
           else if(node_type(j).eq.2) then !! Outflow
         
           end if
           u(i) = zero;v(i) = zero;w(i) = zero  
        end do
        !$omp end parallel do            
        
        if(iproc.eq.0) write(6,*) "Diffusion progress:",100.0d0*dble(it_diff)/dble(nt_diff),"%."!,maxval(u),minval(u)
        
      end do
      
      deallocate(lapu,lapv,lapw)
  
     return
  end subroutine noise_and_diffuse_velocity_field
!! ------------------------------------------------------------------------------------------------  
  subroutine add_random_to_velocity_field(rscal)
     !! Adds some random noise to the velocity field
     real(rkind),intent(in) :: rscal !! scales the random field
     integer(ikind) :: i
     real(rkind) :: rnum,dVi

     !$omp parallel do private(rnum,dVi)
     do i=1,npfb    
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif       
        dVi = one/sqrt(dVi)       
     
        rnum = rand()!call random_number(rnum)
        u(i) = u(i) + rscal*(two*rnum - one)*dVi
        rnum = rand()!call random_number(rnum)
        v(i) = v(i) + rscal*(two*rnum - one)*dVi
        rnum = rand()!call random_number(rnum)
        w(i) = w(i) + rscal*(two*rnum - one)*dVi        
     end do
     !$omp end parallel do      
  
     return
  end subroutine add_random_to_velocity_field 
!! ------------------------------------------------------------------------------------------------    
end module turbulence
