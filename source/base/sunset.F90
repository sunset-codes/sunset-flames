program sunset
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This is the main program of the sunset code. 
  use kind_parameters
  use common_parameter
  use common_vars
  use load_data
  use setup_domain
  use setup_flow
  use output
  use statistics 
  use neighbours
  use labf
  use fd
  use step
#ifdef mp  
  use mpi
#endif  
  implicit none

  integer(ikind) :: n_out,m_out
  real(rkind) :: segment_tstart_main,segment_tend_main

#ifdef mp  
  call MPI_INIT(ierror)
#endif  

  !! Load a little control data, node data, and build MPI decomposition etc
  call load_control_data_LUonly   
  call initial_setup  
  call build_domain
  
  !! Adapt the stencils by reducing h (only if not restarting)
#ifndef restart  
!  call adapt_stencils
  call grow_stencils
#endif  

  !! Shrink the halos to fit, and finalise building the domain
  call refine_and_finalise_domain

  !! Build the neighbour lists
  call find_neighbours
  call order_neighbours

  !! Calculate LABFM and FD weights for all derivatives and filters
  call calc_labf_weights
  if(nb.ne.0) call calc_boundary_weights
  call calc_labf_sums
#ifdef dim3
  call calc_fd_weights
#endif    
  call filter_coefficients   
  
  !! Load data 
  call load_control_data_all
!#ifndef isoT
  call load_chemistry_data
!#endif  
  if(flag_mix_av.eq.1) call load_transport_file

  !! Set initial fields
  call initial_solution
#ifdef restart
  call set_tstep
  dt = half*half*min(dt_cfl,dt_parabolic)
#endif   

  !! Initialise time profiling and output counter...
  n_out = 0;ts_start = omp_get_wtime()
  m_out = 0  ! Stats output counter
        
  !! MAIN TIME LOOP ---------------------------------------------------
  do while (time.le.time_end)

     !! Output, conditionally: at start, subsequently every dt_out
     if(itime.eq.0.or.time.gt.n_out*dt_out) then 
!     if(itime.eq.0.or.mod(itime,1).eq.0)then  !! For debugging, comment line above, and uncomment/modify this line
        n_out = n_out + 1
        call output_layer(n_out)        
        call output_laminar_flame_structure(n_out)
     end if     

     !! Profiling
     segment_tstart_main = omp_get_wtime()

     !! Adjust the time-step
     call set_tstep      !! CFL type stability based dt    
#ifdef react     
     call set_tstep_PID  !! Error estimation based dt (in addition!)
#endif     
 
    
     !! Perform one time step (with Error estimation if reacting flows)
#ifdef react
     call step_rk3_4S_2R_EE     
#else     
     call step_rk3_4S_2R
#endif     

     !! Profiling and write some things to screen
     segment_tend_main = omp_get_wtime()
     segment_time_local(1) = segment_time_local(1) + segment_tend_main - segment_tstart_main
     itime = itime + 1
     call output_to_screen

     !! Call routines to evaluate global statistics and adjust forcing terms if desired
     call statistics_control(m_out)
    
  end do
  !! END MAIN TIME LOOP -----------------------------------------------
  
  !! Deallocate particle properties and neighbour lists
  call deallocate_everything
#ifdef mp    
  call MPI_FINALIZE(ierror)
#endif  
  stop
end program sunset
!! ------------------------------------------------------------------------------------------------
subroutine deallocate_everything
  use kind_parameters
  use common_parameter
  use common_vars
  
  !! Tell the screen we've reached the end
#ifdef mp
  if(iproc.eq.0) write(6,*) "Reached the end of simulation. Cleaning up and stopping!"
#else  
  write(6,*) "Reached the end of simulation. Cleaning up and stopping!"
#endif  
  
  !! Discretisation arrays
  deallocate(rp,s,vol,h)
  
  !! Primary properties
  deallocate(rou,rov,row,ro,roE,Yspec)
  
  !! Secondary properties & transport vars
  deallocate(T,p,cp,visc,Rgas_mix,u,v,w,hrr)
#ifndef isoT
  deallocate(lambda_th)
#endif
  deallocate(roMdiff)
#ifdef sorduf
  deallocate(tdr)
#endif

  !! Neighbours lists and link lists
  deallocate(ij_count,ij_link)
  deallocate(irelation,vrelation)
  deallocate(node_type,internal_list)
  if(allocated(boundary_list)) deallocate(boundary_list)
  
  !! LABFM and FD weightings
  if(allocated(ij_w_grad)) then
     deallocate(ij_w_grad,ij_w_hyp)
     deallocate(ij_w_grad_sum,ij_w_hyp_sum)     
  end if
  if(allocated(ij_wb_grad2)) deallocate(ij_wb_grad2,ij_wb_grad2_sum)
  if(allocated(ij_link_fd)) then
     deallocate(ij_link_fd)
     deallocate(ij_fd_grad,ij_fd_grad2,ij_fd_hyp)
     deallocate(zlayer_index_global)
  end if
  
  !! Transport properties
  if(allocated(molar_mass)) deallocate(molar_mass,one_over_molar_mass)
  if(allocated(one_over_Lewis_number)) deallocate(one_over_Lewis_number)
  
  if(flag_mix_av.eq.1) deallocate(mxav_coef_visc,mxav_coef_lambda,mxav_coef_diff,mxav_coef_tdr)
  
  
  return
end subroutine deallocate_everything
