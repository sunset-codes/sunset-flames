module transport
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !! JRCK               |Jan 2023         |New module to implement mix av transport
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to calculate transport properties.
  
  !! Various options ::
  !!        1) isoT         - isothermal flows, viscosity and molecular diffusivity are set 
  !!                              as constant based on reference values.
  !!        2) not(isoT)    - thermal flows, viscosity, molecular diffusivity and thermal
  !!                          conductivity are dependent on temperature, and if mixture-averaged
  !!                          transport is turned on, also on composition.
  !!
  !! N.B. Currently hard-coded with visc, lambda and D given by 4 polynomial coefficients - reference
  !! for these is given in paper and documentation
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  implicit none
  
  private
  public :: evaluate_transport_properties

  real(rkind) :: segment_tstart,segment_tend

contains
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_transport_properties
     use mirror_boundaries
     !! Uses temperature, cp and density to evaluate thermal conductivity, viscosity and 
     !! molecular diffusivity. For isothermal flows, use reference values.
     integer(ikind) :: ispec,i,j
     real(rkind) :: tmp
     segment_tstart=omp_get_wtime()             
       
#ifndef isoT     
     if(flag_mix_av.eq.0) then
        !! Constant Lewis numbers
        !$omp parallel do private(ispec,tmp)        
        do i=1,npfb
     
           !! Viscosity
           visc(i) = visc_ref*(T(i)/T_ref)**r_temp_dependence
     
           !! Thermal conductivity
           lambda_th(i) = cp(i)*visc(i)/Pr

           !! Molecular diffusivity - actually returning ro*Mdiff
           tmp = visc(i)/Pr
           do ispec=1,nspec
              roMdiff(i,ispec) = tmp*one_over_Lewis_number(ispec)
           end do        
   
        end do
        !$omp end parallel do
     else
        !! For mixture averaged transport, there is no (or very complex) analytic relation
        !! between transport gradients and temperature gradients, so we evaluate the transport
        !! properties at halos and copy to mirrors.

!        !$omp parallel do 
        do i=1,npfb    
           call evaluate_mixture_average_transport_at_node(i)    
        end do
!        !$omp end parallel do     
       
        
        !! Copy to mirrors
!        call mirror_bcs_transport_only
#ifdef mp
     do j=npfb+1,np_nohalo
#else
     do j=npfb+1,np
#endif
        i = irelation(j)
        visc(j) = visc(i)
        lambda_th(j) = lambda_th(i)
        do ispec=1,nspec
           roMdiff(j,ispec) = roMdiff(i,ispec)
           tdr(j,ispec) = tdr(i,ispec)
        end do        
     end do

        
#ifdef mp
        !! Calculate in halos (faster than MPI comms)
!        !$omp parallel do 
        do i=np_nohalo+1,np     
           call evaluate_mixture_average_transport_at_node(i)    
        end do
!        !$omp end parallel do     
#endif   
     
     end if
#else
     !! Isothermal options
     visc(:) = visc_ref
     do i=1,np
        roMdiff(i,:) = visc_ref*one_over_Lewis_number(:) 
     end do
#endif     

     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(8) = segment_time_local(8) + segment_tend - segment_tstart  
  
     return
  end subroutine evaluate_transport_properties
!! ------------------------------------------------------------------------------------------------  
  subroutine evaluate_mixture_average_transport_at_node(inode)
     !! Evaluates the transport properties at an individual node based on mixture averaged
     !! rules for individual species transport properties.
     !!
     !! Optimisations in force:
     !!    a) If we hard-code number of coefficients of polynomials, it is faster to write them out
     !!       explicitly and re-used powers of logT, rather than use a reverse loop.
     !!    b) Powers of visc and lambda required for combination rules are done in log-space.
     !!    c) We store one over species_diffusivity from the polynomials, which replaces a number
     !!       of divisions by multiplications in the combo rule.
     !!    d) Eliminated density from the formulation.
     !!    e) Replace sum of X for jspec!=ispec with one-Y(ispec)
     !!
     integer(ikind),intent(in) :: inode
     integer(ikind) :: ispec,jspec 
     real(rkind) :: logT,roo_molar_mass_mix,logp,tmpro
     real(rkind),dimension(nspec_max) :: species_visc,Xspec,species_lambda
     real(rkind),dimension(nspec_max,nspec_max) :: oospecies_diff,species_tdr
     real(rkind) :: store1,store2,store3
     real(rkind) :: logT2,logT3,Ti,T2,T3
     
     !! Logarithm of temperature ratio
     logT = log(T(inode)/T_ref_mxav)
     
     !! Normalise Ti by reference
     Ti = T(inode)/T_ref_mxav         
     
     !! Powers of logT
     logT2=logT*logT
     logT3=logT2*logT
     
#ifdef sorduf     
     !! Powers of T
     T2 = Ti*Ti
     T3 = Ti*T2
#endif     
     
     !! Store the local density
     tmpro = ro(inode)
     
     !! log of Pressure ratio (relative to reference)
     logp = log(p_ref_mxav/p(inode))
     
     !! Loop over all species and evaluate viscosity, thermal conductivity, diffusivity,
     !! ro/mixture mol mass, and Xspec
     roo_molar_mass_mix = zero
     do ispec=1,nspec

        !! Build Xspec and density over mixture molar mass
        Xspec(ispec) = Yspec(inode,ispec)*one_over_molar_mass(ispec)
        roo_molar_mass_mix = roo_molar_mass_mix + Xspec(ispec)

        !! Evaluate logarithm of species viscosity
        species_visc(ispec) = mxav_coef_visc(ispec,1) + logT*mxav_coef_visc(ispec,2) &
                            + logT2*mxav_coef_visc(ispec,3) + logT3*mxav_coef_visc(ispec,4) 

        !! Evaluate logarithm of species thermal conductivity        
        species_lambda(ispec) = mxav_coef_lambda(ispec,1) + logT*mxav_coef_lambda(ispec,2) &
                              + logT2*mxav_coef_lambda(ispec,3) + logT3*mxav_coef_lambda(ispec,4) 

        !! Evaluate species diffusivity (actually one over species diffusivity) for lower
        !! triangle, then copy across
        do jspec=1,ispec
           store1 = -logp - mxav_coef_diff(ispec,jspec,1) - logT*mxav_coef_diff(ispec,jspec,2) &
                   - logT2*mxav_coef_diff(ispec,jspec,3) - logT3*mxav_coef_diff(ispec,jspec,4) 
           oospecies_diff(ispec,jspec) = exp(store1)
           oospecies_diff(jspec,ispec) = oospecies_diff(ispec,jspec) !! Copy to upper triangle        
        end do
        
#ifdef sorduf        
        !! Evaluate thermal diffusion ratio for lower triangle, then copy accross
        do jspec=1,ispec-1
           species_tdr(ispec,jspec) = mxav_coef_tdr(ispec,jspec,1) + mxav_coef_tdr(ispec,jspec,2)*Ti &
                                    + mxav_coef_tdr(ispec,jspec,3)*T2 + mxav_coef_tdr(ispec,jspec,4)*T3
                                      
           species_tdr(jspec,ispec) = species_tdr(ispec,jspec)
        end do
        species_tdr(ispec,ispec) = zero
#endif        
        
     end do
          
   
     !! Finalise Xspec  
     Xspec(1:nspec) = Xspec(1:nspec)/roo_molar_mass_mix

     !! ---------------------------------------------------------------------------------
     !! Combination rule for viscosity (E&G order 6) ====================================
     store1=zero 
     do ispec=1,nspec
        store2 = exp(six*species_visc(ispec))
        store1 = store1 + Xspec(ispec)*store2
     end do
     visc(inode) = store1**oosix
     
     !! ---------------------------------------------------------------------------------     
     !! Combination rule for thermal conductivity (E&G order 1/4) =======================
     store1=zero 
     do ispec=1,nspec
        store2 = exp(quarter*species_lambda(ispec))
        store1 = store1 + Xspec(ispec)*store2
     end do
     lambda_th(inode) = store1*store1*store1*store1  !! Power of 4 explicit

     !! ---------------------------------------------------------------------------------
     !! Combination rule for molecular diffusivity (Hirschfelder & Curtiss rule) ========
     do ispec = 1,nspec
        !! Initialise with the ispec-ispec contribution subtracted
        store3 = Xspec(ispec) + quitesmall
        store2 = -store3*oospecies_diff(ispec,ispec)
                
        do jspec = 1,nspec
           store3 = Xspec(jspec) + quitesmall
           store2 = store2 + store3*oospecies_diff(ispec,jspec)
        end do
        roMdiff(inode,ispec) = (tmpro-Yspec(inode,ispec))/store2
     end do

#ifdef sorduf     
     !! ---------------------------------------------------------------------------------
     !! Combination rule for thermal diffusion ratio ====================================
     do ispec = 1,nspec
        tdr(inode,ispec) = zero
        do jspec = 1,nspec
           tdr(inode,ispec) = tdr(inode,ispec) + species_tdr(ispec,jspec)*Xspec(jspec)
        end do     
     end do
#endif     



     return
  end subroutine evaluate_mixture_average_transport_at_node
!! ------------------------------------------------------------------------------------------------    
end module transport
