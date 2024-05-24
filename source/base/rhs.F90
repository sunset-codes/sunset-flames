module rhs
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to construct the RHS of all evolution equations
  !! These RHS' are built using calls to routines from the derivatives module, and
  !! calls to specific thermodynamics routines.
  
  !! Although separate subroutines, they must be called in correct order, as they rely
  !! on each other
  
  !! We use the lists internal_list and boundary_list to loop through internal and boundary
  !! nodes respectively. The L arrays for boundaries run from 1 to nb, and so use index j
  !! when within a loop over boundary nodes.
  use kind_parameters
  use common_parameter
  use common_vars
  use derivatives
  use thermodynamics
  use transport
  use characteristic_boundaries
  use chemistry
  implicit none
  
 
  private
  public calc_all_rhs,filter_variables

  !! Allocatable arrays for 1st and 2nd gradients
  real(rkind),dimension(:,:),allocatable :: gradro,gradp  !! Velocity gradients defined in common, as used elsewhere too
  real(rkind),dimension(:,:),allocatable :: gradroE
  real(rkind),dimension(:,:),allocatable :: graddivvel
  real(rkind),dimension(:),allocatable :: lapu,lapv,lapw,lapT
  real(rkind),dimension(:,:),allocatable :: gradT,gradcp
  
  real(rkind),dimension(:),allocatable :: store_diff_E
    
  real(rkind) :: dundn,dutdn,dutdt,dpdn
  real(rkind) :: xn,yn,un,ut
  
  !! Characteristic boundary condition formulation
  real(rkind),dimension(:,:),allocatable :: L  !! The "L" in NSCBC formulation    
  

contains
!! ------------------------------------------------------------------------------------------------
  subroutine calc_all_rhs
     use statistics
     integer(ikind) :: i
     !! Control routine for calculating right hand sides. Does thermodynamic evaluations, finds
     !! gradients, and then calls property-specific RHS routines
          
     !! Some initial allocation of space for boundaries
     if(nb.ne.0) allocate(L(nb,5+nspec)) 

     !! Determine secondary thermodynamic quantities and transport properties
     !! N.B. These are also calculated when setting time step, so not necessary for the first call
     !! to calc_all_rhs in the RK scheme.
     if(iRKstep.ne.1) then   
        call evaluate_temperature_and_pressure
        call evaluate_transport_properties      
     end if

     !! Initialise right hand sides to zero
     rhs_ro=zero;rhs_rou=zero;rhs_rov=zero;rhs_row=zero;rhs_roE=zero;rhs_Yspec=zero
             
     !! Calculate derivatives of primary variables
     allocate(gradro(npfb,dims));gradro=zero  
     allocate(gradu(npfb,dims),gradv(npfb,dims),gradw(npfb,dims));gradw=zero
     call calc_gradient(ro,gradro)     
     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)     
#ifdef dim3
     call calc_gradient(w,gradw)
#endif    
#ifndef isoT
     allocate(gradroE(npfb,dims))
     call calc_gradient(roE,gradroE)

     !! Temperature gradient and Laplacian
     allocate(gradT(npfb,dims),lapT(npfb))
     call calc_gradient(T,gradT)        
     call calc_laplacian_transverse_only_on_bound(T,lapT)                 
#endif         
     
     !! Pressure gradient (method depends on whether isoT or not)
     allocate(gradp(npfb,dims));gradp=zero
#ifndef isoT
     !! Evaluate pressure gradient directly (LABFM...)
     call calc_gradient(p-p_ref,gradp) 


#else
     !! Evaluate the pressure gradient (from density gradient
     !$omp parallel do
     do i=1,npfb  
        gradp(i,:) = csq*gradro(i,:)  !! N.B. not precisely correct for isothermal multispec
     end do
     !$omp end parallel do
#endif                 

     !! Call individual routines to build the RHSs
     !! N.B. second derivatives and derivatives of secondary variables are calculated within
     !! these subroutines
     call calc_rhs_ro
     call calc_rhs_Yspec
     call calc_rhs_rovel
     call calc_rhs_roE

     !! Calculate chemical production rates and add these to rhs of species equation
     call calculate_chemical_production_rate

     !! Evaluate RHS for boundaries
     if(nb.ne.0) then 
        call calc_rhs_nscbc   
        deallocate(sumoverspecies_homega,reaction_rate_bound)
     end if
      
          
     !! If we want to calculate total dissipation rate
     if(.false..and.iRKstep.eq.1) then
        call check_enstrophy
     end if    
     
     !! Clear space no longer required
     deallocate(gradro,gradu,gradv,gradw,gradp)
#ifndef isoT     
     deallocate(gradroE)
     deallocate(gradT,lapT)
#endif     
  
     return
  end subroutine calc_all_rhs
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_ro
     !! Construct the RHS for ro-equation
     integer(ikind) :: i,j
     real(rkind),dimension(dims) :: tmp_vec
     real(rkind) :: tmp_scal
    
     !! Build RHS for internal nodes
     !$omp parallel do private(i,tmp_vec,tmp_scal)
     do j=1,npfb-nb
        i=internal_list(j)
        tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3)= w(i)
        tmp_scal = dot_product(tmp_vec,gradro(i,:))

        rhs_ro(i) = -ro(i)*divvel(i) - tmp_scal   
        
     end do
     !$omp end parallel do

     !! For any boundary nodes make transverse part of rhs
     if(nb.ne.0)then
        !$omp parallel do private(i,tmp_scal,xn,yn,un,ut,dutdt)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0)then  !! in bound norm coords for walls
              xn=rnorm(i,1);yn=rnorm(i,2)
              un = u(i)*xn + v(i)*yn; ut = -u(i)*yn + v(i)*xn  !! Normal and transverse components of velocity           
              dutdt = -yn*gradu(i,2)+xn*gradv(i,2) !! Transverse derivative of transverse velocity...
              

              rhs_ro(i) = - ro(i)*dutdt - ro(i)*gradw(i,3)
              
           else !! In x-y coords for inflow, outflow
              rhs_ro(i) = -v(i)*gradro(i,2) - ro(i)*gradv(i,2) - w(i)*gradro(i,3) - ro(i)*gradw(i,3)
           end if  
        end do
        !$omp end parallel do 
     end if       

     return
  end subroutine calc_rhs_ro
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_rhs_Yspec
     !! Construct the RHS for species Yspec equation
     !! N.B. the variables roVY and divroVY, and their sums, hold *minus* roVY and minus divroVY
     real(rkind),dimension(:),allocatable :: lapYspec,Y_thisspec
     real(rkind),dimension(:,:,:),allocatable :: gradYspec
     real(rkind),dimension(:,:),allocatable :: grad2Yspec
     integer(ikind) :: i,j,ispec
     real(rkind),dimension(dims) :: gradroDY,body_force
     real(rkind) :: tmp_scal,tmpY,divroVY,enthalpy,dcpdT,cpispec,tmpro,tmpT,roDY
     real(rkind) :: q00,q01,q02,q03,q04
     real(rkind),dimension(:,:),allocatable :: speciessum_roVY,speciessum_hgradY
     real(rkind),dimension(:),allocatable :: speciessum_divroVY,speciessum_hY
     real(rkind),dimension(:,:),allocatable :: gradroMdiff
     real(rkind),dimension(:,:),allocatable :: mxav_store1,mxav_store2
     real(rkind),dimension(:),allocatable :: mxav_store3,mxav_store4     
     real(rkind),dimension(:,:),allocatable :: roVY

     !! Allocate space for gradients and stores
     allocate(Y_thisspec(np));Y_thisspec = zero
     allocate(store_diff_E(npfb));store_diff_E = zero
#ifndef isoT     
     allocate(gradcp(npfb,dims));gradcp = zero
#endif     


     allocate(gradYspec(npfb,dims,nspec));gradYspec=zero
     allocate(lapYspec(npfb));lapYspec=zero
     allocate(speciessum_divroVY(npfb));speciessum_divroVY = zero
     allocate(speciessum_roVY(npfb,dims));speciessum_roVY = zero
     allocate(roVY(npfb,dims));roVY=zero
     allocate(speciessum_hY(npfb));speciessum_hY = zero
     allocate(speciessum_hgradY(npfb,dims));speciessum_hgradY = zero
     allocate(gradroMdiff(npfb,dims))
     
     !! For mixture averaged transport, pre-calculate additional diffusion driving forces
     allocate(mxav_store1(npfb,dims),mxav_store2(npfb,dims))
     allocate(mxav_store3(npfb),mxav_store4(npfb))
     if(flag_mix_av.eq.1) then
     
        !! Pre-populate stores 3 and 4 with density and pressure laplacians
        call calc_laplacian_transverse_only_on_bound(ro,mxav_store3)
        call calc_laplacian_transverse_only_on_bound(p-p_ref,mxav_store4)        
        
   
        !$omp parallel do private(tmpro,tmpT,i)
        do j=1,npfb-nb
           i=internal_list(j)
           !! Store inverse of density and temperature
           tmpro = one/ro(i)
           tmpT=one/T(i)        
           
           !! grad(lnT) + grad(lnro)
           mxav_store1(i,:) = tmpT*gradT(i,:) + tmpro*gradro(i,:)
           
           !! grad(p)/(ro*R0*T)
           mxav_store2(i,:) = tmpro*tmpT*gradp(i,:)/Rgas_universal
           
           !! lap(T)/T - grad(lnT)**2 + lap(ro)/ro - grad(lnro).grad(lnro)
           mxav_store3(i) = lapT(i)*tmpT - tmpT*tmpT*dot_product(gradT(i,:),gradT(i,:)) &
                          + mxav_store3(i)*tmpro - tmpro*tmpro*dot_product(gradro(i,:),gradro(i,:))
                          
           !!(1/roR0T)*(lap(p) - gradp.(grad(lnT)+grad(lnro)))               
           mxav_store4(i) = (tmpro*tmpT/Rgas_universal)*(mxav_store4(i) - dot_product(gradp(i,:),mxav_store1(i,:)))


        end do
        !$omp end parallel do 
        
        !! On boundary nodes, we neglect any normal terms
        if(nb.ne.0) then              
           !$omp parallel do private(i)
           do j=1,nb
              i=boundary_list(j)

              !! grad(lnT) + grad(lnro)
              mxav_store1(i,:) = tmpT*gradT(i,:) + tmpro*gradro(i,:)
           
              !! grad(p)/(ro*R0*T)
              mxav_store2(i,:) = tmpro*tmpT*gradp(i,:)/Rgas_universal              
                            
              !! lap(T)/T - grad(lnT)**2 + lap(ro)/ro - grad(lnro).grad(lnro)
              mxav_store3(i) = lapT(i)*tmpT - tmpT*tmpT*dot_product(gradT(i,2:3),gradT(i,2:3)) &
                             + mxav_store3(i)*tmpro - tmpro*tmpro*dot_product(gradro(i,2:3),gradro(i,2:3))
                          
              !!(1/roR0T)*(lap(p) - gradp.(grad(lnT)+grad(lnro)))               
              mxav_store4(i) = (tmpro*tmpT/Rgas_universal)*(mxav_store4(i) &
                             - dot_product(gradp(i,2:3),mxav_store1(i,2:3)))            
           end do
           !$omp end parallel do
        end if
        
     else
        !! These terms are zero for constant Lewis number assumption
        mxav_store1=zero;mxav_store2=zero
        mxav_store3=zero;mxav_store4=zero
     end if
     
     !! Loop over all species
     do ispec=1,nspec
    
        !! Store Y=roY/ro for the species (N.B. this loop is over ALL nodes)
        !$omp parallel do 
        do i=1,np
           Y_thisspec(i) = Yspec(i,ispec)/ro(i)
        end do
        !$omp end parallel do
        
        !! Calculate gradient and Laplacian for Yspec for this species     
        call calc_gradient(Y_thisspec,gradYspec(:,:,ispec))
        call calc_laplacian_transverse_only_on_bound(Y_thisspec,lapYspec)


        !! Diffusivity gradients
        if(flag_mix_av.eq.1) then
           call calc_gradient(roMdiff(:,ispec),gradroMdiff)           
        else
#ifndef isoT
           !$omp parallel do
           do i=1,npfb
              gradroMdiff(i,:) = r_temp_dependence*roMdiff(i,ispec)*gradT(i,:)/T(i)
           end do
           !$omp end parallel do
#else
           gradroMdiff(:,:) = zero
#endif           
        end if
        
        !! Zero diffusive flux
        roVY = zero

segment_tstart=omp_get_wtime()    
      
        !$omp parallel do private(i,tmp_scal,tmpY,divroVY,enthalpy, &
        !$omp dcpdT,cpispec,tmpro,gradroDY,roDY)
        do j=1,npfb-nb
           i=internal_list(j)
           tmpro = one/ro(i) !! tmpro contains 1/ro

           !! Convective term: ro*u.gradY + Y(div.(ro*u))        
           tmp_scal = ro(i)*(u(i)*gradYspec(i,1,ispec) + &
                             v(i)*gradYspec(i,2,ispec) + &
                             w(i)*gradYspec(i,3,ispec) ) - Y_thisspec(i)*rhs_ro(i)             
         
           !! First part of molecular diffusion terms 
           roVY(i,:) = roMdiff(i,ispec)*gradYspec(i,:,ispec)                     
           divroVY = roMdiff(i,ispec)*lapYspec(i) &
                      + dot_product(gradYspec(i,:,ispec),gradroMdiff(i,:))
           
           
           !! roDY and gradroDY - these are used for additional mixture averaged terms
           roDY = roMdiff(i,ispec)*Y_thisspec(i)
           gradroDY = roVY(i,:) + Y_thisspec(i)*gradroMdiff(i,:)
                     
           !! Augment molecular diffusion terms to include mixture averaged driving forces 
           roVY(i,:) = roVY(i,:) + roDY*(mxav_store1(i,:) + mxav_store2(i,:)*molar_mass(ispec))                      
           divroVY = divroVY + dot_product(gradroDY,mxav_store1(i,:)-mxav_store2(i,:)*molar_mass(ispec)) &
                             + roDY*(mxav_store3(i) + mxav_store4(i)*molar_mass(ispec))        
                      
           
           !! Sum roVY and div(roVY) over species
           speciessum_divroVY(i) = speciessum_divroVY(i) + divroVY
           speciessum_roVY(i,:) = speciessum_roVY(i,:) + roVY(i,:)

#ifndef isoT                                 
           !! Evaluate enthalpy, cp and dcp/dT for species ispec
           call evaluate_enthalpy_at_node(T(i),ispec,enthalpy,cpispec,dcpdT)

           !! Add h*div.(roVY) for species ispec to the energy diffusion store
           store_diff_E(i) = store_diff_E(i) + divroVY*enthalpy   
               
           !! Add gradh.roVY for species ispec ( gradh = cp*gradT )
           store_diff_E(i) = store_diff_E(i) + cpispec*dot_product(roVY(i,:),gradT(i,:))
                                               
           !! Add this species contribution to the mixture enthalpy
           speciessum_hY(i) = speciessum_hY(i) + enthalpy*Y_thisspec(i) 
           
           !! Add this species contribution to hgradY
           speciessum_hgradY(i,:) = speciessum_hgradY(i,:) + enthalpy*gradYspec(i,:,ispec)
           
           !! Add this species contrib to gradcp(mix) = YdcpdT*gradT + cp(ispec)gradY
           gradcp(i,:) = gradcp(i,:) + Y_thisspec(i)*dcpdT*gradT(i,:) &
                                     + cpispec*gradYspec(i,:,ispec)
#endif                    
                                          
           !! Add convective and diffusive terms to the RHS
           rhs_Yspec(i,ispec) = -tmp_scal + divroVY 
        end do
        !$omp end parallel do

!! Profiling
segment_tend = omp_get_wtime()
segment_time_local(7) = segment_time_local(7) + segment_tend - segment_tstart

        !! Make L5+ispec and boundary RHS
        if(nb.ne.0)then
           allocate(grad2Yspec(nb,dims));grad2Yspec=zero
           call calc_grad2bound(Y_thisspec,grad2Yspec)      
   
           !$omp parallel do private(i,xn,yn,un,ut,dutdt,tmpY &
           !$omp ,divroVY,enthalpy,tmpro,cpispec,dcpdT,tmp_scal,roDY,gradroDY,q00,q01,q02,q03,q04)
           do j=1,nb
              i=boundary_list(j)
              tmpro = one/ro(i)  !! tmpro contains 1/ro
              
              !! Convective term: ro*u.gradY + Y(div.(ro*u)) with zero boundary normal term
              tmp_scal = ro(i)*(zero*gradYspec(i,1,ispec) + &
                                v(i)*gradYspec(i,2,ispec) + &
                                w(i)*gradYspec(i,3,ispec) ) - Y_thisspec(i)*rhs_ro(i)
                                                               
              !! Molecular diffusion terms 
              roVY(i,:) = roMdiff(i,ispec)*gradYspec(i,:,ispec)                              

              !! NOTE: In all boundary cases we start by constructing only transverse parts of divroVY
              
              !! divroVY for constant Lewis approximation
              divroVY = roMdiff(i,ispec)*lapYspec(i) &
                      + dot_product(gradYspec(i,2:3,ispec),gradroMdiff(i,2:3))  
                                            
              !! roDY and grad(roDY) used in mixture averaged formulation
              roDY = roMdiff(i,ispec)*Y_thisspec(i)
              gradroDY = roVY(i,:) + Y_thisspec(i)*gradroMdiff(i,:)
              
              
              !! Augment to include mixture averaged driving forces
              roVY(i,:) = roVY(i,:) + roDY*(mxav_store1(i,:) + mxav_store2(i,:)*molar_mass(ispec))
              
              !! Add  gradroDY*(grad(lnT)+grad(lnro)-WgradP/roRT) (NB. transverse terms only)
              divroVY = divroVY + dot_product(gradroDY(2:3),mxav_store1(i,2:3)-mxav_store2(i,2:3)*molar_mass(ispec))

              !! Add roDY*(grad(grad(lnT)+grad(lnro)-WgradP/roRT)) (NB. s3 and s4 contain transverse terms only)
              divroVY = divroVY + roDY*(mxav_store3(i) + mxav_store4(i)*molar_mass(ispec))
           
              !! Up to this point, divroVY only contains transverse terms
                   
              !! zero normal components of flux if required
              if(node_type(i).eq.0) then  !! Wall, add modified normal flux derivative to enforce roVY.n=0

                 q00 = Y_thisspec(i)
                 q01 = Y_thisspec(i+1)
                 q02 = Y_thisspec(i+2)
                 q03 = Y_thisspec(i+3)
                 q04 = Y_thisspec(i+4)
                 
                 roVY(i,1) = zero
                 divroVY = divroVY + roMdiff(i,ispec)* &
                                   (-170.0d0*q00 + 216.0d0*q01 - 54.0d0*q02 + 8.0d0*q03)/ &
                                   (36.0d0*s(i)*s(i)*L_char*L_char)
                                                                        
!! DEBUG: how do we add in the mixture averaged parts??!?!
              else if(node_type(i).eq.1.or.node_type(i).eq.2) then !! Inflow or outflow
                 if(znf_mdiff(j)) then

                    !! Setting d(roVY.n)/dn = zero (as in Sutherland & Kennedy and others)
                    !! Set roVY to zero, and don't add any d(roVY.n)/dn term
                    roVY(i,1) = zero                    
                 else
       
                    !! Add d(roVY.n)/dn term
                    divroVY = divroVY + (-25.0d0*roVY(i,1)+48.0d0*roVY(i+1,1)-36.0d0*roVY(i+2,1) &
                                         +16.0d0*roVY(i+3,1)-three*roVY(i+4,1))/(12.0d0*s(i)*L_char) !&
                                      !+ roVY(i,1)/boundary radius of curvature...
                 end if
              end if
                                         
              
              !! Sum roVY and div(roVY) over species
              speciessum_divroVY(i) = speciessum_divroVY(i) + divroVY

              !! sum of species of roVY  
              speciessum_roVY(i,:) = speciessum_roVY(i,:) + roVY(i,:)      

#ifndef isoT
              !! Evaluate enthalpy, cp and dcp/dT for species ispec
              call evaluate_enthalpy_at_node(T(i),ispec,enthalpy,cpispec,dcpdT)

              !! Add h*div.(ro*D*gradY) for species ispec to the energy diffusion store
              store_diff_E(i) = store_diff_E(i) + divroVY*enthalpy
          
              !! Add gradh(ispec).roVY for species ispec (gradh = cp*gradT)
              store_diff_E(i) = store_diff_E(i) + cpispec*dot_product(roVY(i,:),gradT(i,:))

              !! Add this species contribution to the mixture enthalpy
              speciessum_hY(i) = speciessum_hY(i) + enthalpy*Y_thisspec(i)

              !! Add this species contribution to hgradY
              speciessum_hgradY(i,:) = speciessum_hgradY(i,:) + enthalpy*gradYspec(i,:,ispec)
          
              !! Add this species contrib to gradcp(mix) = YdcpdT*gradT + cp(ispec)gradY
              gradcp(i,:) = gradcp(i,:) + Y_thisspec(i)*dcpdT*gradT(i,:) &
                                           + cpispec*gradYspec(i,:,ispec)
#endif                   
                 
              !! Construct RHS (transverse convective and diffusive) (v(i)=w(i)=zero if WALL)
              rhs_Yspec(i,ispec) = -tmp_scal + divroVY

              !! Build the characteristic
              L(j,5+ispec) = u(i)*gradYspec(i,1,ispec)    !! (u(i)=zero if WALL)                         
              
           end do
           !$omp end parallel do 
           deallocate(grad2Yspec)
                       
        end if                         

     end do
    
     segment_tstart = omp_get_wtime()
    
     !! Deallocate any stores no longer required    
     deallocate(lapYspec,Y_thisspec)
     deallocate(gradroMdiff)     
     deallocate(roVY)
     deallocate(mxav_store1,mxav_store2,mxav_store3,mxav_store4)

     !! Run through species again and finalise rhs
     do ispec=1,nspec
     
        !$omp parallel do private(tmpro)
        do i=1,npfb
           tmpro = one/ro(i) !! tmpro contains 1/ro

           !! Add the diffusion correction term to the rhs
           rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) &
                              - tmpro*Yspec(i,ispec)*speciessum_divroVY(i) &  !! Y*sum(divroVY)
                              - dot_product(gradYspec(i,:,ispec),speciessum_roVY(i,:))  !! gradY.sum(roVY)   
                              

        end do
        !$omp end parallel do
     end do
     
     !! Run over all nodes one final time to add terms to energy equation
     !$omp parallel do private(body_force)
     do i=1,npfb
        !! Additional terms for energy equation
#ifndef isoT                                     
        !! Add h*diffusion_correction_term to the energy diffusion store
        store_diff_E(i) = store_diff_E(i) - speciessum_divroVY(i)*speciessum_hY(i)
        
        !! Add sum_over_species(h*gradY).sum_over_species(roVY) 
        store_diff_E(i) = store_diff_E(i) - dot_product(speciessum_hgradY(i,:),speciessum_roVY(i,:))

        !! Add cp(mix)*gradT*sum_over_species(roVY) 
        store_diff_E(i) = store_diff_E(i) - cp(i)*dot_product(gradT(i,:),speciessum_roVY(i,:))

        !! Body force molecular diffusion heating term
        body_force = tmpro*grav + driving_force
        store_diff_E(i) = store_diff_E(i) - dot_product(body_force,speciessum_roVY(i,:))
#endif                                               
     end do
     !$omp end parallel do          

     deallocate(speciessum_divroVY,speciessum_roVY,gradYspec)
     deallocate(speciessum_hY,speciessum_hgradY)

!! Profiling
segment_tend = omp_get_wtime()
segment_time_local(7) = segment_time_local(7) + segment_tend - segment_tstart
          
     return
  end subroutine calc_rhs_Yspec    
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_rovel
     !! Construct the RHS for u-equation
     integer(ikind) :: i,j
     real(rkind),dimension(dims) :: tmp_vec
     real(rkind) :: tmp_scal_u,tmp_scal_v,tmp_scal_w,f_visc_u,f_visc_v,f_visc_w
     real(rkind) :: tmpro,body_force_u,body_force_v,body_force_w
     real(rkind) :: c
     real(rkind),dimension(:,:),allocatable :: grad2uvec,grad2ucross    
     real(rkind),dimension(:,:),allocatable :: gradvisc
     
     !! Gradient of velocity divergence
     allocate(graddivvel(npfb,dims));graddivvel=zero
     call calc_gradient(divvel,graddivvel)
          
     !! Allocate memory for spatial derivatives and stores
     allocate(lapu(npfb),lapv(npfb),lapw(npfb))
     lapu=zero;lapv=zero;lapw=zero

     !! Calculate spatial derivatives
     call calc_laplacian(u,lapu)
     call calc_laplacian(v,lapv)
#ifdef dim3
     call calc_laplacian(w,lapw)          
#endif     
                      
         
     !! Evaluate the viscosity gradient (different methods depending on whether mixture averaged)
#ifndef isoT
     allocate(gradvisc(npfb,dims))
     if(flag_mix_av.eq.1) then
        call calc_gradient(visc,gradvisc)         
     else
        !$omp parallel do
        do i=1,npfb
           gradvisc(i,:) = r_temp_dependence*visc(i)*gradT(i,:)/T(i)
        end do
        !$omp end parallel do
     end if
#endif     
         
     !! Build RHS for internal nodes
     !$omp parallel do private(i,tmp_vec,tmp_scal_u,tmp_scal_v,tmp_scal_w,f_visc_u,f_visc_v,f_visc_w,tmpro &
     !$omp ,body_force_u,body_force_v,body_force_w)
     do j=1,npfb-nb
        i=internal_list(j)
        tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3) = w(i) !! tmp_vec holds (u,v,w) for node i
        tmp_scal_u = ro(i)*dot_product(tmp_vec,gradu(i,:)) - u(i)*rhs_ro(i) !! convective term for u
        tmp_scal_v = ro(i)*dot_product(tmp_vec,gradv(i,:)) - v(i)*rhs_ro(i) !! Convective term for v   
        tmp_scal_w = ro(i)*dot_product(tmp_vec,gradw(i,:)) - w(i)*rhs_ro(i) !! Convective term for w    
                
        !! Viscous term (Lap(U) + (1/3)grad(div.U) formulation - avoids explicit calculation of cross derivs)
#ifndef isoT        
        f_visc_u = visc(i)*(lapu(i) + onethird*graddivvel(i,1))
        f_visc_v = visc(i)*(lapv(i) + onethird*graddivvel(i,2))
        f_visc_w = visc(i)*(lapw(i) + onethird*graddivvel(i,3))
#else
        f_visc_u = visc(i)*(lapu(i))  !! Isothermal flows assume small velocity divergence
        f_visc_v = visc(i)*(lapv(i))
        f_visc_w = visc(i)*(lapw(i))
#endif        
        
        !! Viscous forces due to non-uniform viscosity
#ifndef isoT
        f_visc_u = f_visc_u + gradvisc(i,1)*(fourthirds*gradu(i,1) - twothirds*(gradv(i,2)+gradw(i,3))) &
                            + gradvisc(i,2)*(gradu(i,2)+gradv(i,1)) &
                            + gradvisc(i,3)*(gradu(i,3)+gradw(i,1))
        f_visc_v = f_visc_v + gradvisc(i,1)*(gradu(i,2)+gradv(i,1)) &
                            + gradvisc(i,2)*(fourthirds*gradv(i,2) - twothirds*(gradu(i,1)+gradw(i,3))) &
                            + gradvisc(i,3)*(gradv(i,3)+gradw(i,2))
        f_visc_w = f_visc_w + gradvisc(i,1)*(gradu(i,3)+gradw(i,1)) &
                            + gradvisc(i,2)*(gradv(i,3)+gradw(i,2)) &
                            + gradvisc(i,3)*(fourthirds*gradw(i,3) - twothirds*(gradu(i,1)+gradv(i,2))) 
#endif        
      
        !! Local density 
        tmpro = ro(i)

        !! Body force
        body_force_u = tmpro*grav(1) + driving_force(1)
        body_force_v = tmpro*grav(2) + driving_force(2)
        body_force_w = tmpro*grav(3) + driving_force(3)

        !! Store u.(F_visc + F_body/ro) for use in energy eqn later 
        store_diff_E(i) = store_diff_E(i) + u(i)*(f_visc_u + body_force_u) &
                                          + v(i)*(f_visc_v + body_force_v) &
                                          + w(i)*(f_visc_w + body_force_w)                        
        !! RHS 
        rhs_rou(i) = -tmp_scal_u - gradp(i,1) + body_force_u + f_visc_u 
        rhs_rov(i) = -tmp_scal_v - gradp(i,2) + body_force_v + f_visc_v
#ifdef dim3
        rhs_row(i) = -tmp_scal_w - gradp(i,3) + body_force_w + f_visc_w
#else
        rhs_row(i) = zero
#endif      

     end do
     !$omp end parallel do
         
     !! Make L1,L2,L3,L4,L5 and populate viscous + body force part of rhs' and save transverse
     !! parts of convective terms for later...
     if(nb.ne.0)then
     
        !! Evaluate d2u/dx2,d2v/dx2,d2w/dx2, and d2udxy, d2udxz
        allocate(grad2uvec(nb,3),grad2ucross(nb,2))
        call calc_grad2vecbound(u,v,w,grad2uvec)  
        call calc_grad2crossbound(gradu,grad2ucross)        
   
    
        !$omp parallel do private(i,tmpro,c,xn,yn,un,ut,f_visc_u,f_visc_v,body_force_u,body_force_v &
        !$omp ,dpdn,dundn,dutdn,f_visc_w,tmp_vec,tmp_scal_u,tmp_scal_v,tmp_scal_w)
        do j=1,nb
           i=boundary_list(j)
           tmpro = ro(i)
#ifndef isoT           
           c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
#else
           c=sqrt(csq)
#endif            
           dpdn = gradp(i,1)
           if(node_type(i).eq.0)then !! walls are in bound norm coords
              xn=rnorm(i,1);yn=rnorm(i,2)  !! Bound normals
              un = u(i)*xn + v(i)*yn; ut = -u(i)*yn + v(i)*xn  !! Normal and transverse components of velocity

              dundn = xn*gradu(i,1)+yn*gradv(i,1)
              dutdn = -yn*gradu(i,1)+xn*gradv(i,1)

              L(j,1) = half*(un-c)*(dpdn - tmpro*c*dundn) !! L1 
              L(j,2) = un*(gradro(i,1) - gradp(i,1)/c/c)
              L(j,3) = un*dutdn !! L3 
              L(j,4) = un*(gradw(i,1)*xn + gradw(i,2)*yn) !! L4
              L(j,5) = half*(un+c)*(dpdn + tmpro*c*dundn) !! L5 
              
              !! Should evaluate viscous terms here, but a) rhs_rou,v =0, and b) they aren't needed for
              !! subsequent energy equation 
              rhs_rou(i) = zero
              rhs_rov(i) = zero    
              rhs_row(i) = zero     
              
              !! Don't augment store_diff_E as velocity is zero
!              store_diff_E(i) = store_diff_E(i) + zero 
              
           else    !! In/out is in x-y coord system
              !! Convective terms
              tmp_vec(1) = zero;tmp_vec(2) = v(i);tmp_vec(3) = w(i) !! tmp_vec holds (0,v,w) for node i
              tmp_scal_u = ro(i)*dot_product(tmp_vec,gradu(i,:)) - u(i)*rhs_ro(i) !! convective term for u
              tmp_scal_v = ro(i)*dot_product(tmp_vec,gradv(i,:)) - v(i)*rhs_ro(i) !! Convective term for v   
              tmp_scal_w = ro(i)*dot_product(tmp_vec,gradw(i,:)) - w(i)*rhs_ro(i) !! Convective term for w 


              L(j,1) = half*(u(i)-c)*(dpdn - tmpro*c*gradu(i,1))
              L(j,2) = u(i)*(gradro(i,1) - gradp(i,1)/c/c)
              L(j,3) = u(i)*gradv(i,1)
              L(j,4) = u(i)*gradw(i,1)
              L(j,5) = half*(u(i)+c)*(dpdn + tmpro*c*gradu(i,1))

              !! Build initial stress divergence 
#ifndef isoT              
              f_visc_u = visc(i)*(lapu(i) + onethird*graddivvel(i,1)) 
              f_visc_v = visc(i)*(lapv(i) + onethird*graddivvel(i,2)) 
              f_visc_w = visc(i)*(lapw(i) + onethird*graddivvel(i,3))
#else
              f_visc_u = visc(i)*(lapu(i))  !! Isothermal flows assume small velocity divergence
              f_visc_v = visc(i)*(lapv(i))
              f_visc_w = visc(i)*(lapw(i))
#endif              
#ifndef isoT
              !! non-uniform viscosity terms. 
!              gradvisc(:) = r_temp_dependence*visc(i)*gradT(i,:)/T(i)
              f_visc_u = f_visc_u + gradvisc(i,1)*(fourthirds*gradu(i,1) - twothirds*(gradv(i,2)+gradw(i,3))) &
                                  + gradvisc(i,2)*(gradu(i,2)+gradv(i,1)) &
                                  + gradvisc(i,3)*(gradu(i,3)+gradw(i,1))
              f_visc_v = f_visc_v + gradvisc(i,1)*(gradu(i,2)+gradv(i,1)) &
                                  + gradvisc(i,2)*(fourthirds*gradv(i,2) - twothirds*(gradu(i,1)+gradw(i,3))) &
                                  + gradvisc(i,3)*(gradv(i,3)+gradw(i,2))
              f_visc_w = f_visc_w + gradvisc(i,1)*(gradu(i,3)+gradw(i,1)) &
                                  + gradvisc(i,2)*(gradv(i,3)+gradw(i,2)) &
                                  + gradvisc(i,3)*(fourthirds*gradw(i,3) - twothirds*(gradu(i,1)+gradv(i,2)))     
             
              !! Zero normal flux normal stress divergence
              if(znf_vdiff(j)) then
                 !! dtau_nn/dn = 0
                 f_visc_u = f_visc_u + visc(i)*(twothirds*graddivvel(i,1) - two*grad2uvec(j,1))       
                 f_visc_u = f_visc_u - gradvisc(i,1)*(fourthirds*gradu(i,1) - twothirds*(gradv(i,2)+gradw(i,3)))
              endif

              if(znf_vtdiff(j))then
                 !! dtau_jn/dn = 0 for j=2,3 
!                 f_visc_v = f_visc_v -grad2uvec(j,2) - grad2ucross(j,1)
!                 f_visc_w = f_visc_w -grad2uvec(j,3) - grad2ucross(j,2)

                 !! zero non-uniform viscosity contribution to dtau_jn/dn for j=2,3 
                 f_visc_v = f_visc_v - gradvisc(i,1)*(gradu(i,2)+gradv(i,1)) 
                 f_visc_w = f_visc_w - gradvisc(i,1)*(gradu(i,3)+gradw(i,1)) 
              end if
#endif      
              !! Body force
              body_force_u = tmpro*grav(1) + driving_force(1)
              body_force_v = tmpro*grav(2) + driving_force(2)
              body_force_w = tmpro*grav(3) + driving_force(3)            
                
              !! Viscous dissipation term for energy equation   
              store_diff_E(i) = store_diff_E(i) + u(i)*(f_visc_u + body_force_u) &
                                                + v(i)*(f_visc_v + body_force_v) &
                                                + w(i)*(f_visc_w + body_force_w)

              !! Transverse + visc + source terms only
              rhs_rou(i) = -v(i)*gradu(i,2) - w(i)*gradu(i,3) + f_visc_u + body_force_u  
              rhs_rov(i) = -v(i)*gradv(i,2) - w(i)*gradv(i,3) - gradp(i,2) + f_visc_v + body_force_v
              rhs_row(i) = -v(i)*gradw(i,2) - w(i)*gradw(i,3) - gradp(i,3) + f_visc_w + body_force_w 
                            
           end if
        end do
        !$omp end parallel do 
        deallocate(grad2uvec,grad2ucross)         
     end if
              
         
     !! Deallocate any stores no longer required
     deallocate(lapu,lapv,lapw)
     deallocate(graddivvel) 
#ifndef isoT
     deallocate(gradvisc)
#endif

     return
  end subroutine calc_rhs_rovel
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_rhs_roE
     !! Construct the RHS for energy equation. Do (almost) nothing if isothermal
#ifndef isoT     
     integer(ikind) :: i,j,fla
     real(rkind) :: tmp_conv,tmp_visc,tmp_p,tmp_E
     real(rkind),dimension(:,:),allocatable :: grad2T
     real(rkind) :: lapT_tmp,dutdt,q01,q02
     real(rkind),dimension(:,:),allocatable :: gradlambda  
     real(rkind),dimension(2) :: gradulocal,gradvlocal,gradwlocal   
         
    
     !! Thermal conductivity gradient
     allocate(gradlambda(npfb,dims))
     if(flag_mix_av.eq.1)then
        call calc_gradient(lambda_th,gradlambda)
     else
        !$omp parallel do
        do i=1,npfb
           gradlambda(i,:) = lambda_th(i)*r_temp_dependence*gradT(i,:)/T(i)  + lambda_th(i)*gradcp(i,:)/cp(i)           
        end do
        !$omp end parallel do
     end if
     
     !! Build RHS
     !$omp parallel do private(i,tmp_p,tmp_conv,tmp_visc,tmp_E)
     do j=1,npfb-nb
        i=internal_list(j)

        !! Convection term: -u.grad(roE)
        tmp_conv = -u(i)*gradroE(i,1) - v(i)*gradroE(i,2) - w(i)*gradroE(i,3)

        !! Pressure term: -div.(pu)
        tmp_p = -p(i)*divvel(i) - u(i)*gradp(i,1) - v(i)*gradp(i,2) - w(i)*gradp(i,3) 
        
        !! Energy dilation term: -roE*div.u
        tmp_E = -roE(i)*divvel(i)
     
        !! Viscous heating term: div.(tau u)
        tmp_visc = (fourthirds*gradu(i,1) - twothirds*gradv(i,2) - twothirds*gradw(i,3))*gradu(i,1) &
                 + (fourthirds*gradv(i,2) - twothirds*gradw(i,3) - twothirds*gradu(i,1))*gradv(i,2) &
                 + (fourthirds*gradw(i,3) - twothirds*gradu(i,1) - twothirds*gradv(i,2))*gradw(i,3) &
                 + (gradu(i,2)+gradv(i,1))**two + (gradu(i,3)+gradw(i,1))**two + (gradv(i,3)+gradw(i,2))**two        
        store_diff_E(i) = store_diff_E(i) + visc(i)*tmp_visc

        
        !! Thermal diffusion term: div.(lambda*gradT + heating due to molecular diffusion)
        store_diff_E(i) = store_diff_E(i) + lambda_th(i)*lapT(i)  
        store_diff_E(i) = store_diff_E(i) + dot_product(gradlambda(i,:),gradT(i,:))

                
        !! Build final RHS
        rhs_roE(i) = tmp_conv + tmp_p + tmp_E + store_diff_E(i)

     end do
     !$omp end parallel do
     
     
     !! Populate viscous and transverse parts of rhs
     if(nb.ne.0)then
        allocate(grad2T(nb,3))
        call calc_grad2bound(T,grad2T)          
        !$omp parallel do private(i,tmp_visc,lapT_tmp,gradulocal,gradvlocal,gradwlocal,xn,yn,fla,q01,q02)
        do j=1,nb
           i=boundary_list(j)
              
           !! Viscous energy term: div.(tau u) (N.B. for walls, we need to rotate gradients)
           if(node_type(i).eq.0) then
              xn=rnorm(i,1);yn=rnorm(i,2)           
              gradulocal(1) = xn*gradu(i,1) - yn*gradu(i,2)  !! Convert gradients to x-y-z coordinates
              gradvlocal(1) = xn*gradv(i,1) - yn*gradv(i,2)
              gradwlocal(1) = xn*gradw(i,1) - yn*gradw(i,2)
              gradulocal(2) = yn*gradu(i,1) + xn*gradu(i,2)
              gradvlocal(2) = yn*gradv(i,1) + xn*gradv(i,2)
              gradwlocal(2) = yn*gradw(i,1) + xn*gradw(i,2)           
           
              tmp_visc = (fourthirds*gradulocal(1) - twothirds*gradvlocal(2) - twothirds*gradw(i,3))*gradulocal(1) &
                       + (fourthirds*gradvlocal(2) - twothirds*gradw(i,3) - twothirds*gradulocal(1))*gradvlocal(2) &
                       + (fourthirds*gradw(i,3) - twothirds*gradulocal(1) - twothirds*gradvlocal(2))*gradw(i,3) &
                       + (gradulocal(2)+gradvlocal(1))**two + (gradu(i,3)+gradwlocal(1))**two &
                       + (gradv(i,3)+gradwlocal(2))**two   
           else
              tmp_visc = (fourthirds*gradu(i,1) - twothirds*gradv(i,2) - twothirds*gradw(i,3))*gradu(i,1) &
                       + (fourthirds*gradv(i,2) - twothirds*gradw(i,3) - twothirds*gradu(i,1))*gradv(i,2) &
                       + (fourthirds*gradw(i,3) - twothirds*gradu(i,1) - twothirds*gradv(i,2))*gradw(i,3) &
                       + (gradu(i,2)+gradv(i,1))**two + (gradu(i,3)+gradw(i,1))**two &
                       + (gradv(i,3)+gradw(i,2))**two         
           endif
           store_diff_E(i) = store_diff_E(i) + visc(i)*tmp_visc           

           !! Evaluate temperature Laplacian and zero parts as required
           lapT_tmp = grad2T(j,2) + grad2T(j,3)
           if(node_type(i).eq.0) then
              if(flag_wall_type.eq.0) then !! Adiabatic wall

                 q01 = (-170.0d0*T(i) + 216.0d0*T(i+1) - 54.0d0*T(i+2) + 8.0d0*T(i+3))/(36.0d0*s(i)*s(i)*L_char*L_char)
                 lapT_tmp = lapT_tmp + q01 

                 gradlambda(i,1) = zero
              else   !! Isothermal wall
                 lapT_tmp = lapT_tmp + grad2T(j,1)
              end if            
           else  !! Inflow/outflow
               if(znf_tdiff(j))then
                  gradlambda(i,1) = zero
                  !! Don't add normal term
               else
                  !! Add normal term
                  lapT_tmp = lapT_tmp + grad2T(j,1) 
               endif
           end if                     

           
           !! Add thermal diffusion term: div.(lambda*gradT)              
           store_diff_E(i) = store_diff_E(i) + lambda_th(i)*lapT_tmp + dot_product(gradlambda(i,:),gradT(i,:))
           !! RHS (visc + cond + transverse + source)
           if(node_type(i).eq.0) then          
              !! u,v,w=0 on wall, so velocity gradients are also zero along wall. Only diffusion terms.
              rhs_roE(i) = store_diff_E(i) 
           else                      
              rhs_roE(i) = - v(i)*gradroE(i,2) - roE(i)*gradv(i,2) - w(i)*gradroE(i,3) - roE(i)*gradw(i,3) &
                         - p(i)*gradv(i,2) - v(i)*gradp(i,2) - p(i)*gradw(i,3) - w(i)*gradp(i,3) &
                         + store_diff_E(i) 
           end if                                                               
                                                               
        end do
        !$omp end parallel do
        deallocate(grad2T)
     end if
     deallocate(gradcp)
     deallocate(gradlambda)                    
#else
     rhs_roE(1:npfb) = zero
#endif

     !! Deallocation of spatial derivatives
     deallocate(store_diff_E)
           
     return
  end subroutine calc_rhs_roE
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_nscbc
    !! This routine asks boundaries module to prescribe L as required, then builds the final 
    !! rhs for each equation. It should only be called if nb.ne.0
    integer(ikind) :: i,j,ispec
    real(rkind) :: tmpro,c,tmp_scal,cv,gammagasm1,enthalpy
           
    segment_tstart = omp_get_wtime()           
           
    !! Loop over boundary nodes and specify L as required
    !$omp parallel do private(i)
    do j=1,nb
       i=boundary_list(j)  
       
       !! WALL BOUNDARY 
       if(node_type(i).eq.0) then  
          if(flag_wall_type.eq.1) then    
             call specify_characteristics_isothermal_wall(j,L(j,:))
          else
             call specify_characteristics_adiabatic_wall(j,L(j,:))
          end if
          
       !! INFLOW BOUNDARY
       else if(node_type(i).eq.1) then 
          if(flag_inflow_type.eq.1) then  !! Hard inflow
             call specify_characteristics_hard_inflow(j,L(j,:))       
          else if(flag_inflow_type.eq.2) then !! Pressure-tracking Inflow-outflow
             call specify_characteristics_inflow_outflow(j,L(j,:),gradro(i,:),gradp(i,:),gradu(i,:),gradv(i,:),gradw(i,:))
          else           
             call specify_characteristics_soft_inflow(j,L(j,:),gradro(i,:),gradp(i,:),gradu(i,:),gradv(i,:),gradw(i,:))       
          endif
 
       !! OUTFLOW BOUNDARY 
       else if(node_type(i).eq.2) then   
          call specify_characteristics_outflow(j,L(j,:),gradro(i,:),gradp(i,:),gradu(i,:),gradv(i,:),gradw(i,:))       
       end if          
    end do
    !$omp end parallel do
    
    !! ==================================================================================
    !! Use L to update the rhs on boundary nodes
    !$omp parallel do private(i,tmpro,c,tmp_scal,cv,ispec,gammagasm1,enthalpy)
    do j=1,nb
       i=boundary_list(j)
       tmpro = ro(i)
#ifndef isoT       
       c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
       gammagasm1 = Rgas_mix(i)/(cp(i)-Rgas_mix(i))
#else
       c=sqrt(csq)
#endif       

       !! This quantity appears in rhs_ro and rhs_roE, so save in tmp_scal
       tmp_scal = (c*c*L(j,2) + L(j,5) + L(j,1))/c/c
       
       !! RHS for density logarithm
       rhs_ro(i) = rhs_ro(i) - tmp_scal

       !! Velocity components       
       rhs_rou(i) = rhs_rou(i) - u(i)*tmp_scal - (L(j,5)-L(j,1))/c
       rhs_rov(i) = rhs_rov(i) - v(i)*tmp_scal - ro(i)*L(j,3)
       rhs_row(i) = rhs_row(i) - w(i)*tmp_scal - ro(i)*L(j,4)

       !! Energy
#ifndef isoT
       cv = cp(i) - Rgas_mix(i)
       rhs_roE(i) = rhs_roE(i) - u(i)*(L(j,5)-L(j,1))/c - tmpro*v(i)*L(j,3) -tmpro*w(i)*L(j,4) &
                               - (roE(i)/tmpro - cv*T(i))*tmp_scal - (L(j,5)+L(j,1))/gammagasm1
#endif           
       
       !! Species mass fractions
       do ispec=1,nspec
          rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) - (Yspec(i,ispec)/tmpro)*tmp_scal - ro(i)*L(j,5+ispec) 
#ifndef isoT
          !! Evaluate enthalpy
          call evaluate_enthalpy_only_at_node(T(i),ispec,enthalpy)            

          !! Reduced enthalpy
          enthalpy = enthalpy - cp(i)*T(i)*Rgas_universal*one_over_molar_mass(ispec)/Rgas_mix(i)

          rhs_roE(i) = rhs_roE(i) - tmpro*enthalpy*L(j,5+ispec)
#endif          
       end do
    end do
    !$omp end parallel do
    !! De-allocate L
    deallocate(L)
        
    !! Profiling
    segment_tend = omp_get_wtime()
    segment_time_local(2) = segment_time_local(2) + segment_tend - segment_tstart    

    return  
  end subroutine calc_rhs_nscbc
!! ------------------------------------------------------------------------------------------------  
  subroutine filter_variables
     !! This routine calls the specific filtering routine (within derivatives module) for each
     !! variable - ro,u,v,w,roE,Yspec - and forces certain values on boundaries as required.
     integer(ikind) :: ispec,i
     real(rkind),dimension(:),allocatable :: filter_correction
      
     segment_tstart = omp_get_wtime()
               
     !! Filter density
     call calc_filtered_var(ro)
          
     !! Filter velocity components
     call calc_filtered_var(rou)
     call calc_filtered_var(rov)
#ifdef dim3
     call calc_filtered_var(row)
#endif     
 
     !! Filter energy
#ifndef isoT     
     call calc_filtered_var(roE)       
#endif     

     !! Filter mass fractions
     !! Initialise filter correction term
     allocate(filter_correction(np))
     filter_correction(1:npfb) = zero

     !! Loop over species
     do ispec=1,nspec
     
        !! Filter this species roY
        call calc_filtered_var(Yspec(:,ispec))
 
        !! Add this species contribution to filter correction term
        filter_correction(1:npfb) = filter_correction(1:npfb) + Yspec(1:npfb,ispec)
     end do
    
    
     !! Apply filter correction term to each species
     do ispec=1,nspec
        !$omp parallel do
        do i=1,npfb
           Yspec(i,ispec) = Yspec(i,ispec)*ro(i)/filter_correction(i)
        end do
        !$omp end parallel do
     end do
     
     !! Deallocate space
     deallocate(filter_correction)     
   
     
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(3) = segment_time_local(3) + segment_tend - segment_tstart
     
     return
  end subroutine filter_variables  
!! ------------------------------------------------------------------------------------------------    
end module rhs
