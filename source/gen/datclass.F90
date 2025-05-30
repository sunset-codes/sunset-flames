program datgen
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! DATGEN program to generate node sets for input
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  use kind_parameters
  use common_parameter
  use global_variables 
  implicit none

  !! xbcond/ybcond = 0,1,2,3 for none, periodic or symmetric, no-slip BCs respectively.
  !! btype = 0,1,2,3 for wall, inflow, outflow or periodic/symmetric respectively

  real(rkind) :: x,y

  integer ipart
  integer i,j,icha,nn,ve_model,ii,jj
  double precision h0,r_mag,yl,D_cyl,S_cyl
  integer xbcond_L,xbcond_U,ybcond_L,ybcond_U
  
  double precision :: a0,a1,a2,a3,a4,a5 !! NACA coefficients
  double precision :: temp,tmp2,tmp
  real(rkind),dimension(2) :: tmpN,rn
  integer(ikind) :: iround,nround,nsearch,ipartrow,ipartrowm1,nnear,npdps,idown,iup,ibc
  logical :: keepgoing,keepgoing2,skip
  real(rkind),dimension(:,:),allocatable :: pdp
  real(rkind) :: dist2bound,dxtmp,minpdp,dist2io,dist2iow
  real(rkind),dimension(:),allocatable :: tmpX
  real(rkind),dimension(:),allocatable :: pdp_x,pdp_y,pdp_dist2
  logical,dimension(:),allocatable :: pdp_active
  integer(ikind),dimension(:),allocatable :: pdp_freeindices
  real(rkind),dimension(5) :: pdp_x_new,pdp_y_new
  real(rkind) :: thup,thdown,th_increment
  integer(ikind) :: npdps_new,inew,block_left,block_right,block_new,block_delete,block_end
  double precision, parameter :: r3o2 = 0.5d0*sqrt(3.0d0)
  integer(ikind) :: nbx,nby,nbtot
  

  write(*,*) 'Cases: '
  write(*,*) '  case 1:  Currently empty'
  write(*,*) '  case 2:  Periodic channel'
  write(*,*) '  case 3:  unit torus'
  write(*,*) '  case 4:  Rayleigh-Taylor'
  write(*,*) '  case 5:  Simple flame tube'  
  write(*,*) '  case 6:  Bluff body stabilised flame (from King (2024) CAME)'    
  write(*,*) '  case 7:  Variable resolution demonstration...'
  write(*,*) '  case 8:  Cylinder'  
  write(*,*) '  case 9:  Cylinder array'    
  write(*,*) '  '
  write(*,*) 'Input test case number: '
  read(*,*) itest


  select case (itest) 
!! ------------------------------------------------------------------------------------------------
  case(1) !! EMPTY
!! ------------------------------------------------------------------------------------------------
  case(2) !! Poiseuille flow

     yl=1.0d0
     xl=yl/2.0d0
     dx0=yl/200.0d0
     xbcond_L=1;xbcond_U=1;ybcond_L=0;ybcond_U=0
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 0, 3, 0, 3/)  
     b_node(1,:) = (/-0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/-0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 0;n_blob_coefs=0

     dxmin = dx0/1.0d0
     dx_wall=dxmin;dx_in=1.0d0*dx0;dx_out=dx_in;dx_wallio=dx_in !! dx for solids and in/outs...!! Ratio for scaling far field...  
!! ------------------------------------------------------------------------------------------------
  case(3) !! Kolmogorov flow

     yl=2.0d0*pi
     xl=yl
     dx0=yl/50.0d0
     xbcond_L=1;xbcond_U=1;ybcond_L=1;ybcond_U=1
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 3, 3, 3/)  
     b_node(1,:) = (/-0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/-0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 0;n_blob_coefs=6
!     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,n_blob_coefs),blob_rotation(nb_blobs))
!     blob_centre(1,:)=(/-1.d0,0.d0/); !! Central
!     do i=1,nb_blobs
!        blob_coeffs(i,:)=(/1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/);blob_rotation(i)=0.4d0
!     end do     

     dxmin = dx0/1.0d0
     dx_wall=dxmin;dx_in=1.0d0*dx0;dx_out=dx_in;dx_wallio=dx_in  !! dx for solids and in/outs.
!! ------------------------------------------------------------------------------------------------
case(4) !! Rayleigh-Taylor geometry

     yl=0.25d0!0.0125d0  ! channel width
     xl=1.0d0 ! channel length
     dx0=xl/300.0       !15
     h0 = xl/10.0d0
     xbcond_L=0;xbcond_U=0;ybcond_L=1;ybcond_U=1
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 2, 3, 0/)  
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 1;n_blob_coefs=6
     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,n_blob_coefs),blob_rotation(nb_blobs))
     blob_centre(1,:)=(/0.d0,0.d0/); !! Central
     do i=1,nb_blobs
        blob_coeffs(i,:)=h0*(/1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/);blob_rotation(i)=0.0d0
     end do

     dxmin = dx0/2.0d0
     dx_wall=dxmin;dx_in=1.5d0*dx0;dx_out=dx0*1.5d0;dx_wallio=dx_in  !! dx for solids and in/outs..

!! ------------------------------------------------------------------------------------------------
case(5) !! Inflow/outflow tube for simple flames

     yl=0.055d0!0.0125d0  ! channel width
     xl=1.0d0 ! channel length
     dx0=xl/250.0       !15
     xbcond_L=0;xbcond_U=0;ybcond_L=1;ybcond_U=1
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 2, 3, 1/)  
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 0;n_blob_coefs=6
!     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,n_blob_coefs),blob_rotation(nb_blobs))
!     blob_centre(1,:)=(/0.d0,0.d0/); !! Central
!     do i=1,nb_blobs
!        blob_coeffs(i,:)=h0*(/1.0d0,0.4d0,0.0d0,0.0d0,0.0d0,0.0d0/);blob_rotation(i)=-pi/9.0d0
!     end do

     dxmin = dx0/2.0d0
     dx_wall=dxmin;dx_in=1.0d0*dx0;dx_out=dx0*1.0d0;dx_wallio=dx_in  !! dx for solids and in/outs..

     
!! ------------------------------------------------------------------------------------------------
case(6) !! Bluff body stabilised flame (from CMAME 2024 paper)

     xl=1.0d0 ! channel length
     h0=xl/40.0d0   !cylinder radius
     yl=xl/2.0d0!/10.0d0!(4.0d0/3.0d0)  ! channel width
     dx0=xl/(40.0d0*25.0d0)!25.0       !15
     xbcond_L=0;xbcond_U=0;ybcond_L=1;ybcond_U=1
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 2, 3, 1/)  
     b_node(1,:) = (/ -0.50d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.50d0*xl, 0.5d0*yl /)
     nb_blobs=1
     open(unit=191,file="blob_fcoefs.in")
     read(191,*) n_blob_coefs
     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,n_blob_coefs),blob_rotation(nb_blobs))
     do i=1,n_blob_coefs
        read(191,*) blob_coeffs(1,i)
     end do
     close(191)
!     blob_coeffs(1,:)=0.0d0;blob_coeffs(1,1)=1.0d0
     blob_coeffs(1,:) = blob_coeffs(1,:)*h0;blob_rotation(1)=-0.0d0*pi
!     blob_coeffs(2,:) = blob_coeffs(1,:);blob_rotation(2)=-0.0d0*pi
!     blob_coeffs(3,:) = blob_coeffs(1,:);blob_rotation(3)=-0.0d0*pi    

     blob_centre(1,:)=(/ -0.275d0*xl, 0.0d0*yl/);
!     blob_centre(2,:)=(/ -0.275d0*xl,-0.5d0*yl/);
!     blob_centre(3,:)=(/ -0.275d0*xl, 0.5d0*yl/);     


     dxmin = dx0/2.0d0
     dx_wall=dxmin;dx_in=4.0d0*dx0;dx_out=1.5d0*dx0;dx_wallio=dx_in  !! dx for solids and in/outs...!! 
!! ------------------------------------------------------------------------------------------------
case(7) !! Variable resolution demonstration case

     h0=1.0d0 ! channel length
     xl=60.0d0*h0   !cylinder radius
     yl=xl/10.0d0!/10.0d0!(4.0d0/3.0d0)  ! channel width
     dx0=h0/(20.0d0)!25.0       !15
     xbcond_L=0;xbcond_U=0;ybcond_L=1;ybcond_U=1
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 2, 3, 1/)  
     b_node(1,:) = (/ -0.50d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.50d0*xl, 0.5d0*yl /)
     nb_blobs=1
     open(unit=191,file="blob_fcoefs.in")
     read(191,*) n_blob_coefs
     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,n_blob_coefs),blob_rotation(nb_blobs))
     do i=1,n_blob_coefs
        read(191,*) blob_coeffs(1,i)
     end do
     close(191)
     blob_coeffs(1,:)=0.0d0;blob_coeffs(1,1)=1.0d0
     blob_coeffs(1,:) = blob_coeffs(1,:)*h0;blob_rotation(1)=-0.0d0*pi
!     blob_coeffs(2,:) = blob_coeffs(1,:);blob_rotation(2)=-0.0d0*pi
!     blob_coeffs(3,:) = blob_coeffs(1,:);blob_rotation(3)=-0.0d0*pi    

     blob_centre(1,:)=(/ -0.0d0*xl, 0.0d0*yl/);
!     blob_centre(2,:)=(/ -0.275d0*xl,-0.5d0*yl/);
!     blob_centre(3,:)=(/ -0.275d0*xl, 0.5d0*yl/);     


     dxmin = dx0/3.0d0
     dx_wall=dxmin;dx_in=6.0d0*dx0;dx_out=dx_in;dx_wallio=dx_in  !! dx for solids and in/outs...!!   
!! ------------------------------------------------------------------------------------------------
case(8) !! Arrays of cylinders for lean H2 flame dynamics tests

     D_cyl = 1.0d0;h0 = 0.5d0*D_cyl  !! Cylinder diameter (unity)
     S_cyl = 2.0d0*D_cyl             !! Cylinder spacing (multiples of D_cyl)
     xl = 20.0d0*D_cyl              !! Channel length
     yl = 1.0d0*S_cyl                      !! Channel width 
     dx0 = D_cyl/50                  !! Baseline resolution
     xbcond_L=0;xbcond_U=0;ybcond_L=1;ybcond_U=1
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 2, 3, 1/)  
     b_node(1,:) = (/ -5.0d0*D_cyl, -0.5d0*yl /)
     b_node(2,:) = (/ -5.0d0*D_cyl + xl, -0.5d0*yl /)
     b_node(3,:) = (/ -5.0d0*D_cyl + xl, 0.5d0*yl /)
     b_node(4,:) = (/ -5.0d0*D_cyl, 0.5d0*yl /)
     nb_blobs=1
     open(unit=191,file="blob_fcoefs.in")
     read(191,*) n_blob_coefs
     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,n_blob_coefs),blob_rotation(nb_blobs))
     do i=1,n_blob_coefs
        read(191,*) blob_coeffs(1,i)
     end do
     close(191)
     blob_coeffs(1,:)=0.0d0;blob_coeffs(1,1)=1.0d0
     blob_coeffs(1,:) = blob_coeffs(1,:)*h0;blob_rotation(1)=-0.0d0*pi
 
     blob_centre(1,:)=(/ 0.0d0, 0.0d0/)  !! Row 1
!     blob_centre(2,:)=(/ 0.0d0,-S_cyl/)
!     blob_centre(3,:)=(/ 0.0d0, S_cyl/)
     
!     blob_centre(4,:)=(/ r3o2*S_cyl,0.5d0*S_cyl/)   !! Row 2  
!     blob_centre(5,:)=(/ r3o2*S_cyl,-0.5d0*S_cyl/)     

!     blob_centre(6,:)=(/ 2.0d0*r3o2*S_cyl,0.0d0/)  !! Row 3   
!     blob_centre(7,:)=(/ 2.0d0*r3o2*S_cyl, S_cyl/)     
!     blob_centre(8,:)=(/ 2.0d0*r3o2*S_cyl, -S_cyl/)               

     !! Multiple blobs, copy coefficients and orientations from blob 1
     if(nb_blobs.gt.1)then
        do i=2,nb_blobs
           blob_coeffs(i,:) = blob_coeffs(1,:);blob_rotation(i)=blob_rotation(1)        
        end do
     end if


     dxmin = dx0/2.0d0
     dx_wall=dxmin;dx_in=4.0d0*dx0;dx_out=1.5d0*dx0;dx_wallio=dx_in  !! dx for solids and in/outs...!!     
!! ------------------------------------------------------------------------------------------------
case(9) !! Porous cylinder array

!! JK NOTES TO HB
!! Increase nby to increase cross-stream domain size
!! Increase xl, and adjust the "-5.0d0*D_cyl..." in the definitions of b_node( to make the domain longer
!! For detailed chemistry, likely need slightly finer resolution.

     nbx = 4  !! Number of cylinders along x
     nby = 1   !! Number of cylinders in y    (I think nbx and nby should be even, but it might be okay if not).
     nbtot = nbx*nby + nbx/2

     D_cyl = 1.0d0;h0 = 0.5d0*D_cyl  !! Cylinder diameter (unity)
     S_cyl = (10.0d0/6.0d0)*D_cyl  !1.25           !! Cylinder spacing (multiples of D_cyl)
     yl = dble(max(nby,1))*S_cyl                      !! Channel width 
     xl = 10.0d0*D_cyl !3.0d0*S_cyl + dble(nbx)*r3o2*S_cyl          !! Channel length     
     dx0 = D_cyl/75.0d0                  !! Baseline resolution
     xbcond_L=0;xbcond_U=0;ybcond_L=1;ybcond_U=1
     
!     yl=yl/2.0
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 2, 3, 1/)  
     b_node(1,:) = (/ -5.0d0*D_cyl, -0.5d0*yl /)   !-12
     b_node(2,:) = (/ -5.0d0*D_cyl+xl, -0.5d0*yl /)
     b_node(3,:) = (/ -5.0d0*D_cyl+xl, 0.5d0*yl /)
     b_node(4,:) = (/ -5.0d0*D_cyl, 0.5d0*yl /)
     nb_blobs=nbtot
     if(nb_blobs.ne.0) then
        open(unit=191,file="blob_fcoefs.in")
        read(191,*) n_blob_coefs
        allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,n_blob_coefs),blob_rotation(nb_blobs))
        do i=1,n_blob_coefs
           read(191,*) blob_coeffs(1,i)
        end do
        close(191)
        blob_coeffs(1,:)=0.0d0;blob_coeffs(1,1)=1.0d0
        blob_coeffs(1,:) = blob_coeffs(1,:)*h0;blob_rotation(1)=-0.0d0*pi

        ii = 0 
        do j=1,nbx
 
           if(mod(j,2).eq.0) then !! Even columns
  
              do i=1,nby+1
                 ii = ii+1
                 blob_centre(ii,:) = (/ -(dble(j-1))*r3o2*S_cyl,S_cyl*(dble(nby)/2.0d0 - dble(i)+1.0d0) /)    
              end do   
   
           else !! Odd columns
              do i=1,nby
                 ii = ii+1
                 blob_centre(ii,:) = (/ -dble(j-1)*r3o2*S_cyl,S_cyl*(0.5d0+dble(nby)/2.0d0 - dble(i)) /)
              end do    
           end if
        
        end do

        !! Multiple blobs, copy coefficients and orientations from blob 1
        if(nb_blobs.gt.1)then
           do i=2,nb_blobs
              blob_coeffs(i,:) = blob_coeffs(1,:);blob_rotation(i)=blob_rotation(1)        
           end do
        end if
     else
        n_blob_coefs=6
     end if     

     dxmin = dx0/1.7d0
     dx_wall=dxmin;dx_in=1.0d0*dx0;dx_out=4.0d0*dx0;dx_wallio=dx_in  !! dx for solids and in/outs...!!       
     
!! ------------------------------------------------------------------------------------------------     
end select
!! ------------------------------------------------------------------------------------------------     
     
     !! Create the domain
     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     !! Create the boundary nodes
     call make_boundary_particles
     call make_boundary_blobs               
     ipart = nb   
         
     !! Initialise a line of potential dot points   
     nsearch = ceiling(yb_max-yb_min)/dxmin/2.0d0
     allocate(pdp_x(10*nsearch),pdp_y(10*nsearch))
     npdps = nsearch
     y=yb_min-0.5d0*dxmin
     i=0
     do while (y.lt.yb_max)
        y = y + dxmin/2.0d0;i=i+1
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_x(i) = xb_min + temp
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_y(i) = y 
     end do          
     npdps=i         
     allocate(pdp_dist2(10*nsearch));pdp_dist2=0.0d0   
        
     minpdp = minval(pdp_x(1:npdps))
     do while (minpdp.le.xb_max)  !! Keep going until all PDPs have passed out the end of the domain
        
        !! Pick the left-most pdp
        j=minloc(pdp_x(1:npdps),DIM=1)
        keepgoing = .true.
        x=pdp_x(j);y=pdp_y(j)
        pdp_dist2(1:npdps) = (x-pdp_x(1:npdps))**2.0d0 + (y-pdp_y(1:npdps))**2.0d0
        
        !! How far are we from the boundaries?
        dist2bound = xb_max-xb_min + 100.0d0
        dist2io = xb_max-xb_min + 100.0d0
        dist2iow = xb_max-xb_min + 100.0d0
        i=1
        do while(i.le.nb)
           tmpN(1) = x - xp(i);tmpN(2) = y - yp(i);temp = sqrt(dot_product(tmpN,tmpN))
           if(i.gt.nbio.and.temp.le.dist2bound) dist2bound = temp
           if(i.le.nbio) then
              if(node_type(i).eq.0) then  !! Wall
                 if(temp/dxp(i).le.dist2iow) dist2iow = temp/dxp(i)
              else   !! In-out
                 if(temp/dxp(i).le.dist2io) dist2io = temp/dxp(i)
              end if
           end if
           i=i+1
           if(dist2bound.le.4.25*dxp(i-1)) keepgoing = .false. !! Too close to solid bound, leave it.  !!NEWBC
        end do     

        if(dist2io.le.4.25) keepgoing = .false. !! Too close to io bound, leave it.  !!NEWBC              
        if(dist2iow.le.4.25) keepgoing = .false.
     
        !! Calculate the resolution locally
        call get_resolution(x,y,dist2bound,dx)
                     
        !! Check whether to place a particle here, based on some criteria
        !! --------------------------------------------------------------
           !! Are we within the object!!?!?!
           do i=1,nb_blobs
              temp = sqrt((x-blob_centre(i,1))**2. + (y-blob_centre(i,2))**2.)
              if(x-blob_centre(i,1).ge.0.0d0)then
                 tmp2 = asin((y-blob_centre(i,2))/temp)
              else
                 tmp2 = pi-asin((y-blob_centre(i,2))/temp)
              endif              
              tmp2 = tmp2 - blob_rotation(i)

              r_mag = blob_coeffs(i,1)
              do ibc = 2,n_blob_coefs
                 r_mag = r_mag + blob_coeffs(i,ibc)*cos(dble(ibc-1)*tmp2) 
              end do
                    
              if(temp.le.r_mag)then
                 keepgoing = .false.
              end if
           end do
           
           !! Are we too close to a boundary?
           if(x-0.5d0*dx.le.xb_min.or.x+0.5d0*dx.ge.xb_max.or.y-0.5d0*dx.le.yb_min.or.y+0.5d0*dx.ge.yb_max)then
              keepgoing = .false.
           end if            
                                  

        !! END CRITERIA
        !! --------------------------------------------------------------

        !! Place a particle here
        if(keepgoing) then
           ipart = ipart + 1
           call random_number(temp);temp = temp -0.5d0;xp(ipart) = pdp_x(j) + temp*dxmin*0.0d0
           call random_number(temp);temp = temp -0.5d0;yp(ipart) = pdp_y(j) + temp*dxmin*0.0d0
           dxp(ipart) = dx  
           if(dist2bound.le.4.5d0*dx) node_type(ipart)=998
!           if(dist2bound.le.3.0d0*dx) node_type(ipart)=997
        end if           
        
        !! Deactive all pdps within dx of this pdp
        !! Search down
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i-1;
           if(i.eq.0) then
              idown = 0
              thdown = -0.49999999d0*pi
              keepgoing = .false.
           else if(pdp_dist2(i).ge.dx*dx*(0.8725**2)) then      !! Distance??        
              idown = i
              thdown = atan2((pdp_y(idown)-y),(pdp_x(idown)-x))
              keepgoing = .false.    
 
           end if       
        end do
        !! Search up
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i+1
           if(i.eq.npdps+1) then
              iup = npdps+1
              thup = 0.4999999999d0*pi
              keepgoing = .false.              
           else if(pdp_dist2(i).ge.dx*dx*(0.8725**2)) then !! Distance
              iup = i
              thup = atan2((pdp_y(iup)-y),(pdp_x(iup)-x))              
              keepgoing = .false.
           
           end if
        end do     
        
        !! Temporary store for the new pdps
        th_increment = (thup-thdown)/5.0d0
        inew = 0
        do i=1,5
           temp = y + dx*sin(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)    
           if(temp.ge.yb_min.and.temp.le.yb_max) then
              inew = inew + 1
              pdp_x_new(inew) = x + dx*cos(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)
              pdp_y_new(inew) = temp
           end if            
        end do
 

        block_new = idown + 1
        block_right = block_new + inew
        block_delete = iup - idown - 1
        npdps_new = npdps + inew - block_delete

        !! Shunt indices above pdp of interest
        if(block_delete.gt.inew) then !! Shunt to LEFT
           do i=block_right,npdps_new
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do
        end if
        if(block_delete.lt.inew) then !! Shunt to RIGHT
           do i=npdps_new,block_right,-1
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do

        end if
                  
        !! Insert any new pdps
        if(inew.ne.0) then          
           pdp_x(block_new:block_new+inew-1) = pdp_x_new(1:inew)
           pdp_y(block_new:block_new+inew-1) = pdp_y_new(1:inew)        
        end if       

        npdps = npdps_new
                  
                
        !! How left is the left-most PDP?                 
        minpdp = minval(pdp_x(1:npdps))
                
        write(6,*) ipart,(minval(pdp_x(1:npdps))-xb_min)/(xb_max-xb_min)
     end do  
                         
                                                  
                        
     npfb = ipart
     dx0 = maxval(dxp(1:npfb))
     dxmin = minval(dxp(1:npfb))
     
     write(*,*) 'nb,npfb,nbio= ', nb,npfb,nbio

!! ------------------------------------------------------------------------------------------------
!! Re-order nodes (from left to right)
  
   write(6,*) "About to quicksort"
   call quicksort(xp,1,npfb)
   write(6,*) "Quicksorted nodes ordered increasing x"

!! ------------------------------------------------------------------------------------------------
  !      ** write data out **
  ! 
  
  !! Write to fort file for quick visualisation
!!  write(31,*) npfb
!  do i=1,npfb
!     write(31,*) xp(i),yp(i),dxp(i)
!  end do
!  !! Use Octave/Matlab and run  "A=load('fort.31');scatter(A(:,1),A(:,2),10,A(:,3),'filled');colorbar" 
  !! from this directory for instant visualisation.
  
  !!
  open(13,file='./IPART')
  write(13,*) nb,npfb,dx0,dxmin
  write(13,*) xb_min,xb_max,yb_min,yb_max
  write(13,*) xbcond_L,xbcond_U,ybcond_L,ybcond_U
  do i=1,npfb
     if(node_type(i).ge.0.and.node_type(i).le.2) then
        write(13,*) i,xp(i), yp(i),node_type(i),xnorm(i),ynorm(i),dxp(i)
     else
        write(13,*) i,xp(i), yp(i),node_type(i),0.0d0,0.0d0,dxp(i)
     end if
  end do
  close(13)
  !! we can safely deallocate here
  deallocate(xp, yp,dxp,node_type,xnorm,ynorm)

  deallocate(b_node,b_edge)
  deallocate(b_type)
     ! end boundary

     write(*,*) 'END of DATCLASS'
200  stop
   end program datgen

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
recursive subroutine quicksort(a, first, last)
  implicit none
  double precision  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     call swap_nodes(i,j)
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort
!! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine swap_nodes(i,j)
     use kind_parameters
     use common_parameter
     use global_variables
     implicit none  
     integer :: i,j
     double precision :: tmp
     integer :: itmp

     !! xp is already swapped by sub-routine quicksort     
     tmp = yp(j);yp(j)=yp(i);yp(i)=tmp
     tmp = xnorm(j);xnorm(j)=xnorm(i);xnorm(i)=tmp
     tmp = ynorm(j);ynorm(j)=ynorm(i);ynorm(i)=tmp
     tmp = dxp(j);dxp(j)=dxp(i);dxp(i)=tmp
     itmp = node_type(j);node_type(j)=node_type(i);node_type(i)=itmp                    
     

     return
  end subroutine swap_nodes
!! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine make_boundary_edge_vectors
     use kind_parameters
     use common_parameter
     use global_variables
     implicit none
     integer(ikind) ib,ibp1

     do ib = 1,nb_patches ! loop over all boundary patches
        ibp1 = mod(ib,nb_patches) + 1   
        b_edge(ib,:) = b_node(ibp1,:) - b_node(ib,:)  ! calculate b_edge
     end do

     return 
   end subroutine make_boundary_edge_vectors
!! ------------------------------------------------------------------------------------------------
   subroutine make_boundary_particles
     use kind_parameters
     use common_parameter
     use global_variables 

     implicit none
     integer(ikind) :: ipart,ib,ibm1
     real(rkind) :: x,y,m_be,tmp,tmp2,dx_local,temp
     integer(ikind) :: nround,iround
     real(rkind),dimension(2) :: nrm
   
     ipart=0

     !! we only allocate memory when we need
     allocate(xp(npar), yp(npar))
     allocate(xnorm(npar),ynorm(npar))
     allocate(dxp(npar))
     allocate(node_type(npar));node_type=999     

     !! Wall/inflow/outflow particles
     do ib=1,nb_patches  ! loop over all boundary patches
        ibm1 = mod(ib+nb_patches-2,nb_patches)+1
        if(abs(b_type(ib)).ne.3)then  ! if it is a wall, inflow or outflow patch
           m_be = dsqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
           nrm(1)=-b_edge(ib,2)/m_be;nrm(2)=b_edge(ib,1)/m_be
           if(b_type(ib).eq.1) dxio = dx_in
           if(b_type(ib).eq.2) then
              dxio = dx_out
              call get_resolution(xb_max,zero,1.0d10,dxio)
           end if
           if(b_type(ib).eq.0) dxio = dx_wallio
           tmp = 0.5d0*dxio/m_be
           do while(tmp.lt.1.0-1.0d-10)   ! move along the patch in increments of dx
              ipart = ipart + 1
              x = b_node(ib,1) + tmp*b_edge(ib,1)
              y = b_node(ib,2) + tmp*b_edge(ib,2)
              xp(ipart) = x;yp(ipart) = y
              
           
              !! Hard-coding local dx              
              if(b_type(ib).ne.0) then  !! Inflow or outflow
                 tmp2 = sqrt(0.27d0**2.0d0 + y*y)
                 call get_resolution(x,y,1.0d10,dx_local)     
              else
                 dx_local = dxio
              end if
              
              
              tmp = tmp + dx_local/m_be  ! note, in future we should allow for dx.ne.dy here
              xnorm(ipart) = nrm(1);ynorm(ipart) = nrm(2)
              dxp(ipart)=dx_local
              node_type(ipart) = b_type(ib) !! Set node-type so we know if it's a wall, inflow or outflow...
           end do
        end if
     end do
     nbio = ipart
          
     nb=ipart    
!!
     write(6,*) 'no. of solid boundary particles: ',nb
     return
   end subroutine make_boundary_particles
!! ------------------------------------------------------------------------------------------------
   subroutine get_resolution(x,y,bdist,dx_local)
     use kind_parameters
     use common_parameter
     use global_variables 
     implicit none

     real(rkind),intent(in) :: x,y
     real(rkind),intent(in) :: bdist
     real(rkind),intent(out) :: dx_local
     real(rkind) :: temp,tmp2,r_mag,d2b_local
     real(rkind) :: xhat,yhat

     !! Set dxio to match the end of the domain we're in
     if(x.le.0.0d0) then
        dxio = dx_in
     else
        dxio = dx_out
     end if
     
     !! Copy bdist to local
     d2b_local = bdist


     if(itest.eq.6) then
        xhat = x - blob_centre(1,1)
        yhat = y - blob_centre(1,2)
        !! Stretch high-res region downstream of flameholder        
        if((x-blob_centre(1,1)).gt.0.0d0) then
     
           r_mag = ((xb_max - x)/(xb_max - blob_centre(1,1)))**2.0d0  !! Scale between blob and outlet (0=outlet)

           temp = exp(-(8.0d0*yhat)**4.0d0) !! Blob-side spreading function
           tmp2 = exp(-(4.0d0*yhat)**4.0d0) !! Outflow-side spreading function
           temp = r_mag*temp + (1.0d0-r_mag)*tmp2 !! Linear variation between blob-side and outflow-side
           
           dxio = dx_out + (dx_in - dx_out)*(1.0d0-temp)

        else
           tmp2 = ((blob_centre(1,1)-x)/(blob_centre(1,1)-xb_min))**2.0d0
           d2b_local = d2b_local*(1.5d0*tmp2 + 1.0d0*(1.0d0-tmp2))
           r_mag = sqrt(xhat**2.0d0 + yhat**2.0d0)
           temp = exp(-(8.0d0*r_mag)**4.0d0)
           dxio = dx_out + (dx_in - dx_out)*(1.0d0-temp)
        endif    

     else if(itest.eq.5) then

        !! Over-ride object tests
        temp = 0.05 !! size of refined region
        tmp2 = -0.0d0 !! location of refined region
        tmp2 = x - tmp2 !! Location relative to finest resolution centre
        if(abs(tmp2).le.temp) then
           d2b_local=0.0d0
        else
           d2b_local=min(abs(tmp2-temp),abs(tmp2+temp))
           !   dist2bound=abs(x+0.4)
        endif      
        
     else if(itest.eq.8) then   
        xhat = x - blob_centre(nb_blobs,1)
        yhat = y 
        !! Stretch high-res region downstream of flameholder(s). Note blob nb_blobs is furthest downstream...
        if((x-blob_centre(nb_blobs,1)).gt.0.0d0) then
     
           r_mag = ((xb_max - x)/(xb_max - blob_centre(nb_blobs,1)))**1.0d0  !! Scale between blob and outlet (0=outlet)
           temp = exp(-(0.5d0*yhat)**4.0d0) !! Blob-side spreading function
           tmp2 = exp(-(yhat/3.0d0)**4.0d0) !! Outflow-side spreading function
           temp = r_mag*temp + (1.0d0-r_mag)*tmp2 !! Linear variation between blob-side and outflow-side
           
           dxio = dx_out + (dx_in - dx_out)*(1.0d0-temp)
!           dxio = dx_out*r_mag + (5.0d0/3.0d0)*dx_out*(1.0d0-r_mag)
           dxio = dxio*(r_mag**0.5d0 + (3.0d0/3.0d0)*(1.0d0-r_mag**0.5d0))

        else
           tmp2 = ((blob_centre(1,1)-x)/(blob_centre(1,1)-xb_min))**2.0d0
           d2b_local = d2b_local*(1.5d0*tmp2 + 1.0d0*(1.0d0-tmp2))
           r_mag = sqrt(xhat**2.0d0 + yhat**2.0d0)
           temp = exp(-(0.5d0*r_mag)**4.0d0)
           dxio = dx_out + (dx_in - dx_out)*(1.0d0-temp)
        endif    
        
     else if(itest.eq.9) then   
           if(x.le.0.0d0) then !-1
              dxio = dx_in
           else if(x.le.1.0d0) then !0
              temp = (x+0.0d0)/1.0d0
              dxio = dx_in*(1.0d0-temp) + dx_out*temp
           else
              dxio = dx_out              
           end if
              
           d2b_local = bdist

             
     end if
        
     !! And what is the spacing, based on dist2bound?
     if(d2b_local.le.b0*dx0) then  !! Close - set to dxmin
        dx_local = dxmin                     
     else if(d2b_local.le.b1*dx0)then  !! A bit further out, smoothly vary from dxmin to dx0
        dx_local = 0.5d0*(dx0+dxmin) - 0.5d0*(dx0-dxmin)*cos((d2b_local-b0*dx0)*pi/((b1-b0)*dx0))  
     else if(d2b_local.le.b1*dx0+b2*dxio)then  !! Further still: linearly vary from dx0 to dxio
        dx_local = dx0 + (dxio-dx0)*((d2b_local-b1*dx0)/(b2*dxio))
     else     !! Far out: set to dxio
        dx_local = dxio
     end if  


      
      
      return
   end subroutine get_resolution
!! ------------------------------------------------------------------------------------------------
   subroutine make_boundary_blobs
     use kind_parameters
     use common_parameter
     use global_variables 

     implicit none
     integer(ikind) :: ipart,ib,ibc
     real(rkind) :: x,y,x0,y0,tmp2,r_mag,th,th_sh,dxlocal,th0
     real(rkind) :: a0,a1,a2,a3,a4,a5
   
     ipart=nb

     allocate(Sblob(nb_blobs))
     call evaluate_blob_perimeters


     if(nb_blobs.ne.0)then
        do ib=1,nb_blobs
                  
           !! Estimate approx number of nodes round the blob              
           nb = floor(Sblob(ib)/dx_wall)

           !! Revise node spacing to get integer number round blob
           dxlocal = Sblob(ib)/dble(nb)                      
                      
           a0 = blob_coeffs(ib,1)     
           x0 = blob_centre(ib,1);y0=blob_centre(ib,2)
           th = 0.0d0*0.5d0*dxlocal/a0 !! Initial theta
           do while(th.le.2.0d0*pi-0.0001*dxlocal/a0)
              !! Position node
              th_sh = th - blob_rotation(ib)
                
              !! Evaluate radius
              r_mag = blob_coeffs(ib,1)
              do ibc=2,n_blob_coefs
                 r_mag = r_mag + blob_coeffs(ib,ibc)*cos(dble(ibc-1)*th_sh)
              end do
                 
              !! Relative position
              x = r_mag*cos(th);y = r_mag*sin(th)

              !! If the position is within the domain, create a particle
              if(x+x0.ge.xb_min.and.x+x0.lt.xb_max.and.y+y0.ge.yb_min.and.y+y0.lt.yb_max)then              
                 ipart = ipart + 1 
                 xp(ipart)=x0 + x;yp(ipart)=y0 + y;dxp(ipart)=dx_wall
         
                 !! Calculate normals
                 tmp2 = 0.0
                 do ibc=2,n_blob_coefs
                    tmp2 = tmp2 - dble(ibc-1)*blob_coeffs(ib,ibc)*sin(dble(ibc-1)*th_sh)
                 end do     
                         
                 tmp2 = -atan(tmp2/r_mag)
                 xnorm(ipart) = (x*cos(tmp2) - y*sin(tmp2))/r_mag
                 ynorm(ipart) = (x*sin(tmp2) + y*cos(tmp2))/r_mag
                 node_type(ipart) = 0 !! Identify it as a wall                    
              end if                 
      
              !! Increment angle...
              th  = th + cos(tmp2)*dxlocal/r_mag     
           end do   
                                                                                      
        end do              
     end if

     deallocate(Sblob)
     nb=ipart    
!!
     write(6,*) 'no. of solid boundary particles: ',nb
     return
   end subroutine make_boundary_blobs  
!! --------------------------------------------------------------
   subroutine evaluate_blob_perimeters
     !! Numerically integrate the length of the blob perimeter
     use kind_parameters
     use common_parameter
     use global_variables 

     implicit none
     integer(ikind) :: ipart,ib,ibc
     real(rkind) :: x,y,x0,y0,tmp2,r_mag,th,th0
     real(rkind) :: a0,a1,a2,a3,a4,a5,dstiny
   

     if(nb_blobs.ne.0)then
        do ib=1,nb_blobs
          
           !! Initialise counters etc
           Sblob(ib) = 0.0d0
           a0 = blob_coeffs(ib,1)     
           !! Coordinates of start point
           y0 = 0.0d0
           x0 = a0
           do ibc = 2,n_blob_coefs
              x0 = x0 + blob_coeffs(ib,ibc)
           end do
                           
           th = 0.0d0 !! theta       
           dstiny = 0.005d0*dx_wall       
           do while(th.le.2.0d0*pi-0.01*dstiny/x0)
                
              !! Evaluate radius
              r_mag = blob_coeffs(ib,1)
              do ibc=2,n_blob_coefs
                 r_mag = r_mag + blob_coeffs(ib,ibc)*cos(dble(ibc-1)*th)
              end do
                 
              !! Calculate normals
              tmp2 = 0.0
              do ibc=2,n_blob_coefs
                 tmp2 = tmp2 - dble(ibc-1)*blob_coeffs(ib,ibc)*sin(dble(ibc-1)*th)
              end do                             
              tmp2 = -atan(tmp2/r_mag)
      
              !! Increment angle...
              th  = th + cos(tmp2)*dstiny/r_mag                      
              !! Increment Sblob
              Sblob(ib) = Sblob(ib) + dstiny
           end do                                                                     
        end do              
     end if

     return
   end subroutine evaluate_blob_perimeters  
