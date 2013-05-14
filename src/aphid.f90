
!======================================================================!
!================================================================! MAIN
!======================================================================!

      subroutine aphid_main
      use general

!----------------------------------------------------------------------!
!  Program to solve the full-physics electromagnetic induction problem !
!  in dielectric/conductive media using a A-PHI vector potential      !
!  decomposition on a staggered Yee grid.                              !
!                                                                      !
!  Written by: Chester J Weiss, Virginia Tech Dept of Geosci, 2011     !
!                                                                      !
!                   ---  IMPLEMENTATION NOTES  ---                     !
!                                                                      !
!  1. Mesh input file is identical to that used previously in FDM3D.   !
!  2. Linear system is also solved in the same matrix-free paradigm    !
!  3. Like the global code, TRIFFID, the code is parallelized over     !
!     multiple frequencies. 
!----------------------------------------------------------------------!

      character(20) :: argv,cstr
      integer       :: iargc


! VARIABLES !----------------------------------------------------------!
!  iargc    C-intrinsic that gets the number of arguments from the     !
!           the command line.                                          !
!  argv     Character array which stores a given command line argument !
!  cstr     Dummy array for holding OMP_NUM_THREADS.                   !
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!  Accommodate command line arguments which specify working directory  !
!  different than the present one.                                     !
!                                                                      !
!  If aphi_serial was invoked with command line argument "-w [pwd]",   !
!  then store value of the string 'pwd' in a character array of the    !
!  same name.  The default is pwd = './'                               !
!----------------------------------------------------------------------!

!
!  get number of arguments from the command line
!
      n = iargc() 
!
!  grab the arguments from the command line, one by one...
!
      i=1
 1010 if (i.le.n) then
        call getarg(i,argv)
        if (argv(1:2).eq.'-w') then
          call getarg(i+1,argv)
          read(argv,*) pwd
          i = i + 2
          goto 1010
        end if
      end if

!----------------------------------------------------------------------!
!                      Display the welcome banner                      !
!----------------------------------------------------------------------!

      call welcome_banner

!----------------------------------------------------------------------!
!                  Read in the model input parameters                  !
!----------------------------------------------------------------------!

      call read_input

!----------------------------------------------------------------------!
!           Read in the node coordinates and cell conductivities       !
!----------------------------------------------------------------------!

      call read_mesh
      call dump_sigma

!----------------------------------------------------------------------!
!                          Generate the RHS                            !
!----------------------------------------------------------------------!

      call gen_rhs

!----------------------------------------------------------------------!
!                       Solve the linear system                        !
!----------------------------------------------------------------------!

      call gen_diag

!----------------------------------------------------------------------!
!                       Solve the linear system                        !
!----------------------------------------------------------------------!

      call gen_inprod_x_coefs
      call gen_inprod_y_coefs
      call gen_inprod_z_coefs
      call gen_inprod_p_coefs

      call get_environment_variable("OMP_NUM_THREADS",cstr)
      if (cstr(1:1).eq.'1') call qmr
      if (cstr(1:1).eq.'4') call qmr_mthread

!----------------------------------------------------------------------!
!                        Dump out the results                          !
!----------------------------------------------------------------------!

      call dump_a_phi
      call gen_fields
      call dump_fields

!----------------------------------------------------------------------!
!                                END                                   !
!----------------------------------------------------------------------!

      end subroutine aphid_main


!======================================================================!
!======================================================! WELCOME_BANNER
!======================================================================!

      subroutine welcome_banner
      use general

#ifdef SVN_VERSION_DEF
      character(8) :: svn_version_string
      SVN_VERSION_DEF
#endif

      write(iout,*) ' '

      write(iout,*) '*************************************************',&
                  '*********'
      write(iout,*) '*                                                ',&
                  '        *'
      write(iout,*) '*                          APhiD                 ',&
                  '        *'
      write(iout,*) '*                                                ',&
                  '        *'
      write(iout,*) '*                Cartesian Electromagnetic       ',&
                  '        *'
      write(iout,*) '*                     Induction Code             ',&
                  '        *'
      write(iout,*) '*                          using                 ',&
                  '        *'
      write(iout,*) '*              A-PHI potential Decomposition     ',&
                  '        *'
      write(iout,*) '*                                                ',&
                  '        *'
      write(iout,*) '*                   Chester J Weiss              ',&
                  '        *'
      write(iout,*) '*  Virginia Polytechnic Institute and State ',    &
                  'University   *'
      write(iout,*) '*                                                ',&
                  '        *'
      write(iout,*) '*                    cjweiss@vt.edu              ',&
                  '        *'
      write(iout,*) '*                                                ',&
                  '        *'
      write(iout,*) '*************************************************',&
                  '*********'

      write(iout,*) ' '
#ifdef SVN_VERSION_DEF
      write(iout,*) '!------------------------------------------------',&
                  '--------!'
      write(iout,*) '! SVN Repository Version: ',svn_version_string
      write(iout,*) '!------------------------------------------------',&
                  '--------!'
#endif
      write(iout,*) ' '

      end subroutine welcome_banner


!======================================================================!
!==========================================================! READ_PARAM
!======================================================================!

      subroutine read_input
      use general

!==============================!
! READ THE PRIMARY INPUT FILE  !
!==============================!

      write(iout,*) 'reading model input file '''//trim(pwd)//          &
                                                     '/aphid.input'' '
      open(unit=11,file=trim(pwd)//'/aphid.input',status='unknown')

      read(11,*)
      read(11,*)
      read(11,*)
      read(11,*,err=101) oper_mode
      read(11,*,err=103) qmr_target
      read(11,*,err=104) qmr_max_it
      read(11,*) freq
      read(11,*) csig_flag
      close(11)

      omega = 2.d0 * PI * freq
!
!  error messages for invalid arguments
!
      if (qmr_target.lt.0.) then
        write(iout,*) 'Invalid target residual: value must positive.'
        stop
      end if

      if (qmr_max_it.lt.0) then
        write(iout,*) 'Invalid qmr_max_it: value must positive.'
        stop
      end if

!===========================!
! DUMP OUT WHAT'S BEEN READ !
!===========================!

      write(iout,*) ' '
      write(iout,*) ':: INPUT PARAMETERS '
      write(iout,'(a,i2,a)')' :: Operational mode: ',oper_mode
      write(iout,*) ':: '
      write(iout,*) ':: QMR SOLVER PARAMETERS'
      write(iout,3) ':: target residual reduction: ', qmr_target
      write(iout,*) ':: maximum iterations: ',qmr_max_it
      write(iout,*) ':: '
      write(iout,*) ':: RUNTIME PARAMETERS'
      write(iout,3) ':: frequency [Hz]: ', freq
      if (csig_flag.eq.'R')write(iout,3)':: real conductivity'
      if (csig_flag.eq.'C')write(iout,3)':: complex conductivity'
      if (csig_flag.eq.'D')write(iout,3)':: Cole-Cole dispersion model'
      if (csig_flag.eq.'I')write(iout,3)':: Insulating dielectric model'
      write(iout,3) '::  '
      write(iout,*) ' '

 2    format(1x,":: ",i3,3x,i3,4x,e10.3,4x,e10.3)
 3    format(1x,a,2x,e10.3)
 4    format(1x,a,2x,f8.5)
 10   format(1x,a,a)  
!
!  error messages for absent arguments
!
      return
 101  call input_file_error("operational mode")
 103  call input_file_error("qmr_target")
 104  call input_file_error("qmr_max_it")


      end subroutine read_input


!======================================================================!  
!====================================================! INPUT_FILE_ERROR
!======================================================================!  

      subroutine input_file_error(str)
      character(*) :: str
      write(iout,*) '>> ERROR reading '//trim(str)//' from aphid.input'
      stop
      end subroutine input_file_error

!======================================================================!
!===========================================================! READ_MESH
!======================================================================!

      subroutine read_mesh

!----------------------------------------------------------------------!
!  Routine to 1) read in the node coordinates and cell conductivities  !
!  for the FD mesh, 2) report back the number of degrees of freedom    !
!  in the upcoming linear system, and, 3) determine how many model     !
!  realizations will be computed.                                      !  
!                                                                      !
!  Written by:  C J Weiss, Virginia Tech Dept of Geosciences, 2011     !
!----------------------------------------------------------------------!

      use general
      use mesh

      integer         :: ncel
      real(DP)        :: mindx,mindy,mindz
      real(DP)        :: maxdx,maxdy,maxdz
      real(DP)        :: meanvol,minvol,maxvol,vol
      character(2)    :: calc_type
      real(sp),dimension(:,:,:) :: sig_re,eps_rel,eta,tau,cee

      character(12)   :: cskip
      allocatable     :: sig_re,eps_rel,eta,tau,cee

! Variables !----------------------------------------------------------!
!  ncel:           Number of cells in the mesh.                        !
!  mind*:          Minimum dx, dy and dz values for the mesh.          ![m]
!  maxd*           Maximum dx, dy and dz values for the mesh.          ![m]
!  minvol,meanvol,maxvol,vol: Mesh statistics for the cell volumes.    ![m^3]
!----------------------------------------------------------------------!

!
!  Write an opening message...
!
      write(iout,*)'*********************************'
      write(iout,*)'*  Reading the mesh parameters  *'
      write(iout,*)'*    and conductivity model     *'
      write(iout,*)'*********************************'
      write(iout,*) ' '
      write(iout,*) 'Input file '''//trim(pwd)//'/aphid.mesh'' '
!
!  Open the file containing the mesh conductivity information.
!
      open(unit=11,file=trim(pwd)//'/aphid.mesh',status='unknown')
!
!  Read in the number of nodes in the x, y and z directions.
!
      read(11,*) nx,ny,nz
      read(11,*)
!
!  Allocate storage for the conductivity array and node arrays 
!
      allocate (  xn( nx ),  yn( ny ),  zn( nz ) )
      allocate (  xh(nx-1),  yh(ny-1),  zh(nz-1) )
      allocate (  dx(0:nx),  dy(0:ny),  dz(0:nz) )
      allocate ( dxc(1:nx), dyc(1:ny), dzc(1:nz) )
      allocate (          sig (nx-1,ny-1,nz-1)   )
      allocate (     k2_xgrid (nx-1,ny,nz)       )
      allocate (     k2_ygrid (nx,ny-1,nz)       )
      allocate (     k2_zgrid (nx,ny,nz-1)       )
      allocate (     k2_pgrid (nx,ny,nz)         )
!
!  Read in the locations of the x, y, and z nodes.
!
      read(11,*) (xn(i),i=1,nx)
      read(11,*)
      read(11,*) (yn(i),i=1,ny)
      read(11,*)
      read(11,*) (zn(i),i=1,nz)
!
!  Compute the half-index values of the node coordinates
!
      xh = 0.5d0*( xn(2:nx) + xn(1:nx-1) )
      yh = 0.5d0*( yn(2:ny) + yn(1:ny-1) )
      zh = 0.5d0*( zn(2:nz) + zn(1:nz-1) )
!
!  Compute the nodal spacings.
!
      dx(1:nx-1) = xn(2:nx) - xn(1:nx-1)
      dy(1:ny-1) = yn(2:ny) - yn(1:ny-1)
      dz(1:nz-1) = zn(2:nz) - zn(1:nz-1)
    
      dx(0) = dx(1) ; dx(nx) = dx(nx-1) 
      dy(0) = dy(1) ; dy(ny) = dy(ny-1) 
      dz(0) = dz(1) ; dz(nz) = dz(nz-1) 

      dxc(2:nx-1) = 0.5d0*(dx(1:nx-2)+dx(2:nx-1))
      dyc(2:ny-1) = 0.5d0*(dy(1:ny-2)+dy(2:ny-1))
      dzc(2:nz-1) = 0.5d0*(dz(1:nz-2)+dz(2:nz-1))

      dxc(1) = dxc(2) ; dxc(nx) = dxc(nx-1)
      dyc(1) = dyc(2) ; dyc(ny) = dyc(ny-1)
      dzc(1) = dzc(2) ; dzc(nz) = dzc(nz-1)
!
!  Compute the number of cells in the mesh.
!
      ncel = (nx-1)*(ny-1)*(nz-1)

!----------------------------------------------------------------------!
!  Depending on the value of 'csig_flag' read in the conductivity      !
!  model according to...                                               !
!                                                                      !
!  R:  real-valued conductivity with uniform free-space permittivity   !
!  C:  complex-valued conductivity, frequency independent              !
!  D:  complex-valued conductivity, frequency dependent Cole-Cole      !
!  I:  insulator with relative permittivity
!----------------------------------------------------------------------!

!  for real-valued conductivity models...

      if (csig_flag.eq.'R') then

        allocate (   sig_re(nx-1,ny-1,nz-1)    )

        read(11,*)
        read(11,*) (((sig_re(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)
        sig = sig_re  + II * omega * EPS0

        deallocate (   sig_re   )

      end if

!  for complex-valued conductivity models: sigma,chi

      if (csig_flag.eq.'C') then

        allocate (   sig_re(nx-1,ny-1,nz-1)    )
        allocate (   eps_rel(nx-1,ny-1,nz-1)    )

        read(11,*)
        read(11,*) (((sig_re(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)

        read(11,*)
        read(11,*) (((eps_rel(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)

        sig = sig_re + II * omega * EPS0*eps_rel

        deallocate (   sig_re, eps_rel   )

      end if

!  for Cole-Cole models: sig,eta,tau,cee

      if (csig_flag.eq.'D') then

        allocate (   sig_re(nx-1,ny-1,nz-1)    )
        allocate (   eta(nx-1,ny-1,nz-1)    )
        allocate (   tau(nx-1,ny-1,nz-1)    )
        allocate (   cee(nx-1,ny-1,nz-1)    )

        read(11,*)
        read(11,*) (((sig_re(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)

        read(11,*)
        read(11,*) (((eta(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)

        read(11,*)
        read(11,*) (((tau(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)

        read(11,*)
        read(11,*) (((cee(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)


        sig = sig_re * (1.d0 - eta / (1.d0 + (II*omega*tau)**cee) )

        deallocate (   sig_re,eta,tau,cee  )
      end if

!  for imaginary-valued conductivity models of insulating dielectrics...

      if (csig_flag.eq.'I') then

        allocate (   eps_rel(nx-1,ny-1,nz-1)    )

        read(11,*)
        read(11,*) (((eps_rel(i,j,k),i=1,nx-1),j=1,ny-1),k=1,nz-1)
        sig = II * omega * EPS0 * eps_rel

        deallocate (   eps_rel   )

      end if

      call gen_k2_xgrid
      call gen_k2_ygrid
      call gen_k2_zgrid
      call gen_k2_pgrid
!
!  Close the mesh file.
!
      close(11)
!
!  Compute some mesh statistics and write them out...
!
      mindx = xn(2)-xn(1)
      maxdx = xn(2)-xn(1)
      do i=1,nx-1
        if(xn(i+1)-xn(i).lt.mindx) mindx = xn(i+1)-xn(i)
        if(xn(i+1)-xn(i).gt.maxdx) maxdx = xn(i+1)-xn(i)
      enddo

      mindy = yn(2)-yn(1)
      maxdy = yn(2)-yn(1)
      do i=1,ny-1
        if(yn(i+1)-yn(i).lt.mindy) mindy = yn(i+1)-yn(i)
        if(yn(i+1)-yn(i).gt.maxdy) maxdy = yn(i+1)-yn(i)
      enddo

      mindz = zn(2)-zn(1)
      maxdz = zn(2)-zn(1)
      do i=1,nz-1
        if(zn(i+1)-zn(i).lt.mindz) mindz = zn(i+1)-zn(i)
        if(zn(i+1)-zn(i).gt.maxdz) maxdz = zn(i+1)-zn(i)
      enddo

      minvol = (xn(2)-xn(1))*(yn(2)-yn(1))*(zn(2)-zn(1))
      maxvol = (xn(2)-xn(1))*(yn(2)-yn(1))*(zn(2)-zn(1))
      meanvol = (xn(nx)-xn(1))*(yn(ny)-yn(1))*(zn(nz)-zn(1))/real(ncel)
      do i=1,nx-1
        do j=1,ny-1
          do k=1,nz-1
            vol = (xn(i+1)-xn(i))*(yn(j+1)-yn(j))*(zn(k+1)-zn(k))
            if (vol.lt.minvol) minvol = vol
            if (vol.gt.maxvol) maxvol = vol
          enddo
        enddo
      enddo

      write(iout,*) ' '
      write(iout,*) '====================='
      write(iout,*) '   Mesh properties   '
      write(iout,*) '====================='
      write(iout,*) ' '

      write(iout,98) 'nodes in the x-direction:',nx
      write(iout,98) 'nodes in the y-direction:',ny
      write(iout,98) 'nodes in the z-direction:',nz
      write(iout,*) ' '

      write(iout,99) 'min node position in x[m]:',xn(1)
      write(iout,99) 'max node position in x[m]:',xn(nx)
      write(iout,99) 'min node position in y[m]:',yn(1)
      write(iout,99) 'max node position in y[m]:',yn(ny)
      write(iout,99) 'min node position in z[m]:',zn(1)
      write(iout,99) 'max node position in z[m]:',zn(nz)
      write(iout,*) ' '

      write(iout,99) 'min, mean, max node spacing in x[m]:',            &
                    mindx,(xn(nx)-xn(1))/real(nx-1),maxdx
      write(iout,99) 'min, mean, max node spacing in y[m]:',            &
                   mindy,(yn(ny)-yn(1))/real(ny-1),maxdy
      write(iout,99) 'min, mean, max node spacing in z[m]:',            &
                   mindz,(zn(nz)-zn(1))/real(nz-1),maxdz
      write(iout,99) 'min, mean, max cell volumes [m^3]:  ',            &
                   minvol,meanvol,maxvol

 98   format(1x,a,1x,i3)
 99   format (1x,a,1x,3(e12.5,2x))


      end subroutine read_mesh

!======================================================================!
!=============================================================! GEN_RHS
!======================================================================!

      subroutine gen_rhs
      use general
      use mesh

! Variables !----------------------------------------------------------!
!----------------------------------^-----------------------------------!

      real(dp) :: dxi,dyj,dycj,dzck
      integer  :: i,j,k
      character(2) :: jdir
!
!  allocate storage for the RHS
!
      write(iout,*) ' '
      write(iout,'(a)',advance='NO') ' allocating storage for RHS...'
      allocate ( rhs_x(nx-1,ny,nz) ) 
      allocate ( rhs_y(nx,ny-1,nz) )
      allocate ( rhs_z(nx,ny,nz-1) )
      allocate (  rhs_p(nx,ny,nz)  )
      write(iout,*) ' (done)'
!
!  initialize values for RHS
!
      rhs_x = (0.d0,0.d0)
      rhs_y = (0.d0,0.d0)
      rhs_z = (0.d0,0.d0)
      rhs_p = (0.d0,0.d0)
!
!  oper_mode = 1: elementary x-directed electric dipole source 
!                 with length equal to one edge and the tail-end
!                 sitting at the mesh origin.  Dipole moment is unity.
!
      if (oper_mode.eq.1) then
        do k=1,nz
          if (abs(zn(k)-0.).lt.1.e-4) kp = k
          do j=1,ny
            if (abs(yn(j)- 0.).lt.1.e-4) jp = j
            do i=1,nx
              if (abs(xn(i) - 0.).lt.1.e-4) ip = i
            end do 
          end do 
        end do 
        dxci = 0.5d0*(dx(ip-1)+dx(ip))
        dycj = 0.5d0*(dy(jp-1)+dy(jp))
        dzck = 0.5d0*(dz(kp-1)+dz(kp))

        write(iout,*) ':: Jx dipole source at node',ip,jp,kp
        write(iout,*) (xn(ip)+xn(ip+1))/2.,yn(jp),zn(kp)

        rhs_x(ip,jp,kp)  = 1.d0 * MU0
        rhs_p(ip,jp,kp) =  1.d0 * MU0 / dx(ip)
        rhs_p(ip+1,jp,kp) = -1.d0 * MU0/ dx(ip)

      end if

!
!  oper_mode = 2: x-polarized MT source
!
      if (oper_mode.eq.2) then
        rhs_x(:,:,2) = 1.d0 * MU0

        do i=1,nx-1
          rhs_x(i,:,2) = rhs_x(i,:,2) * dx(i)
        end do

        do j=1,ny
          rhs_x(:,j,2) = rhs_x(:,j,2) * 0.5d0*(dy(j-1)+dy(j))
        end do

        write(iout,*) ':: Jx plane wave source '
      end if
!
!  oper_mode = 3: y-polarized MT source
!
      if (oper_mode.eq.3) then
        rhs_y(:,:,2) = 1.d0 * MU0

        do j=1,ny-1
          rhs_y(:,j,2) = rhs_y(:,j,2) * dy(j)
        end do

        do i=1,nx
          rhs_y(i,:,2) = rhs_y(i,:,2) * 0.5d0*(dx(i-1)+dx(i))
        end do

        write(iout,*) ':: Jy plane wave source '
      end if
 
!  oper_mode = 4: fixed size square loop for testing purposes
!
      if (oper_mode.eq.4) then
        n = 20
        rhs_x(51-n:51+n-1,51+n,51) =  1.d0 * MU0 * dx (51-n:51+n-1)
        rhs_x(51-n:51+n-1,51-n,51) = -1.d0 * MU0 * dx (51-n:51+n-1)
        rhs_y(51-n,51-n:51+n-1,51) =  1.d0 * MU0 * dy (51-n:51+n-1)
        rhs_y(51+n,51-n:51+n-1,51) = -1.d0 * MU0 * dy (51-n:51+n-1)

      end if

!
!  oper_mode = 5: user-specified dipole assemblage, read in from 
!                 file 'aphid.source'
!                 sign convention: index is the tail-end node of the dipole
!
      if (oper_mode.eq.5) then
        isrc_num = 0
        open(unit=11,file=trim(pwd)//'/aphid.source',status='unknown')
        do while (.true.)

          read(11,*,end=1) jdir,ip,jp,kp
          isrc_num = isrc_num + 1

          write(iout,*) jdir,ip,jp,kp

          if ((jdir.eq.'+X').or.(jdir.eq.'-X') .or.                     &
              (jdir.eq.'+Y').or.(jdir.eq.'-Y') .or.                     &
              (jdir.eq.'+Z').or.(jdir.eq.'-Z') ) then

            if (jdir.eq.'+X') then
              if ((ip.lt.1).or.(ip.gt.nx-1) .or.                        &
                           (jp.lt.1).or.(jp.gt.ny  ) .or.               &
                                     (kp.lt.1).or.(kp.gt.nz  ) ) then
                write(iout,*) isrc_num,' error: index out of bounds...' &
                           // ' skipping.'
              else
               
                rhs_x(ip  ,jp,kp) =  rhs_x(ip,jp,kp)  +1.d0*MU0*dx(ip)
                rhs_p(ip  ,jp,kp) =  rhs_p(ip,jp,kp)  +1.d0*MU0 
                rhs_p(ip+1,jp,kp) =  rhs_p(ip+1,jp,kp)-1.d0*MU0 
              end if 
            end if
  
            if (jdir.eq.'-X') then
              if ((ip.lt.2).or.(ip.gt.nx  ) .or.                        &
                           (jp.lt.1).or.(jp.gt.ny  ) .or.               &
                                     (kp.lt.1).or.(kp.gt.nz  ) ) then
                write(iout,*) isrc_num,' error: index out of bounds...' &
                          // ' skipping.'
              else
                rhs_x(ip-1,jp,kp) = rhs_x(ip-1,jp,kp)-1.d0*MU0*dx(ip-1)
                rhs_p(ip-1,jp,kp) = rhs_p(ip-1,jp,kp)-1.d0*MU0 
                rhs_p(ip  ,jp,kp) = rhs_p(ip,jp,kp)  +1.d0*MU0
              end if
            end if
  
            if (jdir.eq.'+Y') then
              if ((ip.lt.1).or.(ip.gt.nx  ) .or.                        &
                           (jp.lt.1).or.(jp.gt.ny-1) .or.               &
                                     (kp.lt.1).or.(kp.gt.nz  ) ) then
                write(iout,*) isrc_num,' error: index out of bounds...' &
                          //  ' skipping.'
              else
                rhs_y(ip,jp  ,kp) =  rhs_y(ip,jp,kp)  +1.d0*MU0*dy(jp)
                rhs_p(ip,jp  ,kp) =  rhs_p(ip,jp,kp)  +1.d0*MU0
                rhs_p(ip,jp+1,kp) =  rhs_p(ip,jp+1,kp)-1.d0*MU0 
              end if 
            end if
  
            if (jdir.eq.'-Y') then
              if ((ip.lt.1).or.(ip.gt.nx  ) .or.                        &
                           (jp.lt.2).or.(jp.gt.ny  ) .or.               &
                                     (kp.lt.1).or.(kp.gt.nz  ) ) then
                write(iout,*) isrc_num,' error: index out of bounds...' &
                           // ' skipping.'
              else
                rhs_y(ip,jp-1,kp) = rhs_y(ip,jp-1,kp)-1.d0*MU0*dy(jp-1)
                rhs_p(ip,jp-1,kp) = rhs_p(ip,jp-1,kp)-1.d0*MU0
                rhs_p(ip,jp  ,kp) = rhs_p(ip,jp  ,kp)+1.d0*MU0 
              end if 
            end if
  
            if (jdir.eq.'+Z') then
              if ((ip.lt.1).or.(ip.gt.nx  ) .or.                        &
                           (jp.lt.1).or.(jp.gt.ny  ) .or.               &
                                     (kp.lt.1).or.(kp.gt.nz-1) ) then
                write(iout,*) isrc_num,' error: index out of bounds...' &
                           // ' skipping.'
              else
                rhs_z(ip,jp,kp)   =  rhs_z(ip,jp,kp)  +1.d0*MU0*dz(kp)
                rhs_p(ip,jp,kp)   =  rhs_p(ip,jp,kp)  +1.d0*MU0
                rhs_p(ip,jp,kp+1) =  rhs_p(ip,jp,kp+1)-1.d0*MU0
              end if 
            end if
  
            if (jdir.eq.'-Z') then
              if ((ip.lt.1).or.(ip.gt.nx  ) .or.                        &
                           (jp.lt.1).or.(jp.gt.ny  ) .or.               &
                                     (kp.lt.2).or.(kp.gt.nz  ) ) then
                write(iout,*) isrc_num,' error: index out of bounds...' &
                           // ' skipping.'
              else
                rhs_z(ip,jp,kp-1) = rhs_z(ip,jp,kp-1)-1.d0*MU0*dz(kp-1)
                rhs_p(ip,jp,kp-1) = rhs_p(ip,jp,kp-1)-1.d0*MU0
                rhs_p(ip,jp  ,kp) = rhs_p(ip,jp,kp)  +1.d0*MU0
              end if 
            end if

          else
            write(iout,*) isrc_num,' error: invalid source direction...'&
                       //  ' skipping'
          end if 



       end do
 1     close(11)

      end if

      end subroutine gen_rhs


!======================================================================!
!==================================================! GEN_INPROD_X_COEFS
!======================================================================!

      subroutine gen_inprod_x_coefs
      use general
      use mesh
      use inprod_x_coefs

      real(dp)     :: dxi,dxcip,dxci
      real(dp)     :: dyj,dyjm,dycj
      real(dp)     :: dzk,dzkm,dzck

      allocate ( coef_2vec(1:nx-1)  , coef_3vec(1:nx-1) )
      allocate ( coef_4vec(1:ny  )  , coef_5vec(1:ny  ) )
      allocate ( coef_6vec(1:nz  )  , coef_7vec(1:nz  ) )
      
      allocate ( ia_vec(1:nx-1)  , ib_vec(1:nx-1)  )
      allocate ( ja_vec(1:ny  )  , jb_vec(1:ny  )  )
      allocate ( ka_vec(1:nz  )  , kb_vec(1:nz  )  )
      
      do k=1,nz
        dzk  = dz(k) 
        dzkm = dz(k-1) 
        dzck = dzc(k)
        coef_6vec(k) = (k/nz - 1)        / dzk   / dzck  
        coef_7vec(k) = ((nz-k)/(nz-1)-1) / dzkm  / dzck
        ka_vec(k) = min(k+1,nz  )   
        kb_vec(k) = max(  1, k-1)
      end do

      do j=1,ny
        dyj  = dy(j)
        dyjm = dy(j-1)
        dycj = dyc(j)
        coef_4vec(j) = (j/ny - 1)        / dyj   / dycj
        coef_5vec(j) = ((ny-j)/(ny-1)-1) / dyjm  / dycj
        ja_vec(j) = min(j+1,ny  )
        jb_vec(j) = max(  1, j-1)
      end do

      do i=1,nx-1
        dxi  = dx(i)
        dxcip = dxc(i+1)
        dxci  = dxc(i)
        coef_2vec(i) = (i/(nx-1)-1)      / dxcip / dxi
        coef_3vec(i) = ((nx-i)/(nx-1)-1) / dxci  / dxi
        ia_vec(i) = min(i+1,nx-1)
        ib_vec(i) = max(  1, i-1)
      end do
 
      end subroutine gen_inprod_x_coefs

!======================================================================!
!============================================================! INPROD_X
!======================================================================!

      subroutine inprod_x(in_x,in_y,in_z,in_p,out_x)
      use general
      use mesh
      use inprod_x_coefs
 
      complex(dp), dimension(nx-1, ny , nz ) :: in_x,out_x
      complex(dp), dimension( nx ,ny-1, nz ) :: in_y
      complex(dp), dimension( nx , ny ,nz-1) :: in_z
      complex(dp), dimension( nx , ny , nz ) :: in_p
      complex(dp), dimension(7)  :: s
 
      real(dp)     :: dxi,dycj,dzck,vol
      complex(dp)  :: k2,kp1,kp2
      integer      :: ia,ib,ja,jb,ka,kb

      complex(dp), dimension(7)  :: coef
      complex(dp), dimension(7)  :: invec

      out_x = 0.d0

!----------------------------------------------------------------------!
!                  Work on the x-directed edges first                  !
!----------------------------------^-----------------------------------!

      do k=1,nz
        dzck = dzc(k)
        coef(6) = coef_6vec(k)
        coef(7) = coef_7vec(k)
        ka = ka_vec(k)
        kb = kb_vec(k)

        do j=1,ny
          dycj = dyc(j)
          coef(4) = coef_4vec(j)
          coef(5) = coef_5vec(j)
          ja = ja_vec(j)
          jb = jb_vec(j)

          do i=1,nx-1
            dxi  = dx(i) 
            coef(2) = coef_2vec(i)
            coef(3) = coef_3vec(i)
            ia = ia_vec(i)
            ib = ib_vec(i)
           
            vol = dxi * dycj * dzck
            k2  = k2_xgrid( i, j,k)
            kp1 = k2_pgrid( i ,j,k)
            kp2 = k2_pgrid(i+1,j,k)

            coef(1) = -sum(coef(2:7)) + k2

            invec(1) = in_x( i , j , k )
            invec(2) = in_x(ia, j , k )
            invec(3) = in_x(ib, j , k )
            invec(4) = in_x( i ,ja, k )
            invec(5) = in_x( i ,jb, k )
            invec(6) = in_x( i , j, ka)
            invec(7) = in_x( i , j, kb)

            out_x(i,j,k) = ( sum( coef*invec ) +                        &
                              ( in_p(i+1,j,k)*(kp2-k2) -                &
                                in_p( i ,j,k)*(kp1-k2)   )/dxi   ) * vol
          
          end do
        end do
      end do

      end subroutine inprod_x


!======================================================================!
!==================================================! GEN_INPROD_Y_COEFS
!======================================================================!

      subroutine gen_inprod_y_coefs
      use general
      use mesh
      use inprod_y_coefs

      real(dp)     :: dxi,dxim,dxci
      real(dp)     :: dyj,dycjp,dycj
      real(dp)     :: dzk,dzkm,dzck

      allocate ( coef_2vec(1:ny-1)  , coef_3vec(1:ny-1) )
      allocate ( coef_4vec(1:nx  )  , coef_5vec(1:nx  ) )
      allocate ( coef_6vec(1:nz  )  , coef_7vec(1:nz  ) )
      
      allocate ( ia_vec(1:nx  )  , ib_vec(1:nx  )  )
      allocate ( ja_vec(1:ny-1)  , jb_vec(1:ny-1)  )
      allocate ( ka_vec(1:nz  )  , kb_vec(1:nz  )  )

      do k=1,nz
        dzk  = dz(k)
        dzkm = dz(k-1)
        dzck = dzc(k)
        ka_vec(k) = min(k+1,nz  )
        kb_vec(k) = max(  1, k-1)
        coef_6vec(k) = (k/nz-1)          / dzk   / dzck
        coef_7vec(k) = ((nz-k)/(nz-1)-1) / dzkm  / dzck
      end do

      do j=1,ny-1
        dyj  = dy(j)
        dycjp = dyc(j+1)
        dycj  = dyc(j)
        ja_vec(j) = min(j+1,ny-1)
        jb_vec(j) = max(  1, j-1)
        coef_2vec(j) = (j/(ny-1)-1)      / dycjp / dyj
        coef_3vec(j) = ((ny-j)/(ny-1)-1) / dycj  / dyj
      end do

      do i=1,nx
        dxi  = dx(i) 
        dxim = dx(i-1) 
        dxci  = dxc(i)
        ia_vec(i) = min(i+1,nx  )
        ib_vec(i) = max(  1, i-1)
        coef_4vec(i) = (i/nx-1)          / dxi   / dxci
        coef_5vec(i) = ((nx-i)/(nx-1)-1) / dxim  / dxci
      end do

      end subroutine gen_inprod_y_coefs


!======================================================================!
!============================================================! INPROD_Y
!======================================================================!

      subroutine inprod_y(in_x,in_y,in_z,in_p,out_y)
      use general
      use mesh
      use inprod_y_coefs
      
      complex(dp), dimension(nx-1, ny , nz ) :: in_x
      complex(dp), dimension( nx ,ny-1, nz ) :: in_y,out_y
      complex(dp), dimension( nx , ny ,nz-1) :: in_z
      complex(dp), dimension( nx , ny , nz ) :: in_p
      complex(dp), dimension(7)  :: s
  
      real(dp)     :: dxci,dyj,dzck,vol
      complex(dp)  :: k2,kp1,kp2
      integer      :: ia,ib,ja,jb,ka,kb

      complex(dp), dimension(7)  :: coef
      complex(dp), dimension(7)  :: invec

      out_y = 0.d0

!----------------------------------------------------------------------!
!                  Work on the y-directed edges next                   !
!----------------------------------^-----------------------------------!

      do k=1,nz
        dzck = dzc(k)
        ka = ka_vec(k)
        kb = kb_vec(k)
        coef(6) = coef_6vec(k)
        coef(7) = coef_7vec(k)

        do j=1,ny-1
          dyj  = dy(j)
          ja = ja_vec(j)
          jb = jb_vec(j)
          coef(2) = coef_2vec(j)
          coef(3) = coef_3vec(j)

          do i=1,nx
            dxci  = dxc(i)
            ia = ia_vec(i)
            ib = ib_vec(i)
            coef(4) = coef_4vec(i)
            coef(5) = coef_5vec(i)

            vol = dxci * dyj * dzck
            k2 =  k2_ygrid(i, j, k)
            kp1 = k2_pgrid(i, j ,k)
            kp2 = k2_pgrid(i,j+1,k)

            coef(1) = -sum(coef(2:7)) + k2

            invec(1) = in_y( i , j , k )
            invec(2) = in_y( i ,ja, k )
            invec(3) = in_y( i ,jb, k )
            invec(4) = in_y(ia, j , k )
            invec(5) = in_y(ib, j , k )
            invec(6) = in_y( i , j ,ka)
            invec(7) = in_y( i , j ,kb)

            out_y(i,j,k) = ( sum( coef*invec ) +                        &
                              ( in_p(i,j+1,k)*(kp2-k2) -                &
                                in_p(i, j ,k)*(kp1-k2)   )/dyj   ) * vol

          end do
        end do
      end do

      end subroutine inprod_y


!======================================================================!
!==================================================! GEN_INPROD_Z_COEFS
!======================================================================!

      subroutine gen_inprod_z_coefs
      use general
      use mesh
      use inprod_z_coefs

      real(dp)     :: dxi,dxim,dxci
      real(dp)     :: dyj,dyjm,dycj
      real(dp)     :: dzk,dzckp,dzck

      allocate ( coef_2vec(1:nz-1)  , coef_3vec(1:nz-1) )
      allocate ( coef_4vec(1:nx  )  , coef_5vec(1:nx  ) )
      allocate ( coef_6vec(1:ny  )  , coef_7vec(1:ny  ) )
      
      allocate ( ia_vec(1:nx  )  , ib_vec(1:nx  )  )
      allocate ( ja_vec(1:ny  )  , jb_vec(1:ny  )  )
      allocate ( ka_vec(1:nz-1)  , kb_vec(1:nz-1)  )

      do i=1,nx
        dxi  = dx(i) 
        dxim = dx(i-1) 
        dxci = dxc(i)
        ia_vec(i) = min(i+1,nx)
        ib_vec(i) = max(  1, i-1)
        coef_4vec(i) = (i/nx-1)          / dxi   / dxci
        coef_5vec(i) = ((nx-i)/(nx-1)-1) / dxim  / dxci
      end do

      do j=1,ny
        dyj  = dy(j)
        dyjm = dy(j-1)
        dycj = dyc(j)
        ja_vec(j) = min(j+1,ny  )
        jb_vec(j) = max(  1, j-1)
        coef_6vec(j) = (j/ny-1)          / dyj   / dycj
        coef_7vec(j) = ((ny-j)/(ny-1)-1) / dyjm  / dycj
      end do

      do k=1,nz-1
        dzk  = dz(k)
        dzckp = dzc(k+1)
        dzck  = dzc(k)
        ka_vec(k) = min(k+1,nz-1)
        kb_vec(k) = max(  1, k-1)
        coef_2vec(k) = (k/(nz-1)-1)      / dzckp / dzk
        coef_3vec(k) = ((nz-k)/(nz-1)-1) / dzck  / dzk
      end do

      end subroutine gen_inprod_z_coefs


!======================================================================!
!============================================================! INPROD_Z
!======================================================================!

      subroutine inprod_z(in_x,in_y,in_z,in_p,out_z)
      use general
      use mesh
      use inprod_z_coefs
      
      complex(dp), dimension(nx-1, ny , nz ) :: in_x
      complex(dp), dimension( nx ,ny-1, nz ) :: in_y
      complex(dp), dimension( nx , ny ,nz-1) :: in_z,out_z
      complex(dp), dimension( nx , ny , nz ) :: in_p
      complex(dp), dimension(7)  :: s
 
      real(dp)     :: dxci,dycj,dzk,vol
      complex(dp)  :: k2,kp1,kp2
      integer      :: ia,ib,ja,jb,ka,kb

      complex(dp), dimension(7)  :: coef
      complex(dp), dimension(7)  :: invec

      out_z = 0.d0
!----------------------------------------------------------------------!
!                  Work on the z-directed edges next                   !
!----------------------------------^-----------------------------------!

      do k=1,nz-1
        dzk  = dz(k)
        ka = ka_vec(k)
        kb = kb_vec(k)
        coef(2) = coef_2vec(k)
        coef(3) = coef_3vec(k)

        do j=1,ny
          dycj  = dyc(j)
          ja = ja_vec(j)
          jb = jb_vec(j)
          coef(6) = coef_6vec(j)
          coef(7) = coef_7vec(j)

          do i=1,nx
            dxci  = dxc(i) 
            ia = ia_vec(i)
            ib = ib_vec(i)
            coef(4) = coef_4vec(i)
            coef(5) = coef_5vec(i)

            vol = dxci * dycj * dzk
            k2  = k2_zgrid(i,j,k)
            kp1 = k2_pgrid(i,j, k )
            kp2 = k2_pgrid(i,j,k+1)

            coef(1) = -sum(coef(2:7)) + k2

            invec(1) = in_z( i , j , k )
            invec(2) = in_z( i , j ,ka)
            invec(3) = in_z( i , j ,kb)
            invec(4) = in_z(ia, j , k )
            invec(5) = in_z(ib, j , k )
            invec(6) = in_z( i ,ja, k )
            invec(7) = in_z( i ,jb, k )

            out_z(i,j,k) = ( sum( coef*invec ) +                        &
                              ( in_p(i,j,k+1)*(kp2-k2) -                &
                                in_p(i,j, k )*(kp1-k2)   )/dzk   ) * vol

          end do
        end do
      end do

      end subroutine inprod_z


!======================================================================!
!==================================================! GEN_INPROD_P_COEFS
!======================================================================!

      subroutine gen_inprod_p_coefs
      use general
      use mesh
      use inprod_p_coefs

      real(dp)     :: dxi,dxim,dxci
      real(dp)     :: dyj,dyjm,dycj
      real(dp)     :: dzk,dzkm,dzck

      allocate ( coef_2vec(2:nx-1)  , coef_3vec(2:nx-1) )
      allocate ( coef_4vec(2:ny-1)  , coef_5vec(2:ny-1) )
      allocate ( coef_6vec(2:nz-1)  , coef_7vec(2:nz-1) )
      
      allocate ( ia_vec(2:nx-1)  , ib_vec(2:nx-1)  , ic_vec(2:nx-1) )
      allocate ( ja_vec(2:ny-1)  , jb_vec(2:ny-1)  , jc_vec(2:ny-1) )
      allocate ( ka_vec(2:nz-1)  , kb_vec(2:nz-1)  , kc_vec(2:nz-1) )

      do i=2,nx-1
        dxi  = dx(i) 
        dxim = dx(i-1) 
        dxci = dxc(i)
        ia_vec(i) = min(i+1,nx  )
        ib_vec(i) = max(  1,i-1 )
        ic_vec(i) = min(i  ,nx-1)
        coef_2vec(i) =  (i/nx-1)          /dxi /dxci
        coef_3vec(i) =  ((nx-i)/(nx-1)-1) /dxim/dxci
      end do

      do j=2,ny-1
        dyj  = dy(j)
        dyjm = dy(j-1)
        dycj = dyc(j)
        ja_vec(j) = min(j+1,ny  )
        jb_vec(j) = max(  1,j-1 )
        jc_vec(j) = min(j  ,ny-1  )
        coef_4vec(j) =  (j/ny-1)          /dyj /dycj
        coef_5vec(j) =  ((ny-j)/(ny-1)-1) /dyjm/dycj
      end do

      do k=2,nz-1
        dzk  = dz(k)
        dzkm = dz(k-1)
        dzck = dzc(k)
        ka_vec(k) = min(k+1,nz  )
        kb_vec(k) = max(  1,k-1 )
        kc_vec(k) = min(k  ,nz-1)
        coef_6vec(k) =  (k/nz-1)          /dzk /dzck
        coef_7vec(k) =  ((nz-k)/(nz-1)-1) /dzkm/dzck
      end do

      end subroutine gen_inprod_p_coefs

!======================================================================!
!============================================================! INPROD_P
!======================================================================!

      subroutine inprod_p(in_x,in_y,in_z,in_p,out_p)
      use general
      use mesh
      use inprod_p_coefs
      
      complex(dp), dimension(nx-1, ny , nz ) :: in_x
      complex(dp), dimension( nx ,ny-1, nz ) :: in_y
      complex(dp), dimension( nx , ny ,nz-1) :: in_z
      complex(dp), dimension( nx , ny , nz ) :: in_p,out_p
      complex(dp), dimension(7)  :: s
  
      real(dp)     :: dxci,dycj,dzck,vol
      complex(dp)  :: k2,kp1,kp2
      integer :: ia,ib,ja,jb,ka,kb,ic,jc,kc

      complex(dp), dimension(7)  :: coef
      complex(dp), dimension(7)  :: invec

      out_p = 0.d0

!----------------------------------------------------------------------!
!                 Finish off by working on the nodes                   !
!----------------------------------^-----------------------------------!

      do k=2,nz-1
        dzck  = dzc(k)
        ka = ka_vec(k)
        kb = kb_vec(k)
        kc = kc_vec(k)

        do j=2,ny-1
          dycj  = dyc(j)
          ja = ja_vec(j)
          jb = jb_vec(j)
          jc = jc_vec(j)

          do i=2,nx-1
            dxci  = dxc(i)
            ia = ia_vec(i)
            ib = ib_vec(i)
            ic = ic_vec(i)

            vol = dxci * dycj * dzck
            k2 = k2_pgrid(i,j,k)

            s(1) =  k2_pgrid(i   ,j   ,k  ) 
            s(2) =  k2_xgrid(ic  ,j   ,k  ) 
            s(3) =  k2_xgrid(ib  ,j   ,k  ) 
            s(4) =  k2_ygrid(i   ,jc  ,k  )
            s(5) =  k2_ygrid(i   ,jb  ,k  )
            s(6) =  k2_zgrid(i   ,j   ,kc )
            s(7) =  k2_zgrid(i   ,j   ,kb )

            coef(2) = coef_2vec(i)*s(2)
            coef(3) = coef_3vec(i)*s(3)
            coef(4) = coef_4vec(j)*s(4)
            coef(5) = coef_5vec(j)*s(5)
            coef(6) = coef_6vec(k)*s(6)
            coef(7) = coef_7vec(k)*s(7)
            coef(1) =  -sum(coef(2:7)) + k2**2

            invec(1) = in_p( i , j , k )
            invec(2) = in_p(ia, j , k ) 
            invec(3) = in_p(ib, j , k ) 
            invec(4) = in_p( i ,ja, k ) 
            invec(5) = in_p( i ,jb, k ) 
            invec(6) = in_p( i , j ,ka) 
            invec(7) = in_p( i , j ,kb) 

            out_p(i,j,k) = sum(coef*invec)*vol

            coef(2) =   (s(2) - s(1) ) / dxci 
            coef(3) =  -(s(3) - s(1) ) / dxci
            coef(4) =   (s(4) - s(1) ) / dycj
            coef(5) =  -(s(5) - s(1) ) / dycj
            coef(6) =   (s(6) - s(1) ) / dzck
            coef(7) =  -(s(7) - s(1) ) / dzck

            invec(2) = in_x(ic, j , k )
            invec(3) = in_x(ib, j , k )
            invec(4) = in_y( i ,jc, k )
            invec(5) = in_y( i ,jb, k )
            invec(6) = in_z( i , j ,kc)
            invec(7) = in_z( i , j ,kb)

            out_p(i,j,k) = out_p(i,j,k) + sum(coef(2:7)*invec(2:7))*vol

          end do
        end do
      end do

      end subroutine inprod_p


!======================================================================!
!========================================================! GEN_K2_XGRID
!======================================================================!

      subroutine gen_k2_xgrid

!
!  function to compute the squared wavenumber i-omega-mu-sigma for the 
!  Voronoi grid centered along x-directed edges.  Conductivity estimate
!  along edge is given by simple volume-weighted harmonic average.
!
      use general
      use mesh

      integer     :: i,j,k
      real(dp)    :: dzk,dzkm,dyj,dyjm
      real(dp),    dimension(4) :: coef
      complex(dp), dimension(4) :: s

      k2_xgrid = 1.e30

      do i=1,nx-1
        do j=2,ny-1
          do k=2,nz-1

            dyj  = dy(j)
            dyjm = dy(j-1)
            dzk  = dz(k)
            dzkm = dz(k-1)
            
      
            coef(1) = dyjm * dzkm
            coef(2) = dyj  * dzkm
            coef(3) = dyj  * dzk 
            coef(4) = dyjm * dzk 
      
            s(1) = sig( i ,j-1,k-1)
            s(2) = sig( i , j ,k-1)
            s(3) = sig( i , j , k )
            s(4) = sig( i ,j-1, k )
      
            k2_xgrid(i,j,k)  =  sum(s*coef)/( (dyjm+dyj)*(dzkm+dzk) )
      
          end do
        end do 
      end do
      k2_xgrid (:,2:ny-1,1 ) = k2_xgrid (:,2:ny-1,2) 
      k2_xgrid (:,2:ny-1,nz) = k2_xgrid (:,2:ny-1,nz-1) 
      k2_xgrid (:,1,:)       = k2_xgrid (:,2,:) 
      k2_xgrid (:,ny,:)      = k2_xgrid (:,ny-1,:) 

      k2_xgrid = k2_xgrid * II * omega * MU0

      end subroutine gen_k2_xgrid


!======================================================================!
!========================================================! GEN_K2_YGRID
!======================================================================!

      subroutine gen_k2_ygrid
!
!  function to compute the squared wavenumber i-omega-mu-sigma for the 
!  Voronoi grid centered along y-directed edges.  Conductivity estimate
!  along edge is given by simple volume-weighted harmonic average.
!
      use general
      use mesh

      integer     :: i,j,k
      real(dp)    :: dzk,dzkm,dxi,dxim
      real(dp),    dimension(4) :: coef
      complex(dp), dimension(4) :: s

      do i=2,nx-1
        do j=1,ny-1
          do k=2,nz-1
      
            dxi  = dx(i) 
            dxim = dx(i-1) 
            dzk  = dz(k)
            dzkm = dz(k-1)
      
            coef(1) = dxim * dzkm 
            coef(2) = dxi  * dzkm 
            coef(3) = dxi  * dzk  
            coef(4) = dxim * dzk  
      
            s(1) = sig(i-1, j ,k-1)
            s(2) = sig( i,  j ,k-1)
            s(3) = sig( i,  j , k )
            s(4) = sig(i-1, j , k )
      
            k2_ygrid(i,j,k) = sum(s*coef)/( (dxim+dxi)*(dzkm+dzk) )

          end do
        end do
      end do

      k2_ygrid (2:nx-1,:,1 ) = k2_ygrid (2:nx-1,:,2) 
      k2_ygrid (2:nx-1,:,nz) = k2_ygrid (2:nx-1,:,nz-1) 
      k2_ygrid (1,:,:)       = k2_ygrid (2,:,:) 
      k2_ygrid (nx,:,:)      = k2_ygrid (nx-1,:,:)

      k2_ygrid = k2_ygrid * II * omega * MU0
 
      end subroutine gen_k2_ygrid


!======================================================================!
!========================================================! GEN_K2_ZGRID
!======================================================================!

      subroutine gen_k2_zgrid
!
!  function to compute the squared wavenumber i-omega-mu-sigma for the 
!  Voronoi grid centered along z-directed edges.  Conductivity estimate
!  along edge is given by simple volume-weighted harmonic average.
!
      use general
      use mesh

      integer     :: i,j,k
      real(dp)    :: dxi,dxim,dyj,dyjm
      real(dp),    dimension(4) :: coef
      complex(dp), dimension(4) :: s
      
      do i=2,nx-1
        do j=2,ny-1
          do k=1,nz-1
      
            dxi  = dx(i) 
            dxim = dx(i-1) 
            dyj  = dy(j)
            dyjm = dy(j-1)
       
            coef(1) = dxim * dyjm
            coef(2) = dxi  * dyjm
            coef(3) = dxi  * dyj 
            coef(4) = dxim * dyj 
      
            s(1) = sig(i-1,j-1,k)
            s(2) = sig( i ,j-1,k)
            s(3) = sig( i , j ,k)
            s(4) = sig(i-1, j, k)
      
            k2_zgrid(i,j,k) =  sum( s * coef ) / ( (dxim+dxi)*(dyjm+dyj) ) 

          end do
        end do
      end do

      k2_zgrid (1,2:ny-1,:) = k2_zgrid (2,2:ny-1,:) 
      k2_zgrid (nx,2:ny-1,:) = k2_zgrid (nx-1,2:ny-1,:) 
      k2_zgrid (:,1,:)       = k2_zgrid (:,2,:) 
      k2_zgrid (:,ny,:)      = k2_zgrid (:,ny-1,:)

      k2_zgrid = k2_zgrid * II * omega * MU0

      end subroutine gen_k2_zgrid


!======================================================================!
!========================================================! GEN_K2_PGRID
!======================================================================!

      subroutine gen_k2_pgrid
!
!  function to compute the squared wavenumber i-omega-mu-sigma for the 
!  Voronoi grid centered along primary grid nodes.  Conductivity
!  estimate at node is given by simple volume-weighted harmonic average.
!
      use general
      use mesh

      integer     :: i,j,k
      real(dp)    :: dxi,dxim,dyj,dyjm,dzk,dzkm
      real(dp),    dimension(8) :: coef
      complex(dp), dimension(8) :: s

      do i=2,nx-1
        do j=2,ny-1
          do k=2,nz-1
      
            dxi  = dx(i) 
            dxim = dx(i-1) 
            dyj  = dy(j)
            dyjm = dy(j-1)
            dzk  = dz(k)
            dzkm = dz(k-1)
       
            coef(1) = dxim * dyjm * dzk
            coef(2) = dxi  * dyjm * dzk
            coef(3) = dxi  * dyj  * dzk
            coef(4) = dxim * dyj  * dzk
            coef(5) = dxim * dyjm * dzkm
            coef(6) = dxi  * dyjm * dzkm
            coef(7) = dxi  * dyj  * dzkm
            coef(8) = dxim * dyj  * dzkm
      
            s(1) = sig(i-1,j-1,k)
            s(2) = sig( i ,j-1,k)
            s(3) = sig( i , j ,k)
            s(4) = sig(i-1, j, k)
            s(5) = sig(i-1,j-1,k-1)
            s(6) = sig( i ,j-1,k-1)
            s(7) = sig( i , j ,k-1)
            s(8) = sig(i-1, j, k-1)

            k2_pgrid(i,j,k)= sum(s*coef)/((dxim+dxi)*(dyjm+dyj)*(dzkm+dzk)) 
          end do
        end do
      end do

      k2_pgrid (1  ,2:ny-1,2:nz-1)  = k2_pgrid(2,    2:ny-1,2:nz-1)
      k2_pgrid (nx ,2:ny-1,2:nz-1)  = k2_pgrid(nx-1, 2:ny-1,2:nz-1)

      k2_pgrid (2:nx-1, 1  ,2:nz-1) = k2_pgrid(2:nx-1, 2,   2:nz-1)
      k2_pgrid (2:nx-1, ny ,2:nz-1) = k2_pgrid(2:nx-1, ny-1,2:nz-1)

      k2_pgrid (2:nx-1, 2:ny-1, 1 ) = k2_pgrid(2:nx-1,2:ny-1, 2   )
      k2_pgrid (2:nx-1, 2:ny-1, nz) = k2_pgrid(2:nx-1,2:ny-1, nz-1)

      k2_pgrid( 1, 1,2:nz-1) = k2_pgrid(   2,   2,2:nz-1)
      k2_pgrid( 1,ny,2:nz-1) = k2_pgrid(   2,ny-1,2:nz-1)
      k2_pgrid(nx, 1,2:nz-1) = k2_pgrid(nx-1,   2,2:nz-1)
      k2_pgrid(nx,ny,2:nz-1) = k2_pgrid(nx-1,ny-1,2:nz-1)

      k2_pgrid(2:nx-1, 1, 1) = k2_pgrid(2:nx-1,   2,   2)
      k2_pgrid(2:nx-1,ny, 1) = k2_pgrid(2:nx-1,ny-1,   2)
      k2_pgrid(2:nx-1, 1,nz) = k2_pgrid(2:nx-1,   2,nz-1)
      k2_pgrid(2:nx-1,ny,nz) = k2_pgrid(2:nx-1,ny-1,nz-1)

      k2_pgrid( 1,2:ny-1, 1) = k2_pgrid(   2,2:ny-1,   2)
      k2_pgrid(nx,2:ny-1, 1) = k2_pgrid(nx-1,2:ny-1,   2)
      k2_pgrid( 1,2:ny-1,nz) = k2_pgrid(   2,2:ny-1,nz-1)
      k2_pgrid(nx,2:ny-1,nz) = k2_pgrid(nx-1,2:ny-1,nz-1)

      k2_pgrid( 1, 1,1) = k2_pgrid(   2,   2,2)
      k2_pgrid(nx, 1,1) = k2_pgrid(nx-1,   2,2)
      k2_pgrid( 1,ny,1) = k2_pgrid(   2,ny-1,2)
      k2_pgrid(nx,ny,1) = k2_pgrid(nx-1,ny-1,2)

      k2_pgrid( 1, 1,nz) = k2_pgrid(   2,   2,nz-1)
      k2_pgrid(nx, 1,nz) = k2_pgrid(nx-1,   2,nz-1)
      k2_pgrid( 1,ny,nz) = k2_pgrid(   2,ny-1,nz-1)
      k2_pgrid(nx,ny,nz) = k2_pgrid(nx-1,ny-1,nz-1)

      k2_pgrid = k2_pgrid * II * omega * MU0

      end subroutine gen_k2_pgrid


!======================================================================!
!============================================================! GEN_DIAG
!======================================================================!

      subroutine gen_diag
      use general
      use mesh

      real(dp)     :: dxi,dxim,dxip,dxci,dxcip
      real(dp)     :: dyj,dyjm,dyjp,dycj,dycjp
      real(dp)     :: dzk,dzkm,dzkp,dzck,dzckp
      real(dp)     :: vol
      complex(dp)  :: k2

      complex(dp),    dimension(7)  :: coef
      complex(dp), dimension(7)  :: s
!
!  Allocate storage for the diagonal elements of the coefficient mtx
!
      allocate (diag_x (nx-1, ny , nz ))
      allocate (diag_y ( nx ,ny-1, nz ))
      allocate (diag_z ( nx , ny ,nz-1))
      allocate (diag_p ( nx , ny , nz ))
!
!  Initialize the diagonals with the value unity
!     
      diag_x = 1.d0
      diag_y = 1.d0
      diag_z = 1.d0
      diag_p = 1.d0
!
!  x-edges
!
      do k=1,nz
        dzk  = dz(k) 
        dzkm = dz(k-1) 
        dzck = 0.5d0*(dzk+dzkm)

        coef(6) = (k/nz - 1)        / dzk   / dzck 
        coef(7) = ((nz-k)/(nz-1)-1) / dzkm  / dzck

        do j=1,ny
          dyj  = dy(j)
          dyjm = dy(j-1)
          dycj = 0.5d0*(dyj+dyjm)

          coef(4) = (j/ny - 1)        / dyj   / dycj
          coef(5) = ((ny-j)/(ny-1)-1) / dyjm  / dycj

          do i=1,nx-1
            dxip = dx(i+1) 
            dxi  = dx(i) 
            dxim = dx(i-1) 
            dxcip = 0.5d0*(dxi+dxip)
            dxci  = 0.5d0*(dxi+dxim)
           
            vol = dxi * dycj * dzck
            k2  = k2_xgrid( i, j,k)

            coef(2) = (i/(nx-1)-1)      / dxcip / dxi
            coef(3) = ((nx-i)/(nx-1)-1) / dxci  / dxi

            diag_x(i,j,k) = (-sum(coef(2:7)) + k2)*vol

          end do
        end do
      end do
!
!  y-edges
!

      do k=1,nz
        dzk  = dz(k)
        dzkm = dz(k-1)
        dzck = 0.5d0*(dzk+dzkm)

        coef(6) = (k/nz-1)          / dzk   / dzck
        coef(7) = ((nz-k)/(nz-1)-1) / dzkm  / dzck

        do j=1,ny-1
          dyjp = dy(j+1)
          dyj  = dy(j)
          dyjm = dy(j-1)
          dycjp = 0.5d0*(dyj+dyjp)
          dycj  = 0.5d0*(dyj+dyjm)

          coef(2) = (j/(ny-1)-1)      / dycjp / dyj
          coef(3) = ((ny-j)/(ny-1)-1) / dycj  / dyj

          do i=1,nx
            dxi  = dx(i) 
            dxim = dx(i-1) 
            dxci  = 0.5d0*(dxi+dxim)

            vol = dxci * dyj * dzck
            k2 =  k2_ygrid(i, j, k)

            coef(4) = (i/nx-1)          / dxi   / dxci
            coef(5) = ((nx-i)/(nx-1)-1) / dxim  / dxci

            diag_y(i,j,k) = (-sum(coef(2:7)) + k2)*vol

          end do
        end do
      end do
!
!  z-edges
!

      do k=1,nz-1
        dzkp = dz(k+1)
        dzk  = dz(k)
        dzkm = dz(k-1)
        dzckp = 0.5d0*(dzk+dzkp)
        dzck  = 0.5d0*(dzk+dzkm)

        coef(2) = (k/(nz-1)-1)      / dzckp / dzk
        coef(3) = ((nz-k)/(nz-1)-1) / dzck  / dzk

        do j=1,ny
          dyj  = dy(j)
          dyjm = dy(j-1)
          dycj  = 0.5d0*(dyj+dyjm)
          coef(6) = (j/ny-1)          / dyj   / dycj
          coef(7) = ((ny-j)/(ny-1)-1) / dyjm  / dycj

          do i=1,nx
            dxi  = dx(i) 
            dxim = dx(i-1) 
            dxci  = 0.5d0*(dxi+dxim)

            vol = dxci * dycj * dzk
            k2  = k2_zgrid(i,j,k)

            coef(4) = (i/nx-1)          / dxi   / dxci
            coef(5) = ((nx-i)/(nx-1)-1) / dxim  / dxci

            diag_z(i,j,k) = (-sum(coef(2:7)) + k2) *vol

          end do
        end do
      end do
!
!  p-nodes
!

      do k=2,nz-1
        dzk  = dz(k)
        dzkm = dz(k-1)
        dzck  = 0.5d0*(dzk+dzkm)

        do j=2,ny-1
          dyj  = dy(j)
          dyjm = dy(j-1)
          dycj  = 0.5d0*(dyj+dyjm)

          do i=2,nx-1
            dxi  = dx(i) 
            dxim = dx(i-1) 
            dxci  = 0.5d0*(dxi+dxim)

            vol = dxci * dycj * dzck
            k2 = k2_pgrid(i,j,k)


            s(1) =  k2_pgrid(i  ,j  ,k  ) 
            s(2) =  k2_xgrid(min(i,nx-1)  ,j  ,k  ) 
            s(3) =  k2_xgrid(max(i-1,1),j  ,k  ) 
            s(4) =  k2_ygrid(i  ,min(j,ny-1)  ,k  )
            s(5) =  k2_ygrid(i  ,max(j-1,1),k  )
            s(6) =  k2_zgrid(i  ,j  ,min(k,nz-1)  )
            s(7) =  k2_zgrid(i  ,j  ,max(k-1,1))

            coef(2) =  (i/nx-1)          * s(2)/dxi /dxci
            coef(3) =  ((nx-i)/(nx-1)-1) * s(3)/dxim/dxci
            coef(4) =  (j/ny-1)          * s(4)/dyj /dycj
            coef(5) =  ((ny-j)/(ny-1)-1) * s(5)/dyjm/dycj
            coef(6) =  (k/nz-1)          * s(6)/dzk /dzck
            coef(7) =  ((nz-k)/(nz-1)-1) * s(7)/dzkm/dzck

            diag_p(i,j,k) =  (-sum(coef(2:7)) + k2**2)*vol


          end do
        end do
      end do

      end subroutine gen_diag

!======================================================================!
!==========================================================! DUMP_A_PHI
!======================================================================!

      subroutine dump_a_phi
      use general
      use mesh

      write(iout,*) 'writing out potentials A and Phi...'

      write(11,'(a)') '# vtk DataFile Version 3.0'
      write(11,'(a)') 'vtk output'
      write(11,'(a)') 'ASCII'
      write(11,'(a)') 'DATASET RECTILINEAR_GRID'
      write(11,'(a,3(i3,1x))') 'DIMENSIONS ',nx-2,ny-2,nz-2

      write(11,'(a,i8,a)') 'X_COORDINATES ',nx-2,' float'
      write(11,1) xn(2:nx-1)

      write(11,'(a,i8,a)') 'Y_COORDINATES ',ny-2,' float'
      write(11,1) yn(2:ny-1)

      write(11,'(a,i8,a)') 'Z_COORDINATES ',nz-2,' float'
      write(11,1) zn(2:nz-1)

      write(11,'(a,i8,a)') 'POINT_DATA ' ,(nx-2)*(ny-2)*(nz-2)

      write(11,'(a,i8,a)') 'SCALARS Ax_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) 0.5 * real(rhs_x(1:nx-2,2:ny-1,2:nz-1) +              &
                             rhs_x(2:nx-1,2:ny-1,2:nz-1)   )

      write(11,'(a,i8,a)') 'SCALARS Ax_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) 0.5 * imag(rhs_x(1:nx-2,2:ny-1,2:nz-1) +              &
                             rhs_x(2:nx-1,2:ny-1,2:nz-1)   )

      write(11,'(a,i8,a)') 'SCALARS Ay_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) 0.5 * real(rhs_y(2:nx-1,1:ny-2,2:nz-1) +              &
                             rhs_y(2:nx-1,2:ny-1,2:nz-1)   )

      write(11,'(a,i8,a)') 'SCALARS Ay_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) 0.5 * imag(rhs_y(2:nx-1,1:ny-2,2:nz-1) +              &
                             rhs_y(2:nx-1,2:ny-1,2:nz-1)   )

      write(11,'(a,i8,a)') 'SCALARS Az_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) 0.5 * real(rhs_z(2:nx-1,2:ny-1,1:nz-2) +              &
                             rhs_z(2:nx-1,2:ny-1,2:nz-1)   )

      write(11,'(a,i8,a)') 'SCALARS Az_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) 0.5 * imag(rhs_z(2:nx-1,2:ny-1,1:nz-2) +              &
                             rhs_z(2:nx-1,2:ny-1,2:nz-1)   )

      write(11,'(a,i8,a)') 'SCALARS phi_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( rhs_p(2:nx-1,2:ny-1,2:nz-1) )

      write(11,'(a,i8,a)') 'SCALARS phi_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( rhs_p(2:nx-1,2:ny-1,2:nz-1) )
  
      close(11)

 1    format(1000(e14.7,1x))

      end subroutine dump_a_phi


!======================================================================!
!==========================================================! DUMP_SIGMA
!======================================================================!

      subroutine dump_sigma
      use general
      use mesh

      write(iout,*) 'Writing file: sigma.vtk'
      open(unit=11,file=trim(pwd)//'/sigma.vtk',status='unknown') 

      write(11,'(a)') '# vtk DataFile Version 3.0'
      write(11,'(a)') 'vtk output'
      write(11,'(a)') 'ASCII'
      write(11,'(a)') 'DATASET RECTILINEAR_GRID'
      write(11,'(a,3(i3,1x))') 'DIMENSIONS ',nx,ny,nz

      write(11,'(a,i8,a)') 'X_COORDINATES ',nx,' float'
      write(11,1) xn

      write(11,'(a,i8,a)') 'Y_COORDINATES ',ny,' float'
      write(11,1) yn

      write(11,'(a,i8,a)') 'Z_COORDINATES ',nz,' float'
      write(11,1) zn

      write(11,'(a,i14,a)') 'CELL_DATA ' ,(nx-1)*(ny-1)*(nz-1)

      write(11,'(a,i8,a)') 'SCALARS sig_real float 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real(sig)

      write(11,'(a,i8,a)') 'SCALARS sig_imag float 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag(sig)

      close(11)

 1    format(1000(e14.7,1x))
      end subroutine dump_sigma

!======================================================================!
!=========================================================! DUMP_FIELDS
!======================================================================!

      subroutine dump_fields
      use general
      use mesh

      write(iout,*) 'Writing file: Bx.vtk'
      open(unit=11,file=trim(pwd)//'/Bx.vtk',status='unknown') 

      write(11,'(a)') '# vtk DataFile Version 3.0'
      write(11,'(a)') 'vtk output'
      write(11,'(a)') 'ASCII'
      write(11,'(a)') 'DATASET RECTILINEAR_GRID'
      write(11,'(a,3(i3,1x))') 'DIMENSIONS ',nx,ny-1,nz-1

      write(11,'(a,i8,a)') 'X_COORDINATES ',nx,' float'
      write(11,1) xn

      write(11,'(a,i8,a)') 'Y_COORDINATES ',ny-1,' float'
      write(11,1) 0.5d0*( yn(1:ny-1) + yn(2:ny) )

      write(11,'(a,i8,a)') 'Z_COORDINATES ',nz-1,' float'
      write(11,1) 0.5d0*( zn(1:nz-1) + zn(2:nz) )

      write(11,'(a,i14,a)') 'POINT_DATA ' ,nx*(ny-1)*(nz-1)

      write(11,'(a,i8,a)') 'SCALARS Bx_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real(Bx)

      write(11,'(a,i8,a)') 'SCALARS Bx_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag(Bx)

      close(11)


      write(iout,*) 'Writing file: By.vtk'
      open(unit=11,file=trim(pwd)//'/By.vtk',status='unknown') 

      write(11,'(a)') '# vtk DataFile Version 3.0'
      write(11,'(a)') 'vtk output'
      write(11,'(a)') 'ASCII'
      write(11,'(a)') 'DATASET RECTILINEAR_GRID'
      write(11,'(a,3(i3,1x))') 'DIMENSIONS ',nx-1,ny,nz-1

      write(11,'(a,i8,a)') 'X_COORDINATES ',nx-1,' float'
      write(11,1) 0.5d0*( xn(1:nx-1) + xn(2:nx) )

      write(11,'(a,i8,a)') 'Y_COORDINATES ',ny,' float'
      write(11,1) yn

      write(11,'(a,i8,a)') 'Z_COORDINATES ',nz-1,' float'
      write(11,1) 0.5d0*( zn(1:nz-1) + zn(2:nz) )

      write(11,'(a,i14,a)') 'POINT_DATA ' ,(nx-1)*ny*(nz-1)

      write(11,'(a,i8,a)') 'SCALARS By_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real(By)

      write(11,'(a,i8,a)') 'SCALARS By_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag(By)

      close(11)

      write(iout,*) 'Writing file: Bz.vtk'
      open(unit=11,file=trim(pwd)//'/Bz.vtk',status='unknown') 

      write(11,'(a)') '# vtk DataFile Version 3.0'
      write(11,'(a)') 'vtk output'
      write(11,'(a)') 'ASCII'
      write(11,'(a)') 'DATASET RECTILINEAR_GRID'
      write(11,'(a,3(i3,1x))') 'DIMENSIONS ',nx-1,ny-1,nz

      write(11,'(a,i8,a)') 'X_COORDINATES ',nx-1,' float'
      write(11,1) 0.5d0*( xn(1:nx-1) + xn(2:nx) )

      write(11,'(a,i8,a)') 'Y_COORDINATES ',ny-1,' float'
      write(11,1) 0.5d0*( yn(1:ny-1) + yn(2:ny) )

      write(11,'(a,i8,a)') 'Z_COORDINATES ',nz,' float'
      write(11,1) zn

      write(11,'(a,i14,a)') 'POINT_DATA ' ,(nx-1)*(ny-1)*nz

      write(11,'(a,i8,a)') 'SCALARS Bz_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real(Bz)

      write(11,'(a,i8,a)') 'SCALARS Bz_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag(Bz)

      close(11)


      write(iout,*) 'Writing file: Ex.vtk'
      open(unit=11,file=trim(pwd)//'/Ex.vtk',status='unknown') 

      write(11,'(a)') '# vtk DataFile Version 3.0'
      write(11,'(a)') 'vtk output'
      write(11,'(a)') 'ASCII'
      write(11,'(a)') 'DATASET RECTILINEAR_GRID'
      write(11,'(a,3(i3,1x))') 'DIMENSIONS ',nx-1,ny,nz

      write(11,'(a,i8,a)') 'X_COORDINATES ',nx-1,' float'
      write(11,1) 0.5d0*( xn(1:nx-1) + xn(2:nx) )

      write(11,'(a,i8,a)') 'Y_COORDINATES ',ny,' float'
      write(11,1) yn

      write(11,'(a,i8,a)') 'Z_COORDINATES ',nz,' float'
      write(11,1) zn

      write(11,'(a,i14,a)') 'POINT_DATA ' ,(nx-1)*ny*nz

      write(11,'(a,i8,a)') 'SCALARS Ex_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real(Ex)

      write(11,'(a,i8,a)') 'SCALARS Ex_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag(Ex)

      close(11)

      write(iout,*) 'Writing file: Ey.vtk'
      open(unit=11,file=trim(pwd)//'/Ey.vtk',status='unknown') 

      write(11,'(a)') '# vtk DataFile Version 3.0'
      write(11,'(a)') 'vtk output'
      write(11,'(a)') 'ASCII'
      write(11,'(a)') 'DATASET RECTILINEAR_GRID'
      write(11,'(a,3(i3,1x))') 'DIMENSIONS ',nx,ny-1,nz

      write(11,'(a,i8,a)') 'X_COORDINATES ',nx,' float'
      write(11,1) xn

      write(11,'(a,i8,a)') 'Y_COORDINATES ',ny-1,' float'
      write(11,1) 0.5d0*( yn(1:ny-1) + yn(2:ny) )

      write(11,'(a,i8,a)') 'Z_COORDINATES ',nz,' float'
      write(11,1) zn

      write(11,'(a,i8,a)') 'POINT_DATA ' ,nx*(ny-1)*nz

      write(11,'(a,i8,a)') 'SCALARS Ey_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real(Ey)

      write(11,'(a,i8,a)') 'SCALARS Ey_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag(Ey)

      close(11)

      write(iout,*) 'Writing file: Ez.vtk'
      open(unit=11,file=trim(pwd)//'/Ez.vtk',status='unknown') 

      write(11,'(a)') '# vtk DataFile Version 3.0'
      write(11,'(a)') 'vtk output'
      write(11,'(a)') 'ASCII'
      write(11,'(a)') 'DATASET RECTILINEAR_GRID'
      write(11,'(a,3(i3,1x))') 'DIMENSIONS ',nx,ny,nz-1

      write(11,'(a,i8,a)') 'X_COORDINATES ',nx,' float'
      write(11,1) xn

      write(11,'(a,i8,a)') 'Y_COORDINATES ',ny,' float'
      write(11,1) yn

      write(11,'(a,i8,a)') 'Z_COORDINATES ',nz-1,' float'
      write(11,1) 0.5d0*( zn(1:nz-1) + zn(2:nz) )

      write(11,'(a,i8,a)') 'POINT_DATA ' ,nx*ny*(nz-1)

      write(11,'(a,i8,a)') 'SCALARS Ez_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real(Ez)

      write(11,'(a,i8,a)') 'SCALARS Ez_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag(Ez)

      close(11)
!
!  note: simple *unweighted* interpolation onto grid nodes.
!
      write(iout,*) 'Writing file: EM_grid.vtk'
      open(unit=11,file=trim(pwd)//'/EM_grid.vtk',status='unknown') 

      write(11,'(a)') '# vtk DataFile Version 3.0'
      write(11,'(a)') 'vtk output'
      write(11,'(a)') 'ASCII'
      write(11,'(a)') 'DATASET RECTILINEAR_GRID'
      write(11,'(a,3(i3,1x))') 'DIMENSIONS ',nx-2,ny-2,nz-2


      write(11,'(a,i8,a)') 'X_COORDINATES ',nx-2,' float'
      write(11,1) xn(2:nx-1)

      write(11,'(a,i8,a)') 'Y_COORDINATES ',ny-2,' float'
      write(11,1) yn(2:ny-1)

      write(11,'(a,i8,a)') 'Z_COORDINATES ',nz-2,' float'
      write(11,1) zn(2:nz-1)

      write(11,'(a,i8,a)') 'POINT_DATA ' ,(nx-2)*(ny-2)*(nz-2)

      write(11,'(a,i8,a)') 'SCALARS Ex_grid_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( ( Ex(1:nx-2,2:ny-1,2:nz-1) +                    &
                          Ex(2:nx-1,2:ny-1,2:nz-1)   )/2.)

      write(11,'(a,i8,a)') 'SCALARS Ex_grid_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( ( Ex(1:nx-2,2:ny-1,2:nz-1) +                    &
                          Ex(2:nx-1,2:ny-1,2:nz-1)   )/2.)

      write(11,'(a,i8,a)') 'SCALARS Ey_grid_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( ( Ey(2:nx-1,1:ny-2,2:nz-1) +                    &
                          Ey(2:nx-1,2:ny-1,2:nz-1)   )/2.)

      write(11,'(a,i8,a)') 'SCALARS Ey_grid_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( ( Ey(2:nx-1,1:ny-2,2:nz-1) +                    &
                          Ey(2:nx-1,2:ny-1,2:nz-1)   )/2.)

      write(11,'(a,i8,a)') 'SCALARS Ez_grid_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( ( Ez(2:nx-1,2:ny-1,1:nz-2) +                    &
                          Ez(2:nx-1,2:ny-1,2:nz-1)   )/2.)

      write(11,'(a,i8,a)') 'SCALARS Ez_grid_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( ( Ez(2:nx-1,2:ny-1,1:nz-2) +                    &
                          Ez(2:nx-1,2:ny-1,2:nz-1)   )/2.)

      write(11,'(a,i8,a)') 'SCALARS Jx_grid_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( (                                               &
            Ex(1:nx-2,2:ny-1,2:nz-1)*k2_xgrid(1:nx-2,2:ny-1,2:nz-1) +   &
            Ex(2:nx-1,2:ny-1,2:nz-1)*k2_xgrid(2:nx-1,2:ny-1,2:nz-1) ) / &
                                              (2. * II * omega * MU0) )

      write(11,'(a,i8,a)') 'SCALARS Jx_grid_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( (                                               &
            Ex(1:nx-2,2:ny-1,2:nz-1)*k2_xgrid(1:nx-2,2:ny-1,2:nz-1) +   &
            Ex(2:nx-1,2:ny-1,2:nz-1)*k2_xgrid(2:nx-1,2:ny-1,2:nz-1) )/  &
                                              (2. * II * omega * MU0) )

      write(11,'(a,i8,a)') 'SCALARS Jy_grid_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( (                                               &
            Ey(2:nx-1,1:ny-2,2:nz-1)*k2_ygrid(2:nx-1,1:ny-2,2:nz-1) +   &
            Ey(2:nx-1,2:ny-1,2:nz-1)*k2_ygrid(2:nx-1,2:ny-1,2:nz-1) )/  &
                                              (2. * II * omega * MU0) )

      write(11,'(a,i8,a)') 'SCALARS Jy_grid_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( (                                               &
            Ey(2:nx-1,1:ny-2,2:nz-1)*k2_ygrid(2:nx-1,1:ny-2,2:nz-1) +   &
            Ey(2:nx-1,2:ny-1,2:nz-1)*k2_ygrid(2:nx-1,2:ny-1,2:nz-1) )/  &
                                              (2. * II * omega * MU0) )

      write(11,'(a,i8,a)') 'SCALARS Jz_grid_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( (                                               &
            Ez(2:nx-1,2:ny-1,2:nz-1)*k2_zgrid(2:nx-1,2:ny-1,2:nz-1) +   &
            Ez(2:nx-1,2:ny-1,1:nz-2)*k2_zgrid(2:nx-1,2:ny-1,1:nz-2) )/  &
                                              (2. * II * omega * MU0) )

      write(11,'(a,i8,a)') 'SCALARS Jz_grid_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( (                                               &
            Ez(2:nx-1,2:ny-1,2:nz-1)*k2_zgrid(2:nx-1,2:ny-1,2:nz-1) +   &
            Ez(2:nx-1,2:ny-1,1:nz-2)*k2_zgrid(2:nx-1,2:ny-1,1:nz-2) )/  &
                                              (2. * II * omega * MU0) )

      write(11,'(a,i8,a)') 'SCALARS Bx_grid_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( (                                               &
               Bx(2:nx-1,1:ny-2,1:nz-2) + Bx(2:nx-1,2:ny-1,1:nz-2) +   &
               Bx(2:nx-1,1:ny-2,2:nz-1) + Bx(2:nx-1,2:ny-1,2:nz-1) )/4.)

      write(11,'(a,i8,a)') 'SCALARS Bx_grid_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( (                                               &
               Bx(2:nx-1,1:ny-2,1:nz-2) + Bx(2:nx-1,2:ny-1,1:nz-2) +   &
               Bx(2:nx-1,1:ny-2,2:nz-1) + Bx(2:nx-1,2:ny-1,2:nz-1) )/4.)

      write(11,'(a,i8,a)') 'SCALARS By_grid_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( (                                               &
               By(1:nx-2,2:ny-1,1:nz-2) + By(2:nx-1,2:ny-1,1:nz-2) +   &
               By(2:nx-1,2:ny-1,1:nz-2) + By(2:nx-1,2:ny-1,2:nz-1) )/4.)

      write(11,'(a,i8,a)') 'SCALARS By_grid_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( (                                               &
               By(1:nx-2,2:ny-1,1:nz-2) + By(2:nx-1,2:ny-1,1:nz-2) +    &
               By(2:nx-1,2:ny-1,1:nz-2) + By(2:nx-1,2:ny-1,2:nz-1) )/4.)

      write(11,'(a,i8,a)') 'SCALARS Bz_grid_real double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) real( (                                               &
               Bz(1:nx-2,1:ny-2,2:nz-1) + Bz(1:nx-2,2:ny-1,2:nz-1) +    &
               Bz(2:nx-1,1:ny-2,2:nz-1) + Bz(2:nx-1,2:ny-1,2:nz-1) )/4.)

      write(11,'(a,i8,a)') 'SCALARS Bz_grid_imag double 1'
      write(11,'(a,i8,a)') 'LOOKUP_TABLE default'
      write(11,1) imag( (                                               &
               Bz(1:nx-2,1:ny-2,2:nz-1) + Bz(1:nx-2,2:ny-1,2:nz-1) +    &
               Bz(2:nx-1,1:ny-2,2:nz-1) + Bz(2:nx-1,2:ny-1,2:nz-1) )/4.)


      close(11)

 1    format(1000(e14.7,1x))

      end subroutine dump_fields


!======================================================================!
!==========================================================! GEN_FIELDS
!======================================================================!

      subroutine gen_fields
      use general
      use mesh

      allocate ( Bx(  nx  , ny-1 , nz-1 ) )
      allocate ( By( nx-1 ,  ny  , nz-1 ) )
      allocate ( Bz( nx-1 , ny-1 ,  nz  ) )

      allocate ( Ex( nx-1 ,  ny  ,  nz  ) )
      allocate ( Ey(  nx  , ny-1 ,  nz  ) )
      allocate ( Ez(  nx  ,  ny  , nz-1 ) )
!
!  Compute B at the centroid of primary cell faces: B = curl(A)
!
      Bx = 0.d0
      By = 0.d0
      Bz = 0.d0

      do k=1,nz-1
        do j=1,ny-1
          do i=1,nx
            Bx(i,j,k) = rhs_z(i,j+1,k)/dy(j) - rhs_y(i,j,k+1)/dz(k)     &
                       -rhs_z(i, j ,k)/dy(j) + rhs_y(i,j, k )/dz(k)  
          end do
        end do
      end do

      do k=1,nz-1
        do j=1,ny
          do i=1,nx-1
            By(i,j,k) = rhs_x(i,j,k+1)/dz(k) - rhs_z(i+1,j,k)/dx(i)     &
                       -rhs_x(i,j, k )/dz(k) + rhs_z( i ,j,k)/dx(i)  
          end do
        end do
      end do

      do k=1,nz
        do j=1,ny-1
          do i=1,nx-1
            Bz(i,j,k) = rhs_y(i+1,j,k)/dx(i) - rhs_x(i,j+1,k)/dy(j)     &
                       -rhs_y( i ,j,k)/dx(i) + rhs_x(i, j ,k)/dy(j)  
          end do
        end do
      end do

!
!  compute E at primary grid edges: E = -i*omega(A-grad(Phi))
!
      Ex = 0.d0
      Ey = 0.d0
      Ez = 0.d0

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            Ex(i,j,k) = -II*omega*( rhs_x(i, j ,k)      -               &
                             (rhs_p(i+1,j,k)-rhs_p(i,j,k))/dx(i) )
          end do
        end do
      end do

      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            Ey(i,j,k) = -II*omega*( rhs_y(i, j ,k)      -               &
                             (rhs_p(i,j+1,k)-rhs_p(i,j,k))/dy(j) )
          end do
        end do
      end do

      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            Ez(i,j,k) = -II*omega*( rhs_z(i, j ,k)      -               &
                             (rhs_p(i,j,k+1)-rhs_p(i,j,k))/dz(k) )
          end do
        end do
      end do

      end subroutine gen_fields
