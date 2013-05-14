!----------------------------------------------------------------------!
!-----------------------------------------------------------------! QMR
!----------------------------------------------------------------------!

      subroutine qmr

!----------------------------------------------------------------------!
!                                                                      !
!  Routine to solve a complex symmetric system of linear equations     !
!  using the Quasi-Minimal-Residual (QMR) method.  The specific        !
!  algorithm used here is the coupled two-term recurrence scheme       !
!  without 'look-ahead' outlined in Freund and Nachtigal, SIAM J Sci.  !
!  Comput., v15, 313-337 (1994), algorithm 8.1.  Routine implements a  !
!  diagonal preconditioning matrix into the Lanczos iterates.          !
!  Multiplication by the preconditioner of the right-hand-side vector  !
!  upon entering or multiplication of the solution upon exiting the    !
!  routine are NOT necessary.  Nor is it necessary to modify the       !
!  routine 'inprod' which computes the action of the FD matrix upon an !
!  arbitrary vector.                                                   !
!                                                                      !
!  Written by: C J Weiss, Virginia Tech Dept of Geosciences, 2012-13   !
!                                                                      !
!----------------------------------------------------------------------!


      use general
      use mesh
      include 'mpif.h'

      complex(dp), dimension(:,:,:) :: inv_diag_x
      complex(dp), dimension(:,:,:) :: tmp_x,r_x,s_x
      complex(dp), dimension(:,:,:) :: p_x,d_x,vt_x,v_x

      complex(dp), dimension(:,:,:) :: inv_diag_y
      complex(dp), dimension(:,:,:) :: tmp_y,r_y,s_y
      complex(dp), dimension(:,:,:) :: p_y,d_y,vt_y,v_y

      complex(dp), dimension(:,:,:) :: inv_diag_z
      complex(dp), dimension(:,:,:) :: tmp_z,r_z,s_z
      complex(dp), dimension(:,:,:) :: p_z,d_z,vt_z,v_z

      complex(dp), dimension(:,:,:) :: inv_diag_p
      complex(dp), dimension(:,:,:) :: tmp_p,r_p,s_p
      complex(dp), dimension(:,:,:) :: p_p,d_p,vt_p,v_p

      allocatable  :: inv_diag_x
      allocatable  :: tmp_x,r_x,s_x
      allocatable  :: p_x,d_x,vt_x,v_x

      allocatable  :: inv_diag_y
      allocatable  :: tmp_y,r_y,s_y
      allocatable  :: p_y,d_y,vt_y,v_y

      allocatable  :: inv_diag_z
      allocatable  :: tmp_z,r_z,s_z
      allocatable  :: p_z,d_z,vt_z,v_z

      allocatable  :: inv_diag_p
      allocatable  :: tmp_p,r_p,s_p
      allocatable  :: p_p,d_p,vt_p,v_p

      integer      :: ndf
      complex(DP)  :: eps,del,beta,eta,k1
      real(DP)     :: k2,rho_n,rho_np
      real(DP)     :: c_nm,c_n,theta_n,theta_nm
      integer      :: iter,nunk
      real(DP)     :: t1,t2
      real(DP)     :: rnorm,rhsnorm
      logical      :: converge_flag

! Variables !----------------------------------------------------------!
!                                                                      !
!  rhs:            Upon entering, the right hand side for the linear   !
!                    system;  once inside the main loop, the current   !
!                    solution estimate.                                !
!  diag:           Diagonal of coefficient matrix.                     !
!  inv_diag:       Inverse of the diagonal of coefficient matrix, used !
!                    here as a Jacobi-style preconditioner.            !
!  nx,ny,nz:       Number of mesh nodes in the x, y and z directions.  !
!  xn,yn,zn:       Node coordinates for the x, y and z directions.     !
!  sig:            Array of mesh cell conductivities.                  ![S/m]
!  freq:           Frequency of the source.                            ![Hz]
!  ndf:            Number of degrees of freedom in the linear system.  !
!  d:              Direction vector for solution updates.              !
!  vt:             Unnormalized Lanczos vector.                        !
!  v,p:            The two 'Lanczos vectors'.                          !
!  tmp:            Action of the coeficient matrix on a given vector.  !
!  r:              Residual vector.                                    !
!  s:              Supplementary vector used to compute residual via   !
!                    vector updates instead of matrix-vector products  !
!                    at each interation.                               !
!  k1,k2:          Convenient constants for complex double and real    !
!                  double precision, respectively.                     !
!  eps:                 \                                              !
!  del:                  |                                             !
!  beta:                 |  See Freund's paper...                      !
!  eta:                  |                                             !
!  theta_nm,theta_n:    /                                              !
!  rho_n,rho_np:   L2 norms of the 'v' Lanczos vector for interations  !
!                    n and n+1.                                        !
!  c_nm,cn:        Elements of the Given's rotation matrix used in the !
!                    QR decomposition of the coefficient matrix for    !
!                    iterations n and n+1.                             !
!  iter:           Iteration counter.                                  !
!  dtime,time:     System call for elapsed time in QMR routine.        !
!  t1,t2:          Starting and ending times for QMR routine.          !
!  rnorm:          L2 norm of the residual vector.                     !
!  rhsnorm:        L2 norm of the RHS vector.                          !
!                                                                      !
!----------------------------------------------------------------------!

!
!  Allocate storage for the local arrays used in the QMR routine.
!

      allocate (inv_diag_x(nx-1,  ny  ,  nz ))
      allocate (       p_x(nx-1,  ny  ,  nz ))
      allocate (       d_x(nx-1,  ny  ,  nz ))
      allocate (      vt_x(nx-1,  ny  ,  nz ))
      allocate (       v_x(nx-1,  ny  ,  nz ))
      allocate (     tmp_x(nx-1,  ny  ,  nz ))
      allocate (       r_x(nx-1,  ny  ,  nz ))
      allocate (       s_x(nx-1,  ny  ,  nz ))

      allocate (inv_diag_y( nx , ny-1 ,  nz ))
      allocate (       p_y( nx , ny-1 ,  nz ))
      allocate (       d_y( nx , ny-1 ,  nz ))
      allocate (      vt_y( nx , ny-1 ,  nz ))
      allocate (       v_y( nx , ny-1 ,  nz ))
      allocate (     tmp_y( nx , ny-1 ,  nz ))
      allocate (       r_y( nx , ny-1 ,  nz ))
      allocate (       s_y( nx , ny-1 ,  nz ))

      allocate (inv_diag_z( nx ,  ny  , nz-1))
      allocate (       p_z( nx ,  ny  , nz-1))
      allocate (       d_z( nx ,  ny  , nz-1))
      allocate (      vt_z( nx ,  ny  , nz-1))
      allocate (       v_z( nx ,  ny  , nz-1))
      allocate (     tmp_z( nx ,  ny  , nz-1))
      allocate (       r_z( nx ,  ny  , nz-1))
      allocate (       s_z( nx ,  ny  , nz-1))

      allocate (inv_diag_p( nx ,  ny  ,  nz ))
      allocate (       p_p( nx ,  ny  ,  nz ))
      allocate (       d_p( nx ,  ny  ,  nz ))
      allocate (      vt_p( nx ,  ny  ,  nz ))
      allocate (       v_p( nx ,  ny  ,  nz ))
      allocate (     tmp_p( nx ,  ny  ,  nz ))
      allocate (       r_p( nx ,  ny  ,  nz ))
      allocate (       s_p( nx ,  ny  ,  nz ))

!
!  Compute the number of unknowns in the linear system. 
!  This parameter is computed only as a courtesey to the 
!  user and is not used explicitly in the main QMR algorithm.
!
      nunk = (nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1) + nx*ny*nz
!
!  Initialize time counter...
!
      t1 = mpi_wtime(ierr)
!
!  Open the output file.
!
      open(unit=11,file=trim(pwd)//'/qmr.iter',status='unknown')
!
!  Compute the norm of the RHS vector
!
      rhsnorm = sum(rhs_x*dconjg(rhs_x)) + sum(rhs_y*dconjg(rhs_y)) +   &
                sum(rhs_z*dconjg(rhs_z)) + sum(rhs_p*dconjg(rhs_p))   

      rhsnorm = dsqrt(rhsnorm)
!
!  Compute the inverse of the diagonal.  Do this once since a complex
!  double presion multiply of the inverse is FAR less expensive than
!  computing the quotient. 
!
      inv_diag_x = 1.d0/diag_x
      inv_diag_y = 1.d0/diag_y
      inv_diag_z = 1.d0/diag_z
      inv_diag_p = 1.d0/diag_p
!
!  Politely introduce yourself.
!
      write(iout,*) ' '
      write(iout,*) '****************************************'
      write(iout,*) '*  Entering single-thread QMR Routine  *'
      write(iout,*) '****************************************'
      write(iout,*) ' '
 
      write(iout,*) '---------------------------------'
      write(iout,*) '  2-term coupled recurrence QMR  '
      write(iout,*) '    with Jacobi preconditioner   '
      write(iout,*) '================================='
      write(iout,*) ' '                                  
      write(iout,'(a,1x,e12.5)') ' Target residual norm: ',qmr_target*rhsnorm
      write(iout,'(a,1x,e12.5)') ' RHS norm:             ',rhsnorm
      write(iout,*) 'Maximum number of iterations: ',qmr_max_it
      write(iout,*) ' '
      write(iout,*) 'A list of iterate numbers and corresponding'
      write(iout,*) 'residual norms is located in the file: qmr.iter'
      write(iout,*) ' '
 
      write(iout,99) 'Number of unknowns:  ', nunk
      write(iout,*) '  '
 99   format(1x,a21,1x,i8)                           

!
!  Initialize some vectors
!
      p_x = (0.d0,0.d0)
      p_y = (0.d0,0.d0)
      p_z = (0.d0,0.d0)
      p_p = (0.d0,0.d0)

      d_x = (0.d0,0.d0)
      d_y = (0.d0,0.d0)
      d_z = (0.d0,0.d0)
      d_p = (0.d0,0.d0)
!
!  Squirrel away the RHS vector into the vector 's'
!  since 's' isn't used until the main loop.  Store in
!  the vector 'rhs' the initial guess at the solution.
!
      s_x = rhs_x
      s_y = rhs_y
      s_z = rhs_z
      s_p = rhs_p

      rhs_x = (0.d0,0.d0)
      rhs_y = (0.d0,0.d0)
      rhs_z = (0.d0,0.d0)
      rhs_p = (0.d0,0.d0)

!----------------------------! ATTENTION !-----------------------------!
!                                                                      !
!  From this point forward in the algorithm, the vector 'rhs' will     !
!  contain the current solution.                                       !
!                                                                      !
!----------------------------------------------------------------------!       
 
!
!  Compute the initial residual vector.
!
      call inprod_x(rhs_x,rhs_y,rhs_z,rhs_p,tmp_x)
      call inprod_y(rhs_x,rhs_y,rhs_z,rhs_p,tmp_y)
      call inprod_z(rhs_x,rhs_y,rhs_z,rhs_p,tmp_z)
      call inprod_p(rhs_x,rhs_y,rhs_z,rhs_p,tmp_p)

      vt_x = s_x - tmp_x
      vt_y = s_y - tmp_y
      vt_z = s_z - tmp_z
      vt_p = s_p - tmp_p

      r_x = s_x - tmp_x
      r_y = s_y - tmp_y
      r_z = s_z - tmp_z
      r_p = s_p - tmp_p

      rho_np = sum(r_x*dconjg(r_x)) + sum(r_y*dconjg(r_y)) +            &
               sum(r_z*dconjg(r_z)) + sum(r_p*dconjg(r_p))

      rho_np = dsqrt(rho_np)
 
      write(iout,'(a,1x,e14.5)') ' Initial residual: ',rho_np
      write(iout,*) ' '

      v_x = vt_x/rho_np
      v_y = vt_y/rho_np
      v_z = vt_z/rho_np
      v_p = vt_p/rho_np
!
!  Now, set the vector 's' to its proper starting value...
!
      s_x = (0.d0,0.d0)
      s_y = (0.d0,0.d0)
      s_z = (0.d0,0.d0)
      s_p = (0.d0,0.d0)
!
!  ... and initialize these chaps as well.
!
      c_n      =  1.d0
      eps      = ( 1.d0, 0.d0)
      theta_n  =  0.d0
      eta      = (-1.d0, 0.d0)

!----------------------------------------------------------------------!
!                       Start main iterative loop                      !
!----------------------------------------------------------------------!             
 
      write(iout,*) 'Now computing the QMR iterates...'
      write(iout,*) ' '

      converge_flag = .false.
      iter = 1

100   continue

      c_nm     = c_n
      rho_n    = rho_np
      theta_nm = theta_n
      del = (0.d0,0.d0)
!
! step 1
!
      tmp_x = v_x*inv_diag_x
      tmp_y = v_y*inv_diag_y
      tmp_z = v_z*inv_diag_z
      tmp_p = v_p*inv_diag_p

      del = sum(v_x*tmp_x)+sum(v_y*tmp_y)+sum(v_z*tmp_z)+sum(v_p*tmp_p) 

      k1 = rho_n * del / eps
!
! step 2
!
      p_x = tmp_x - p_x*k1
      p_y = tmp_y - p_y*k1
      p_z = tmp_z - p_z*k1
      p_p = tmp_p - p_p*k1
!
! step 3
!
      call inprod_x(p_x,p_y,p_z,p_p,tmp_x)
      call inprod_y(p_x,p_y,p_z,p_p,tmp_y)
      call inprod_z(p_x,p_y,p_z,p_p,tmp_z)
      call inprod_p(p_x,p_y,p_z,p_p,tmp_p)

      eps = sum(p_x*tmp_x)+sum(p_y*tmp_y)+sum(p_z*tmp_z)+sum(p_p*tmp_p) 
      beta = eps / del
 
      vt_x = tmp_x - v_x*beta
      vt_y = tmp_y - v_y*beta
      vt_z = tmp_z - v_z*beta
      vt_p = tmp_p - v_p*beta

      rho_np = sqrt(  sum(vt_x*conjg(vt_x)) + sum(vt_y*conjg(vt_y))     &
                    + sum(vt_z*conjg(vt_z)) + sum(vt_p*conjg(vt_p))  )
!
! step 4
!
      theta_n = rho_np*rho_np/c_nm/(dble(beta)**2+dimag(beta)**2)
      c_n = 1.d0/(1.d0 + theta_n)
      eta = -eta*rho_n*c_n/(beta*c_nm)
      k2 = theta_nm*c_n

      d_x   =   p_x*eta + d_x*k2
      rhs_x = rhs_x     + d_x
      s_x   = tmp_x*eta + s_x*k2
      r_x   =   r_x     - s_x
      v_x   = vt_x/rho_np

      d_y   =   p_y*eta + d_y*k2
      rhs_y = rhs_y     + d_y
      s_y   = tmp_y*eta + s_y*k2
      r_y   =   r_y     - s_y
      v_y   = vt_y/rho_np

      d_z   =   p_z*eta + d_z*k2
      rhs_z = rhs_z     + d_z
      s_z   = tmp_z*eta + s_z*k2
      r_z   =   r_z     - s_z
      v_z   = vt_z/rho_np

      d_p   =   p_p*eta + d_p*k2
      rhs_p = rhs_p     + d_p
      s_p   = tmp_p*eta + s_p*k2
      r_p   =   r_p     - s_p
      v_p   = vt_p/rho_np

      rnorm = sqrt(  sum(r_x*conjg(r_x)) + sum(r_y*conjg(r_y))          &
                   + sum(r_z*conjg(r_z)) + sum(r_p*conjg(r_p))  )
!
! step 5
!
      if (mod(iter,10).eq.0) then
        t2 = mpi_wtime(ierr)
        write(iout,98) iter, qmr_max_it,rnorm,t2-t1
      endif

      write(11,*) iter,rnorm,t2-t1

!
! Check to see if the solution has converged...
!
      if (rnorm.lt.qmr_target*rhsnorm) converge_flag = .true.
      if ((.not.converge_flag).and.(iter.le.qmr_max_it)) iter = iter + 1


      if ((.not.converge_flag).and.(iter.le.qmr_max_it)) goto 100

          
      t2 = mpi_wtime(ierr)
      if (converge_flag) then 
        write(iout,*) 'Solution converged in ',iter,' iterations'
      else
        write(iout,*) '********** WARNING ************'
        write(iout,*) '  Maximum iterations exceeded  '
        write(iout,*) '*******************************'

      end if

      write(iout,'(a,1x,e11.4)') '   Residual:',rnorm
      write(iout,*) ' '
      write(iout,*) 'Time spent in QMR [s]: ',t2 - t1
      write(iout,*) 'Seconds per iteration:    ',(t2-t1)/real(iter)
      write(iout,*) ' '
      write(iout,*) '  '
      close(11)
!
!  Deallocate storage for the local arrays used in the QMR routine.
! 
      deallocate (inv_diag_x,p_x,d_x,vt_x,v_x)
      deallocate (tmp_x,r_x,s_x)
      deallocate (inv_diag_y,p_y,d_y,vt_y,v_y)
      deallocate (tmp_y,r_y,s_y)
      deallocate (inv_diag_z,p_z,d_z,vt_z,v_z)
      deallocate (tmp_z,r_z,s_z)
      deallocate (inv_diag_p,p_p,d_p,vt_p,v_p)
      deallocate (tmp_p,r_p,s_p)
!
!  Exit gracefully
!
      write(iout,*) ' '
      write(iout,*) '*************************'
      write(iout,*) '*  Exiting QMR Routine  *'
      write(iout,*) '*************************'
      write(iout,*) ' '
 
 
 98   format('it: ',i6,'/',i6,'    r: ',e13.6,'    t:',e13.6,1x,f8.4)
 101  format(10(f8.5,1x))

      end


!----------------------------------------------------------------------!
!---------------------------------------------------------! QMR_MTHREAD
!----------------------------------------------------------------------!

      subroutine qmr_mthread

!----------------------------------------------------------------------!
!                                                                      !
!  Routine to solve a complex symmetric system of linear equations     !
!  using the Quasi-Minimal-Residual (QMR) method.  The specific        !
!  algorithm used here is the coupled two-term recurrence scheme       !
!  without 'look-ahead' outlined in Freund and Nachtigal, SIAM J Sci.  !
!  Comput., v15, 313-337 (1994), algorithm 8.1.  Routine implements a  !
!  diagonal preconditioning matrix into the Lanczos iterates.          !
!  Multiplication by the preconditioner of the right-hand-side vector  !
!  upon entering or multiplication of the solution upon exiting the    !
!  routine are NOT necessary.  Nor is it necessary to modify the       !
!  routine 'inprod' which computes the action of the FD matrix upon an !
!  arbitrary vector.                                                   !
!                                                                      !
!  Written by: C J Weiss, Virginia Tech Dept of Geosciences, 2012-13   !
!                                                                      !
!----------------------------------------------------------------------!


      use general
      use mesh
      use omp_lib

      complex(dp), dimension(:,:,:) :: inv_diag_x
      complex(dp), dimension(:,:,:) :: tmp_x,r_x,s_x
      complex(dp), dimension(:,:,:) :: p_x,d_x,vt_x,v_x

      complex(dp), dimension(:,:,:) :: inv_diag_y
      complex(dp), dimension(:,:,:) :: tmp_y,r_y,s_y
      complex(dp), dimension(:,:,:) :: p_y,d_y,vt_y,v_y

      complex(dp), dimension(:,:,:) :: inv_diag_z
      complex(dp), dimension(:,:,:) :: tmp_z,r_z,s_z
      complex(dp), dimension(:,:,:) :: p_z,d_z,vt_z,v_z

      complex(dp), dimension(:,:,:) :: inv_diag_p
      complex(dp), dimension(:,:,:) :: tmp_p,r_p,s_p
      complex(dp), dimension(:,:,:) :: p_p,d_p,vt_p,v_p

      allocatable  :: inv_diag_x
      allocatable  :: tmp_x,r_x,s_x
      allocatable  :: p_x,d_x,vt_x,v_x

      allocatable  :: inv_diag_y
      allocatable  :: tmp_y,r_y,s_y
      allocatable  :: p_y,d_y,vt_y,v_y

      allocatable  :: inv_diag_z
      allocatable  :: tmp_z,r_z,s_z
      allocatable  :: p_z,d_z,vt_z,v_z

      allocatable  :: inv_diag_p
      allocatable  :: tmp_p,r_p,s_p
      allocatable  :: p_p,d_p,vt_p,v_p

      integer      :: ndf
      complex(DP)  :: eps,del,beta,eta,k1
      real(DP)     :: k2,rho_n,rho_np
      real(DP)     :: c_nm,c_n,theta_n,theta_nm
      integer      :: iter,nunk
      real(DP)     :: t1,t2
      real(DP)     :: t11,t22,t(6)
      real(DP)     :: rnorm,rhsnorm
      real(DP)     :: psum_real(0:3)
      complex(DP)  :: psum_complex(0:3)
      complex(DP)  :: psum
      logical      :: converge_flag

! Variables !----------------------------------------------------------!
!                                                                      !
!  rhs:            Upon entering, the right hand side for the linear   !
!                    system;  once inside the main loop, the current   !
!                    solution estimate.                                !
!  diag:           Diagonal of coefficient matrix.                     !
!  inv_diag:       Inverse of the diagonal of coefficient matrix, used !
!                    here as a Jacobi-style preconditioner.            !
!  nx,ny,nz:       Number of mesh nodes in the x, y and z directions.  !
!  xn,yn,zn:       Node coordinates for the x, y and z directions.     !
!  sig:            Array of mesh cell conductivities.                  ![S/m]
!  freq:           Frequency of the source.                            ![Hz]
!  ndf:            Number of degrees of freedom in the linear system.  !
!  d:              Direction vector for solution updates.              !
!  vt:             Unnormalized Lanczos vector.                        !
!  v,p:            The two 'Lanczos vectors'.                          !
!  tmp:            Action of the coeficient matrix on a given vector.  !
!  r:              Residual vector.                                    !
!  s:              Supplementary vector used to compute residual via   !
!                    vector updates instead of matrix-vector products  !
!                    at each interation.                               !
!  k1,k2:          Convenient constants for complex double and real    !
!                  double precision, respectively.                     !
!  eps:                 \                                              !
!  del:                  |                                             !
!  beta:                 |  See Freund's paper...                      !
!  eta:                  |                                             !
!  theta_nm,theta_n:    /                                              !
!  rho_n,rho_np:   L2 norms of the 'v' Lanczos vector for interations  !
!                    n and n+1.                                        !
!  c_nm,cn:        Elements of the Given's rotation matrix used in the !
!                    QR decomposition of the coefficient matrix for    !
!                    iterations n and n+1.                             !
!  iter:           Iteration counter.                                  !
!  dtime,time:     System call for elapsed time in QMR routine.        !
!  t1,t2:          Starting and ending times for QMR routine.          !
!  rnorm:          L2 norm of the residual vector.                     !
!  rhsnorm:        L2 norm of the RHS vector.                          !
!  nunk:           Number of unknowns, less than ndf since the matrix  !
!                    contains several rows with only unity on the      !
!                    diagonal.                                         !
!----------------------------------------------------------------------!

!
!  Allocate storage for the local arrays used in the QMR routine.
!

      allocate (inv_diag_x(nx-1,  ny  ,  nz ))
      allocate (       p_x(nx-1,  ny  ,  nz ))
      allocate (       d_x(nx-1,  ny  ,  nz ))
      allocate (      vt_x(nx-1,  ny  ,  nz ))
      allocate (       v_x(nx-1,  ny  ,  nz ))
      allocate (     tmp_x(nx-1,  ny  ,  nz ))
      allocate (       r_x(nx-1,  ny  ,  nz ))
      allocate (       s_x(nx-1,  ny  ,  nz ))

      allocate (inv_diag_y( nx , ny-1 ,  nz ))
      allocate (       p_y( nx , ny-1 ,  nz ))
      allocate (       d_y( nx , ny-1 ,  nz ))
      allocate (      vt_y( nx , ny-1 ,  nz ))
      allocate (       v_y( nx , ny-1 ,  nz ))
      allocate (     tmp_y( nx , ny-1 ,  nz ))
      allocate (       r_y( nx , ny-1 ,  nz ))
      allocate (       s_y( nx , ny-1 ,  nz ))

      allocate (inv_diag_z( nx ,  ny  , nz-1))
      allocate (       p_z( nx ,  ny  , nz-1))
      allocate (       d_z( nx ,  ny  , nz-1))
      allocate (      vt_z( nx ,  ny  , nz-1))
      allocate (       v_z( nx ,  ny  , nz-1))
      allocate (     tmp_z( nx ,  ny  , nz-1))
      allocate (       r_z( nx ,  ny  , nz-1))
      allocate (       s_z( nx ,  ny  , nz-1))

      allocate (inv_diag_p( nx ,  ny  ,  nz ))
      allocate (       p_p( nx ,  ny  ,  nz ))
      allocate (       d_p( nx ,  ny  ,  nz ))
      allocate (      vt_p( nx ,  ny  ,  nz ))
      allocate (       v_p( nx ,  ny  ,  nz ))
      allocate (     tmp_p( nx ,  ny  ,  nz ))
      allocate (       r_p( nx ,  ny  ,  nz ))
      allocate (       s_p( nx ,  ny  ,  nz ))

!
!  Compute the number of unknowns since some of the rows of the
!  coefficient matrix contain only a 1 on the diagonal.  These
!  rows correspond to the E-field components which lie on the
!  outermost boundary of the mesh and are thus subject to a
!  Dirichlet condition.  Note: this parameter is computed only
!  as a courtesey to the user and is not used explicitly in
!  the main QMR algorithm.
!
      nunk = (nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1) + nx*ny*nz
!
!  Initialize time counter...
!
      t1 = omp_get_wtime()
!
!  Open the output file.
!
      open(unit=11,file=trim(pwd)//'/qmr.iter',status='unknown')
!
!  Compute the norm of the RHS vector
!
      rhsnorm = sum(rhs_x*dconjg(rhs_x)) + sum(rhs_y*dconjg(rhs_y)) +   &
                sum(rhs_z*dconjg(rhs_z)) + sum(rhs_p*dconjg(rhs_p))   


      rhsnorm = dsqrt(rhsnorm)
!
!  Compute the inverse of the diagonal.  Do this once since a complex
!  double presion multiply of the inverse is FAR less expensive than
!  computing the quotient.  (6 flops versus  11 or 14? flops)
!
      inv_diag_x = 1.d0/diag_x
      inv_diag_y = 1.d0/diag_y
      inv_diag_z = 1.d0/diag_z
      inv_diag_p = 1.d0/diag_p
!
!  Politely introduce yourself.
!
      write(iout,*) ' '
      write(iout,*) '*****************************************'
      write(iout,*) '*  Entering multi-threaded QMR Routine  *'
      write(iout,*) '*****************************************'
      write(iout,*) ' '
 
      write(iout,*) '---------------------------------'
      write(iout,*) '  2-term coupled recurrence QMR  '
      write(iout,*) '    with Jacobi preconditioner   '
      write(iout,*) '================================='
      write(iout,*) ' '                                  
      write(iout,'(a,1x,e12.5)') ' Target residual norm: ',qmr_target*rhsnorm
      write(iout,'(a,1x,e12.5)') ' RHS norm:             ',rhsnorm
      write(iout,*) 'Maximum number of iterations: ',qmr_max_it
      write(iout,*) ' '
      write(iout,*) 'A list of iterate numbers and corresponding'
      write(iout,*) 'residual norms is located in the file: qmr.iter'
      write(iout,*) ' '
 
      write(iout,99) 'Number of unknowns:  ', nunk
      write(iout,*) '  '
 99   format(1x,a21,1x,i8)                           

!
!  Initialize some vectors
!
      p_x = (0.d0,0.d0)
      p_y = (0.d0,0.d0)
      p_z = (0.d0,0.d0)
      p_p = (0.d0,0.d0)

      d_x = (0.d0,0.d0)
      d_y = (0.d0,0.d0)
      d_z = (0.d0,0.d0)
      d_p = (0.d0,0.d0)
!
!  Squirrel away the RHS vector into the vector 's'
!  since 's' isn't used until the main loop.  Store in
!  the vector 'rhs' the initial guess at the solution.
!
      s_x = rhs_x
      s_y = rhs_y
      s_z = rhs_z
      s_p = rhs_p

      rhs_x = (0.d0,0.d0)
      rhs_y = (0.d0,0.d0)
      rhs_z = (0.d0,0.d0)
      rhs_p = (0.d0,0.d0)

!----------------------------! ATTENTION !-----------------------------!
!                                                                      !
!  From this point forward in the algorithm, the vector 'rhs' will     !
!  contain the current solution.                                       !
!                                                                      !
!----------------------------------------------------------------------!       
 
!
!  Compute the initial residual vector.
!
      call inprod_x(rhs_x,rhs_y,rhs_z,rhs_p,tmp_x)
      call inprod_y(rhs_x,rhs_y,rhs_z,rhs_p,tmp_y)
      call inprod_z(rhs_x,rhs_y,rhs_z,rhs_p,tmp_z)
      call inprod_p(rhs_x,rhs_y,rhs_z,rhs_p,tmp_p)

      vt_x = s_x - tmp_x
      vt_y = s_y - tmp_y
      vt_z = s_z - tmp_z
      vt_p = s_p - tmp_p

      r_x = s_x - tmp_x
      r_y = s_y - tmp_y
      r_z = s_z - tmp_z
      r_p = s_p - tmp_p

      rho_np = sum(r_x*dconjg(r_x)) + sum(r_y*dconjg(r_y)) +            &
               sum(r_z*dconjg(r_z)) + sum(r_p*dconjg(r_p))

      rho_np = dsqrt(rho_np)
 
      write(iout,'(a,1x,e11.5)') ' Initial residual: ',rho_np
      write(iout,*) ' '

      v_x = vt_x/rho_np
      v_y = vt_y/rho_np
      v_z = vt_z/rho_np
      v_p = vt_p/rho_np
!
!  Now, set the vector 's' to its proper starting value...
!
      s_x = (0.d0,0.d0)
      s_y = (0.d0,0.d0)
      s_z = (0.d0,0.d0)
      s_p = (0.d0,0.d0)
!
!  ... and initialize these chaps as well.
!
      c_n      =  1.d0
      eps      = ( 1.d0, 0.d0)
      theta_n  =  0.d0
      eta      = (-1.d0, 0.d0)

!----------------------------------------------------------------------!
!                       Start main iterative loop                      !
!----------------------------------------------------------------------!             
 
      write(iout,*) 'Now computing the QMR iterates...'
      write(iout,*) ' '

      converge_flag = .false.
      iter = 1

!$omp parallel private (my_pid,psum) num_threads(4)
      my_pid = omp_get_thread_num()

100   continue

!$omp barrier
      if (my_pid.eq.0) then 
        c_nm     = c_n
        rho_n    = rho_np
        theta_nm = theta_n
        del = (0.d0,0.d0)
      end if
!$omp barrier

!
! step 1
!
      if (my_pid.eq.0) then
        tmp_x = v_x*inv_diag_x
        psum = sum(v_x*tmp_x) 
      end if

      if (my_pid.eq.1) then
        tmp_y = v_y*inv_diag_y
        psum = sum(v_y*tmp_y) 
      end if

      if (my_pid.eq.2) then
        tmp_z = v_z*inv_diag_z
        psum = sum(v_z*tmp_z) 
      end if

      if (my_pid.eq.3) then
        tmp_p = v_p*inv_diag_p
        psum = sum(v_p*tmp_p) 
      end if

!omp critical
      del = del + psum
!omp end critical

!$omp barrier
      if (my_pid.eq.0) k1 = rho_n * del / eps
!$omp barrier
 
!
! step 2
!
      if (my_pid.eq.0) p_x = tmp_x - p_x*k1
      if (my_pid.eq.1) p_y = tmp_y - p_y*k1
      if (my_pid.eq.2) p_z = tmp_z - p_z*k1
      if (my_pid.eq.3) p_p = tmp_p - p_p*k1
!$omp barrier


!
! step 3
!
      if (my_pid.eq.0) then
        call inprod_x(p_x,p_y,p_z,p_p,tmp_x)
        psum_complex(my_pid) =  sum(p_x*tmp_x) 
      end if

      if (my_pid.eq.1) then
        call inprod_y(p_x,p_y,p_z,p_p,tmp_y)
        psum_complex(my_pid) =  sum(p_y*tmp_y) 
      end if

      if (my_pid.eq.2) then
        call inprod_z(p_x,p_y,p_z,p_p,tmp_z)
        psum_complex(my_pid) =  sum(p_z*tmp_z) 
      end if

      if (my_pid.eq.3) then
        call inprod_p(p_x,p_y,p_z,p_p,tmp_p)
        psum_complex(my_pid) =  sum(p_p*tmp_p) 
      end if
!$omp barrier

      eps  =  sum(psum_complex)
      beta = eps / del
!$omp barrier
 
      if (my_pid.eq.0) then
        vt_x = tmp_x - v_x*beta
        psum_real(my_pid) =  sum(vt_x*dconjg(vt_x)) 
      end if 

      if (my_pid.eq.1) then
        vt_y = tmp_y - v_y*beta
        psum_real(my_pid) =  sum(vt_y*dconjg(vt_y)) 
      end if 

      if (my_pid.eq.2) then
        vt_z = tmp_z - v_z*beta
        psum_real(my_pid) =  sum(vt_z*dconjg(vt_z)) 
      end if 

      if (my_pid.eq.3) then
        vt_p = tmp_p - v_p*beta
        psum_real(my_pid) =  sum(vt_p*dconjg(vt_p)) 
      end if 
!$omp barrier
      rho_np = dsqrt(sum(psum_real))
!
! step 4
!
      if (my_pid.eq.0) then
        theta_n = rho_np*rho_np/c_nm/(dble(beta)**2+dimag(beta)**2)
        c_n = 1.d0/(1.d0 + theta_n)
        eta = -eta*rho_n*c_n/(beta*c_nm)
        k2 = theta_nm*c_n
      end if
!$omp barrier

      if (my_pid.eq.0) then
        d_x   =   p_x*eta + d_x*k2
        rhs_x = rhs_x     + d_x
        s_x   = tmp_x*eta + s_x*k2
        r_x   =   r_x     - s_x

        v_x   = vt_x/rho_np
        psum_real(my_pid) = sum(r_x*conjg(r_x)) 
      end if 

      if (my_pid.eq.1) then
        d_y   =   p_y*eta + d_y*k2
        rhs_y = rhs_y     + d_y
        s_y   = tmp_y*eta + s_y*k2
        r_y   =   r_y     - s_y

        v_y   = vt_y/rho_np
        psum_real(my_pid) = sum(r_y*conjg(r_y)) 
      end if 

      if (my_pid.eq.2) then
        d_z   =   p_z*eta + d_z*k2
        rhs_z = rhs_z     + d_z
        s_z   = tmp_z*eta + s_z*k2
        r_z   =   r_z     - s_z

        v_z   = vt_z/rho_np
        psum_real(my_pid) = sum(r_z*conjg(r_z)) 
      end if 

      if (my_pid.eq.3) then
        d_p   =   p_p*eta + d_p*k2
        rhs_p = rhs_p     + d_p
        s_p   = tmp_p*eta + s_p*k2
        r_p   =   r_p     - s_p

        v_p   = vt_p/rho_np
        psum_real(my_pid) = sum(r_p*conjg(r_p)) 
      end if 

!$omp barrier
      rnorm = sqrt(sum(psum_real))
!$omp barrier

!
! step 5
!

      if (my_pid.eq.0) then

        if (mod(iter,10).eq.0) then
          t2 = omp_get_wtime()
          write(iout,98) iter, qmr_max_it,rnorm,t2-t1
        endif

        write(11,*) iter,rnorm,t2-t1

!
! Check to see if the solution has converged...
!
        if (rnorm.lt.qmr_target*rhsnorm) converge_flag = .true.
        if ((.not.converge_flag).and.(iter.le.qmr_max_it)) iter = iter + 1

      end if 

!$omp barrier

      if ((.not.converge_flag).and.(iter.le.qmr_max_it)) goto 100

!$omp end parallel        
          
        t2 = omp_get_wtime()
        if (converge_flag) then 
          write(iout,*) 'Solution converged in ',iter,' iterations'
        else
          write(iout,*) '********** WARNING ************'
          write(iout,*) '  Maximum iterations exceeded  '
          write(iout,*) '*******************************'

        end if

        write(iout,'(a,1x,e11.4)') '   Residual:',rnorm
        write(iout,*) ' '
        write(iout,*) 'Time spent in QMR [s]: ',t2 - t1
        write(iout,*) 'Seconds per iteration:    ',(t2-t1)/real(iter)
        write(iout,*) ' '
        write(iout,99) 'Number of unknowns:  ', nunk
        write(iout,*) '  '
        close(11)
!
!  Deallocate storage for the local arrays used in the QMR routine.
! 
        deallocate (inv_diag_x,p_x,d_x,vt_x,v_x)
        deallocate (tmp_x,r_x,s_x)
        deallocate (inv_diag_y,p_y,d_y,vt_y,v_y)
        deallocate (tmp_y,r_y,s_y)
        deallocate (inv_diag_z,p_z,d_z,vt_z,v_z)
        deallocate (tmp_z,r_z,s_z)
        deallocate (inv_diag_p,p_p,d_p,vt_p,v_p)
        deallocate (tmp_p,r_p,s_p)
!
!  Exit gracefully
!
        write(iout,*) ' '
        write(iout,*) '*************************'
        write(iout,*) '*  Exiting QMR Routine  *'
        write(iout,*) '*************************'
        write(iout,*) ' '
 
 
 98   format('it: ',i6,'/',i6,'    r: ',e13.6,'    t:',e13.6,1x,f8.4)
 101  format(10(f8.5,1x))

      end
