!======================================================================!  
!==============================================================! GENERAL
!======================================================================!  

      module general

      integer,     parameter :: DP = 8
      integer,     parameter :: SP = 4
      real(DP),    parameter :: PI = 3.141592653589793d0
      real(DP),    parameter :: MU0 = PI * 4.d-7
      complex(DP), parameter :: II = (0.d0,1.d0)
      real(DP),    parameter :: EPS0 = 8.8541878176d-12
      real(DP),    parameter :: SP_TOL = 1.d-37

      character(1)  :: sigma_dump_flag,csig_flag
      character(12) :: output_field
      character(30) :: pwd = "./"
      integer       :: my_rank,num_procs,ierr
      integer       :: iout = 6

      integer       :: oper_mode
      real(dp)      :: freq,omega
      real(dp)      :: zeta
      real(dp)      :: xs,ys,zs,thetas,phis
       
      real(dp)      :: qmr_target
      integer       :: qmr_max_it

      integer       :: MAX_RECUR

      end module general

!======================================================================!  
!================================================================! MESH
!======================================================================!  

      module mesh
      use general

      integer      :: nx,ny,nz

      real(dp),    dimension(:)      :: xn,yn,zn
      real(dp),    dimension(:)      :: xh,yh,zh
      real(dp),    dimension(:)      :: dx,dy,dz
      real(dp),    dimension(:)      :: dxc,dyc,dzc
      complex(dp), dimension(:,:,:)  :: sig
      complex(dp), dimension(:,:,:)  :: k2_xgrid,k2_ygrid
      complex(dp), dimension(:,:,:)  :: k2_zgrid,k2_pgrid
      complex(dp), dimension(:,:,:)  :: rhs_x,rhs_y,rhs_z,rhs_p
      complex(dp), dimension(:,:,:)  :: diag_x,diag_y,diag_z,diag_p
      complex(dp), dimension(:,:,:)  :: cpl_x,cpl_y,cpl_z,cpl_p

      complex(dp), dimension(:,:,:)  :: Bx,By,Bz
      complex(dp), dimension(:,:,:)  :: Ex,Ey,Ez

      allocatable  :: xn,yn,zn,xh,yh,zh,dx,dy,dz
      allocatable  :: dxc,dyc,dzc
      allocatable  :: sig
      allocatable  :: k2_xgrid,k2_ygrid
      allocatable  :: k2_zgrid,k2_pgrid
      allocatable  :: rhs_x,rhs_y,rhs_z,rhs_p
      allocatable  :: diag_x,diag_y,diag_z,diag_p
      allocatable  :: Bx,By,Bz,Ex,Ey,Ez
      allocatable  :: cpl_x,cpl_y,cpl_z,cpl_p

! Variables !----------------------------------------------------------!
!                                                                      !
!  nx,ny,nz:       Number of nodes in the x, y and z directions.       !
!  xn,yn,zn:       Node coordinates in the x, y and z directions.      ![m]
!  dx,dy,dz:       Node spacing in the x, y and z direction.           ![m]
!  sig:            Cell conductivities.                                ![S/m]
!                                                                      !
!----------------------------------------------------------------------!

      end module mesh


!======================================================================!  
!======================================================! INPROD_X_COEFS
!======================================================================!  

      module inprod_x_coefs
      use general

      real(dp), dimension(:)  :: coef_2vec,coef_3vec
      real(dp), dimension(:)  :: coef_4vec,coef_5vec
      real(dp), dimension(:)  :: coef_6vec,coef_7vec

      integer, dimension(:)   :: ia_vec,ib_vec
      integer, dimension(:)   :: ja_vec,jb_vec
      integer, dimension(:)   :: ka_vec,kb_vec

      allocatable  :: coef_2vec,coef_3vec
      allocatable  :: coef_4vec,coef_5vec
      allocatable  :: coef_6vec,coef_7vec

      allocatable  :: ia_vec,ib_vec
      allocatable  :: ja_vec,jb_vec
      allocatable  :: ka_vec,kb_vec

      end module inprod_x_coefs


!======================================================================!  
!======================================================! INPROD_Y_COEFS
!======================================================================!  

      module inprod_y_coefs
      use general

      real(dp), dimension(:)  :: coef_2vec,coef_3vec
      real(dp), dimension(:)  :: coef_4vec,coef_5vec
      real(dp), dimension(:)  :: coef_6vec,coef_7vec

      integer, dimension(:)   :: ia_vec,ib_vec
      integer, dimension(:)   :: ja_vec,jb_vec
      integer, dimension(:)   :: ka_vec,kb_vec

      allocatable  :: coef_2vec,coef_3vec
      allocatable  :: coef_4vec,coef_5vec
      allocatable  :: coef_6vec,coef_7vec

      allocatable  :: ia_vec,ib_vec
      allocatable  :: ja_vec,jb_vec
      allocatable  :: ka_vec,kb_vec

      end module inprod_y_coefs


!======================================================================!  
!======================================================! INPROD_Z_COEFS
!======================================================================!  

      module inprod_z_coefs
      use general

      real(dp), dimension(:)  :: coef_2vec,coef_3vec
      real(dp), dimension(:)  :: coef_4vec,coef_5vec
      real(dp), dimension(:)  :: coef_6vec,coef_7vec

      integer, dimension(:)   :: ia_vec,ib_vec
      integer, dimension(:)   :: ja_vec,jb_vec
      integer, dimension(:)   :: ka_vec,kb_vec

      allocatable  :: coef_2vec,coef_3vec
      allocatable  :: coef_4vec,coef_5vec
      allocatable  :: coef_6vec,coef_7vec

      allocatable  :: ia_vec,ib_vec
      allocatable  :: ja_vec,jb_vec
      allocatable  :: ka_vec,kb_vec

      end module inprod_z_coefs


!======================================================================!  
!======================================================! INPROD_P_COEFS
!======================================================================!  

      module inprod_p_coefs
      use general

      real(dp), dimension(:)  :: coef_2vec,coef_3vec
      real(dp), dimension(:)  :: coef_4vec,coef_5vec
      real(dp), dimension(:)  :: coef_6vec,coef_7vec

      integer, dimension(:)   :: ia_vec,ib_vec,ic_vec
      integer, dimension(:)   :: ja_vec,jb_vec,jc_vec
      integer, dimension(:)   :: ka_vec,kb_vec,kc_vec

      allocatable  :: coef_2vec,coef_3vec
      allocatable  :: coef_4vec,coef_5vec
      allocatable  :: coef_6vec,coef_7vec

      allocatable  :: ia_vec,ib_vec,ic_vec
      allocatable  :: ja_vec,jb_vec,jc_vec
      allocatable  :: ka_vec,kb_vec,kc_vec

      end module inprod_p_coefs
