      program main
!
!  Program to generate aphid meshes for the T=1,10s MT examples
!  from Computers & Geosciences article.
!
!  compile and exectute, piping output to aphid.mesh
!
!  Written by C. J Weiss, Dept of Geosciences, Virginia Tech, 2012
!
      integer, parameter :: NAIR_CELLS = 51
      integer, parameter :: NEAR_CELLS = 69

      real(4), parameter :: AIR_MIN = -100.e3
      real(4), parameter :: AIR_MAX = -1.e1
      real(4), parameter :: EAR_MIN =  1.e1
      real(4), parameter :: EAR_MAX =  50.e3
      real(4), parameter :: XMIN = -100.e3
      real(4), parameter :: XMAX =  100.e3
      real(4), parameter :: YMIN = -100.e3
      real(4), parameter :: YMAX =  100.e3

      integer, parameter :: NV = NAIR_CELLS + NEAR_CELLS + 1
      integer, parameter :: NH = 50

      real(4)    :: x(NH),y(NH),z(NV)
      real(4) :: sig((NH-1)*(NH-1)*(NV-1))
      real(4) :: eps((NH-1)*(NH-1)*(NV-1))

!=======================!
!                       !
      x(1)  = XMIN
      x(NH) = XMAX

      y(1)  = YMIN
      y(NH) = YMAX

      z(1)  = log10(-AIR_MIN)
      z(NV) = log10( EAR_MAX)
!                       !
!=======================!

      dx=abs(x(NH)-x(1))/real(NH-1)
      dy=abs(y(NH)-y(1))/real(NH-1)
!
      do i=2,NH
        x(i) = x(i-1) + dx
        y(i) = y(i-1) + dy
      enddo

      dz_air = (log10(-AIR_MIN) - log10(-AIR_MAX)) / real(NAIR_CELLS-1)
      dz_ear = (log10( EAR_MAX) - log10( EAR_MIN)) / real(NEAR_CELLS-1)

      z(NAIR_CELLS+1) = 1.e30
      do i=2,NAIR_CELLS
        z(i) = z(i-1) - dz_air
      end do
      do i=NV-1,NV-NEAR_CELLS+1,-1
        z(i) = z(i+1) - dz_ear
      end do
      z = 10**z
      z(NAIR_CELLS+1) = 0. 
      z(1:NAIR_CELLS+1) = -z(1:NAIR_CELLS+1) 

      dz_air = -AIR_MIN / real(NAIR_CELLS)
      z(1) = AIR_MIN
      do i=2,NAIR_CELLS+1
        z(i) = z(i-1) + dz_air
      end do

      dz_norm = 1.0d3
      dz_small = 0.05d3
      do i=NAIR_CELLS+2,NAIR_CELLS+NEAR_CELLS+1
        if ((i.ge.NAIR_CELLS+2+19).and.(i.lt.NAIR_CELLS+41)) then 
          z(i) = z(i-1) + dz_small
        else
          z(i) = z(i-1) + dz_norm
        end if
 
      end do
!
!  Populate the cells of the grid with some conductivity value.
!
      ind = 0 
      do k=1,NV-1
        z1 = z(k)
        z2 = z(k+1)
        zc = (z1+z2)/2.
 
        do j=1,NH-1
          y1 = y(j)
          y2 = y(j+1)
          yc = (y1+y2)/2.

          do i=1,NH-1
            x1 = x(i)
            x2 = x(i+1)
            xc = (x1+x2)/2.

            ind = ind + 1 

            sig(ind) = 0.01
            if (zc.lt.20.e3)  sig(ind) =   0.10 
            if (zc.lt.19.e3)  sig(ind) =   0.001
            if (zc.lt.  0.e3) sig(ind) =   1.e-8
            
          enddo
        enddo
      enddo 
      eps = 1.
!
!  Write it out!
! 
      write(6,*) NH,NH,NV
      write(6,*) ' '
      write(6,1) x
      write(6,*) ' '
      write(6,1) y
      write(6,*) ' '
      write(6,1) z
      write(6,*) ' REAL CONDUCTIVITY'
      write(6,1) sig
      write(6,*) ' RELATIVE ELECTRIC PERMITTIVITY'
      write(6,1) eps

 1    format(6(e12.4,1x))
      end program main
      
