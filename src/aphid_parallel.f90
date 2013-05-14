      program main
      use general
      use omp_lib
      include 'mpif.h'

      integer      :: num_freq,num_source
      character(2) :: source_type
      character(7) :: tmp_dir
      character(80):: cstr

      real(8), dimension(:) :: freq_array
      integer, dimension(:) :: mode_array
      allocatable           :: freq_array,mode_array

!
!  test that multi-threading option is set to a valid number of threads,
!  either 1 or 4.
!
      call get_environment_variable("OMP_NUM_THREADS",cstr)
      write(6,*) cstr
      if ( ( (cstr(1:1).ne.'1').and.(cstr(1:1).ne.'4') )                &
                                         .or. (cstr(2:2).ne.' ') ) then
        write(6,*) 'Error: OMP_NUM_THREADS = ',cstr(1:5)
        write(6,*) 'Value must be 1 or 4.  Please reset and try again.'
        stop
      end if
!
!  initialize MPI
!
      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,num_procs,ierr)

      iout = 1001 + my_rank

      allocate ( freq_array(num_procs) , mode_array(num_procs) )
!
!  Input file management is handled by the process 0.
!
      if (my_rank.eq.0) then
!
!  open the frequency array file and load its contents into the 
!  first 'num_freq' indices of 'freq_array'
!
        open(unit=11,file='aphid_mpi.freq',status='unknown')
        read(11,*) 
        read(11,*) 
        read(11,*) 
        read(11,*) num_freq
        do ifreq=1,num_freq
          read(11,*) freq_array(ifreq)
        end do
        close(11)
!
!  open the input file to determin whether this is an MT problem
!  or a set of discrete, dipole arrays.
!
        open(unit=12,file='aphid_mpi.input',status='unknown')
        read(12,*) 
        read(12,*) 
        read(12,*) 
        read(12,*) source_type
        close(12)
!
!  If this is an MT problem, then there are two forward calcs required
!  for each frequency if impedances are to be computed at the end.  Load
!  the set of frequencies block-wise into the 'freq_array'.  The plane
!  wave polarization is specified in 'mode_array'.  Set the first 
!  'num_freq' values for x-polarization and the second block of values 
!  for y-polarization.  Therefore, the x/y polarization pair for 
!  frequency i(>0) is handled by processes i-1 and i-1+num_freq, or 
!  equivalently, stored in directories 'tmp-i' and 'tmp-j' where 
!  j=i+num_freq.
!
        if (source_type.eq.'MT') then
          num_source = 2

          ia =    num_freq+1
          ib =  2*num_freq
          freq_array(ia:ib) = freq_array(1:num_freq)

          mode_array(1 :num_freq) = 2  ! x-polarization
          mode_array(ia:ib) = 3        ! y-polarization
!
!  If this is a dipole array problem, load the frequencies block-wise, 
!  as done before for the MT problem.  In this case, however, the 
!  'operation_mode' is '5' for all calculations.
!
        else
          open(unit=13,file='aphid_mpi.source',status='unknown')
          read(13,*) 
          read(13,*) 
          read(13,*) 
          read(13,*) num_source
          close(13)

          do i=2,num_source
            ia = (i-1)*num_freq+1
            ib =     i*num_freq
            freq_array(ia:ib) = freq_array(1:num_freq)
          end do 
          mode_array = 5

        end if 
        
      end if
!
!  broadcast across processes the arrays containing frequency and 
!  operation mode info.  For a given process, only one value from 
!  each array is used, but they're small and this is easier than
!  parsing 'em out one-by-one from process 0 to the right recipient.
!
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call mpi_bcast(freq_array,num_procs,                              &
                             MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(mode_array,num_procs,                              &
                                      MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(num_freq,1,      MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(num_source,1,    MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(source_type,2, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!
!  create an array of working directories for parallel execution.
!
      if (num_procs.le.9) then
        num_digits = 1 
      else if (num_procs.le.99) then
        num_digits = 2 
      else if (num_procs.le.999) then
        num_digits = 3 
      else
        write(iout,*) 'error: number of MPI processes > 1000.'
        stop
      end  if

      itmp = my_rank+1
      tmp_dir = "tmp-   "
      do n=1,num_digits
        ifac = 10**(num_digits - n)
        tmp_dir(n+4:n+4) = char(48+itmp/ifac)
        itmp = itmp - (itmp/ifac)*ifac
      end do
      call system("mkdir "//tmp_dir) 
!
!  load each working directory with the appropriate input file.
!
      pwd = tmp_dir
      open(unit=11,file='aphid_mpi.input',status='unknown')
      open(unit=12,file=trim(pwd)//"//aphid.input",status='unknown')
      read(11,'(a)') cstr; write(12,'(a)',ADVANCE='YES') cstr
      read(11,'(a)') 
      cstr = 'APhiD input file generated from ../aphid_mpi.input'
      write(12,'(a)',ADVANCE='YES') cstr
      read(11,'(a)') cstr; write(12,'(a)',ADVANCE='YES') cstr
      read(11,*)
      write(12,'(i4,2x,a)') mode_array(my_rank+1),' !! operational mode'
      read(11,'(a)') cstr; write(12,'(a)',ADVANCE='YES') cstr
      read(11,'(a)') cstr; write(12,'(a)',ADVANCE='YES') cstr
      write(12,'(e14.7,2x,a)')freq_array(my_rank+1),' !! frequency [Hz]'
      do while (.true.)
        read(11,'(a)',end=1) cstr; write(12,'(a)',ADVANCE='YES') cstr
      end do
 1    close(11)
      close(12)
!
!  if required, load each working directory with the appropriate dipole 
!  array info.
!
      if (source_type.eq.'CS') then
        open(unit=11,file='aphid_mpi.source',status='unknown')
        open(unit=12,file=trim(pwd)//"//aphid.source",status='unknown')

        read(11,*)
        read(11,*)
        read(11,*)
        read(11,*)
        read(11,*)

      
        do i=1,my_rank/num_freq + 1
          read(11,*)
          read(11,*) 
          read(11,*) 
          read(11,*) ndip
          do j=1,ndip
            read(11,'(a)') cstr 
            if (my_rank/num_freq+1.eq.i)                                &
                                      write(12,'(a)',ADVANCE='YES') cstr
          end do
          read(11,*) 
        end do
        close(11)
        close(12)
      end if 
!
!  copy mesh file into each working directory
!
      call system("cp aphid.mesh "//tmp_dir)
!
!  launch the multi-threaded jobs in each working directory
!
      call aphid_main
      write(6,*) 'Process complete: ',my_rank 
!
!  shut down parallel processes and exit program
!
      call mpi_finalize(ierr)

      end program main
