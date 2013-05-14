      program main

      character(80) :: cstr
!
!  checking to see if the job is 
!  multi-threaded (OMP_NUM_THREADS=4) or single-threaded
!  (OMP_NUM_THREADS=1), and halting execution otherwise.
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
! launch the serial job.  note that the user can specify
! a working directory 'pwd' on the command line by including the flag
! '-wd pwd' at job initiation.
!
      call aphid_main

      end program main
