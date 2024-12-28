      subroutine getline( unitno, line )
c  Copyright (C) 2018, Jerry Sellwood
      implicit none
c This routine is calls readdat on the master node only.  Readdat
c   is separated from getline in order to minimize duplication of
c   code in the single processor and mpi versions
c This is the single processor version
c
c calling arguments
      integer unitno
      character*(*) line
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
c external
      integer lnblnk
c
c local variables
      integer i, ifail
c
c check for nonsense arguments
      if( ( unitno .le. 0 ) .or. ( len( line ) .le. 0 ) )then
        print *, 'unitno =', unitno, '   maxline =', len( line )
        call crash( 'GETLINE', 'nonsense calling arguments' )
      end if
c
      if( numprocs .gt. 1 )call crash( 'GETLINE', 'mpi version needed' )
c must be master node
      call readdat( unitno, line, ifail )
c check error flag
      if( ifail .ne. 0 )then
        i = lnblnk( line )
        print *, line( 1:i )
        if( ifail .eq. 1 )then
          call crash( 'GETLINE', 'line not long enough to detab' )
        else if( ifail .eq. 2 )then
          call crash( 'GETLINE', 'End of file encountered' )
        else
          call crash( 'GETLINE', 'problem reading line' )
        end if
      end if
      return
      end
