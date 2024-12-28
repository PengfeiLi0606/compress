      subroutine getset
c  Copyright (C) 2015, Jerry Sellwood
      implicit none
c general start up routine that sets run number and field method type:
c   sets flags for single or parallel processing
c   reads in the run number
c   verifies that the appropriate .dat file is present and opens it
c   verifies that .dat file contains the input run number
c   reads method type from the .dat file and calls codeset
c   opens .lis file with access = 'append' if possible
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
      integer iunit( 6 + 2 * mgrids )
      equivalence ( ni, iunit( 1 ) )
c
c local variables
      character line*80
      character*3 code
      integer i, j, istate
c
c set flags for single processor operation
      master = .true.
      myid = 0
      numprocs = 1
      parallel = .false.
c set all logical unit nos to nonsense values
      do i = 1, 6 + 2 * mgrids
        iunit( i ) = -1
      end do
c read run number from stdin
c      if( master )call gtintg( 'Enter run no', irun )
      irun = 100
c open input file - single processor version
c      call opnfil( ni, 'dat', 'formatted', 'old', 'seq', i )
c the above line is commented and replaced by the following 4 lines --PLi
      if( ni .le. 0 )call new_unit( ni ) 
      open( ni, file = 'input.dat',
     +          status = 'old', form = 'formatted', iostat = istate )
      if( istate .ne. 0 )call crash( 'GETSET', '.dat file not found' )
c check run number
      i =100
c cc Following lines are commented --PLi
c      call getline( ni, line )
c      read( line( 11:80 ), * )i
c      if( i .ne. irun )then
c        if( master )then
c          print *, 'Run number input was', irun
c          print *, 'Number read from data file was', i
c        end if
c        call crash( 'GETSET', 'Incorrect data file for run' )
c      end if
c open output file - single processor version
c      call opnfil( no, 'lis', 'formatted', 'old', 'append', i )
c the above line is commented and replaced by the following 3 lines--PLi
      if( no .le. 0 )call new_unit( no )
      open(no, file='run.log',
     +         status='old', form='formatted', iostat=istate)
      if( istate .ne. 0 )then
        print *, 'Old run.log file not found - opening a new file'
c        call opnfil( no, 'lis', 'formatted', 'new', 'seq', i )
        open(no, file='run.log',
     +           status='new', form='formatted', iostat=istate)
        if( istate.ne.0)call crash( 'GETSET', 'Error opening .lis file')
      end if
c read and interpret field determination method
c cc Following lines are commented --PLi
c      call getline( ni, line )
c      i = 11
c      do while ( line( i:i ) .eq. ' ' )
c        i = i + 1
c      end do
c      read( line( i:i+2 ), '( a )' )code
c ccccccc set code=non, PLi 2024/10/05 ccccccc
      code = 'NON' 
c ccccccccccccccccccccccc
      if( code .eq. 'MLT' )code = 'TWO'
      hybrid = .false.
      twogrd = .false.
      if( ( code .eq. 'HYB' ) .or. ( code .eq. 'TWO' ) )then
        if( code .eq. 'HYB' )hybrid = .true.
        if( code .eq. 'TWO' )twogrd = .true.
        read( line( 21:80 ), * )ngrid
        call getline( ni, line )
        i = 11
        do while ( line( i:i ) .eq. ' ' )
          i = i + 1
        end do
        read( line( i:i+2 ), '( a )' )code
      else
        ngrid = 1
      end if
      do j = 1, ngrid
        call lowercase( code )
        if( master )print
     +        '( '' Field determination method selected: '', a3 )', code
        i = -1
        if( code .eq. 'non' )i = 0
        if( code .eq. 'p2d' )i = 1
        if( code .eq. 'c3d' )i = 2
        if( code .eq. 'sf2' )i = 3
        if( code .eq. 'p3a' )i = 4
        if( code .eq. 'p3d' )i = 5
        if( code .eq. 's3d' )i = 6
        if( code .eq. 'sf3' )i = 7
        if( code .eq. 'c2d' )i = 8
        if( code .eq. 'scf' )i = 9
        if( code .eq. 'dr3' )i = 10
        if( code .eq. 'bht' )i = 11
        if( i .lt. 0 )call crash( 'GETSET', 'Unrecognized method' )
        igrid( j ) = i
        if( j .lt. ngrid )then
          call getline( ni, line )
          i = 11
          do while ( line( i:i ) .eq. ' ' )
            i = i + 1
          end do
          read( line( i:i+2 ), '( a )' )code
        end if
      end do
      call codeset
c report selected option
      if( master )then
        if( ngrid .eq. 1 )then
          write( no,
     +        '( '' Field determination method selected: '', a3 )' )code
        else if( hybrid )then
          write( no, '( '' Hybrid methods selected: '', 3i2 )' )
     +                                      ( igrid( i ), i = 1, ngrid )
        else
          write( no, '( '' Multiple grids selected: '', 3i2 )' )
     +                                      ( igrid( i ), i = 1, ngrid )
        end if
      end if
      return
      end
