      real*8 function gmdtab( r )
      use data_array
      implicit none
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      integer lnblnk
      real*8 splint2
c
c local variables
      character*80 line
      integer i, ifail, iuse
      real a
      include 'inc/pi.f'
      save iuse
      data iuse / 0 /
c
      if( iuse .eq. 0 )then
c read in data
        call getdat
c        call gtreal( 'Enter desired disk M/L', mtold )
        print *, 'Disk M/L is', mtold
c convert to surface mass density
        do i = 1, nvd
          sdsd( i ) = mtold * sdsd( i ) * 1.e-4
        end do
c convert to mass enclosed assuming a linear variation between table values
        msd2( 1 ) = 0
        do i = 2, nvd
          msd2( i ) = msd2( i - 1 ) +
     +     pi * ( rdtab( i ) - rdtab( i - 1 ) ) *
     +  ( rdtab( i - 1 ) * ( 5. * sdsd( i )     - 2. * sdsd( i - 1 ) ) +
     +    rdtab( i )     * ( 5. * sdsd( i - 1 ) - 2. * sdsd( i ) ) ) / 3.
        end do
        rdmax = rdtab( nvd )
        mdtot = msd2( nvd )
        print '( ''stellar disk mass to '', f6.3, '' kpc ='', f8.4, '//
     +                      ' '' x 10^10 solar masses'' )', rdmax, mdtot
c create array in internal units and rescale masses
        allocate ( rdtabs( nvd ) )
        do i = 1, nvd
          rdtabs( i ) = rdtab( i ) * rscale( 1 )
          msd2( i ) = msd2( i ) * cmpmas( 1 )
        end do
c set spline coeffs
        gmdtab = splint2( rdtabs, msd2, nvd, rdtabs( 2 ),
     +                    sdlam, sdc, .true. )
c adjust veld for adopted M/L
        do i = 1, nvd
          veld( i ) = sqrt( mtold ) * veld( i )
        end do
c initialize disk velocity spline
        a = splint2( rdtab( 2 ), veld( 2 ), nvd - 1, rdtab( 2 ),
     +               vdlam, vdc, .true. )
        iused = 1
        if( lgas )then
c allow for He and convert to internal units
          do i = 1, nvg
            sdsg( i ) = sdsg( i ) * 1.33e-4
          end do
c convert to mass enclosed assuming a linear variation between table values
          msg2( 1 ) = 0
          do i = 2, nvg
            msg2( i ) = msg2( i - 1 ) +
     +        pi * ( rdtab( i ) - rdtab( i - 1 ) ) *
     +  ( rdtab( i - 1 ) * ( 5. * sdsg( i )     - 2. * sdsg( i - 1 ) ) +
     +   rdtab( i )     * ( 5. * sdsg( i - 1 ) - 2. * sdsg( i ) ) ) / 3.
          end do
          print '( ''gas disk mass to '', f6.3, '' kpc ='', f8.4, '//
     +         ' '' x 10^10 solar masses'' )', rdtab( nvg ), msg2( nvg )
c convert to internal units
          do i = 1, nvg
            msg2( i ) = msg2( i ) * cmpmas( 1 )
          end do
c set spline coeffs
          gmdtab = gmdtab + splint2( rdtabs, msg2, nvg, rdtabs( 2 ),
     +                               sglam, sgc, .true. )
c adjust velg for Helium
          do i = 1, nvg
            velg( i ) = sqrt( 1.33 ) * velg( i ) ! M/L for gas
          end do
c initialize gas velocity spline
          a = splint2( rdtab( 2 ), velg( 2 ), nvg - 1, rdtab( 2 ),
     +                 vglam, vgc, .true. )
        else
          print *, 'No gas data found'
        end if
c bulge density and velocity data
        if( lbulge )then
c          call gtreal( 'Enter desired bulge M/L', mtolb )
          print *, 'Bulge M/L is', mtolb
c allow M/L and convert to internal units
          do i = 1, nvb
            sdsb( i ) = mtolb * sdsb( i ) * 1.e-4  ! M/L for bulge
          end do
c convert to mass enclosed assuming a linear variation between table values
          msb2( 1 ) = 0
          do i = 2, nvb
            msb2( i ) = msb2( i - 1 ) +
     +       pi * ( rdtab( i ) - rdtab( i - 1 ) ) *
     +  ( rdtab( i - 1 ) * ( 5. * sdsb( i )     - 2. * sdsb( i - 1 ) ) +
     +   rdtab( i )     * ( 5. * sdsb( i - 1 ) - 2. * sdsb( i ) ) ) / 3.
          end do
          rdmax = max( rdtab( nvb ), rdmax )
          print '( ''bulge mass to '', f6.3, '' kpc ='', f8.4, '//
     +         ' '' x 10^10 solar masses'' )', rdtab( nvb ), msb2( nvb )
c convert bulge mass to internal units
          do i = 1, nvb
            msb2( i ) = msb2( i ) * cmpmas( 1 )
          end do
c set spline coeffs
          gmdtab = gmdtab + splint2( rdtabs, msb2, nvb, rdtabs( 2 ),
     +                               sblam, sbc, .true. )
c adjust velb for bulge M/L
          do i = 1, nvb
            velb( i ) = sqrt( mtolb ) * velb( i )
          end do
c initialize bulge velocity spline
          a = splint2( rdtab( 2 ), velb( 2 ), nvb - 1, rdtab( 2 ),
     +                 vblam, vbc, .true. )
        else
          print *, 'No bulge data found'
        end if
        iuse = 1
      end if
c get interior disk mass at scaled radius
      if( r .le. 0. )then
        gmdtab = 0
      else if( r .gt. rdtabs( nvd ) )then
        gmdtab = msd2( nvd )
      else if( r .lt. rdtabs( 2 ) )then
        gmdtab = msd2( 2 ) * ( r / rdtabs( 2 ) )**2
      else
c interpolate from table
        gmdtab = splint2( rdtabs, msd2, nvd, r, sdlam, sdc, .false. )
      end if
      if( lgas )then
c add interior gas mass at scaled radius
        if( r .gt. rdtabs( nvg ) )then
          gmdtab = gmdtab + msg2( nvg )
        else if( r .gt. rdtabs( 2 ) )then
c interpolate from table
          gmdtab = gmdtab +
     +              splint2( rdtabs, msg2, nvg, r, sglam, sgc, .false. )
        else if( r .gt. 0. )then
          gmdtab = gmdtab + msg2( 2 ) * ( r / rdtabs( 2 ) )**2
        end if
      end if
      if( lbulge )then
c add interior bulge mass at scaled radius
        if( r .gt. rdtabs( nvb ) )then
          gmdtab = gmdtab + msb2( nvb )
        else if( r .gt. rdtabs( 2 ) )then
c interpolate from table
          gmdtab = gmdtab +
     +              splint2( rdtabs, msb2, nvb, r, sblam, sbc, .false. )
        else if( r .gt. 0. )then
          gmdtab = gmdtab + msb2( 2 ) * ( r / rdtabs( 2 ) )**2
        end if
      end if
      return
      end
c ############################################################################
      subroutine getdat
      use data_array
      implicit none
c
c scratch array
      real, allocatable :: data( :,: )
c
      character extnam*6, line*80
      integer i, ifail, k, n, lnblnk, n2, nbulge, ngas
c
      if( iused .ne. 0 )call crash( 'GETDAT', 'Repeat call' )
! JAS: extnam is declared as just 6 characters
c attempt to open a _total.dat file for the bulge or gas data
      open( 1, file = 'data/MassModeling/'//rootname( 1:nch )//
     +     '.dat', 
     +     form = 'formatted', status = 'old', iostat = ifail )
      if( ifail .ne. 0 )call crash( 'GETDAT', 'Data file not found' )
c _total.dat file was found
      print '( ''Reading file '', a )', rootname( 1:nch )//'.dat'
c count header lines
      n = 0
      line( 1: 1 ) = '#'
      do while ( line( 1:1 ) .eq. '#' )
        read( 1, '( a )' )line
        n = n + 1
      end do
      n = n - 1
      print *, n, ' header lines skipped'
c count lines of data
      k = 0
      ifail = 0
      do while( ifail .eq. 0 )
        k = k + 1
        read( 1, '( a )', iostat = ifail )line
      end do
c ############ naoise removed this next line since it missed last lint of data #######
c      k = k - 1
      print *, k, ' lines of data found'
c allocate scratch array
      allocate ( data( 7, k ) )
c restart file - skipping header
      rewind 1
      do i = 1, n
        read( 1, * )
      end do
c extract raw data
      ngas = 0
      nbulge = 0
      do n = 1, k
        read( 1, '( a )' )line
        read( line, * )( data( i, n ), i = 1, 7 )
c count number of non-zero Vbul
        if( data( 4, n ) .gt. 0.0)Nbulge = Nbulge + 1
c count number of non-zero Vgas
        if( data( 2, n ) .ne. 0.0)Ngas = Ngas + 1
      end do
      close( 1 )
      lbulge = nbulge .gt. 0
      lgas = ngas .gt. 0
c allocate space
      nvd = k
      nvb = nbulge
      nvg = ngas
      allocate ( rdtab( nvd ), veld( nvd ) )
      allocate ( vdlam( nvd + 4 ), vdc( nvd + 4 ) )
      allocate ( sdsd( nvd ), msd2( nvd ) )
      allocate ( sdlam( nvd + 4 ), sdc( nvd + 4 ) )
      if( lbulge )then
        allocate ( velb( nvb ) )
        allocate ( vblam( nvb + 4 ), vbc( nvb + 4 ) )
        allocate ( sdsb( nvb ), msb2( nvb ) )
        allocate ( sblam( nvb + 4 ), sbc( nvb + 4 ) )
      end if
      if( lgas )then
        allocate ( velg( nvg ) )
        allocate ( vglam( nvg + 4 ), vgc( nvg + 4 ) )
        allocate ( sdsg( nvg ), msg2( nvg ) )
        allocate ( sglam( nvg + 4 ), sgc( nvg + 4 ) )
      end if
c save to appropriate arrays
      do n = 1, k
!        read( line, * )r, Vgas, Vdisk, Vbul, SBgas, SBdisk, SBbul ! HighRes
        rdtab( n ) = data( 1, n )
        veld( n ) = data( 3, n )
        sdsd( n ) = data( 6, n )
        if( lgas .and. ( n .le. ngas ) )then
          velg( n ) = data( 2, n )
          sdsg( n ) = data( 5, n )
        end if
        if( lbulge .and. ( n .le. nbulge ) )then
          velb( n ) = data( 4, n )
          sdsb( n ) = data( 7, n )
        end if
      end do
      deallocate ( data )
      return
      end
