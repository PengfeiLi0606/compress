      include 'data_array.f'
c      
      program compress
c program to implement Young's (1980) algorithm for compression of a
c  spherical system with a known DF by the addition of new matter in
c  the center.  This version is for an NFW halo having an initial 
c  isotropic DF that is compressed by baryonic infall to create a disk.
c  The mass profile of the disk has to be supplied in tabular form
      use aarrays
      use data_array
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/lunits.f'
c
      include 'inc/model.f'
c
      include 'fit_quality.f'
c externals
      logical gtlogl
      real*8 cofa, distfn, E0tab, gmassh, gmdtab, gmext, Phihal, Phinx
      real*8 rhohal, Phitot
c
c local variables
      character*20 args( 5 )
      integer i, j, k, nargs, status
      logical convrg, conserve, cuspy, line_arg
      real a, x, y
      real*8 c, diffm, diffo, E, errmas, r, totmas, extmas
      include 'inc/pi.f'
c
c check whether there are command line parameters
      nargs = command_argument_count( )
      line_arg = nargs .gt. 0
      if( line_arg )then
c read and interpret command line parameters
        if( nargs .gt. 5 )print *, 'Ignoring extra command arguments'
        if( nargs .lt. 4 )
     +       call crash( 'COMPRESS', 'Too few calling arguments' )
        do i = 1, nargs
          call get_command_argument( i, args( i ), j, status )
          if( status .ne. 0 )then
            print *, 'error reading command parameter #', i, args( i )
            call crash( 'COMPRESS', 'command line error' )
          end if
        end do
        rootname = trim( args( 1 ) )
c set halo mass and scale (relative to disk)
        
        read( args( 2 ), *, err = 99 )hmass
        read( args( 3 ), *, err = 99 )hrad
        read( args( 4 ), *, err = 99 )mtold
c bulge M/L is optional
        if( nargs .eq. 5 )then
          read( args( 5 ), *, err = 99 )mtolb
        else
          mtolb = 0
        end if
c      else
c no command line arguments -> enter parameters from input stream
c        call gtreal( 'Enter halo M_* in 10^10 MSun', hmass )
c        call gtreal( 'Enter halo scale radius in kpc', hrad )
c        call gtchar( 'Enter galaxy rootname', rootname )
c        call gtreal( 'Enter desired disk M/L', mtold )
c        call gtreal( 'Enter desired bulge M/L', mtolb )
      end if
      print *, 'Parameters read for galaxy ' // trim( rootname )
      print '( a, f10.4  )', 'halo mass ', hmass
      print '( a, f10.4  )', 'halo scale', hrad
      print '( a, f10.4  )', 'disk M/L  ', mtold
      print '( a, f10.4  )', 'bulge M/L ', mtolb
      nch = lnblnk( rootname )
      iused = 0
c
      lnag = .false.
      no = 6
      master = .true.
c set initial model
      ni = 1
      if( ni .ne. 0 )call getset
      call msetup
c determine type of compress
      ncomp = 0
      nrigid = 0
      extmas = 0
      cuspy = .false.
      r2cusp = .false.
      snglph = .true.
      i = 0
      do icmp = 1, ncmp
        if( idftyp( icmp ) .gt. 1 )then
          rigidp( icmp ) = .false.
          ncomp = ncomp + 1
          if( ncomp .gt. 2
     +       )call crash( 'COMPRESS', 'More than 2 compressible halos' )
c remember properties of uncompressed halo(s)
          iccmp( ncomp ) = icmp
          j = imtyp( icmp )
          ihalo0( ncomp ) = j
          cuspy = cuspy .and. ( ( j .eq. 4 ) .or. ( j .eq. 6 ) .or.
     +                        ( j .eq. 16 ) )
          cuspy = cuspy .or. ( ( j .eq. 9 ) .or. ( j .eq. 20 ) .or.
     +                       ( j .eq. 25 ) ) .or. ( j .eq. 32 )
          if( j .eq. 9 )r2cusp = .true.
          idfn0( ncomp ) = idftyp( icmp )
          dfm0( ncomp ) = dfcns( 1, icmp )
          dfb0( ncomp ) = dfcns( 2, icmp )
          cc0( ncomp ) = dfcns( 3, icmp )
          r = rtrunc( icmp )
          Emaxo( ncomp ) = Phihal( r )
        else
          nrigid = nrigid + 1
          rigidp( icmp ) = .true.
          if( nrigid .gt. 2
     +         )call crash( 'COMPRESS', 'More than 2 rigid components' )
c remember compressing mass component(s)
          ircmp( nrigid ) = icmp
          extmas = extmas + cmpmas( icmp )
c set disk thickness and functional form - Gaussian with z0=rscale/10
          if( disc( icmp ) )then
            i = i + 1
            if( i .gt. 1 )call crash( 'COMPRESS', 'more than 1 disk' )
            z0init( icmp ) = .1 * rscale( icmp )
            iztyp( icmp ) = 2
c set up table for gmext
            r = .5 * rtrunc( icmp )
            x = gmext( r )
          end if
        end if
      end do
      r1cusp = cuspy
      r0core = .not. cuspy
c allocate arrays
      allocate ( arad( mradad ) )
      allocate ( phin( mradad ) )
      allocate ( gmn( mradad, ncomp ) )
      allocate ( rhon( mradad, ncomp ) )
      allocate ( plamda( mradad + 4 ) )
      allocate ( phinc( mradad + 4 ) )
      allocate ( glamda( mradad + 4, ncomp ) )
      allocate ( gmnc( mradad + 4, ncomp ) )
      allocate ( rlamda( mradad + 4, ncomp ) )
      allocate ( rhonc( mradad + 4, ncomp ) )
c set revised energy bounds for this component only
      icmp = iccmp( 1 )
      x = rtrunc( icmp )
      call cutoff( x )
c build tables of E_0(J1,j2 ) and df( E, Lz)
      icomp = 1
      icmp = iccmp( icomp )
      E = E0tab( 0.d0, 0.d0 )
      E = max( E, -1.d0 )
      E = distfn( E, 0.d0 )
c set form of compressing mass - type 6 = exp disc; type = 21 for tabulated
      imtyp( 1 ) = 21
      ctype( 1 ) = 'MTAB'
      cmpmas( 1 ) = 1. / hmass
      rscale( 1 ) = 1. / hrad
c initialize the disk mass profile
      if( rdmax .le. 0. )c = gmdtab( 0.d0 )
c compute concentration and v200
      a = 3e6 / ( 8. * pi * 1.30 ) * ( hmass / hrad**3 )
      a = sqrt( a )
      c = cofa( dble( a ) )
      cNFW = c
      v200 = 0.7 * c * hrad
      m200 = hmass * ( log( 1. + c ) - c / ( 1. + c ) )
      print *, 'M_* halo =', hmass
      print *, 'r_s =', hrad
      print *, 'concentration and a( c )', cNFW, a
      print *, 'v200 and  m200', v200, m200
      a = cNFW
      call cutoff( a )
c initialize the main tables
      nradad = min( 501, mradad )
      do i = 1, nradad
        r = dble( i - 1 ) / dble( 2 * ( nradad - 1 ) )
        if( r0core )then
          r = rmax * r / ( 1. - r )
          r = max( r, 5.d-4 )
        else if( r1cusp .or. r2cusp )then
          r = rmax * ( r / ( 1. - r ) )**2
          r = max( r, 1.d-8 )
        else
          call crash( 'MAIN', 'Unknown rule for abscissae' )
        end if
        arad( i ) = r
        Phin( i ) = Phihal( r )
        gmn( i, 1 ) = gmassh( r )
        rhon( i, 1 ) = rhohal( r )
        if( r1cusp .or. r2cusp )rhon( i, 1 ) = log10( rhon( i, 1 ) )
      end do
      totmas = gmn( nradad, 1 )
c compute truncated mass
      call renew_rho( diffm )
c change halo type
      imtyp( 2 ) = 23
      ctype( 2 ) = 'ADIA'
      idftyp( 2 ) = 25
      cdft( 2 ) = 'COMP'
      snglph = .false.
      cmprssd = .true.
c compute initial potential
      call renew_Phi( diffm )
      print *, 'Notional and truncated total masses', sngl( totmas ),
     +                                          sngl( gmn( nradad, 1 ) )
      totmas = gmn( nradad, 1 )
      conserve = .false. !totmas .gt. cmpmas( 1 )
      if( conserve )cmpmas( 2 ) = ( totmas - cmpmas( 1 ) ) / totmas
c
c iterate
      iterad = 1
      convrg = .false.
      diffo = 1.d10
      do while ( .not. convrg )
        r = rmax
        print *, 'rcut = ', sngl (r)
        call cutoff( sngl( r ) )
        print *, ' after cutoff', sngl( emine ), sngl( Phimax )
c set flag to recompute tables
        iusead = 0
        call renew_rho( diffm )
        call renew_Phi( diffm )
        print *, 'Iteration', iterad, ' Max pot diff', sngl( diffm )
        if( cmpmas( icmp )
     +           .lt. 0. )call crash( 'COMPRESS', '-ve halo mass' )
        errmas = totmas - gmn( nradad, 1 )
        if( conserve )errmas = errmas - cmpmas( 1 )
        errmas = errmas / totmas
        if( diffm .lt. diffo )then
          diffo = diffm
c            convrg = diffm .lt. 1.d-5
          convrg = diffm .lt. 1.d-3
        else
          print *, 'convergence is not making good progress'
          convrg = .true. !gtlogl( 'Is the present error level acceptable' )
        end if
c        print *, 'Fractional error in total mass =', sngl( errmas )
        iterad = iterad + 1
      end do
      errmas = gmn( nradad, 1 )
      if( conserve )errmas = errmas + cmpmas( 1 )
      print *, 'Original and current total masses', sngl( totmas ),
     +                                              sngl( errmas )
c output results
      call contrast
      stop
c$$$C ###### these lines should not be needed ############
c$$$      close( ni )
c$$$      close( no )
c$$$      deallocate ( arad )
c$$$      deallocate ( phin )
c$$$      deallocate ( gmn )
c$$$      deallocate ( rhon )
c$$$      deallocate ( plamda )
c$$$      deallocate ( phinc )
c$$$      deallocate ( glamda )
c$$$      deallocate ( gmnc )
c$$$      deallocate ( rlamda )
c$$$      deallocate ( rhonc )
c ###### P. Li ############
      stop
   99 call crash( 'COMPRESS', 'Error reading command line parameters' )
      end
c ##############################################################################
      subroutine contrast
      use data_array
      implicit none
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
      include 'fit_quality.f'
c externals
      character dtype*4, htype*4
      external vbaryn, vbulge, vdisk, vgas, vtot, vhalo
      real roundup
      real*8 gmassh, splint2, vtot, vdisc, vbaryn, vbulge, vdisk, vgas
      real*8 vhalo, gmdtab, rhohal
c
c local arrays
      integer m
      parameter ( m = 2000 )
      real v( 9, m )
      real Rad( m ), V_obs( m ), eV_obs( m ), V_tot
c
c local variables
      character type*4, outs*120
      integer i, id, ifail, ih, j, n, ncols, nobs
      real a, b, m2save, mhalo, rs, vm, vobs, verr
      real*8 c, E, rkpc, rcomp, sr
      character*20 hmass_str, hrad_str, ml_str
      integer nhmass_str, nhrad_str, nml_str
c
c compute G10^10Msun/1kpc in km/s and set velocity scale
      vscale = 6.668 * 1.999e4 / 3.086
      vscale = sqrt( vscale * rscale( 1 ) / cmpmas( 1 ) )
c read observed velocities
      open( 3, file = 'data/RotationCurve/'//rootname( 1:nch )//
     +     '.dat',
     +     form = 'formatted', status = 'old' )
      do i = 1, 3
        read( 3, * )
      end do
      read( 3, *, iostat = ifail )a, vobs, verr
      i = 0
      b = 1000
      vm = 0
      do while ( ifail .eq. 0 )
        i = i + 1
        if( i .gt. m )call crash( 'CONTRAST',
     +                                'local array too small for data' )
        v( 1, i ) = a
        v( 2, i ) = vobs
        v( 3, i ) = verr
        Rad( i ) = a
        V_obs( i ) = vobs
        eV_obs( i ) = verr
        b = min( a, b )
        vm = max( vm, vobs )
        read( 3, *, iostat = ifail )a, vobs, verr
      end do
      close( 3 )
      nobs = i
c restore parameters and flags for the uncompressed halo
      imtyp( icmp ) = ihalo0( 1 )
      type = htype( ihalo0( 1 ) )
      read( type, '( a4 )' )ctype( icmp )
      m2save = cmpmas( 2 )
      cmpmas( 2 ) = 1
c set up table of radii
      rdmax = max( rdmax, Rad( nobs ) )
      rkpc = 0
      n = 0
      do while ( rkpc .lt. 1.2 * rdmax )
        if( rkpc .lt. 0.99 )then
          rkpc = rkpc + .05
        else
          rkpc = rkpc + .25
        end if
        n = n + 1
        if( n .gt. m )call crash( 'CONTRAST', 'Array too small' )
        v( 1, n ) = rkpc
c set value of V_uh
        v( 2, n ) = vhalo( rkpc * rscale( 1 ) ) * vscale
      end do
c restore parameters and flags for the compressed halo
      imtyp( icmp ) = 23
      type = 'ADIA'
      read( type, '( a4 )' )ctype( icmp )
      cmpmas( 2 ) = m2save
c create header for .out file
      open( 2, file = 'output/'//rootname( 1:nch )//'.out', 
     +      form = 'formatted',
     +      status = 'unknown' )
      write( 2, '( ''# Disk parameters'' )' )
      write( 2, '( ''# M_d/e10M_sun'', f10.4 )' )mdtot
      write( 2, '( ''#    Disk M/L '', f10.4 )' )mtold
      if( lbulge )then
        write( 2, '( ''# Bulge included'' )' )
        write( 2,
     +        '( ''# M_g/e10M_sun'', f10.4 )' )msb2( nvb ) / cmpmas( 1 )
        write( 2, '( ''#   Bulge M/L '', f10.4 )' )mtolb
      end if
      if( lgas )then
        write( 2, '( ''# Gas included'' )' )
        write( 2,
     +        '( ''# M_g/e10M_sun'', f10.4 )' )msg2( nvg ) / cmpmas( 1 )
      end if
      write( 2, '( ''# NFW halo parameters'' )' )
      write( 2, '( ''#  M_*/e10M_sun'', f10.4 )' )mhalo
      write( 2, '( ''#     r_s/1~kpc'', f10.4 )' )rs
      write( 2, '( ''# concentration'', f10.4 )' )cNFW
      write( 2, '( ''#   v_200/1km/s'', f10.4 )' )v200
      write( 2, '( ''# M_200/e10Msun'', f10.4 )' )m200
      write( 2, '( ''#+++++++++++++++'' )' )
c compose column headers
      ncols = 6
      if( lgas .or. lbulge )then
        ncols = 7
        if( lgas )ncols = ncols + 1
        if( lbulge )ncols = ncols + 1
      end if
      if( ncols .eq. 6 )then
        outs = '#   Rkpc      V_uh       V_b      V_ch     V_tot'
     +     // '     r2rho'
c               12345678901234567890123456789012345678901234567890
      else
        if( lgas .and. lbulge )then
          outs = '#   Rkpc      V_uh       V_b      V_ch     V_tot'
     +       // '     r2rho       V_disk     V_gas   V_bulge'
c                 12345678901234567890123456789012345678901234567890
        else if( lgas )then
          outs = '#   Rkpc      V_uh       V_b      V_ch     V_tot'
     +       // '     r2rho       V_disk     V_gas'
c                 12345678901234567890123456789012345678901234567890
        else
          outs = '#   Rkpc      V_uh       V_b      V_ch     V_tot'
     +       // '     r2rho       V_disk   V_bulge'
        end if
      end if
      write( 2, '( a )' )trim( outs )
      write( 2, '( ''#'' )' )
c tabulate values for the .out file
      do i = 1, n
        rkpc = v( 1, i )
        sr = rkpc * rscale( 1 )
        v( 3, i ) = vbaryn( rkpc )
        v( 4, i ) = vhalo( sr ) * vscale
        v( 5, i ) = vtot( sr ) * vscale
        v( 6, i ) = sr**2 * rhohal( sr ) * rscale( 1 ) / cmpmas( 1 )
        if( lgas .or. lbulge )then
          j = 7
          v( j, i ) = vdisk( rkpc )
          if( lgas )then
            j = j + 1
            v( j, i ) = vgas( rkpc )
          end if
          if( lbulge )then
            j = j + 1
            v( j, i ) = vbulge( rkpc )
          end if
        end if
        write( 2, '( 5f10.3, (e12.3), 3f10.3 )' ),
     +                                       ( v( j, i ), j = 1, ncols )
      end do
      close( 2 )
c Calculate chi^2 and d chi^2 
      j = 1
      chi_sq = 0
c      dchi_sq = 0
      do i = 1, nobs
        do while ( v( 1, j ) .lt. Rad( i ) )
          j = j + 1
        end do
        if( ( v( 1, j - 1 ) .le. Rad( i ) ) .and.
     +      ( v( 1, j ) .ge. Rad( i ) ) )then
          V_tot = ( v( 5, j - 1 ) * ( v( 1, j ) - Rad( i ) ) +
     +              v( 5, j ) * ( Rad( i ) - v( 1, j - 1) ) ) /
     +                                     ( v( 1, j ) - v( 1, j - 1 ) )
          chi_sq = chi_sq + ( ( V_obs( i ) - V_tot ) / eV_obs( i ) )**2
c          dchi_sq = dchi_sq + ( V_obs( i ) - V_tot ) / eV_obs( i )**2
        else
          call crash( 'CONTRAST', 'Error in interpolation')
        end if
      end do
      print *, 'Reduced chi^2 = ', chi_sq / nobs
c      print *, 'hmass = ',hmass,int(hmass*10000),int(hmass*10000)/1.0D4
      if(hmass<10)then
          write(hmass_str, '(f6.4)') int(hmass*10000)/1.0D4
      else if(hmass<100) then
          write(hmass_str, '(f7.4)') int(hmass*10000)/1.0D4
      else if(hmass<1000) then
          write(hmass_str, '(f8.4)') int(hmass*10000)/1.0D4
      else if(hmass<10000) then
          write(hmass_str, '(f9.4)') int(hmass*10000)/1.0D4
      else if(hmass<100000) then
          write(hmass_str, '(f10.4)') int(hmass*10000)/1.0D4
      end if

      if(hrad<10)then
          write(hrad_str, '(f6.4)') int(hrad*10000)/1.0D4
      else if(hrad<100) then
          write(hrad_str, '(f7.4)') int(hrad*10000)/1.0D4
      else if(hrad<1000) then
          write(hrad_str, '(f8.4)') int(hrad*10000)/1.0D4
      else if(hrad<10000) then
          write(hrad_str, '(f9.4)') int(hrad*10000)/1.0D4
      end if

      write(ml_str, '(f4.2,a1,f4.2)') mtold, '_', mtolb
      nml_str = lnblnk(ml_str)
      nhmass_str = lnblnk(hmass_str)
      nhrad_str = lnblnk(hrad_str)
      open( 2, file = 'output/'//rootname(1:nch)//'_'//ml_str(1:nml_str)
     + //'_'//hmass_str(1:nhmass_str)//'_'//hrad_str(1:nhrad_str)
     + //'.fitquality', 
     + form ='formatted', status = 'unknown' )
      write( 2, * )chi_sq / nobs, nobs
      print *, '\n', '######## End of this iteration #########', '\n'
      close( 2 )
      return
      end
c ##############################################################################
      real*8 function vtot( r )
      use data_array 
      implicit none
c
      real*8 r, frhalo, vbaryn
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      vtot = -r * frhalo( r ) + ( vbaryn( r / rscale(1) ) / vscale )**2
      vtot = sqrt( vtot )
      return
      end
c ##############################################################################
      real*8 function vbaryn( r )
      use data_array
      implicit none
c
      real*8 r
c
      real*8 vbulge, vdisk, vgas
c
      real*8 v2
c
      vbaryn = vdisk( r )
      if( lgas .or. lbulge )then
        v2 = vbaryn**2
        if( lgas )v2 = v2 + vgas( r )**2
        if( lbulge )v2 = v2 + vbulge( r )**2
        vbaryn = sqrt( v2 )
      end if
      return
      end
c ##############################################################################
      real*8 function vdisk( r )
      use data_array
      implicit none
c
      real*8 r, splint2
c
      if( iused .le. 0 )call crash( 'VDISK',
     +                                    'called before data is read' )
c linear rise inside innermost point
      if( r .lt. rdtab( 2 ) )then
        vdisk = veld( 2 ) * r / rdtab( 2 )
      else if( r .gt. rdtab( nvd ) )then
        vdisk = veld( nvd )
      else
c use spline
        vdisk = splint2( rdtab( 2 ), veld( 2 ), nvd - 1,
     +                   r, vdlam, vdc, .false. )
      end if
      return
      end
c ##############################################################################
      real*8 function vbulge( r )
      use data_array
      implicit none
c
      real*8 r, splint2
c
      if( iused .le. 0 )call crash( 'VBULGE',
     +                                    'called before data is read' )
c linear rise inside innermost point
      if( r .lt. rdtab( 2 ) )then
        vbulge = velb( 2 ) * r / rdtab( 2 )
      else if( r .gt. rdtab( nvb ) )then
        vbulge = velb( nvb )
      else
c use spline
        vbulge = splint2( rdtab( 2 ), velb( 2 ), nvb - 1,
     +                    r, vblam, vbc, .false. )
      end if
      return
      end
c ##############################################################################
      real*8 function vgas( r )
      use data_array
      implicit none
c
      real*8 r, splint2
c
      if( iused .le. 0 )call crash( 'VGAS',
     +                                    'called before data is read' )
c linear rise inside innermost point
      if( r .lt. rdtab( 2 ) )then
        vgas = velg( 2 ) * r / rdtab( 2 )
      else if( r .gt. rdtab( nvg ) )then
        vgas = velg( nvg )
      else
c use spline
        vgas = splint2( rdtab( 2 ), velg( 2 ), nvg - 1,
     +                   r, vglam, vgc, .false. )
      end if
      return
      end
c ##############################################################################
      subroutine renew_Phi( difm )
      use aarrays
      implicit none
c
c calling argument
      real*8 difm
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 gmnx, Phinx, splint2
c
c local variables
      integer i
      real*8 r, s
c
c compute new mass array
      do i = 1, nradad
        r = arad( i )
        gmn( i, 1 ) = gmnx( r )
      end do
c fit new spline
      r = arad( 1 )
      s = splint2( arad, gmn, nradad, r, glamda, gmnc, .true. )
c compute new potential array
      difm = 0
      do i = 1, nradad
        r = arad( i )
        s = Phinx( r )
c        difm = max( difm, abs( s - Phin( i ) ) )
        if( abs( s - Phin( i ) ) .gt. difm )then
          difm = abs( s - Phin( i ) )
        end if
        Phin( i ) = s
      end do
c fit new spline
      r = arad( 1 )
      s = splint2( arad, Phin, nradad, r, plamda, Phinc, .true. )
      return
      end
c ##############################################################################
      real*8 function gmnx( r )
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
c external
      real*8 rhohal
c
c local arrays
      integer npts
      parameter ( npts = 64 )
      real*8 absc( npts ), wght( npts )
c
c local variables
      integer i, n
      real*8 rp
      include 'inc/pi.f'
c
      gmnx = 0
      if( r .gt. 0.d0 )then
c Gauss-Legendre quadrature over interior mass
        n = npts
        if( r .lt. 1.d0 )n = r * real( npts )
        n = max( n, 2 )
        n = min( n, npts )
        call GLqtab( 0.d0, r, npts, wght, absc )
        do i = 1, npts
          rp = absc( i )
          gmnx = gmnx + wght( i ) * rp**2 * rhohal( rp )
        end do
        gmnx = 4. * pi * gmnx
      end if
      return
      end
c ##############################################################################
      real*8 function Phinx( r )
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
      external frhalo
      real*8 gmassh, phidtab, quad_gnr
c
c local variables
      integer i, ifail, n
      real*8 abserr, epsr, rp
c
      Phinx = 0
      if( r .lt. rmax )then
c integrate central attraction
        epsr = 1.d-30
        ifail = 1
        Phinx = quad_gnr( frhalo, r, rmax, epsr, epsr, ifail )
        if( ifail .gt. 5 )then
          print *, 'ifail =', ifail, ' from quad_gnr'
          call crash( 'PHINX', 'Quadrature error' )
        end if
      end if
c add the contribution to infinity assuming density cuts off at rmax
      rp = max( r, rmax )
      Phinx = Phinx - gmassh( rmax ) / rp
c add the disk contribution
      Phinx = Phinx + phidtab( r )
      return
      end
c ##############################################################################
      subroutine renew_rho( difm )
      use aarrays
      implicit none
c
c calling argument
      real*8 difm
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c externals
      real*8 rhoint, splint2
c
c local variables
      integer i
      real*8 r, rho, rhol
c
c update local density
      difm = 0
      do i = 1, nradad - 1
        r = arad( i )
        rho = rhoint( r )
c compute max change
        rhol = rhon( i, 1 )
        if( r1cusp .or. r2cusp )then
          if( rhol .gt. -8.d0 )then
            rhol = 10.**rhol
          else
            rhol = 0
          end if
        end if
        difm = max( difm, abs( rhol - rho ) )
        if( rho .lt. 0.d0 )then
          print *, 'Negative halo density!'
          print *, i, sngl( arad( i ) ), sngl( rho ), sngl( rhol )
          if( rho .lt. -1.d-4 )read *, rhol
          rho = 0
        end if
c save new value
        rhon( i, 1 ) = rho
        if( r1cusp .or. r2cusp )then
          if( rho .gt. 0.d0 )then
            rhon( i, 1 ) = log10( rho )
          else
            rhon( i, 1 ) = -8
          end if
        end if
      end do
      if( r1cusp .or. r2cusp )then
        rhon( nradad, 1 ) = -8
      else
        rhon( nradad, 1 ) = 0
      end if
c fit new spline
      r = arad( 1 )
      rho = splint2( arad, rhon, nradad, r, rlamda, rhonc, .true. )
      return
      end
c ##############################################################################
      real*8 function rhoint( r )
      implicit none
c Performs a double integration over all allowed velocities
c
c calling argument
      real*8 r
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
      real*8 E, rr, pot
      common / intrho / E, rr, pot
c
c externals
      external rhoifn
      real*8 distfn, Phitot, quad_Pat
c
c local arrays
      integer npts
      parameter ( npts = 32 )
      real*8 absc( npts ), wght( npts )
c
c local variables
      integer i, ifail, n, nlim
      real*8 epsr, f, lmax, rel, u
      include 'inc/pi.f'
c
c set integration range
      rr = r
      pot = Phitot( rr )
c Gauss-Legendre quadrature over energies
      ifail = 0
      call GLqtab( pot, Phimax, npts, wght, absc )
      rhoint = 0
      do i = 1, npts
        E = absc( i )
        if( r .gt. 0.d0 )then
c adaptive integration over L
          lmax = rr * sqrt( 2. * ( E - pot ) )
          nlim = 0
          epsr = 1.e-6
          ifail = 1
          f = quad_Pat( rhoifn, 0.d0, lmax, epsr, ifail )
c ifail = 1 means requested precision not achieved
          if( ifail .gt. 1 )then
            print *, 'IFAIL =', ifail, ' from quad_Pat'
            call crash( 'RHOINT', 'Quadrature error' )
          end if
          rhoint = rhoint + wght( i ) * f
        else
          u = 2 * ( E - pot )
          if( u .lt. 0. )call crash( 'RHOINT', 'E out of range' )
          u = sqrt( u )
          rhoint = rhoint + wght( i ) * u * distfn( E, 0.d0 )
        end if
      end do
      rhoint = 4. * pi * rhoint
      return
      end
c ##############################################################################
      real*8 function rhoifn( L )
      implicit none
c
c calling argument
      real*8 L
c
c common blocks
c
      real*8 E, rr, pot
      common / intrho / E, rr, pot
c
      include 'inc/params.f'
c
      include 'inc/model.f'
c
c external
      real*8 distfn
c
c local variable
      real*8 u
c
      u = 2. * ( E - pot ) - ( L / rr )**2
      if( u .le. 0. )then
        if( u .lt. -1.d-5 )then
          print *, E, pot, L, rr, u
          call crash( 'RHOIFN', 'Argument out of range' )
        end if
        rhoifn = 0
      else
        u = sqrt( u )
        rhoifn = L * distfn( E, L ) / ( rr**2 * u )
      end if
      return
      end
c ##############################################################################
      real*8 function cofa( a )
c returns the value of c for a given A(c) for NFW halos
c
c calling argument
      real*8 a
c
c external
      real*8 aofc
c
c local array
      real*8 w( 17 )
c
c local variables
      integer ifail, ind, ir
      real*8 c1, c2, res, tol
c
      ir = 0
      ind = 1
      ifail = 0
      tol = 1.e-8
      c1 = 1.e-3
      c2 = 100.
      do while ( ind .ne. 0 )
        call fndzro( c1, c2, res, tol, ir, w, ind, ifail )
        res = a - aofc( c1 )
      end do
      cofa = c1
      return
      end
c ##############################################################################
      real*8 function aofc( c )
c returns the function A(c) needed for the velocity scales of NFW halos
c
c calling argument
      real*8 c
c
      aofc = sqrt( c**3 / ( log( 1. + c ) - c / ( 1. + c ) ) )
      return
      end

C      include 'frtot.f'
c      include 'getset.f'
c      include 'gmdtab.f'
c      include 'phidtab.f'
c      include 'phitot.f'
c      include 'rlims.f'
