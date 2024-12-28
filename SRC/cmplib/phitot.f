      real*8 function phitot( r )
      use aarrays
      implicit none
c returns potential in the mid-plane due to both disc and halo mass compnts
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
      external frtot
      real*8 frtot, phidsc, phihal, phispl, phidtab, splint2
c, phitab
c
c local arrays
      integer npts
      parameter ( npts = 32 )
      real*8 absc( npts ), wght( npts )
c
c local variables
      integer i, jcmp, k
      real*8 gm, haloc, hrad, r1, r2
c
c preserve icmp
      jcmp = icmp
      Phitot = 0
c case needed only at start of compress
      if( snglph )then
        Phitot = Phihal( r )
      else if( numcal )then
c numerical potential, value in mid-plane only for spheroidal systems
        do icmp = 1, ncmp
          if( sphrod( icmp ) )phitot = phispl( r, 0.d0 )
        end do
      else if( fixed )then
c fixed rotation curve model
        icmp = 0
        do i = 1, ncmp
          if( .not. disc( i ) )then
            haloc = dfcns( 3, i )
            hrad = rscale( i )
            icmp = i
          end if
        end do
        if( icmp .eq. 0 )call crash( 'PHITOT',
     +            'No halo cmp found for a fixed rotation curve model' )
c Fall-Efstathiou model
        if( ctype( icmp ) .eq. 'FE  ' )then
          phitot = phitot + .5 * haloc**2 * log( ( r**2 + hrad**2 )
     +         / ( rmax**2 + hrad**2 ) )
c power law rotation curve - choose zero at r = 1
        else if( ctype( icmp ) .eq. 'POWE' )then
          r1 = r
          if( haloc .lt. .01 )r1 = max( r, 1.d-12 )
          if( abs( haloc ) .lt. .01 )then
            phitot = log( r1 )
          else
            phitot = ( r1**( 2. * haloc ) - 1. ) / ( 2. * haloc )
          end if
        else
c assume halo cuts off sharply at rmax
          phitot = rmax * frtot( rmax )
          if( r .ge. rmax )then
            phitot = phitot * rmax / r
          else
c add integral of radial force from r to rmax
            r2 = r
            call GLqtab( r2, rmax, npts, wght, absc )
            do i = 1, npts
              r2 = absc( i )
              phitot = phitot + wght( i ) * frtot( r2 )
            end do
          end if
        end if
c adiabatically compressed halo - common potential for all components
      else if( cmprssd )then
        if( r .gt. arad( nradad ) )then
          Phitot = Phin( nradad )
c rigid components may have non-zero density at these radii
          do irigid = 1, nrigid
            icmp = ircmp( irigid )
            Phitot = Phitot - Phidtab( arad( nradad ) ) + Phidtab( r )
          end do
c compressible components may have zero density at these radii
          gm = 0
          do i = 1, ncomp
            gm = gm + gmn( nradad, i )
          end do
          Phitot = Phitot + gm * ( 1. / arad( nradad ) - 1. / r )
        else if( r .lt. arad( 1 ) )then
c linear extrapolation to zero radius
          r1 = ( Phin( 2 ) - Phin( 1 ) ) / ( arad( 2 ) - arad( 1 ) )
          Phitot = Phin( 1 ) + r1 * ( r - arad( 1 ) )
        else
c spline interpolation - assume knots and coeffs are already set
          Phitot = splint2( arad, Phin, nradad, r, Plamda, Phinc,
     +                      .false. )
        end if
      else
        do icmp = 1, ncmp
          if( cmpmas( icmp ) .gt. 0. )then
            if( disc( icmp ) )then
c disc component
              if( ctype( icmp ) .eq. 'MTAB' )then
                phitot = phitot + phidtab( r )
              else
                phitot = phitot + phidsc( r )
              end if
            else
              if( ctype( icmp ) .eq. 'POWC' )then
c power law with a core
                call crash( 'PHITOT',
     +                      'Power-law option not programmed' )
c              if( abs( dfcns( 3, icmp ) ) .lt. .01 )then
c                phitot = .5 * ( log( 1. + ( r / hrad )**2 ) - 5. ) / hrad
c              else
c                phitot = ( 1. + ( r / hrad )**2 )**( -.5 * dfcns( 3, icmp ) )
c                if( dfcns( 3, icmp ) .lt. 0.d0 )
c     +          phitot = ( rmax / hrad )**( -dfcns( 3, icmp ) ) - phitot
c                phitot = -phitot / hrad
c              end if
              else
c halo component
                phitot = phitot + phihal( r )
c adiabatically compressed halo
                if( ctype( icmp ) .eq. 'ADIA' )then
                  if( ctype( 1 ) .eq. 'EXP ' )then
c exponential sphere
                    r2 = r / rscale( 1 )
                    if( r2 .gt. 1.d-4 )then
                      phitot = phitot -
     +                             cmpmas( 1 ) * ( 1. - exp( -r2 ) ) / r
                    else
c expansion for small arguments
                      phitot = phitot - cmpmas( 1 ) *
     +         ( 1. - .5 * r2 + r2**2 / 6. - r2**3 / 24. ) / rscale( 1 )
                    end if
                  else if( ctype( 1 ) .eq. 'MTAB' )then
                    r2 = phitot
                    phitot = phitot + phidtab( r )
                  else
c softened point mass
c                   phitot = phitot -
c     +                 cmpmas( 1 ) / sqrt( r**2 + rscale( 1 )**2 )
                    call crash( 'PHITOT',
     +                           'need formula for external potential' )
                  end if
                end if
              end if
            end if
          end if
        end do
      end if
c restore icmp
      icmp = jcmp
      return
      end
