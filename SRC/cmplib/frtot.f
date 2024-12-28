      real*8 function frtot( r )
      implicit none
c returns total radial force in the mid-plane for mass model
c
c calling argument
      real*8 r
c
c common block
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c externals
      external phitot
      real*8 deriv2, frdisc, frhalo, gmdtab, vcirc
c
c local variables
      integer j, jcmp, jp, kp
      real*8 r2
c
      frtot = 0.
      if( r .eq. 0. )return
c arbitrary model - numerical potential only
      if( numcal )then
        frtot = -deriv2( r, phitot )
      else if( fixed )then
c fixed rotation curve models
        frtot = -vcirc( r )**2 / r
      else
        jcmp = icmp
c adiabatically compressed model
        if( cmprssd )then
          j = icomp
          do icomp = 1, ncomp
            icmp = iccmp( icomp )
c live halo component
            if( cmpmas( icmp ) .gt. 0. )frtot = frtot + frhalo( r )
          end do
          icomp = j
          j = irigid
          do irigid = 1, nrigid
            icmp = ircmp( irigid )
c rigid adiabatically compressing masses are assumed spherical
            if( cmpmas( icmp ) .gt. 0. )then
              if( rigidp( icmp ) .or. ( .not. disc( icmp ) ) )then
c mean spherical attraction from the disk
                r2 = r / rscale( icmp )
                if( ctype( icmp ) .eq. 'EXP ' )then
c exponential sphere
                  if( r2 .gt. 1.d-4 )then
                    frtot = frtot - cmpmas( icmp ) *
     +                          ( 1. - ( 1. + r2 ) * exp( -r2 ) ) / r**2
                  else
c expansion for small arguments
                    frtot = frtot - cmpmas( icmp ) *
     +                 ( .5 + r2 / 3. + r2**2 / 8. ) / rscale( icmp )**2
                  end if
                else if( ctype( icmp ) .eq. 'MTAB' )then
c disk M(r) from a table
                  if( r .gt. 0.d0 )frtot = frtot - gmdtab( r ) / r**2
                else
c softened point mass
c                 frtot = frtot - 1.d-3 * r / ( r * r + 1.d-4 )**1.5
                  call crash( 'FRTOT',
     +                           'need formula for external potential' )
                end if
c               frtot = frtot + frext( r )
              else
c theoretical attraction in the disk mid-plane
                frtot = frtot + frdisc( r )
              end if
            end if
          end do
          irigid = j
        else
c regular model
          if( snglph )then
            jp = icmp
            kp = icmp
          else
            jp = 1
            kp = ncmp
          end if
          do icmp = jp, kp
c       do icmp = 1, ncmp
            if( cmpmas( icmp ) .gt. 0. )then
              if( disc( icmp ) )then
                frtot = frtot + frdisc( r )
              else
                frtot = frtot + frhalo( r )
              end if
            end if
          end do
        end if
c restore icmp
        icmp = jcmp
      end if
      return
      end
