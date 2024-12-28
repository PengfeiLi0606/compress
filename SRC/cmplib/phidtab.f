      real*8 function phidtab( r )
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
c external
      real*8 gmdtab
c
c local variables
      integer i, npts
      real*8 r2
      parameter ( npts = 64 )
      real*8 absc( npts ), wght( npts )
c
c disk mass is zero beyond rdmax
      r2 = rdmax / rscale( 1 )
      if( r .lt. r2 )then
c integral of radial force from r to rdmax
        call GLqtab( r, r2, npts, wght, absc )
        phidtab = 0
        do i = 1, npts
          r2 = absc( i )
          phidtab = phidtab - wght( i ) * gmdtab( r2 ) / r2**2
        end do
c add value at rdmax
        phidtab = phidtab - gmdtab( r2 ) / r2
      else
c field point at large radius
        phidtab = -gmdtab( r ) / r
      end if
      return
      end
