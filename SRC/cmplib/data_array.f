      module data_array 
c mostly for input data arrays
      implicit none
c
c "common" variables
      character*40 rootname
      integer iused, nch, nvb, nvd, nvg
      logical lgas, lbulge
c ###### naoise added kappa as shape parameter for einasto in next line ############
      real mtolb, mtold, rd, cNFW, v200, m200!, kappa
      real rdmax, mdtot, vscale
      real*8 hmass, hrad
      save iused, nch, nvb, nvd, nvg, lgas, lbulge!, kappa
      save mtolb, mtold, rd, cNFW, v200, m200, rdmax, mdtot, vscale
      save hmass, hrad
c
c stellar disk
      real*8, allocatable, save :: rdtab( : ), veld( : )
      real*8, allocatable, save :: vdlam( : ), vdc( : )
      real*8, allocatable, save :: rdtabs( : )
c
      real*8, allocatable, save :: sdsd( : ), msd2( : )
      real*8, allocatable, save :: sdlam( : ), sdc( : )
c
c stellar bulge
      real*8, allocatable, save :: velb( : )
      real*8, allocatable, save :: vblam( : ), vbc( : )
c
      real*8, allocatable, save :: sdsb( : ), msb2( : )
      real*8, allocatable, save :: sblam( : ), sbc( : )
c
c gas disk
      real*8, allocatable, save :: velg( : )
      real*8, allocatable, save :: vglam( : ), vgc( : )
c
      real*8, allocatable, save :: sdsg( : ), msg2( : )
      real*8, allocatable, save :: sglam( : ), sgc( : )
c
      end module data_array
