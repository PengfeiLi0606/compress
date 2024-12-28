      real*8 function dfwght( E, Lz )
c function to weight DF when unequal particle masses are requested
c  value of this function is the inverse of the desired particle weights
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c
c calling arguments
      real*8 E, Lz
c
c common blocks
c
      include 'inc/params.f'
c
      include 'inc/admin.f'
c
      include 'inc/model.f'
c
c external
      logical gtlogl
c
c local array
      logical firstc( mcmp )
      save firstc
c
      character*2 cn
c
      data firstc / mcmp * .true. /
c
      if( firstc( icmp ) )then
        if( uqmass )then
          write( cn, '( i2 )' )icmp
          if( master )print *, 'This default version of DFwght assumes'
     +     // ' all particles in component ' // cn // ' have equal mass'
          if( .not. gtlogl( 'Is this what you want' )
     +                             )call crash( 'DFWGHT', 'User abort' )
          uqmass = .false.
        end if
        firstc( icmp ) = .false.
      end if
c all particles in this component will have equal masses
      dfwght = 1
c example: weight as inverse of angular momentum with a small offset
c      dfwght = 1. / ( .01 + Lz )
      return
      end
