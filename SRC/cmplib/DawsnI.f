      real*8 function DawsnI( x, ifail )
      implicit none
c returs the value of Dawson's integral for the argument x
c   follows the algorithm in the NAG source code
c   the 2nd ifail is unused
c
c calling arguments
      integer ifail
      real*8 x
c
c local variables
      real*8 t, u, xmax, xmin, xmxinv, y
c
      data xmin,xmax / 3.1d-9, 3.2d+8 /
      data xmxinv / 2.24d+307 /
c     xmxinv = 0.5/minreal  (rounded down)
c
      ifail = 0
      u = abs(x)
      if( x .eq. 0.d0 )then
        DawsnI = 0.d0
      else if( u .lt. xmin )then
        DawsnI = x
      else
c small arguments
        if( u .le. 4.d0 )then
          t = x * x * 0.125d0 - 1.d0
c Chebyshev expansion evaluated as y(t) - precision 17e.18
          y = ( ( ( ( ( ( ( ( ( ( ( ( ( ( ( 4.73110787583905890d-10
     +   * t - 1.84156612203899571d-9 ) * t + 3.39447752478095416d-9 )
     +   * t - 1.19747678960170492d-8 ) * t + 5.26701186204891318d-8 )
     +   * t - 1.76391854013499939d-7 ) * t + 5.45102274319414278d-7 )
     +   * t - 1.68329667633696962d-6 ) * t + 5.02166959936255662d-6 )
     +   * t - 1.42799530892351998d-5 ) * t + 3.88215258394090029d-5 )
     +   * t - 1.00839699013042906d-4 ) * t + 2.49599593876339565d-4 )
     +   * t - 5.87289824114168026d-4 ) * t + 1.31037931970268204d-3 )
     +   * t - 2.76505243767561823d-3 ) * t + 5.50158555493458225d-3
          y = ( ( ( ( ( ( ( ( ( ( ( ( ( y * t - 1.02887338764824911d-2 )
     +   * t + 1.80231969369657062d-2 )
     +   * t - 2.94652501039055316d-2 ) * t + 4.47867518277388141d-2 )
     +   * t - 6.30546446294856549d-2 ) * t + 8.19485953308367479d-2 )
     +   * t - 9.80825026698305700d-2 ) * t + 1.08086253592738758d-1 )
     +   * t - 1.10084144800233882d-1 ) * t + 1.04702304155811420d-1 )
     +   * t - 9.47947597433254029d-2 ) * t + 8.39163150531250491d-2 )
     +   * t - 7.45921286034793947d-2 ) * t + 6.75990739533505166d-2
          DawsnI = x * y
        else
c large arguments
          if( u .lt. xmax )then
            t = 32.d0 / ( x * x ) - 1.d0
c Chebyshev expansion evaluated as y(t) -- precision 17e.18
            y = ( ( ( ( ( ( ( ( ( ( ( ( ( ( ( 2.74025015401477308d-10
     +   * t - 1.92777411862777037d-10 ) * t - 2.30082253585696739d-9 )
     +   * t + 1.68458106112580151d-9 ) * t + 8.95290316416285986d-9 )
     +   * t - 6.92502020367956994d-9 ) * t - 2.16297541013440026d-8 )
     +   * t + 1.80075901126820663d-8 ) * t + 3.70364318218408061d-8 )
     +   * t - 3.37996293801057122d-8 ) * t - 4.97301908488584026d-8 )
     +   * t + 4.97624015903312482d-8 ) * t + 5.90832390144876555d-8 )
     +   * t - 6.03608014702259991d-8 ) * t - 7.23521212676254469d-8 )
     +   * t + 5.77122948307144822d-8 ) * t + 9.90333958007913327d-8
            y = ( ( ( ( ( ( ( ( ( ( ( ( y * t - 2.05663666463271863d-8 )
     +   * t - 1.26392804533869908d-7 )
     +   * t - 9.09755613786389276d-8 ) * t + 3.64246653663605826d-8 )
     +   * t + 1.64902654547901545d-7 ) * t + 2.97061775695647310d-7 )
     +   * t + 5.96822092643564129d-7 ) * t + 1.73792937740632787d-6 )
     +   * t + 7.50767774865078295d-6 ) * t + 4.79259370751892485d-5 )
     +   * t + 4.76832273761536912d-4 ) * t +
     +         8.64607314481504465d-3 ) * t + 5.08210986449041748d-1
            DawsnI = y / x
          else
c very large arguments
            if( u .gt. xmxinv )then
              DawsnI = 0.d0
            else
              DawsnI = 0.5d0 / x
            end if
          end if
        end if
      end if
      return
      end
