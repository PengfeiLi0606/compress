      subroutine readdat( unitno, line, ifail )
c Created: 12 Apr 1992 by David Earn
c Extensively modified by Jerry Sellwood and included in this package
c   with permission
      implicit none
c This routine is called only from getline and only on the master node.
c   It is separated from getline in order to minimize duplication of
c   code in the mpi version
c
c            : Return a string (line) of length len(line) containing
c either the ascii end of file character or the next input line from the
c text input stream with logical unit number unitno.  All alphabetic
c characters are converted to uppercase.  All tabs are converted to
c single spaces (to avoid undesirable behaviour: e.g., apparently blank
c lines containing tabs are treated as strings otherwise).  Blank lines
c and comments are ignored.  A comment begins with a sharp sign (#) and
c ends with a newline (return) as in unix.  In addition, if any line
c begins with the words `end of file' then getline returns line = eof.
c This makes it possible to store alternate data cards at the end of a
c data file without having to comment them out.  Getline should normally
c be used in constructions such as the following.
c
c       call getline( ni, line )
c       do while ( line .ne. eof )
c           read( line, ... ) ...
c           ...
c           call getline( ni, line )
c       end do
c
c Note: (1) statement labels are not required in the calling module in
c order to deal with end of file, (2) string comparisons against line
c make it possible to use obvious ending codes rather than nonsense
c values of numbers being input, (3) line can be picked apart in steps
c by reading it several times, and (4) the character variable eof must
c be set to char(4) in the calling module.
c
c calling arguments
      integer ifail, unitno
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
      integer i, io, j, k, l, m, maxline
      character eof, tab
c
      ifail = 0
c eof is ascii end of file (or end of transmission)
      eof = char( 4 )
c tab character
      tab = char( 9 )
c make line blank to begin with
      maxline = len( line )
      do i = 1, maxline
        line( i:i ) = ' '
      end do
c read line(s) from the formatted input file until a meaningful line is found
      do while ( line .eq. ' ' )
        read( unitno, '( a )', iostat = io )line
        if( io .eq. 0 )then
c remove comments
          i = index( line, '#' )
          if( i .gt. 0 )then
            do j = i, maxline
              line( j:j ) = ' '
            end do
          end if
c search for and replace tabs with a suitable number of blanks
c new code added by JAS
          j = 0
          do while ( j .lt. maxline )
            j = j + 1
            if( line( j:j ) .eq. tab )then
c this assumes modulo 10 for tabs
              k = 10 * ( 1 + j / 10 )
              l = lnblnk( line )
c move subsequent non-blank characters out along the line
              if( l .gt. j )then
c number of blank characters to insert to replace tab
                m = k - j
c abandon detab if line is not long enough
                if( m + l .gt. maxline )then
                  ifail = 1
                  return
                end if
                do i = l, j + 1, -1
                  line( i+m:i+m ) = line( i:i )
                end do
              end if
c blank out tab, together with added space
              do i = j, k
                line( i:i ) = ' '
              end do
              j = k
            end if
          end do
        else
c read failed
          ifail = 2
          line = eof
        end if
      end do
c convert all alphabetic characters to uppercase
      call uppercase( line )
c watch out for end of file indicator
      if( ( maxline .ge. 11 ) .and.
     +    ( line .eq. 'END OF FILE' ) )line = eof
      return
      end
