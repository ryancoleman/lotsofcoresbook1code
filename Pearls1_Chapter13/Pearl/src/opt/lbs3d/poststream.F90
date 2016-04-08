!-------------------------------------------------------------------------------
! Subroutine : PostStream
! Revision   : 1.2 (2013/09/10)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Update f in nodes at processor boundaries.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
! Copyright 2008 Carlos Rosales Fernandez, David S. Whyte and IHPC (A*STAR).
!
! This file is part of MP-LABS.
!
! MP-LABS is free software: you can redistribute it and/or modify it under the
! terms of the GNU GPL version 3 or (at your option) any later version.
!
! MP-LABS is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! MP-LABS, in the file COPYING.txt. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------

 SUBROUTINE PostStream

!  Common Variables
 USE Domain
 USE FluidParams, ONLY : phi
 USE LBMParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k

!$OMP PARALLEL

! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     f2( NX + NXG*( j + NYG*k ) + nxt ) = f2( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     f1( 1 + NXG*( j + NYG*k ) + nxt ) = f1( NX+1 + NYG*( j + NZG*k ) + nxt )
   END DO
 END DO

! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     f4( i + NXG*( NY + NYG*k ) + nxt ) = f4( i + NXG*( NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     f3( i + NXG*( 1 + NYG*k ) + nxt ) = f3( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
 END DO

! Z direction
!$OMP DO PRIVATE(i)
 DO j = 1, NY
    DO i = 1, NX
     f6( i + NXG*( j + NYG*NZ ) + nxt ) = f6( i + NXG*( j ) + nxt )
   END DO
   DO i = 1, NX
     f5( i + NXG*( j + NYG ) + nxt ) = f5( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
 END DO

!
! At this stage all values of f needed for the next step have been updated
!


!-------- Exchange data for phi and inward g components ------------------------
! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     phi( NX+1 + NXG*( j + NYG*k ) ) = f0( 1 + NXG*( j + NYG*k )       ) &
                                     + f1( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f2( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f3( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f4( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f5( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f6( 1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g2( NX + NXG*( j + NYG*k ) + nxt ) = g2( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g8( NX + NXG*( j + NYG*k ) + nxt ) = g8( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g10( NX + NXG*( j + NYG*k ) + nxt ) = g10( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g12( NX + NXG*( j + NYG*k ) + nxt ) = g12( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g14( NX + NXG*( j + NYG*k ) + nxt ) = g14( NXG*( j + NYG*k ) + nxt )
   END DO
 END DO

!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     phi( NXG*( j + NYG*k ) ) = f0( NX + NXG*( j + NYG*k )       ) &
                              + f1( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f2( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f3( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f4( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f5( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f6( NX + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g1( 1 + NXG*( j + NYG*k ) + nxt ) = g1( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g7( 1 + NXG*( j + NYG*k ) + nxt ) = g7( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g9( 1 + NXG*( j + NYG*k ) + nxt ) = g9( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g11( 1 + NXG*( j + NYG*k ) + nxt ) = g11( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     g13( 1 + NXG*( j + NYG*k ) + nxt ) = g13( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
 END DO

! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     phi( i + NXG*( NY+1 + NYG*k ) ) = f0( i + NXG*( 1 + NYG*k )       ) &
                                     + f1( i + NXG*( 1 + NYG*k ) + nxt ) &
                                     + f2( i + NXG*( 1 + NYG*k ) + nxt ) &
                                     + f3( i + NXG*( 1 + NYG*k ) + nxt ) &
                                     + f4( i + NXG*( 1 + NYG*k ) + nxt ) &
                                     + f5( i + NXG*( 1 + NYG*k ) + nxt ) &
                                     + f6( i + NXG*( 1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     g4( i + NXG*( NY + NYG*k ) + nxt ) = g4( i + NXG*( NYG*k) + nxt )
   END DO
   DO i = 1, NX
     g8( i + NXG*( NY + NYG*k ) + nxt ) = g8( i + NXG*( NYG*k) + nxt )
   END DO
   DO i = 1, NX
     g9( i + NXG*( NY + NYG*k ) + nxt ) = g9( i + NXG*( NYG*k) + nxt )
   END DO
   DO i = 1, NX
     g16( i + NXG*( NY + NYG*k ) + nxt ) = g16( i + NXG*( NYG*k) + nxt )
   END DO
   DO i = 1, NX
     g18( i + NXG*( NY + NYG*k ) + nxt ) = g18( i + NXG*( NYG*k) + nxt )
   END DO
 END DO

!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     phi( i + NXG*( NYG*k ) ) = f0( i + NXG*( NY + NYG*k )       ) &
                              + f1( i + NXG*( NY + NYG*k ) + nxt ) &
                              + f2( i + NXG*( NY + NYG*k ) + nxt ) &
                              + f3( i + NXG*( NY + NYG*k ) + nxt ) &
                              + f4( i + NXG*( NY + NYG*k ) + nxt ) &
                              + f5( i + NXG*( NY + NYG*k ) + nxt ) &
                              + f6( i + NXG*( NY + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     g3( i + NXG*( 1 + NYG*k ) + nxt ) = g3( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     g7( i + NXG*( 1 + NYG*k ) + nxt ) = g7( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     g10( i + NXG*( 1 + NYG*k ) + nxt ) = g10( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     g15( i + NXG*( 1 + NYG*k ) + nxt ) = g15( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     g17( i + NXG*( 1 + NYG*k ) + nxt ) = g17( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
 END DO


! Z direction
!$OMP DO PRIVATE(i)
 DO j = 1, NY
   DO i = 1, NX
     phi( i + NXG*( j + NYG*(NZ+1) ) ) = f0( i + NXG*( j + NYG )       ) &
                                       + f1( i + NXG*( j + NYG ) + nxt ) &
                                       + f2( i + NXG*( j + NYG ) + nxt ) &
                                       + f3( i + NXG*( j + NYG ) + nxt ) &
                                       + f4( i + NXG*( j + NYG ) + nxt ) &
                                       + f5( i + NXG*( j + NYG ) + nxt ) &
                                       + f6( i + NXG*( j + NYG ) + nxt )
   END DO
   DO i = 1, NX
     g6( i + NXG*( j + NYG*NZ ) + nxt ) = g6( i + NXG*j + nxt )
   END DO
   DO i = 1, NX
     g12( i + NXG*( j + NYG*NZ ) + nxt ) = g12( i + NXG*j + nxt )
   END DO
   DO i = 1, NX
     g13( i + NXG*( j + NYG*NZ ) + nxt ) = g13( i + NXG*j + nxt )
   END DO
   DO i = 1, NX
     g16( i + NXG*( j + NYG*NZ ) + nxt ) = g16( i + NXG*j + nxt )
   END DO
   DO i = 1, NX
     g17( i + NXG*( j + NYG*NZ ) + nxt ) = g17( i + NXG*j + nxt )
   END DO
 END DO

!$OMP DO PRIVATE(i)
 DO j = 1, NY
   DO i = 1, NX
     phi( i + NXG*j ) = f0( i + NXG*( j + NYG*NZ )       ) &
                      + f1( i + NXG*( j + NYG*NZ ) + nxt ) &
                      + f2( i + NXG*( j + NYG*NZ ) + nxt ) &
                      + f3( i + NXG*( j + NYG*NZ ) + nxt ) &
                      + f4( i + NXG*( j + NYG*NZ ) + nxt ) &
                      + f5( i + NXG*( j + NYG*NZ ) + nxt ) &
                      + f6( i + NXG*( j + NYG*NZ ) + nxt )
   END DO
   DO i = 1, NX
     g5( i + NXG*( j + NYG ) + nxt ) = g5( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
   DO i = 1, NX
     g11( i + NXG*( j + NYG ) + nxt ) = g11( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
   DO i = 1, NX
     g14( i + NXG*( j + NYG ) + nxt ) = g14( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
   DO i = 1, NX
     g15( i + NXG*( j + NYG ) + nxt ) = g15( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
   DO i = 1, NX
     g18( i + NXG*( j + NYG ) + nxt ) = g18( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
 END DO

! Edges along X
!$OMP DO
 DO i = 1, NX
   phi( i ) = f0( i + NXG*( NY + NYG*NZ )       ) &
            + f1( i + NXG*( NY + NYG*NZ ) + nxt ) &
            + f2( i + NXG*( NY + NYG*NZ ) + nxt ) &
            + f3( i + NXG*( NY + NYG*NZ ) + nxt ) &
            + f4( i + NXG*( NY + NYG*NZ ) + nxt ) &
            + f5( i + NXG*( NY + NYG*NZ ) + nxt ) &
            + f6( i + NXG*( NY + NYG*NZ ) + nxt )
   g15( i + NXG*( 1 + NYG ) + nxt ) = g15( i + NXG*( NY+1 + NYG*(NZ+1) ) + nxt )
   phi( i + NXG*(NY+1) ) = f0( i + NXG*( 1 + NYG*NZ )       ) &
                         + f1( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                         + f2( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                         + f3( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                         + f4( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                         + f5( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                         + f6( i + NXG*( 1 + NYG*NZ ) + nxt )
   g18( i + NXG*( NY + NYG ) + nxt ) = g18( i + NXG*NYG*(NZ+1) + nxt )
   phi( i + NXG*( NY+1 + NYG*(NZ+1) ) ) = f0( i + NXG*( 1 + NYG )       ) &
                                        + f1( i + NXG*( 1 + NYG ) + nxt ) &
                                        + f2( i + NXG*( 1 + NYG ) + nxt ) &
                                        + f3( i + NXG*( 1 + NYG ) + nxt ) &
                                        + f4( i + NXG*( 1 + NYG ) + nxt ) &
                                        + f5( i + NXG*( 1 + NYG ) + nxt ) &
                                        + f6( i + NXG*( 1 + NYG ) + nxt )
   g16( i + NXG*( NY + NYG*NZ ) + nxt ) = g16( i + nxt )
   phi( i + NXG*NYG*(NZ+1) ) = f0( i + NXG*( NY + NYG )       ) &
                             + f1( i + NXG*( NY + NYG ) + nxt ) &
                             + f2( i + NXG*( NY + NYG ) + nxt ) &
                             + f3( i + NXG*( NY + NYG ) + nxt ) &
                             + f4( i + NXG*( NY + NYG ) + nxt ) &
                             + f5( i + NXG*( NY + NYG ) + nxt ) &
                             + f6( i + NXG*( NY + NYG ) + nxt )
   g17( i + NXG*( 1 + NYG*NZ ) + nxt ) = g17( i + NXG*( NY+1 ) + nxt )
 END DO


! Edges along Y
!$OMP DO
 DO j = 1, NY
   phi( NXG*j ) = f0( NX + NXG*( j + NYG*NZ )       ) &
                + f1( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f2( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f3( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f4( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f5( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f6( NX + NXG*( j + NYG*NZ ) + nxt )
   g11( 1 + NXG*( j + NYG ) + nxt ) = g11( NX+1 + NXG*( j + NYG*(NZ+1) ) + nxt )
   phi( NX+1 + NXG*j ) = f0( 1 + NXG*( j + NYG*NZ )       ) &
                       + f1( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                       + f2( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                       + f3( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                       + f4( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                       + f5( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                       + f6( 1 + NXG*( j + NYG*NZ ) + nxt )
   g14( NX + NXG*( j + NYG ) + nxt ) = g14( NXG*( j + NYG*(NZ+1) ) + nxt )
   phi( NX + 1 + NXG*(j + NYG*(NZ+1) ) ) = f0( 1 + NXG*( j + NYG )       ) &
                                         + f1( 1 + NXG*( j + NYG ) + nxt ) &
                                         + f2( 1 + NXG*( j + NYG ) + nxt ) &
                                         + f3( 1 + NXG*( j + NYG ) + nxt ) &
                                         + f4( 1 + NXG*( j + NYG ) + nxt ) &
                                         + f5( 1 + NXG*( j + NYG ) + nxt ) &
                                         + f6( 1 + NXG*( j + NYG ) + nxt )
   g12( NX + NXG*( j + NYG*NZ ) + nxt ) = g12( NXG*j + nxt )
   phi( NXG*( j + NYG*(NZ+1) ) ) = f0( NX + NXG*( j + NYG )       ) &
                                 + f1( NX + NXG*( j + NYG ) + nxt ) &
                                 + f2( NX + NXG*( j + NYG ) + nxt ) &
                                 + f3( NX + NXG*( j + NYG ) + nxt ) &
                                 + f4( NX + NXG*( j + NYG ) + nxt ) &
                                 + f5( NX + NXG*( j + NYG ) + nxt ) &
                                 + f6( NX + NXG*( j + NYG ) + nxt )
   g13( 1 + NXG*( j + NYG*NZ ) + nxt ) = g13( NX + 1 + NXG*j + nxt ) 
 END DO

! Edges along Z
!$OMP DO
 DO k = 1, NZ
   phi( NXG*NYG*k ) = f0( NX + NXG*( NY + NYG*k )       ) &
                    + f1( NX + NXG*( NY + NYG*k ) + nxt ) &
                    + f2( NX + NXG*( NY + NYG*k ) + nxt ) &
                    + f3( NX + NXG*( NY + NYG*k ) + nxt ) &
                    + f4( NX + NXG*( NY + NYG*k ) + nxt ) &
                    + f5( NX + NXG*( NY + NYG*k ) + nxt ) &
                    + f6( NX + NXG*( NY + NYG*k ) + nxt )
   g7( 1 + NXG*( 1 + NYG*k ) + nxt ) = g7( NX+1 + NXG*( NY+1 + NYG*k ) + nxt )
   phi( NX+1 + NXG*NYG*k ) = f0( 1 + NXG*( NY + NYG*k )       ) &
                           + f1( 1 + NXG*( NY + NYG*k ) + nxt ) &
                           + f2( 1 + NXG*( NY + NYG*k ) + nxt ) &
                           + f3( 1 + NXG*( NY + NYG*k ) + nxt ) &
                           + f4( 1 + NXG*( NY + NYG*k ) + nxt ) &
                           + f5( 1 + NXG*( NY + NYG*k ) + nxt ) &
                           + f6( 1 + NXG*( NY + NYG*k ) + nxt )
   g10( NX + NXG*( 1 + NYG*k ) + nxt ) = g10( NXG*( NY+1 + NYG*k ) + nxt )
   phi( NX+1 + NXG*( NY+1 + NYG*k ) ) = f0( 1 + NXG*( 1 +NYG*k )       ) &
                                      + f1( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                                      + f2( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                                      + f3( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                                      + f4( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                                      + f5( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                                      + f6( 1 + NXG*( 1 +NYG*k ) + nxt )
   g8( NX + NXG*( NY + NYG*k ) + nxt ) = g8( NXG*NYG*k + nxt )
   phi( NXG*( NY+1 + NYG*k ) ) = f0( NX + NXG*( 1 + NYG*k )       ) &
                               + f1( NX + NXG*( 1 + NYG*k ) + nxt ) &
                               + f2( NX + NXG*( 1 + NYG*k ) + nxt ) &
                               + f3( NX + NXG*( 1 + NYG*k ) + nxt ) &
                               + f4( NX + NXG*( 1 + NYG*k ) + nxt ) &
                               + f5( NX + NXG*( 1 + NYG*k ) + nxt ) &
                               + f6( NX + NXG*( 1 + NYG*k ) + nxt )
    g9( 1 + NXG*( NY + NYG*k ) + nxt ) = g9( NX+1 + NXG*NYG*k + nxt )
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE PostStream

