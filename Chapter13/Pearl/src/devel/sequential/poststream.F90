!-------------------------------------------------------------------------------
! Subroutine : PostStream
! Revision   : 1.3 (2013/11/12)
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


! X direction
 DO k = 1, NZ
   DO j = 1, NY
     f(NX,j,k,2,nxt ) = f(0,j,k,2,nxt )
   END DO
   DO j = 1, NY
     f(1,j,k,1,nxt ) = f(NX+1,j,k,1,nxt )
   END DO
 END DO

! Y direction
 DO k = 1, NZ
   DO i = 1, NX
     f(i,NY,k,4,nxt ) = f(i,0,k,4,nxt )
   END DO
   DO i = 1, NX
     f(i,1,k,3,nxt ) = f(i,NY+1,k,3,nxt )
   END DO
 END DO

! Z direction
 DO j = 1, NY
    DO i = 1, NX
     f(i,j,NZ,6,nxt ) = f(i,j,0,6,nxt )
   END DO
   DO i = 1, NX
     f(i,j,1,5,nxt ) = f(i,j,NZ+1,5,nxt )
   END DO
 END DO

!
! At this stage all values of f needed for the next step have been updated
!


!-------- Exchange data for phi and inward g components ------------------------
! X direction
 DO k = 1, NZ
   DO j = 1, NY
     phi( NX+1,j,k ) = f( 1,j,k,0,nxt ) &
                                     + f( 1,j,k,1,nxt ) &
                                     + f( 1,j,k,2,nxt ) &
                                     + f( 1,j,k,3,nxt ) &
                                     + f( 1,j,k,4,nxt ) &
                                     + f( 1,j,k,5,nxt ) &
                                     + f( 1,j,k,6,nxt )
   END DO
   DO j = 1, NY
     g( NX,j,k,2,nxt ) = g( 0,j,k,2,nxt )
   END DO
   DO j = 1, NY
     g( NX,j,k,8,nxt ) = g( 0,j,k,8,nxt )
   END DO
   DO j = 1, NY
     g( NX,j,k,10,nxt ) = g( 0,j,k,10,nxt )
   END DO
   DO j = 1, NY
     g( NX,j,k,12,nxt ) = g( 0,j,k,12,nxt )
   END DO
   DO j = 1, NY
     g( NX,j,k,14,nxt ) = g( 0,j,k,14,nxt )
   END DO
 END DO

 DO k = 1, NZ
   DO j = 1, NY
     phi( 0,j,k ) = f( NX,j,k,0,nxt ) &
                              + f( NX,j,k,1,nxt ) &
                              + f( NX,j,k,2,nxt ) &
                              + f( NX,j,k,3,nxt ) &
                              + f( NX,j,k,4,nxt ) &
                              + f( NX,j,k,5,nxt ) &
                              + f( NX,j,k,6,nxt )
   END DO
   DO j = 1, NY
     g( 1,j,k,1,nxt ) = g( NX+1,j,k,1,nxt )
   END DO
   DO j = 1, NY
     g( 1,j,k,7,nxt ) = g( NX+1,j,k,7,nxt )
   END DO
   DO j = 1, NY
     g( 1,j,k,9,nxt ) = g( NX+1,j,k,9,nxt )
   END DO
   DO j = 1, NY
     g( 1,j,k,11,nxt ) = g( NX+1,j,k,11,nxt )
   END DO
   DO j = 1, NY
     g( 1,j,k,13,nxt ) = g( NX+1,j,k,13,nxt )
   END DO
 END DO

! Y direction
 DO k = 1, NZ
   DO i = 1, NX
     phi( i,NY+1,k ) = f( i,1,k,0,nxt ) &
                                     + f( i,1,k,1,nxt ) &
                                     + f( i,1,k,2,nxt ) &
                                     + f( i,1,k,3,nxt ) &
                                     + f( i,1,k,4,nxt ) &
                                     + f( i,1,k,5,nxt ) &
                                     + f( i,1,k,6,nxt )
   END DO
   DO i = 1, NX
     g( i,NY,k,4,nxt ) = g( i,0,k,4,nxt )
   END DO
   DO i = 1, NX
     g( i,NY,k,8,nxt ) = g( i,0,k,8,nxt )
   END DO
   DO i = 1, NX
     g( i,NY,k,9,nxt ) = g( i,0,k,9,nxt )
   END DO
   DO i = 1, NX
     g( i,NY,k,16,nxt ) = g( i,0,k,16,nxt )
   END DO
   DO i = 1, NX
     g( i,NY,k,18,nxt ) = g( i,0,k,18,nxt )
   END DO
 END DO

 DO k = 1, NZ
   DO i = 1, NX
     phi( i,0,k ) = f( i,NY,k,0,nxt ) &
                              + f( i,NY,k,1,nxt ) &
                              + f( i,NY,k,2,nxt ) &
                              + f( i,NY,k,3,nxt ) &
                              + f( i,NY,k,4,nxt ) &
                              + f( i,NY,k,5,nxt ) &
                              + f( i,NY,k,6,nxt )
   END DO
   DO i = 1, NX
     g( i,1,k,3,nxt ) = g( i,NY+1,k,3,nxt )
   END DO
   DO i = 1, NX
     g( i,1,k,7,nxt ) = g( i,NY+1,k,7,nxt )
   END DO
   DO i = 1, NX
     g( i,1,k,10,nxt ) = g( i,NY+1,k,10,nxt )
   END DO
   DO i = 1, NX
     g( i,1,k,15,nxt ) = g( i,NY+1,k,15,nxt )
   END DO
   DO i = 1, NX
     g( i,1,k,17,nxt ) = g( i,NY+1,k,17,nxt )
   END DO
 END DO


! Z direction
 DO j = 1, NY
   DO i = 1, NX
     phi( i,j,NZ+1 ) = f( i,j,1,0,nxt ) &
                                       + f( i,j,1,1,nxt ) &
                                       + f( i,j,1,2,nxt ) &
                                       + f( i,j,1,3,nxt ) &
                                       + f( i,j,1,4,nxt ) &
                                       + f( i,j,1,5,nxt ) &
                                       + f( i,j,1,6,nxt )
   END DO
   DO i = 1, NX
     g( i,j,NZ,6,nxt ) = g( i,j,0,6,nxt )
   END DO
   DO i = 1, NX
     g( i,j,NZ,12,nxt ) = g( i,j,0,12,nxt )
   END DO
   DO i = 1, NX
     g( i,j,NZ,13,nxt ) = g( i,j,0,13,nxt )
   END DO
   DO i = 1, NX
     g( i,j,NZ,16,nxt ) = g( i,j,0,16,nxt )
   END DO
   DO i = 1, NX
     g( i,j,NZ,17,nxt ) = g( i,j,0,17,nxt )
   END DO
 END DO

 DO j = 1, NY
   DO i = 1, NX
     phi( i,j,0 ) = f( i,j,NZ,0,nxt ) &
                      + f( i,j,NZ,1,nxt ) &
                      + f( i,j,NZ,2,nxt ) &
                      + f( i,j,NZ,3,nxt ) &
                      + f( i,j,NZ,4,nxt ) &
                      + f( i,j,NZ,5,nxt ) &
                      + f( i,j,NZ,6,nxt )
   END DO
   DO i = 1, NX
     g( i,j,1,5,nxt ) = g( i,j,NZ+1,5,nxt )
   END DO
   DO i = 1, NX
     g( i,j,1,11,nxt ) = g( i,j,NZ+1,11,nxt )
   END DO
   DO i = 1, NX
     g( i,j,1,14,nxt ) = g( i,j,NZ+1,14,nxt )
   END DO
   DO i = 1, NX
     g( i,j,1,15,nxt ) = g( i,j,NZ+1,15,nxt )
   END DO
   DO i = 1, NX
     g( i,j,1,18,nxt ) = g( i,j,NZ+1,18,nxt )
   END DO
 END DO

! Edges along X
 DO i = 1, NX
   phi( i,0,0 ) = f( i,NY,NZ,0,nxt ) &
            + f( i,NY,NZ,1,nxt ) &
            + f( i,NY,NZ,2,nxt ) &
            + f( i,NY,NZ,3,nxt ) &
            + f( i,NY,NZ,4,nxt ) &
            + f( i,NY,NZ,5,nxt ) &
            + f( i,NY,NZ,6,nxt )
   g( i,1,1,15,nxt ) = g( i,NY+1,NZ+1,15,nxt )
   phi( i,NY+1,0 ) = f( i,1,NZ,0,nxt ) &
                         + f( i,1,NZ,1,nxt ) &
                         + f( i,1,NZ,2,nxt ) &
                         + f( i,1,NZ,3,nxt ) &
                         + f( i,1,NZ,4,nxt ) &
                         + f( i,1,NZ,5,nxt ) &
                         + f( i,1,NZ,6,nxt )
   g( i,NY,1,18,nxt ) = g( i,0,NZ+1,18,nxt )
   phi( i,NY+1,NZ+1 ) = f( i,1,1,0,nxt ) &
                                        + f( i,1,1,1,nxt ) &
                                        + f( i,1,1,2,nxt ) &
                                        + f( i,1,1,3,nxt ) &
                                        + f( i,1,1,4,nxt ) &
                                        + f( i,1,1,5,nxt ) &
                                        + f( i,1,1,6,nxt )
   g( i,NY,NZ,16,nxt ) = g( i,0,0,16,nxt )
   phi( i,0,NZ+1 ) = f( i,NY,1,0,nxt ) &
                             + f( i,NY,1,1,nxt ) &
                             + f( i,NY,1,2,nxt ) &
                             + f( i,NY,1,3,nxt ) &
                             + f( i,NY,1,4,nxt ) &
                             + f( i,NY,1,5,nxt ) &
                             + f( i,NY,1,6,nxt )
   g( i,1,NZ,17,nxt ) = g( i,NY+1,0,17,nxt )
 END DO


! Edges along Y
 DO j = 1, NY
   phi( 0,j,0 ) = f( NX,j,NZ,0,nxt ) &
                + f( NX,j,NZ,1,nxt ) &
                + f( NX,j,NZ,2,nxt ) &
                + f( NX,j,NZ,3,nxt ) &
                + f( NX,j,NZ,4,nxt ) &
                + f( NX,j,NZ,5,nxt ) &
                + f( NX,j,NZ,6,nxt )
   g( 1,j,1,11,nxt ) = g( NX+1,j,NZ+1,11,nxt )
   phi( NX+1,j,0 ) = f( 1,j,NZ,0,nxt ) &
                       + f( 1,j,NZ,1,nxt ) &
                       + f( 1,j,NZ,2,nxt ) &
                       + f( 1,j,NZ,3,nxt ) &
                       + f( 1,j,NZ,4,nxt ) &
                       + f( 1,j,NZ,5,nxt ) &
                       + f( 1,j,NZ,6,nxt )
   g( NX,j,1,14,nxt ) = g( 0,j,NZ+1,14,nxt )
   phi( NX+1,j,NZ+1 ) = f( 1,j,1,0,nxt ) &
                                         + f( 1,j,1,1,nxt ) &
                                         + f( 1,j,1,2,nxt ) &
                                         + f( 1,j,1,3,nxt ) &
                                         + f( 1,j,1,4,nxt ) &
                                         + f( 1,j,1,5,nxt ) &
                                         + f( 1,j,1,6,nxt )
   g( NX,j,NZ,12,nxt ) = g( 0,j,0,12,nxt )
   phi( 0,j,NZ+1 ) = f( NX,j,1,0,nxt ) &
                                 + f( NX,j,1,1,nxt ) &
                                 + f( NX,j,1,2,nxt ) &
                                 + f( NX,j,1,3,nxt ) &
                                 + f( NX,j,1,4,nxt ) &
                                 + f( NX,j,1,5,nxt ) &
                                 + f( NX,j,1,6,nxt )
   g( 1,j,NZ,13,nxt ) = g( NX+1,j,0,13,nxt ) 
 END DO

! Edges along Z
 DO k = 1, NZ
   phi( 0,0,k ) = f( NX,NY,k,0,nxt ) &
                    + f( NX,NY,k,1,nxt ) &
                    + f( NX,NY,k,2,nxt ) &
                    + f( NX,NY,k,3,nxt ) &
                    + f( NX,NY,k,4,nxt ) &
                    + f( NX,NY,k,5,nxt ) &
                    + f( NX,NY,k,6,nxt )
   g( 1,1,k,7,nxt ) = g( NX+1,NY+1,k,7,nxt )
   phi( NX+1,0,k ) = f( 1,NY,k,0,nxt ) &
                           + f( 1,NY,k,1,nxt ) &
                           + f( 1,NY,k,2,nxt ) &
                           + f( 1,NY,k,3,nxt ) &
                           + f( 1,NY,k,4,nxt ) &
                           + f( 1,NY,k,5,nxt ) &
                           + f( 1,NY,k,6,nxt )
   g( NX,1,k,10,nxt ) = g( 0,NY+1,k,10,nxt )
   phi( NX+1,NY+1,k ) = f( 1,1,k,0,nxt ) &
                                      + f( 1,1,k,1,nxt ) &
                                      + f( 1,1,k,2,nxt ) &
                                      + f( 1,1,k,3,nxt ) &
                                      + f( 1,1,k,4,nxt ) &
                                      + f( 1,1,k,5,nxt ) &
                                      + f( 1,1,k,6,nxt )
   g( NX,NY,k,8,nxt ) = g( 0,0,k,8,nxt )
   phi( 0,NY+1,k ) = f( NX,1,k,0,nxt ) &
                               + f( NX,1,k,1,nxt ) &
                               + f( NX,1,k,2,nxt ) &
                               + f( NX,1,k,3,nxt ) &
                               + f( NX,1,k,4,nxt ) &
                               + f( NX,1,k,5,nxt ) &
                               + f( NX,1,k,6,nxt )
    g( 1,NY,k,9,nxt ) = g( NX+1,0,k,9,nxt )
 END DO

 RETURN
 END SUBROUTINE PostStream

