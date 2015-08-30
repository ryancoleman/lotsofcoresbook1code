!-------------------------------------------------------------------------------
! Subroutine : UpdatePhi
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Update f in nodes at processor boundaries.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
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

 SUBROUTINE UpdatePhi

!  Common Variables
 USE Domain
 USE FluidParams, ONLY : phi
 USE LBMParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, offset

 offset = NX*NZ

!$OMP PARALLEL

! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = NY_MIC+1, NY
     f2( NX + NXG*( j + NYG*k ) + nxt ) = f2( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     f1( 1 + NXG*( j + NYG*k ) + nxt ) = f1( NX+1 + NYG*( j + NZG*k ) + nxt )
   END DO
 END DO

! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     f4( i + NXG*( NY + NYG*k ) + nxt )   = buff_mic( i + NX*(k-1))
   END DO
   DO i = 1, NX
     f3( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 6*offset  )
   END DO
 END DO

! Z direction
!$OMP DO PRIVATE(i)
 DO j = NY+1, NY
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
! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ

   DO i = 1, NX
     buff_cpu_phi( i + NX*(k-1) ) = f0( i + NXG*( NY_MIC+1 + NYG*k )       ) &
                             + f1( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                             + f2( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                             + f3( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                             + f4( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                             + f5( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                             + f6( i + NXG*( NY_MIC+1 + NYG*k ) + nxt )
   END DO

   DO i = 1, NX
     buff_cpu_phi( i + NX*(k-1) + offset ) = f0( i + NXG*( NY + NYG*k )       ) &
                                      + f1( i + NXG*( NY + NYG*k ) + nxt ) &
                                      + f2( i + NXG*( NY + NYG*k ) + nxt ) &
                                      + f3( i + NXG*( NY + NYG*k ) + nxt ) &
                                      + f4( i + NXG*( NY + NYG*k ) + nxt ) &
                                      + f5( i + NXG*( NY + NYG*k ) + nxt ) &
                                      + f6( i + NXG*( NY + NYG*k ) + nxt )
   END DO

 END DO




! Edges along X
!$OMP DO
 DO i = 1, NX

   buff_cpu_phi_edge_x( i )   = f0( i + NXG*( NY_MIC+1 + NYG )           ) &
                             + f1( i + NXG*( NY_MIC+1 + NYG ) + nxt ) &
                             + f2( i + NXG*( NY_MIC+1 + NYG ) + nxt ) &
                             + f3( i + NXG*( NY_MIC+1 + NYG ) + nxt ) &
                             + f4( i + NXG*( NY_MIC+1 + NYG ) + nxt ) &
                             + f5( i + NXG*( NY_MIC+1 + NYG ) + nxt ) &
                             + f6( i + NXG*( NY_MIC+1 + NYG ) + nxt )

   buff_cpu_phi_edge_x( i + NX ) = f0( i + NXG*( NY + NYG )           ) &
                                + f1( i + NXG*( NY + NYG ) + nxt ) &
                                + f2( i + NXG*( NY + NYG ) + nxt ) &
                                + f3( i + NXG*( NY + NYG ) + nxt ) &
                                + f4( i + NXG*( NY + NYG ) + nxt ) &
                                + f5( i + NXG*( NY + NYG ) + nxt ) &
                                + f6( i + NXG*( NY + NYG ) + nxt )

   buff_cpu_phi_edge_x( i + 2*NX ) = f0( i + NXG*( NY_MIC+1 + NYG*NZ )           ) &
                                  + f1( i + NXG*( NY_MIC+1 + NYG*NZ ) + nxt ) &
                                  + f2( i + NXG*( NY_MIC+1 + NYG*NZ ) + nxt ) &
                                  + f3( i + NXG*( NY_MIC+1 + NYG*NZ ) + nxt ) &
                                  + f4( i + NXG*( NY_MIC+1 + NYG*NZ ) + nxt ) &
                                  + f5( i + NXG*( NY_MIC+1 + NYG*NZ ) + nxt ) &
                                  + f6( i + NXG*( NY_MIC+1 + NYG*NZ ) + nxt )

   buff_cpu_phi_edge_x( i + 3*NX ) = f0( i + NXG*( NY + NYG*NZ )           ) &
                                  + f1( i + NXG*( NY + NYG*NZ ) + nxt ) &
                                  + f2( i + NXG*( NY + NYG*NZ ) + nxt ) &
                                  + f3( i + NXG*( NY + NYG*NZ ) + nxt ) &
                                  + f4( i + NXG*( NY + NYG*NZ ) + nxt ) &
                                  + f5( i + NXG*( NY + NYG*NZ ) + nxt ) &
                                  + f6( i + NXG*( NY + NYG*NZ ) + nxt )


 END DO

! Edges along Z
!$OMP DO
 DO k = 1, NZ

   buff_cpu_phi_edge_z( k ) = f0( 1 + NXG*( NY_MIC+1 + NYG*k )           ) &
                           + f1( 1 + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                           + f2( 1 + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                           + f3( 1 + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                           + f4( 1 + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                           + f5( 1 + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                           + f6( 1 + NXG*( NY_MIC+1 + NYG*k ) + nxt )

   buff_cpu_phi_edge_z( k + NZ ) = f0( NX + NXG*( NY_MIC+1 + NYG*k )           ) &
                                + f1( NX + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                                + f2( NX + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                                + f3( NX + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                                + f4( NX + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                                + f5( NX + NXG*( NY_MIC+1 + NYG*k ) + nxt ) &
                                + f6( NX + NXG*( NY_MIC+1 + NYG*k ) + nxt )


   buff_cpu_phi_edge_z( k + 2*NZ ) = f0( 1 + NXG*( NY + NYG*k )           ) &
                                  + f1( 1 + NXG*( NY + NYG*k ) + nxt ) &
                                  + f2( 1 + NXG*( NY + NYG*k ) + nxt ) &
                                  + f3( 1 + NXG*( NY + NYG*k ) + nxt ) &
                                  + f4( 1 + NXG*( NY + NYG*k ) + nxt ) &
                                  + f5( 1 + NXG*( NY + NYG*k ) + nxt ) &
                                  + f6( 1 + NXG*( NY + NYG*k ) + nxt )

   buff_cpu_phi_edge_z( k + 3*NZ ) = f0( NX + NXG*( NY + NYG*k )           ) &
                                  + f1( NX + NXG*( NY + NYG*k ) + nxt ) &
                                  + f2( NX + NXG*( NY + NYG*k ) + nxt ) &
                                  + f3( NX + NXG*( NY + NYG*k ) + nxt ) &
                                  + f4( NX + NXG*( NY + NYG*k ) + nxt ) &
                                  + f5( NX + NXG*( NY + NYG*k ) + nxt ) &
                                  + f6( NX + NXG*( NY + NYG*k ) + nxt )
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE UpdatePhi

