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
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, offset

 offset = NX*NZ

!$OMP PARALLEL

!-------- Exchange data for phi and inward g components ------------------------
! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = NY_MIC+1, NY
     phi( NX+1 + NXG*( j + NYG*k ) ) = f0( 1 + NXG*( j + NYG*k )       ) &
                                     + f1( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f2( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f3( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f4( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f5( 1 + NXG*( j + NYG*k ) + nxt ) &
                                     + f6( 1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g2( NX + NXG*( j + NYG*k ) + nxt ) = g2( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g8( NX + NXG*( j + NYG*k ) + nxt ) = g8( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g10( NX + NXG*( j + NYG*k ) + nxt ) = g10( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g12( NX + NXG*( j + NYG*k ) + nxt ) = g12( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g14( NX + NXG*( j + NYG*k ) + nxt ) = g14( NXG*( j + NYG*k ) + nxt )
   END DO
 END DO

!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = NY_MIC+1, NY
     phi( NXG*( j + NYG*k ) ) = f0( NX + NXG*( j + NYG*k )       ) &
                              + f1( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f2( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f3( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f4( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f5( NX + NXG*( j + NYG*k ) + nxt ) &
                              + f6( NX + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g1( 1 + NXG*( j + NYG*k ) + nxt ) = g1( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g7( 1 + NXG*( j + NYG*k ) + nxt ) = g7( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g9( 1 + NXG*( j + NYG*k ) + nxt ) = g9( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g11( 1 + NXG*( j + NYG*k ) + nxt ) = g11( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = NY_MIC+1, NY
     g13( 1 + NXG*( j + NYG*k ) + nxt ) = g13( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
 END DO

! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ

   ! From j=1 on the MIC
   DO i = 1, NX
     phi( i + NXG*( NY+1 + NYG*k ) ) = buff_mic_phi( i + NX*(k-1) )
   END DO

   ! From j = 0 on the MIC
   DO i = 1, NX
     g4( i + NXG*( NY + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + offset )
   END DO
   DO i = 1, NX
     g8( i + NXG*( NY + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 2*offset )
   END DO
   DO i = 1, NX
     g9( i + NXG*( NY + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 3*offset )
   END DO
   DO i = 1, NX
     g16( i + NXG*( NY + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 4*offset )
   END DO
   DO i = 1, NX
     g18( i + NXG*( NY + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 5*offset )
   END DO

 END DO

!$OMP DO PRIVATE(i)
 DO k = 1, NZ

   ! From j = NY_MIC on the MIC
   DO i = 1, NX
     phi( i + NXG*( NY_MIC + NYG*k ) ) = buff_mic_phi( i + NX*(k-1) + offset )
   END DO

   ! From j = NY_MIC+1 on the MIC
   DO i = 1, NX
     g3( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 7*offset )
   END DO
   DO i = 1, NX
     g7( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 8*offset )
   END DO
   DO i = 1, NX
     g10( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 9*offset )
   END DO
   DO i = 1, NX
     g15( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 10*offset )
   END DO
   DO i = 1, NX
     g17( i + NXG*( NY_MIC+1 + NYG*k ) + nxt ) = buff_mic( i + NX*(k-1) + 11*offset )
   END DO
 END DO


! Z direction
!$OMP DO PRIVATE(i)
 DO j = NY_MIC+1, NY
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
 DO j = NY_MIC+1, NY
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
   phi( i + NXG*NY_MIC )                  = buff_mic_phi_edge_x( i + 3*NX ) ! From (NY_MIC,NZ) on the MIC
   phi( i + NXG*(NY+1) )                  = buff_mic_phi_edge_x( i + 2*NX ) ! From (1,NZ)      on the MIC
   phi( i + NXG*( NY_MIC + NYG*(NZ+1) ) ) = buff_mic_phi_edge_x( i +   NX ) ! From (NY_MIC,1)  on the MIC
   phi( i + NXG*( (NY+1) + NYG*(NZ+1) ) ) = buff_mic_phi_edge_x( i        ) ! From (1,1)       on the MIC


   g15( i + NXG*( NY_MIC+1 + NYG ) + nxt )    = buff_mic_edge_x(i       ) ! From (NY_MIC+1,NZ+1) on the MIC
   g16( i + NXG*( NY + NYG*NZ ) + nxt )       = buff_mic_edge_x(i +   NX) ! From (0,0)           on the MIC
   g17( i + NXG*( NY_MIC+1 + NYG*NZ ) + nxt ) = buff_mic_edge_x(i + 2*NX) ! From (NY_MIC+1,0)    on the MIC
   g18( i + NXG*( NY + NYG ) + nxt )          = buff_mic_edge_x(i + 3*NX) ! From (0,NZ+1)        on the MIC
 END DO


! Edges along Y
!$OMP DO
 DO j = NY_MIC+1, NY
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
   phi(        NXG*( NY_MIC + NYG*k ) ) = buff_mic_phi_edge_z( k + 3*NZ ) ! From (NX,NY_MIC) on the MIC
   phi( NX+1 + NXG*( NY_MIC + NYG*k ) ) = buff_mic_phi_edge_z( k + 2*NZ ) ! From (1,NY_MIC)  on the MIC
   phi(        NXG*( NY+1   + NYG*k ) ) = buff_mic_phi_edge_z( k +   NZ ) ! From (NX,1)      on the MIC
   phi( NX+1 + NXG*( NY+1   + NYG*k ) ) = buff_mic_phi_edge_z( k        ) ! From (1,1)       on the MIC

   g7( 1   + NXG*( NY_MIC+1 + NYG*k ) + nxt ) = buff_mic_edge_z(k       ) ! From (NX+1,NY_MIC+1) on the MIC
   g8( NX  + NXG*( NY       + NYG*k ) + nxt ) = buff_mic_edge_z(k +   NZ) ! From (0,0)           on the MIC
   g9( 1   + NXG*( NY       + NYG*k ) + nxt ) = buff_mic_edge_z(k + 2*NZ) ! From (NX+1,0)        on the MIC
   g10( NX + NXG*( NY_MIC+1 + NYG*k ) + nxt ) = buff_mic_edge_z(k + 3*NZ) ! From (0,NY_MIC+1)    on the MIC

 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE PostStream

