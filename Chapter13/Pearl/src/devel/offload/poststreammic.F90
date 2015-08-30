!-------------------------------------------------------------------------------
! Subroutine : PostStreamMIC
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

 SUBROUTINE PostStreamMIC
!DIR$ ATTRIBUTES OFFLOAD:mic :: PostStreamMIC

!  Common Variables
 USE Domain
 USE FluidParams, ONLY : phi_mic
 USE LBMParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, offset

  offset = NX*NZ

!$OMP PARALLEL

!-------- Exchange data for phi_mic and inward g components ------------------------
! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY_MIC
     phi_mic( NX+1 + NXG*( j + NYG_MIC*k ) ) = f0_mic( 1 + NXG*( j + NYG_MIC*k )       ) &
                                             + f1_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                             + f2_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                             + f3_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                             + f4_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                             + f5_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                             + f6_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g2_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g2_mic( NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g8_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g8_mic( NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g10_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g10_mic( NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g12_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g12_mic( NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g14_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g14_mic( NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
 END DO

!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY_MIC
     phi_mic( NXG*( j + NYG_MIC*k ) ) = f0_mic( NX + NXG*( j + NYG_MIC*k )       ) &
                                      + f1_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                      + f2_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                      + f3_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                      + f4_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                      + f5_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) &
                                      + f6_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g1_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g1_mic( NX+1 + NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g7_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g7_mic( NX+1 + NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g9_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g9_mic( NX+1 + NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g11_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g11_mic( NX+1 + NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY_MIC
     g13_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic ) = g13_mic( NX+1 + NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
 END DO

! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     phi_mic( i + NXG*( NY_MIC+1 + NYG_MIC*k ) ) = buff_cpu_phi( i + NX*(k-1) )
   END DO
   DO i = 1, NX
     g4_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + offset )
   END DO
   DO i = 1, NX
     g8_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + 2*offset )
   END DO
   DO i = 1, NX
     g9_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + 3*offset )
   END DO
   DO i = 1, NX
     g16_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + 4*offset )
   END DO
   DO i = 1, NX
     g18_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + 5*offset )
   END DO
 END DO

!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     phi_mic( i + NXG*( NYG_MIC*k ) ) = buff_cpu_phi( i + NX*(k-1) + offset )
   END DO
   DO i = 1, NX
     g3_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + 7*offset )
   END DO
   DO i = 1, NX
     g7_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + 8*offset )
   END DO
   DO i = 1, NX
     g10_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + 9*offset )
   END DO
   DO i = 1, NX
     g15_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + 10*offset )
   END DO
   DO i = 1, NX
     g17_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) + 11*offset )
   END DO
 END DO


! Z direction
!$OMP DO PRIVATE(i)
 DO j = 1, NY_MIC
   DO i = 1, NX
     phi_mic( i + NXG*( j + NYG_MIC*(NZ+1) ) ) = f0_mic( i + NXG*( j + NYG_MIC )       ) &
                                               + f1_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f2_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f3_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f4_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f5_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f6_mic( i + NXG*( j + NYG_MIC ) + nxt_mic )
   END DO
   DO i = 1, NX
     g6_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) = g6_mic( i + NXG*j + nxt_mic )
   END DO
   DO i = 1, NX
     g12_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) = g12_mic( i + NXG*j + nxt_mic )
   END DO
   DO i = 1, NX
     g13_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) = g13_mic( i + NXG*j + nxt_mic )
   END DO
   DO i = 1, NX
     g16_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) = g16_mic( i + NXG*j + nxt_mic )
   END DO
   DO i = 1, NX
     g17_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) = g17_mic( i + NXG*j + nxt_mic )
   END DO
 END DO

!$OMP DO PRIVATE(i)
 DO j = 1, NY_MIC
   DO i = 1, NX
     phi_mic( i + NXG*j ) = f0_mic( i + NXG*( j + NYG_MIC*NZ )       ) &
                          + f1_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                          + f2_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                          + f3_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                          + f4_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                          + f5_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                          + f6_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic )
   END DO
   DO i = 1, NX
     g5_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) = g5_mic( i + NXG*( j + NYG_MIC*(NZ+1) ) + nxt_mic )
   END DO
   DO i = 1, NX
     g11_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) = g11_mic( i + NXG*( j + NYG_MIC*(NZ+1) ) + nxt_mic )
   END DO
   DO i = 1, NX
     g14_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) = g14_mic( i + NXG*( j + NYG_MIC*(NZ+1) ) + nxt_mic )
   END DO
   DO i = 1, NX
     g15_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) = g15_mic( i + NXG*( j + NYG_MIC*(NZ+1) ) + nxt_mic )
   END DO
   DO i = 1, NX
     g18_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) = g18_mic( i + NXG*( j + NYG_MIC*(NZ+1) ) + nxt_mic )
   END DO
 END DO

! Edges along X
!$OMP DO
 DO i = 1, NX
   phi_mic( i + NXG*(NY_MIC+1 + NYG_MIC*(NZ+1) ) ) = buff_cpu_phi_edge_x( i )        ! From (NY_MIC+1,1) on CPU
   phi_mic( i + NXG*(           NYG_MIC*(NZ+1) ) ) = buff_cpu_phi_edge_x( i + NX )   ! From (NY,1) on CPU
   phi_mic( i + NXG*(NY_MIC+1)                   ) = buff_cpu_phi_edge_x( i + 2*NX ) ! From (NY_MIC+1,NZ) on CPU
   phi_mic( i                                    ) = buff_cpu_phi_edge_x( i + 3*NX ) ! From (NY,NZ) on CPU

   g15_mic( i + NXG*( 1      + NYG_MIC    ) + nxt_mic ) = buff_cpu_edge_x( i )        ! From (NY+1,NZ+1) on the CPU
   g16_mic( i + NXG*( NY_MIC + NYG_MIC*NZ ) + nxt_mic ) = buff_cpu_edge_x( i + NX )   ! From (NY_MIC,0) on the CPU
   g17_mic( i + NXG*( 1      + NYG_MIC*NZ ) + nxt_mic ) = buff_cpu_edge_x( i + 2*NX ) ! From (NY+1,0) on the CPU
   g18_mic( i + NXG*( NY_MIC + NYG_MIC    ) + nxt_mic ) = buff_cpu_edge_x( i + 3*NX ) ! From (NY_MIC,NZ+1) on the CPU
 END DO


! Edges along Y
!$OMP DO
 DO j = 1, NY_MIC
   phi_mic( NXG*j ) = f0_mic( NX + NXG*( j + NYG_MIC*NZ )       ) &
                    + f1_mic( NX + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                    + f2_mic( NX + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                    + f3_mic( NX + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                    + f4_mic( NX + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                    + f5_mic( NX + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                    + f6_mic( NX + NXG*( j + NYG_MIC*NZ ) + nxt_mic )
   g11_mic( 1 + NXG*( j + NYG_MIC ) + nxt_mic ) = g11_mic( NX+1 + NXG*( j + NYG_MIC*(NZ+1) ) + nxt_mic )
   phi_mic( NX+1 + NXG*j ) = f0_mic( 1 + NXG*( j + NYG_MIC*NZ )       ) &
                           + f1_mic( 1 + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                           + f2_mic( 1 + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                           + f3_mic( 1 + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                           + f4_mic( 1 + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                           + f5_mic( 1 + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) &
                           + f6_mic( 1 + NXG*( j + NYG_MIC*NZ ) + nxt_mic )
   g14_mic( NX + NXG*( j + NYG_MIC ) + nxt_mic ) = g14_mic( NXG*( j + NYG_MIC*(NZ+1) ) + nxt_mic )
   phi_mic( NX+1 + NXG*(j + NYG_MIC*(NZ+1) ) ) = f0_mic( 1 + NXG*( j + NYG_MIC )       ) &
                                               + f1_mic( 1 + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f2_mic( 1 + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f3_mic( 1 + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f4_mic( 1 + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f5_mic( 1 + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                               + f6_mic( 1 + NXG*( j + NYG_MIC ) + nxt_mic )
   g12_mic( NX + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) = g12_mic( NXG*j + nxt_mic )
   phi_mic( NXG*( j + NYG_MIC*(NZ+1) ) ) = f0_mic( NX + NXG*( j + NYG_MIC )       ) &
                                         + f1_mic( NX + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                         + f2_mic( NX + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                         + f3_mic( NX + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                         + f4_mic( NX + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                         + f5_mic( NX + NXG*( j + NYG_MIC ) + nxt_mic ) &
                                         + f6_mic( NX + NXG*( j + NYG_MIC ) + nxt_mic )
   g13_mic( 1 + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) = g13_mic( NX+1 + NXG*j + nxt_mic ) 
 END DO

! Edges along Z
!$OMP DO
 DO k = 1, NZ
   phi_mic( NX+1 + NXG*( NY_MIC+1 + NYG_MIC*k ) ) = buff_cpu_phi_edge_z( k )         ! From (1,NY_MIC+1) on the CPU
   phi_mic(        NXG*( NY_MIC+1 + NYG_MIC*k ) ) = buff_cpu_phi_edge_z( k + NZ )    ! From (NX,NY_MIC+1) on the CPU
   phi_mic( NX+1 + NXG*(            NYG_MIC*k ) ) = buff_cpu_phi_edge_z( k + 2*NZ )  ! From (1,NY) on the CPU
   phi_mic(        NXG*(            NYG_MIC*k ) ) = buff_cpu_phi_edge_z( k + 3*NZ )  ! From (NX,NY) on the CPU


   g7_mic( 1   + NXG*( 1      + NYG_MIC*k ) + nxt_mic ) = buff_cpu_edge_z( k )        ! From (NX+1,NY+1) on the CPU
   g8_mic( NX  + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) = buff_cpu_edge_z( k + NZ )   ! From (0,NY_MIC) on the CPU
   g9_mic( 1   + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) = buff_cpu_edge_z( k + 2*NZ ) ! From (NX+1,NY_MIC) on the CPU
   g10_mic( NX + NXG*( 1      + NYG_MIC*k ) + nxt_mic ) = buff_cpu_edge_z( k + 3*NZ ) ! From (0,NY+1) on the CPU
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE PostStreamMIC

