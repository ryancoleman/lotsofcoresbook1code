!-------------------------------------------------------------------------------
! Subroutine : UpdatePhiMIC
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

 SUBROUTINE UpdatePhiMIC
!DIR$ ATTRIBUTES OFFLOAD:mic :: UpdatePhiMIC

!  Common Variables
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, offset

  offset = NX*NZ

!$OMP PARALLEL

! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY_MIC
     f2_mic( NX + NXG*( j + NYG_MIC*k ) + nxt_mic ) = f2_mic( NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
   DO j = 1, NY
     f1_mic( 1 + NXG*( j + NYG_MIC*k ) + nxt_mic )  = f1_mic( NX+1 + NXG*( j + NYG_MIC*k ) + nxt_mic )
   END DO
 END DO

! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     f4_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) = buff_cpu( i + NX*(k-1) )
   END DO
   DO i = 1, NX
     f3_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic )      = buff_cpu( i + NX*(k-1) + 6*offset )
   END DO
 END DO

! Z direction
!$OMP DO PRIVATE(i)
 DO j = 1, NY_MIC
    DO i = 1, NX
     f6_mic( i + NXG*( j + NYG_MIC*NZ ) + nxt_mic ) = f6_mic( i + NXG*( j ) + nxt_mic )
   END DO
   DO i = 1, NX
     f5_mic( i + NXG*( j + NYG_MIC ) + nxt_mic ) = f5_mic( i + NXG*( j + NYG_MIC*(NZ+1) ) + nxt_mic )
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
     buff_mic_phi(i+NX*(k-1)) = f0_mic( i + NXG*( 1 + NYG_MIC*k )       ) &
                          + f1_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                          + f2_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                          + f3_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                          + f4_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                          + f5_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                          + f6_mic( i + NXG*( 1 + NYG_MIC*k ) + nxt_mic )
   END DO

   DO i = 1, NX
     buff_mic_phi(i+NX*(k-1) + offset) = f0_mic( i + NXG*( NY_MIC + NYG_MIC*k )           ) &
                                   + f1_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                   + f2_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                   + f3_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                   + f4_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                   + f5_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                   + f6_mic( i + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic )
   END DO

 END DO


! Edges along X
!$OMP DO
 DO i = 1, NX

   buff_mic_phi_edge_x( i )   = f0_mic( i + NXG*( 1 + NYG_MIC )           ) &
                             + f1_mic( i + NXG*( 1 + NYG_MIC ) + nxt_mic ) &
                             + f2_mic( i + NXG*( 1 + NYG_MIC ) + nxt_mic ) &
                             + f3_mic( i + NXG*( 1 + NYG_MIC ) + nxt_mic ) &
                             + f4_mic( i + NXG*( 1 + NYG_MIC ) + nxt_mic ) &
                             + f5_mic( i + NXG*( 1 + NYG_MIC ) + nxt_mic ) &
                             + f6_mic( i + NXG*( 1 + NYG_MIC ) + nxt_mic )

   buff_mic_phi_edge_x( i + NX ) = f0_mic( i + NXG*( NY_MIC + NYG_MIC )           ) &
                                + f1_mic( i + NXG*( NY_MIC + NYG_MIC ) + nxt_mic ) &
                                + f2_mic( i + NXG*( NY_MIC + NYG_MIC ) + nxt_mic ) &
                                + f3_mic( i + NXG*( NY_MIC + NYG_MIC ) + nxt_mic ) &
                                + f4_mic( i + NXG*( NY_MIC + NYG_MIC ) + nxt_mic ) &
                                + f5_mic( i + NXG*( NY_MIC + NYG_MIC ) + nxt_mic ) &
                                + f6_mic( i + NXG*( NY_MIC + NYG_MIC ) + nxt_mic )

   buff_mic_phi_edge_x( i + 2*NX ) = f0_mic( i + NXG*( 1 + NYG_MIC*NZ )           ) &
                                  + f1_mic( i + NXG*( 1 + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f2_mic( i + NXG*( 1 + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f3_mic( i + NXG*( 1 + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f4_mic( i + NXG*( 1 + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f5_mic( i + NXG*( 1 + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f6_mic( i + NXG*( 1 + NYG_MIC*NZ ) + nxt_mic )

   buff_mic_phi_edge_x( i + 3*NX ) = f0_mic( i + NXG*( NY_MIC + NYG_MIC*NZ )           ) &
                                  + f1_mic( i + NXG*( NY_MIC + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f2_mic( i + NXG*( NY_MIC + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f3_mic( i + NXG*( NY_MIC + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f4_mic( i + NXG*( NY_MIC + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f5_mic( i + NXG*( NY_MIC + NYG_MIC*NZ ) + nxt_mic ) &
                                  + f6_mic( i + NXG*( NY_MIC + NYG_MIC*NZ ) + nxt_mic )


 END DO

! Edges along Z
!$OMP DO
 DO k = 1, NZ

   buff_mic_phi_edge_z( k ) = f0_mic( 1 + NXG*( 1 + NYG_MIC*k )           ) &
                           + f1_mic( 1 + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                           + f2_mic( 1 + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                           + f3_mic( 1 + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                           + f4_mic( 1 + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                           + f5_mic( 1 + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                           + f6_mic( 1 + NXG*( 1 + NYG_MIC*k ) + nxt_mic )

   buff_mic_phi_edge_z( k + NZ ) = f0_mic( NX + NXG*( 1 + NYG_MIC*k )           ) &
                                + f1_mic( NX + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                                + f2_mic( NX + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                                + f3_mic( NX + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                                + f4_mic( NX + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                                + f5_mic( NX + NXG*( 1 + NYG_MIC*k ) + nxt_mic ) &
                                + f6_mic( NX + NXG*( 1 + NYG_MIC*k ) + nxt_mic )


   buff_mic_phi_edge_z( k + 2*NZ ) = f0_mic( 1 + NXG*( NY_MIC + NYG_MIC*k )           ) &
                                  + f1_mic( 1 + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f2_mic( 1 + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f3_mic( 1 + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f4_mic( 1 + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f5_mic( 1 + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f6_mic( 1 + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic )

   buff_mic_phi_edge_z( k + 3*NZ ) = f0_mic( NX + NXG*( NY_MIC + NYG_MIC*k )           ) &
                                  + f1_mic( NX + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f2_mic( NX + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f3_mic( NX + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f4_mic( NX + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f5_mic( NX + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic ) &
                                  + f6_mic( NX + NXG*( NY_MIC + NYG_MIC*k ) + nxt_mic )
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE UpdatePhiMIC

