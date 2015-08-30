!-------------------------------------------------------------------------------
! Subroutine : MemAlloc
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Allocate (FLAG == 1) or deallocate (FLAG != 1) memory for common arrays.
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

 SUBROUTINE MemAlloc(FLAG)

!  Common Variables
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

 !   Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN
! Host arrays
 ALLOCATE( phi( 0:NG ) )
 ALLOCATE( lapPhi( 0:NG ) )
 ALLOCATE( gradPhiX( 0:NG ) )
 ALLOCATE( gradPhiY( 0:NG ) )
 ALLOCATE( gradPhiZ( 0:NG ) )

 ALLOCATE( f0( 0:NG ) )
 ALLOCATE( f1( 0:2*NG ) )
 ALLOCATE( f2( 0:2*NG ) )
 ALLOCATE( f3( 0:2*NG ) )
 ALLOCATE( f4( 0:2*NG ) )
 ALLOCATE( f5( 0:2*NG ) )
 ALLOCATE( f6( 0:2*NG ) )

 ALLOCATE( g0( 0:NG ) )
 ALLOCATE( g1( 0:2*NG ) )
 ALLOCATE( g2( 0:2*NG ) )
 ALLOCATE( g3( 0:2*NG ) )
 ALLOCATE( g4( 0:2*NG ) )
 ALLOCATE( g5( 0:2*NG ) )
 ALLOCATE( g6( 0:2*NG ) )
 ALLOCATE( g7( 0:2*NG ) )
 ALLOCATE( g8( 0:2*NG ) )
 ALLOCATE( g9( 0:2*NG ) )
 ALLOCATE( g10( 0:2*NG ) )
 ALLOCATE( g11( 0:2*NG ) )
 ALLOCATE( g12( 0:2*NG ) )
 ALLOCATE( g13( 0:2*NG ) )
 ALLOCATE( g14( 0:2*NG ) )
 ALLOCATE( g15( 0:2*NG ) )
 ALLOCATE( g16( 0:2*NG ) )
 ALLOCATE( g17( 0:2*NG ) )
 ALLOCATE( g18( 0:2*NG ) )

! MIC arrays
 ALLOCATE( phi_mic( 0:NG_MIC ) )
 ALLOCATE( lapPhi_mic( 0:NG_MIC ) )
 ALLOCATE( gradPhiX_mic( 0:NG_MIC ) )
 ALLOCATE( gradPhiY_mic( 0:NG_MIC ) )
 ALLOCATE( gradPhiZ_mic( 0:NG_MIC ) )

 ALLOCATE( f0_mic( 0:NG_MIC ) )
 ALLOCATE( f1_mic( 0:2*NG_MIC ) )
 ALLOCATE( f2_mic( 0:2*NG_MIC ) )
 ALLOCATE( f3_mic( 0:2*NG_MIC ) )
 ALLOCATE( f4_mic( 0:2*NG_MIC ) )
 ALLOCATE( f5_mic( 0:2*NG_MIC ) )
 ALLOCATE( f6_mic( 0:2*NG_MIC ) )

 ALLOCATE( g0_mic( 0:NG_MIC ) )
 ALLOCATE( g1_mic( 0:2*NG_MIC ) )
 ALLOCATE( g2_mic( 0:2*NG_MIC ) )
 ALLOCATE( g3_mic( 0:2*NG_MIC ) )
 ALLOCATE( g4_mic( 0:2*NG_MIC ) )
 ALLOCATE( g5_mic( 0:2*NG_MIC ) )
 ALLOCATE( g6_mic( 0:2*NG_MIC ) )
 ALLOCATE( g7_mic( 0:2*NG_MIC ) )
 ALLOCATE( g8_mic( 0:2*NG_MIC ) )
 ALLOCATE( g9_mic( 0:2*NG_MIC ) )
 ALLOCATE( g10_mic( 0:2*NG_MIC ) )
 ALLOCATE( g11_mic( 0:2*NG_MIC ) )
 ALLOCATE( g12_mic( 0:2*NG_MIC ) )
 ALLOCATE( g13_mic( 0:2*NG_MIC ) )
 ALLOCATE( g14_mic( 0:2*NG_MIC ) )
 ALLOCATE( g15_mic( 0:2*NG_MIC ) )
 ALLOCATE( g16_mic( 0:2*NG_MIC ) )
 ALLOCATE( g17_mic( 0:2*NG_MIC ) )
 ALLOCATE( g18_mic( 0:2*NG_MIC ) )

! Allocate memory for exchange buffers
 ALLOCATE( f_buff_mic(1:NXG*NZG*2) )
 ALLOCATE( f_buff_cpu(1:NXG*NZG*2) )

 ALLOCATE( buff_mic(1:NXG*NZG*12) )
 ALLOCATE( buff_cpu(1:NXG*NZG*12) )

 ALLOCATE( buff_mic_edge_x(1:NXG*4) )
 ALLOCATE( buff_mic_edge_z(1:NZG*4) )
 ALLOCATE( buff_cpu_edge_x(1:NXG*4) )
 ALLOCATE( buff_cpu_edge_z(1:NZG*4) )

 ALLOCATE( buff_mic_phi(1:NXG*NZG*2) )
 ALLOCATE( buff_cpu_phi(1:NXG*NZG*2) )

 ALLOCATE( buff_mic_phi_edge_x(1:NXG*4) )
 ALLOCATE( buff_mic_phi_edge_z(1:NZG*4) )
 ALLOCATE( buff_cpu_phi_edge_x(1:NXG*4) )
 ALLOCATE( buff_cpu_phi_edge_z(1:NZG*4) )

 ELSE

   DEALLOCATE( phi )
   DEALLOCATE( lapPhi )
   DEALLOCATE( gradPhiX )
   DEALLOCATE( gradPhiY )
   DEALLOCATE( gradPhiZ )
   DEALLOCATE( f0 )
   DEALLOCATE( f1 )
   DEALLOCATE( f2 )
   DEALLOCATE( f3 )
   DEALLOCATE( f4 )
   DEALLOCATE( f5 )
   DEALLOCATE( f6 )
   DEALLOCATE( g0 )
   DEALLOCATE( g1 )
   DEALLOCATE( g2 )
   DEALLOCATE( g3 )
   DEALLOCATE( g4 )
   DEALLOCATE( g5 )
   DEALLOCATE( g6 )
   DEALLOCATE( g7 )
   DEALLOCATE( g8 )
   DEALLOCATE( g9 )
   DEALLOCATE( g10 )
   DEALLOCATE( g11 )
   DEALLOCATE( g12 )
   DEALLOCATE( g13 )
   DEALLOCATE( g14 )
   DEALLOCATE( g15 )
   DEALLOCATE( g16 )
   DEALLOCATE( g17 )
   DEALLOCATE( g18 )

   DEALLOCATE( phi_mic )
   DEALLOCATE( lapPhi_mic )
   DEALLOCATE( gradPhiX_mic )
   DEALLOCATE( gradPhiY_mic )
   DEALLOCATE( gradPhiZ_mic )
   DEALLOCATE( f0_mic )
   DEALLOCATE( f1_mic )
   DEALLOCATE( f2_mic )
   DEALLOCATE( f3_mic )
   DEALLOCATE( f4_mic )
   DEALLOCATE( f5_mic )
   DEALLOCATE( f6_mic )
   DEALLOCATE( g0_mic )
   DEALLOCATE( g1_mic )
   DEALLOCATE( g2_mic )
   DEALLOCATE( g3_mic )
   DEALLOCATE( g4_mic )
   DEALLOCATE( g5_mic )
   DEALLOCATE( g6_mic )
   DEALLOCATE( g7_mic )
   DEALLOCATE( g8_mic )
   DEALLOCATE( g9_mic )
   DEALLOCATE( g10_mic )
   DEALLOCATE( g11_mic )
   DEALLOCATE( g12_mic )
   DEALLOCATE( g13_mic )
   DEALLOCATE( g14_mic )
   DEALLOCATE( g15_mic )
   DEALLOCATE( g16_mic )
   DEALLOCATE( g17_mic )
   DEALLOCATE( g18_mic )

   DEALLOCATE( f_buff_mic )
   DEALLOCATE( f_buff_cpu )
   DEALLOCATE( buff_mic )
   DEALLOCATE( buff_cpu )
   DEALLOCATE( buff_mic_edge_x )
   DEALLOCATE( buff_mic_edge_z )
   DEALLOCATE( buff_cpu_edge_x )
   DEALLOCATE( buff_cpu_edge_z )
   DEALLOCATE( buff_mic_phi )
   DEALLOCATE( buff_cpu_phi )
   DEALLOCATE( buff_mic_phi_edge_x )
   DEALLOCATE( buff_mic_phi_edge_z )
   DEALLOCATE( buff_cpu_phi_edge_x )
   DEALLOCATE( buff_cpu_phi_edge_z )

 END IF

 RETURN
 END SUBROUTINE MemAlloc
