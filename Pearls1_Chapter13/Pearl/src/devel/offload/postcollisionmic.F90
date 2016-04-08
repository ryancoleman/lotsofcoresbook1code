!-------------------------------------------------------------------------------
! Subroutine : PostCollisionMIC
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Update g in nodes at processor boundaries and outward f in ghosts.
! This is faster than using a neighbor array in collision.
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

 SUBROUTINE PostCollisionMIC
!DIR$ ATTRIBUTES OFFLOAD:mic :: PostCollisionMIC

!  Common Variables
 USE Domain
 USE LBMParams
 IMPLICIT NONE

!  Local variables
 INTEGER :: i, j, k

!$OMP PARALLEL

! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY_MIC
     f1_mic( (NX+1) + NXG*( j + NYG_MIC*k ) + now_mic) = f1_mic( 1 + NXG*( j + NYG_MIC*k ) + now_mic)
   END DO
   DO j = 1, NY_MIC
     f2_mic(          NXG*( j + NYG_MIC*k ) + now_mic) = f2_mic( NX + NXG*( j + NYG_MIC*k ) + now_mic)
   END DO
 END DO

! This comes from the CPU
! We copy the f4 from the cpu exchange buffer into the bottom ghost
! and the f3 into the top ghost 
! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     f3_mic( i + NXG*( NY_MIC+1 + NYG_MIC*k ) + now_mic) = f_buff_cpu( i + NX*(k-1) ) ! From j=NY_MIC+1 on CPU
   END DO
   DO i = 1, NX
     f4_mic( i + NXG*(            NYG_MIC*k ) + now_mic) = f_buff_cpu( i + NX*(k-1) + NX*NZ )    ! From j=NY on CPU
   END DO
 END DO

! Z direction
!$OMP DO PRIVATE(i)
 DO j = 1, NY_MIC
   DO i = 1, NX
     f5_mic( i + NXG*( j + NYG_MIC*(NZ+1) ) + now_mic) = f5_mic( i + NXG*( j + NYG_MIC ) + now_mic)
   END DO
   DO i = 1, NX
     f6_mic( i + NXG*( j                  ) + now_mic) = f6_mic( i + NXG*( j + NYG_MIC*NZ ) + now_mic)
   END DO
 END DO

!
! At this stage we are ready to do the relaxation step in stream.f90
!

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE PostCollisionMIC

