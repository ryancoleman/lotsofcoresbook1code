!-------------------------------------------------------------------------------
! Subroutine : Stream
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

 SUBROUTINE StreamMIC
!DIR$ ATTRIBUTES OFFLOAD:mic :: StreamMIC

! Common Variables
 USE Domain
 USE FluidParams
 USE LBMParams
 USE OMP_LIB
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, m

! Here we MUST avoid ghost nodes
!$OMP PARALLEL DO PRIVATE(i,j,m)
 DO k = 1, NZ
   DO j = 1, NY_MIC
!DIR$ IVDEP
     DO i = 1, NX

       m = i + NXG*( j + NYG_MIC*k )

! Relax and stream the order paramter distribution f
       f1_mic( m + 1                 + nxt_mic ) = eta*f1_mic(m+now_mic) + eta2*f1_mic(m+1+now_mic)
       f2_mic( m - 1                 + nxt_mic ) = eta*f2_mic(m+now_mic) + eta2*f2_mic(m-1+now_mic)
       f3_mic( m     + NXG           + nxt_mic ) = eta*f3_mic(m+now_mic) + eta2*f3_mic(m+NXG+now_mic)
       f4_mic( m     - NXG           + nxt_mic ) = eta*f4_mic(m+now_mic) + eta2*f4_mic(m-NXG+now_mic)
       f5_mic( m           + NXG*NYG_MIC + nxt_mic ) = eta*f5_mic(m+now_mic) + eta2*f5_mic(m+NXG*NYG_MIC+now_mic)
       f6_mic( m           - NXG*NYG_MIC + nxt_mic ) = eta*f6_mic(m+now_mic) + eta2*f6_mic(m-NXG*NYG_MIC+now_mic)

     END DO
   END DO
 END DO

 RETURN
 END SUBROUTINE StreamMIC
