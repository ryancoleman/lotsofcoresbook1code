!-------------------------------------------------------------------------------
! Subroutine : Stream
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Relaxation and streaming steps for the distribution function f
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

 SUBROUTINE Stream

! Common Variables
 USE Domain,      ONLY : NX, NXG, NY, NYG, NZ, NZG, now, nxt
 USE FluidParams, ONLY : eta, eta2
 USE LBMParams
 USE OMP_LIB
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, m

! Here we MUST avoid ghost nodes
!$OMP PARALLEL DO PRIVATE(i,j,m)
 DO k = 1, NZ
   DO j = 1, NY
!DIR$ IVDEP
     DO i = 1, NX

       m = i + NXG*( j + NYG*k )

! Relax and stream the order paramter distribution f
       f1( m + 1                 + nxt ) = eta*f1(m+now) + eta2*f1(m+1+now)
       f2( m - 1                 + nxt ) = eta*f2(m+now) + eta2*f2(m-1+now)
       f3( m     + NXG           + nxt ) = eta*f3(m+now) + eta2*f3(m+NXG+now)
       f4( m     - NXG           + nxt ) = eta*f4(m+now) + eta2*f4(m-NXG+now)
       f5( m           + NXG*NYG + nxt ) = eta*f5(m+now) + eta2*f5(m+NXG*NYG+now)
       f6( m           - NXG*NYG + nxt ) = eta*f6(m+now) + eta2*f6(m-NXG*NYG+now)

     END DO
   END DO
 END DO

 RETURN
 END SUBROUTINE Stream
