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
 DO k = 1, NZ
   DO j = 1, NY
     DO i = 1, NX

! Relax and stream the order paramter distribution f
       f( i+1,j,k,1, nxt ) = eta*f(i,j,k,1,now) + eta2*f(i+1,j,k,1,now)
       f( i-1,j,k,2, nxt ) = eta*f(i,j,k,2,now) + eta2*f(i-1,j,k,2,now)
       f( i,j+1,k,3, nxt ) = eta*f(i,j,k,3,now) + eta2*f(i,j+1,k,3,now)
       f( i,j-1,k,4, nxt ) = eta*f(i,j,k,4,now) + eta2*f(i,j-1,k,4,now)
       f( i,j,k+1,5, nxt ) = eta*f(i,j,k,5,now) + eta2*f(i,j,k+1,5,now)
       f( i,j,k-1,6, nxt ) = eta*f(i,j,k,6,now) + eta2*f(i,j,k-1,6,now)

     END DO
   END DO
 END DO

 RETURN
 END SUBROUTINE Stream
