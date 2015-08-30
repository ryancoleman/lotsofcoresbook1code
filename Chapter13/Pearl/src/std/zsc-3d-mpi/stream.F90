!-------------------------------------------------------------------------------
! Subroutine : Stream
! Revision   : 1.0 (2008-06-15)
! Author     : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @file
!! Relaxation and streaming steps for the distribution function f
!> @details
!! Relaxation and streaming step for distribution function f in the parallel
!! D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM.

!-------------------------------------------------------------------------------
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
 USE Domain,      ONLY : ni, now, nxt, xl, xu, yl, yu, zl, zu
 USE FluidParams, ONLY : eta, eta2
 USE LBMParams,   ONLY : f
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, ie, iw, jn, js, kt, kb


 DO k = zl, zu
   DO j = yl, yu
     DO i = xl, xu

! Identify neighbours
       ie = ni(i,j,k,1)
       iw = ni(i,j,k,2)
       jn = ni(i,j,k,3)
       js = ni(i,j,k,4)
       kt = ni(i,j,k,5)
       kb = ni(i,j,k,6)

! Relax and stream the order paramter distribution f
       f(ie,j ,k ,1, nxt) = eta*f(i,j,k,1,now) + eta2*f(ie,j ,k ,1,now)
       f(iw,j ,k ,2, nxt) = eta*f(i,j,k,2,now) + eta2*f(iw,j ,k ,2,now)
       f(i ,jn,k ,3, nxt) = eta*f(i,j,k,3,now) + eta2*f(i ,jn,k ,3,now)
       f(i ,js,k ,4, nxt) = eta*f(i,j,k,4,now) + eta2*f(i ,js,k ,4,now)
       f(i ,j ,kt,5, nxt) = eta*f(i,j,k,5,now) + eta2*f(i ,j ,kt,5,now)
       f(i ,j ,kb,6, nxt) = eta*f(i,j,k,6,now) + eta2*f(i ,j ,kb,6,now)

     END DO
   END DO
 END DO

 RETURN
 END SUBROUTINE Stream
