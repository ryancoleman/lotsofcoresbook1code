!-------------------------------------------------------------------------------
! Subroutine : Stream
! Revision   : 1.0 (2008-06-15)
! Author     : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @file
!! Relaxation and streaming steps for the distribution function f
!> @details
!! Relaxation and streaming step for distribution function f in the parallel
!! D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.

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
 USE Domain,      ONLY : ni, now, nxt, xl, xu, yl, yu
 USE LBMParams,   ONLY : f
 USE FluidParams, ONLY : eta, eta2
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js

 DO j = yl, yu
   DO i = xl, xu

     ie = ni(i,j,1)
     jn = ni(i,j,2)
     iw = ni(i,j,3)
     js = ni(i,j,4)

     f(ie, j,1,nxt) = eta*f(i,j,1,now) + eta2*f(ie,j,1,now)
     f( i,jn,2,nxt) = eta*f(i,j,2,now) + eta2*f(i,jn,2,now)
     f(iw, j,3,nxt) = eta*f(i,j,3,now) + eta2*f(iw,j,3,now)
     f( i,js,4,nxt) = eta*f(i,j,4,now) + eta2*f(i,js,4,now)

   END DO
 END DO

 RETURN
 END SUBROUTINE Stream
