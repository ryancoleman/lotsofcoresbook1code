!-------------------------------------------------------------------------------
! Subroutine : Stream
! Revision   : 1.0 (2008-06-15)
! Author     : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @file
!! Relaxation and streaming steps for the distribution function f
!> @details
!! Relaxation and streaming step for distribution function f in the dual grid
!! parallel D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.

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
 USE Domain,      ONLY : ni_f, xl_f, xu_f, yl_f, yu_f
 USE LBMParams,   ONLY : f, fcol
 USE FluidParams, ONLY : eta, eta2
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js

! Relaxation and streaming steps
 DO j = yl_f, yu_f
   DO i = xl_f, xu_f
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

     f(ie, j,1) = eta*fcol(i,j,1) + eta2*fcol(ie,j,1)
     f( i,jn,2) = eta*fcol(i,j,2) + eta2*fcol(i,jn,2)
     f(iw, j,3) = eta*fcol(i,j,3) + eta2*fcol(iw,j,3)
     f( i,js,4) = eta*fcol(i,j,4) + eta2*fcol(i,js,4)
   END DO
 END DO

 RETURN
 END SUBROUTINE Stream
