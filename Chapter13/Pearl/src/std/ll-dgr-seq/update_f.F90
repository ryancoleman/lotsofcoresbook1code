!-------------------------------------------------------------------------------
! Subroutine : UpdateF
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update velocity and pressure in the order parameter grid.
!> @details
!! Transfer the updated velocity and pressure values from the pressure/momentum
!! grid to the finer order parameter grid for the dual grid D2Q9 Lee-Lin
!! multiphase LBM. Bilinear interpolation is used for the information transfer.

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

 SUBROUTINE UpdateF

! Common Variables
 USE Domain,      ONLY : ni, ni_f, xmax, xmin, ymax, ymin
 USE FluidParams, ONLY : p, p_f, u, u_f
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j
 INTEGER :: xs1, xs2, xt1, xt2, xt3, ys1, ys2, yt1, yt2, yt3


 DO j = ymin, ymax-1
   DO i = xmin, xmax-1

! Coordinates of source nodes for the interpolation (xs,ys)
     xs1 = i
     ys1 = j
     xs2 = ni(i,j,1)
     ys2 = ni(i,j,2)

! Coordinates of target nodes for the interpolation (xt,yt)
     xt1 = 2*i - 1
     yt1 = 2*j - 1
     xt2 = ni_f(xt1,yt1,1)
     yt2 = ni_f(xt1,yt1,2)
     xt3 = ni_f(xt2,yt2,1)
     yt3 = ni_f(xt2,yt2,2)
 
! Target 1
     u_f(xt1,yt1,1) = u(xs1,ys1,1)
     u_f(xt1,yt1,2) = u(xs1,ys1,2)
     p_f(xt1,yt1)   = p(xs1,ys1)

! Target 2
     u_f(xt2,yt1,1) = 0.5D0*( u(xs1,ys1,1) + u(xs2,ys1,1) )
     u_f(xt2,yt1,2) = 0.5D0*( u(xs1,ys1,2) + u(xs2,ys1,2) )
     p_f(xt2,yt1)   = 0.5D0*( p(xs1,ys1) + p(xs2,ys1) )

! Target 3
     u_f(xt3,yt1,1) = u(xs2,ys1,1)
     u_f(xt3,yt1,2) = u(xs2,ys1,2)
     p_f(xt3,yt1)   = p(xs2,ys1)

! Target 4
    u_f(xt1,yt2,1) = 0.5D0*( u(xs1,ys1,1) + u(xs1,ys2,1) )
    u_f(xt1,yt2,2) = 0.5D0*( u(xs1,ys1,2) + u(xs1,ys2,2) )
    p_f(xt1,yt2)   = 0.5D0*( p(xs1,ys1) + p(xs1,ys2) )

! Target 5
    u_f(xt2,yt2,1) = 0.25D0*( u(xs1,ys1,1) + u(xs2,ys1,1) + u(xs1,ys2,1) + u(xs2,ys2,1) )
    u_f(xt2,yt2,2) = 0.25D0*( u(xs1,ys1,2) + u(xs2,ys1,2) + u(xs1,ys2,2) + u(xs2,ys2,2) )
    p_f(xt2,yt2)   = 0.25D0*( p(xs1,ys1)   + p(xs2,ys1)   + p(xs1,ys2)   + p(xs2,ys2) )

! Target 6
    u_f(xt3,yt2,1) = 0.5D0*( u(xs2,ys1,1) + u(xs2,ys2,1) )
    u_f(xt3,yt2,2) = 0.5D0*( u(xs2,ys1,2) + u(xs2,ys2,2) )
    p_f(xt3,yt2)   = 0.5D0*( p(xs2,ys1) + p(xs2,ys2) )

! Target 7
    u_f(xt1,yt3,1) = u(xs1,ys2,1)
    u_f(xt1,yt3,2) = u(xs1,ys2,2)
    p_f(xt1,yt3)   = p(xs1,ys2)

! Target 8
    u_f(xt2,yt3,1) = 0.5D0*( u(xs1,ys2,1) + u(xs2,ys2,1) )
    u_f(xt2,yt3,2) = 0.5D0*( u(xs1,ys2,2) + u(xs2,ys2,2) )
    p_f(xt2,yt3)   = 0.5D0*( p(xs1,ys2) + p(xs2,ys2) )

! Target 9
    u_f(xt3,yt3,1) = u(xs2,ys2,1)
    u_f(xt3,yt3,2) = u(xs2,ys2,2)
    p_f(xt3,yt3)   = p(xs2,ys2)

   END DO
 END DO

 RETURN
 END SUBROUTINE UpdateF
