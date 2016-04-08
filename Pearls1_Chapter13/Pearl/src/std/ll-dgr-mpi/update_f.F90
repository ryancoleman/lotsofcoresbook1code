!-------------------------------------------------------------------------------
! Subroutine : UpdateF
! Revision   : 1.1 (2008-11-10)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update velocity and pressure in the order parameter grid.
!> @details
!! Transfer the updated velocity and pressure values from the pressure/momentum
!! grid to the finer order parameter grid for the parallel dual grid D2Q9
!! Lee-Lin multiphase LBM. Bilinear interpolation is used for the information
!! transfer.
!!
!! Special care is taken of order parameter grid nodes that lay between pressure
!! grid nodes at the domain edges. These points use the values of the velocity
!! and presure in the ghost layers of the pressure grid for the interpolation,
!! which must have been updated before this function is called.

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
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, xs1, xs2, xt1, xt2, xt3, ys1, ys2, yt1, yt2, yt3

!--------- Interpolation of velocity and pressure on the order parameter mesh --
 DO j = yl, yu-1
   DO i = xl, xu-1

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
     p_f(xt2,yt1)   = 0.5D0*( p(xs1,ys1)   + p(xs2,ys1) )

! Target 3
     u_f(xt3,yt1,1) = u(xs2,ys1,1)
     u_f(xt3,yt1,2) = u(xs2,ys1,2)
     p_f(xt3,yt1)   = p(xs2,ys1)

! Target 4
    u_f(xt1,yt2,1) = 0.5D0*( u(xs1,ys1,1) + u(xs1,ys2,1) )
    u_f(xt1,yt2,2) = 0.5D0*( u(xs1,ys1,2) + u(xs1,ys2,2) )
    p_f(xt1,yt2)   = 0.5D0*( p(xs1,ys1)   + p(xs1,ys2) )

! Target 5
    u_f(xt2,yt2,1) = 0.25D0*( u(xs1,ys1,1) + u(xs2,ys1,1) + u(xs1,ys2,1) + u(xs2,ys2,1) )
    u_f(xt2,yt2,2) = 0.25D0*( u(xs1,ys1,2) + u(xs2,ys1,2) + u(xs1,ys2,2) + u(xs2,ys2,2) )
    p_f(xt2,yt2)   = 0.25D0*( p(xs1,ys1)   + p(xs2,ys1)   + p(xs1,ys2)   + p(xs2,ys2) )

! Target 6
    u_f(xt3,yt2,1) = 0.5D0*( u(xs2,ys1,1) + u(xs2,ys2,1) )
    u_f(xt3,yt2,2) = 0.5D0*( u(xs2,ys1,2) + u(xs2,ys2,2) )
    p_f(xt3,yt2)   = 0.5D0*( p(xs2,ys1)   + p(xs2,ys2) )

! Target 7
    u_f(xt1,yt3,1) = u(xs1,ys2,1)
    u_f(xt1,yt3,2) = u(xs1,ys2,2)
    p_f(xt1,yt3)   = p(xs1,ys2)

! Target 8
    u_f(xt2,yt3,1) = 0.5D0*( u(xs1,ys2,1) + u(xs2,ys2,1) )
    u_f(xt2,yt3,2) = 0.5D0*( u(xs1,ys2,2) + u(xs2,ys2,2) )
    p_f(xt2,yt3)   = 0.5D0*( p(xs1,ys2)   + p(xs2,ys2) )

! Target 9
    u_f(xt3,yt3,1) = u(xs2,ys2,1)
    u_f(xt3,yt3,2) = u(xs2,ys2,2)
    p_f(xt3,yt3)   = p(xs2,ys2)

   END DO
 END DO

! Update leftmost f values which are outside this processor's momentum grid
 DO j = yl, yu-1

! Coordinates of source nodes for the interpolation (xs,ys)
   xs1 = xlg
   ys1 = j
   xs2 = xl
   ys2 = ni(xlg,j,2)

! Coordinates of target nodes for the interpolation (xt,yt)
   xt1 = xl_f
   yt1 = 2*j - 1
   yt2 = ni_f(xt1,yt1,2)
   yt3 = ni_f(xt1,yt2,2)

! Target 2
   u_f(xt1,yt1,1) = 0.5D0*( u(xs1,ys1,1) + u(xs2,ys1,1) )
   u_f(xt1,yt1,2) = 0.5D0*( u(xs1,ys1,2) + u(xs2,ys1,2) )
   p_f(xt1,yt1)   = 0.5D0*( p(xs1,ys1)   + p(xs2,ys1) )

! Target 5
   u_f(xt1,yt2,1) = 0.25D0*( u(xs1,ys1,1) + u(xs2,ys1,1) + u(xs1,ys2,1) + u(xs2,ys2,1) )
   u_f(xt1,yt2,2) = 0.25D0*( u(xs1,ys1,2) + u(xs2,ys1,2) + u(xs1,ys2,2) + u(xs2,ys2,2) )
   p_f(xt1,yt2)   = 0.25D0*( p(xs1,ys1)   + p(xs2,ys1)   + p(xs1,ys2)   + p(xs2,ys2) )

! Target 8
   u_f(xt1,yt3,1) = 0.5D0*( u(xs1,ys2,1) + u(xs2,ys2,1) )
   u_f(xt1,yt3,2) = 0.5D0*( u(xs1,ys2,2) + u(xs2,ys2,2) )
   p_f(xt1,yt3)   = 0.5D0*( p(xs1,ys2)   + p(xs2,ys2) )
 END DO

! Update lowermost f values which are outside this processor's momentum grid
DO i = xl, xu-1

! Coordinates of source nodes for the interpolation (xs,ys)
   xs1 = i
   ys1 = ylg
   xs2 = ni(i,ylg,1)
   ys2 = yl

! Coordinates of target nodes for the interpolation (xt,yt)
   xt1 = 2*i - 1
   yt1 = yl_f
   xt2 = ni_f(xt1,yt1,1)
   xt3 = ni_f(xt2,yt1,1)

! Target 4
   u_f(xt1,yt1,1) = 0.5D0*( u(xs1,ys1,1) + u(xs1,ys2,1) )
   u_f(xt1,yt1,2) = 0.5D0*( u(xs1,ys1,2) + u(xs1,ys2,2) )
   p_f(xt1,yt1)   = 0.5D0*( p(xs1,ys1)   + p(xs1,ys2) )

! Target 5
   u_f(xt2,yt1,1) = 0.25D0*( u(xs1,ys1,1) + u(xs2,ys1,1) + u(xs1,ys2,1) + u(xs2,ys2,1) )
   u_f(xt2,yt1,2) = 0.25D0*( u(xs1,ys1,2) + u(xs2,ys1,2) + u(xs1,ys2,2) + u(xs2,ys2,2) )
   p_f(xt2,yt1)   = 0.25D0*( p(xs1,ys1)   + p(xs2,ys1)   + p(xs1,ys2)   + p(xs2,ys2) )

! Target 6
   u_f(xt3,yt1,1) = 0.5D0*( u(xs2,ys1,1) + u(xs2,ys2,1) )
   u_f(xt3,yt1,2) = 0.5D0*( u(xs2,ys1,2) + u(xs2,ys2,2) )
   p_f(xt3,yt1)   = 0.5D0*( p(xs2,ys1)   + p(xs2,ys2) )
 END DO

! Update southwest f value outside this processor's momentum grid
! Coordinates of source nodes for the interpolation (xs,ys)
 xs1 = xlg
 ys1 = ylg
 xs2 = xl
 ys2 = yl

! Coordinates of target nodes for the interpolation (xt,yt)
 xt1 = xl_f
 yt1 = yl_f

! Target 5
 u_f(xt1,yt1,1) = 0.25D0*( u(xs1,ys1,1) + u(xs2,ys1,1) + u(xs1,ys2,1) + u(xs2,ys2,1) )
 u_f(xt1,yt1,2) = 0.25D0*( u(xs1,ys1,2) + u(xs2,ys1,2) + u(xs1,ys2,2) + u(xs2,ys2,2) )
 p_f(xt1,yt1)   = 0.25D0*( p(xs1,ys1)   + p(xs2,ys1)   + p(xs1,ys2)   + p(xs2,ys2) )

 RETURN
 END SUBROUTINE updateF
