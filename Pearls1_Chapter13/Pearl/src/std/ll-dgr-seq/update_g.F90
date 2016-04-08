!-------------------------------------------------------------------------------
! Subroutine : UpdateG
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update density and its differentials in the pressure/momentum grid.
!> @details
!! Copy the updated density and differential terms from the order parameter grid
!! to the overlapping nodes in the coarser pressure/momentum grid for the dual
!! grid D2Q9 Lee-Lin multiphase LBM.

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

 SUBROUTINE UpdateG

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain,    ONLY : ni_f, xmax, xmax_f, xmin, xmin_f, ymax, ymax_f, ymin, ymin_f
 USE FluidParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js
 INTEGER :: xs, ys
 REAL(KIND = DBL) :: gradRhoSqX_f, gradRhoSqY_f
 REAL(KIND = DBL) :: gradRhoXX2_f, gradRhoXY2_f, gradRhoYX2_f, gradRhoYY2_f


! Square of the density gradient (components and full)
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f
     gradRhoXX_f(i,j) = gradRhoX_f(i,j)*gradRhoX_f(i,j)
     gradRhoXY_f(i,j) = gradRhoX_f(i,j)*gradRhoY_f(i,j)
     gradRhoYY_f(i,j) = gradRhoY_f(i,j)*gradRhoY_f(i,j)
     gradRhoSq_f(i,j) = gradRhoXX_f(i,j) + gradRhoYY_f(i,j)
   END DO
 END DO

!--------- Interpolate the density in the pressure/momentum mesh ---------------
 DO j = ymin, ymax
   DO i = xmin, xmax

! Relative positions of the source and target nodes
     xs = 2*i - 1
     ys = 2*j - 1

! Identify first neighbors
     ie = ni_f(xs,ys,1)
     jn = ni_f(xs,ys,2)
     iw = ni_f(xs,ys,3)
     js = ni_f(xs,ys,4)

! Gradient of the square of the density gradient
     gradRhoSqX_f = 4.D0*( gradRhoSq_f(ie,ys) - gradRhoSq_f(iw,ys) ) &
                  + gradRhoSq_f(ie,jn) - gradRhoSq_f(iw,js)          &
                  + gradRhoSq_f(ie,js) - gradRhoSq_f(iw,jn)

     gradRhoSqY_f = 4.D0*( gradRhoSq_f(xs,jn) - gradRhoSq_f(xs,js) ) &
                  + gradRhoSq_f(ie,jn) - gradRhoSq_f(iw,js)          &
                  + gradRhoSq_f(iw,jn) - gradRhoSq_f(ie,js)

! Second derivatives of rho
     gradRhoXX2_f = 4.D0*( gradRhoXX_f(ie,ys) - gradRhoXX_f(iw,ys) ) &
                  + gradRhoXX_f(ie,jn) - gradRhoXX_f(iw,js)          &
                  + gradRhoXX_f(ie,js) - gradRhoXX_f(iw,jn)

     gradRhoXY2_f = 4.D0*( gradRhoXY_f(xs,jn) - gradRhoXY_f(xs,js) ) &
                  + gradRhoXY_f(ie,jn) - gradRhoXY_f(iw,js)          &
                  + gradRhoXY_f(iw,jn) - gradRhoXY_f(ie,js)

     gradRhoYX2_f = 4.D0*( gradRhoXY_f(ie,ys) - gradRhoXY_f(iw,ys) ) &
                  + gradRhoXY_f(ie,jn) - gradRhoXY_f(iw,js)          &
                  + gradRhoXY_f(ie,js) - gradRhoXY_f(iw,jn)

     gradRhoYY2_f = 4.D0*( gradRhoYY_f(xs,jn) - gradRhoYY_f(xs,js) ) &
                  + gradRhoYY_f(ie,jn) - gradRhoYY_f(iw,js)          &
                  + gradRhoYY_f(iw,jn) - gradRhoYY_f(ie,js)

! Arrays necessary in preStream and postStream
     GX(i,j) = kappaEf*( gradRhoSqX_f - gradRhoXX2_f - gradRhoXY2_f )
     GY(i,j) = kappaEf*( gradRhoSqY_f - gradRhoYX2_f - gradRhoYY2_f )

! Copy the value of rho from the overlapping f-mesh nodes
     rho(i,j) = rho_f(xs,ys)

! Copy values of the gradient of rho from the overlapping f-mesh nodes
     gradRhoX(i,j) = 2.D0*gradRhoX_f(xs,ys)
     gradRhoY(i,j) = 2.D0*gradRhoY_f(xs,ys)

   END DO
 END DO

 RETURN
 END SUBROUTINE UpdateG
