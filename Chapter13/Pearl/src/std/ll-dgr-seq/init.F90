!-------------------------------------------------------------------------------
! Subroutine : Init
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Initialize all variables and arrays and save initialized data to file.
!> @details
!! Initialization step for the dual grid D2Q9 Lee-Lin multiphase LBM. Smeared
!! interface initialized using equilibrium order parameter function for each
!! drop defined in the input (in the range [R-IntWidth,R+IntWidth]). The
!! distribution functions f and g are initialized to their equilibrium values
!! for zero velocity. Biliner interpolation used for the transfer of the density
!! values from the order parameter grid to the pressure grid.

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

 SUBROUTINE Init

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, m, ie, iw, jn, js, xs, ys
 REAL(KIND = DBL) :: R, lapRho_f, gradRhoSqX_f, gradRhoSqY_f
 REAL(KIND = DBL) :: gradRhoXX2_f, gradRhoXY2_f, gradRhoYX2_f, gradRhoYY2_f


! Initialize counters
 iStep = 0
 tCall = 1
 eps   = 1.0D0

! Define vtk limits
 NX   = xmax - xmin + 1
 NY   = ymax - ymin + 1
 NX_f = xmax_f - xmin_f + 1
 NY_f = ymax_f - ymin_f + 1


!---------Near neighbors in the order parameter mesh ---------------------------
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f
      ni_f(i,j,1) = i + 1
      ni_f(i,j,2) = j + 1
      ni_f(i,j,3) = i - 1
      ni_f(i,j,4) = j - 1
   END DO
 END DO
! Fix near neighbours at edges
 ni_f(xmin_f,:,3) = xmax_f
 ni_f(xmax_f,:,1) = xmin_f
 ni_f(:,ymin_f,4) = ymax_f
 ni_f(:,ymax_f,2) = ymin_f

!---------Near neighbors in the momentum mesh ----------------------------------
 DO j = ymin, ymax
   DO i = xmin, xmax
     ni(i,j,1) = i + 1
     ni(i,j,2) = j + 1
     ni(i,j,3) = i - 1
     ni(i,j,4) = j - 1
   END DO
 END DO
! Fix near neighbours at edges
 ni(xmin,:,3) = xmax
 ni(xmax,:,1) = xmin
 ni(:,ymin,4) = ymax
 ni(:,ymax,2) = ymin


!---------- Initialization for order parameter grid arrays ---------------------
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

! Initialize density rho
     rho_f(i,j) = rhoH
     DO m = 1, nBubbles
       R =  DSQRT( ( DBLE(i)-bubbles(m,1) )**2 + ( DBLE(j)-bubbles(m,2) )**2 )
       IF ( R <= (DBLE(bubbles(m,3)) + IntWidth) ) THEN
          rho_f(i,j) = rhoStar - 0.5*(rhoH - rhoL)*TANH( 2.D0*( DBLE(bubbles(m,3)) - R )/IntWidth )
       END IF
     END DO

! Initialize velocity and pressure
     u_f(i,j,1) = 0.D0
     u_f(i,j,2) = 0.D0
     p_f(i,j)   = 1.D0

! Initialize the order parameter distribution function
     f(i,j,0  ) = W0*rho_f(i,j)
     f(i,j,1:4) = W1*rho_f(i,j)
     f(i,j,5:8) = W2*rho_f(i,j)
   END DO
 END DO

!---------- Initialization for the momentum/pressure grid arrays ---------------
 DO j = ymin, ymax
   DO i = xmin, xmax

! Copy the value of rho from the overlapping f-mesh nodes
     xs = 2*i - 1
     ys = 2*j - 1
     rho(i,j) = rho_f(xs,ys)

! Define velocity and pressure
     u(i,j,1) = 0.D0
     u(i,j,2) = 0.D0
     p(i,j)   = 1.D0

! Equilibrium values of Gamma, and the phase and pressure distributions
     g(i,j,0  ) = W0C*p(i,j)
     g(i,j,1:4) = W1C*p(i,j)
     g(i,j,5:8) = W2C*p(i,j)
   END DO
 END DO


!---------- Differential terms in order parameter grid (stress formulation) ----
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

! Identify first neighbors
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

! Gradient of the density rho
     gradRhoX_f(i,j) = ( 4.D0*( rho_f(ie,j ) - rho_f(iw,j ) ) + rho_f(ie,jn) &
                     - rho_f(iw,jn) + rho_f(ie,js) - rho_f(iw,js) )*inv12

     gradRhoY_f(i,j) = ( 4.D0*( rho_f(i ,jn) - rho_f(i ,js) ) + rho_f(ie,jn) &
                     - rho_f(ie,js) + rho_f(iw,jn) - rho_f(iw,js) )*inv12

! Square of the density gradient (components and full)
     gradRhoXX_f(i,j) = gradRhoX_f(i,j)*gradRhoX_f(i,j)
     gradRhoXY_f(i,j) = gradRhoX_f(i,j)*gradRhoY_f(i,j)
     gradRhoYY_f(i,j) = gradRhoY_f(i,j)*gradRhoY_f(i,j)
     gradRhoSq_f(i,j) = gradRhoXX_f(i,j) + gradRhoYY_f(i,j)

! Laplacian of the density rho
     lapRho_f = 4.D0*( rho_f(ie,j ) + rho_f(iw,j ) + rho_f(i ,jn) + rho_f(i ,js) ) &
              +  rho_f(ie,jn) + rho_f(ie,js) + rho_f(iw,jn) + rho_f(iw,js)         &
              - 20.D0*rho_f(i,j)

! Define Psi
     psi_f(i,j) = 4.D0*beta*( rho_f(i,j) - rhoStar )*( rho_f(i,j) - rhoL )*( rho_f(i,j) - rhoH ) &
                - kappa_6*lapRho_f

   END DO
 END DO


!---------- Interpolate differential terms in the momentum mesh ----------------
! 1/12 factor embedded in kappaEf
 DO j = ymin, ymax
   DO i = xmin, xmax

! Relative positions of the source nodes (overlapping nodes in f and g grids)
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

! Copy and rescale values of the gradient of rho from the overlapping f-mesh nodes
     gradRhoX(i,j) = 2.D0*gradRhoX_f(xs,ys)
     gradRhoY(i,j) = 2.D0*gradRhoY_f(xs,ys)

   END DO
 END DO

! Save initialized data to file
 CALL Stats
 CALL VtkPlane
 
 RETURN
 END SUBROUTINE Init

