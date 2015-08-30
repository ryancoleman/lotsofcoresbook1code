!-------------------------------------------------------------------------------
! Subroutine : Init
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Initialize all variables and arrays and save initialized data to file.
!> @details
!! Initialization step for the parallel dual grid D2Q9 Lee-Lin multiphase LBM.
!! Smeared interface initialized using equilibrium order parameter function for
!! each drop defined in the input (in the range [R-IntWidth,R+IntWidth]). The
!! distribution functions f and g are initialized to their equilibrium values
!! for zero velocity.

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
!--------------------------------------------------------------------------------

 SUBROUTINE Init

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MPIParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, m, ie, iw, jn, js, xs, ys
 REAL(KIND = DBL) :: R, lapRho_f
 REAL(KIND = DBL) :: gradRhoSqX_f, gradRhoSqY_f
 REAL(KIND = DBL) :: gradRhoXX2_f, gradRhoXY2_f, gradRhoYX2_f, gradRhoYY2_f

! Initialize counters
 iStep = 0
 tCall = 1
 eps   = 1.0D0

! Set vtk limits
 NX = xu - xl + 1
 NY = yu - yl + 1
 NX_f = xu_f - xl_f + 1
 NY_f = yu_f - yl_f + 1

!---------Near neighbors in the order parameter mesh ---------------------------
 DO j = ylg2_f, yug2_f
   DO i = xlg2_f, xug2_f
      ni_f(i,j,1) = i + 1
      ni_f(i,j,2) = j + 1
      ni_f(i,j,3) = i - 1
      ni_f(i,j,4) = j - 1
   END DO
 END DO

!---------Near neighbors in the momentum mesh ----------------------------------
 DO j = ylg, yug
   DO i = xlg, xug
     ni(i,j,1) = i + 1
     ni(i,j,2) = j + 1
     ni(i,j,3) = i - 1
     ni(i,j,4) = j - 1
   END DO
 END DO


! Initialize density for order parameter (rhoH assigned to continuous phase)
 DO j = ylg3_f, yug3_f
   DO i = xlg3_f, xug3_f

     rho_f(i,j) = rhoH

     DO m = 1, nBubbles
       R =  DSQRT( ( DBLE(i)-bubbles(m,1) )**2 + ( DBLE(j)-bubbles(m,2) )**2 )
       IF ( R <= (DBLE(bubbles(m,3)) + IntWidth) ) THEN
          rho_f(i,j) = rhoStar - 0.5*(rhoH - rhoL)*TANH( 2.D0*( DBLE(bubbles(m,3)) - R )/IntWidth )
       END IF
     END DO

   END DO
 END DO

! Initialize velocity, pressure and distribution function for order parameter
 DO j = ylg_f, yug_f
   DO i = xlg_f, xug_f

     u_f(i,j,1) = 0.D0
     u_f(i,j,2) = 0.D0
     p_f(i,j)   = 1.D0

     f(i,j,0  ) = W0*rho_f(i,j)
     f(i,j,1:4) = W1*rho_f(i,j)
     f(i,j,5:8) = W2*rho_f(i,j)

   END DO
 END DO

! Initialization for velocity, pressure, density and distribution function
! in the momentum grid
 DO j = ylg, yug
   DO i = xlg, xug

     xs = 2*i - 1
     ys = 2*j - 1
     rho(i,j) = rho_f(xs,ys)

     u(i,j,1) = 0.D0
     u(i,j,2) = 0.D0
     p(i,j)   = 1.D0

     g(i,j,0  ) = W0C
     g(i,j,1:4) = W1C
     g(i,j,5:8) = W2C

   END DO
 END DO


!-------- Differential terms for stress formulation of interfacial force -------
 DO j = ylg2_f, yug2_f
   DO i = xlg2_f, xug2_f 

! Identify first neighbors
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

! Laplacian of the density rho
     lapRho_f = 4.D0*( rho_f(ie,j ) + rho_f(iw,j ) + rho_f(i ,jn)            &
              + rho_f(i ,js) ) +  rho_f(ie,jn) + rho_f(ie,js) + rho_f(iw,jn) &
              + rho_f(iw,js) - 20.D0*rho_f(i,j)

! Define Psi
     psi_f(i,j) = beta4*( rho_f(i,j) - rhoStar )*( rho_f(i,j) - rhoL )*( rho_f(i,j) - rhoH ) &
                - kappa_6*lapRho_f
   END DO
 END DO

!-------- Differential terms for potential formulation of interfacial force ----
 DO j = ylg_f, yug_f
   DO i = xlg_f, xug_f

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

   END DO
 END DO

!-------- Interpolate differential terms in the momentum mesh ------------------
 DO j = yl, yu
   DO i = xl, xu

! Define the relative positions of the source nodes
     xs = 2*i - 1
     ys = 2*j - 1

! Identify first neighbors
     ie = ni_f(xs,ys,1)
     jn = ni_f(xs,ys,2)
     iw = ni_f(xs,ys,3)
     js = ni_f(xs,ys,4)

! Gradient of the square of the density gradient
     gradRhoSqX_f = ( 4.D0*( gradRhoSq_f(ie,ys) - gradRhoSq_f(iw,ys) ) &
                  + gradRhoSq_f(ie,jn) - gradRhoSq_f(iw,js)            &
                  + gradRhoSq_f(ie,js) - gradRhoSq_f(iw,jn) )*inv12

     gradRhoSqY_f = ( 4.D0*( gradRhoSq_f(xs,jn) - gradRhoSq_f(xs,js) ) &
                  + gradRhoSq_f(ie,jn) - gradRhoSq_f(iw,js)            &
                  + gradRhoSq_f(iw,jn) - gradRhoSq_f(ie,js) )*inv12

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

! Copy values of the gradient of rho from the overlapping f-mesh nodes
     gradRhoX(i,j) = 2.D0*gradRhoX_f(xs,ys)
     gradRhoY(i,j) = 2.D0*gradRhoY_f(xs,ys)

   END DO
 END DO

! Make sure all ghost layers are available in the momentum grid
 CALL MPI_UpdateRhoG

! Save initial data to file
 CALL Stats
 CALL VtkPlane

 RETURN
 END SUBROUTINE Init

