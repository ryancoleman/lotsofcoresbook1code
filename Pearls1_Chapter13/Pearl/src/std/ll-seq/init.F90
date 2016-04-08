!-------------------------------------------------------------------------------
! Subroutine : Init
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Initialize all variables and arrays and save initialized data to file.
!> @details
!! Initialization step for the D2Q9 Lee-Lin multiphase LBM. Smeared
!! interface initialized using equilibrium order parameter function for each
!! drop defined in the input (in the range [R-IntWidth,R+IntWidth]). The
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
!-------------------------------------------------------------------------------

 SUBROUTINE Init

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain,    ONLY : inv12, iStep, ni, NX, NY, tCall, xmax, xmin, ymax, ymin
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, m, ie, iw, jn, js
 REAL(KIND = DBL) :: R, rhon, lapRho

! Initialize counters
 iStep = 0
 tCall = 1
 eps   = 1.D0

! Set vtk limits
 NX = xmax - xmin + 1
 NY = ymax - ymin + 1

!-------- Define near neighbours and initialize the density rho ----------------
 DO j = ymin, ymax
   DO i = xmin, xmax

! Initialize density
      rho(i,j) = rhoH
      DO m = 1, nBubbles
        R =  DSQRT( ( DBLE(i)-bubbles(m,1) )**2 + ( DBLE(j)-bubbles(m,2) )**2 )
        IF (R <= (DBLE(bubbles(m,3)) + IntWidth)) THEN
          rho(i,j) = rhoStar - 0.5*(rhoH - rhoL)*TANH(2.D0*(DBLE(bubbles(m,3)) - R)/IntWidth)
        END IF
      END DO

! Generate neighbors array
      ni(i,j,1) = i + 1
      ni(i,j,2) = j + 1
      ni(i,j,3) = i - 1
      ni(i,j,4) = j - 1
   END DO
 END DO

! Fix neighbours at edges
 ni(xmin,:,3) = xmax
 ni(xmax,:,1) = xmin
 ni(:,ymin,4) = ymax
 ni(:,ymax,2) = ymin


!---------- Set equilibrium distributions, pressure and velocity ---------------
 DO j = ymin, ymax
   DO i = xmin, xmax

! Define velocity and pressure
     rhon = rho(i,j)
     u(i,j,1) = 0.D0
     u(i,j,2) = 0.D0
     p(i,j)   = 1.D0

! Equilibrium values of Gamma, and the phase and pressure distributions
     f(i,j,0)   = W0*rhon
     f(i,j,1:4) = W1*rhon
     f(i,j,5:8) = W2*rhon

     g(i,j,0)   = W0C
     g(i,j,1:4) = W1C
     g(i,j,5:8) = W2C

   END DO
 END DO

!--------- Differential terms for the stress form of the interfacial force -----
 DO j = ymin, ymax
   DO i = xmin, xmax

! Identify neighbors
     ie = ni(i,j,1)
     jn = ni(i,j,2)
     iw = ni(i,j,3)
     js = ni(i,j,4)

! Gradient of the density rho
     gradRhoX(i,j) = ( 4.D0*( rho(ie,j ) - rho(iw,j ) ) + rho(ie,jn) &
                   - rho(iw,js) + rho(ie,js) - rho(iw,jn) )*inv12

     gradRhoY(i,j) = ( 4.D0*( rho(i ,jn) - rho(i ,js) ) + rho(ie,jn) &
                   - rho(iw,js) + rho(iw,jn) - rho(ie,js) )*inv12

! Square of the density gradient (components and full)
     gradRhoXX(i,j) = gradRhoX(i,j)*gradRhoX(i,j)
     gradRhoXY(i,j) = gradRhoX(i,j)*gradRhoY(i,j)
     gradRhoYY(i,j) = gradRhoY(i,j)*gradRhoY(i,j)
     gradRhoSq(i,j) = gradRhoXX(i,j) + gradRhoYY(i,j)

! Laplacian of the density rho
     lapRho = 4.D0*( rho(ie,j ) + rho(iw,j ) + rho(i ,jn) + rho(i ,js) )       &
            + rho(ie,jn) + rho(ie,js) + rho(iw,jn) + rho(iw,js) - 20.D0*rho(i,j)

! Define the chemical potential Psi
     psi(i,j) = beta4*( rho(i,j) - rhoStar )*( rho(i,j) - rhoL )*( rho(i,j) - rhoH ) &
              - kappa_6*lapRho
   END DO
 END DO

! Save initialized data to file
 CALL Stats
 CALL VtkPlane

 RETURN
 END SUBROUTINE Init

