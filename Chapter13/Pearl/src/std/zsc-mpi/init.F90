!-------------------------------------------------------------------------------
! Subroutine  : Init
! Revision    : 1.0 (2008-06-15)
! Author      : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @file
!! Initialize all variables and arrays and save initialized data to file.
!> @details
!! Initialization step for the parallel D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.
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
!-------------------------------------------------------------------------------

 SUBROUTINE Init

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, m, ie, iw, jn, js
 REAL(KIND = DBL) :: R, rhon, phin, muPhin, lapPhi
 REAL(KIND = DBL) :: Af0, Af1
 REAL(KIND = DBL) :: Ag0, Ag1, Eg1A, Eg2A


! Initialize counters
 now   = 0
 nxt   = 1
 iStep = 0
 tCall = 1
 eps   = 1.D0

! Domain size for vtk output
 NX = xu - xl + 1
 NY = yu - yl + 1

! Set initial average density
 rhon = 0.5D0*(rhoH + rhoL)

!--------- Initialize neighbours and the order parameter phi -------------------
 DO j = ylg, yug
   DO i = xlg, xug

      phi(i,j) = -phistar
      DO m = 1, nBubbles
        R =  DSQRT((DBLE(i)-bubbles(m,1))**2 + (DBLE(j)-bubbles(m,2))**2)
        IF ( R <= ( DBLE(bubbles(m,3)) + IntWidth ) ) THEN
          phi(i,j) = phistar*TANH( 2.D0*( DBLE(bubbles(m,3)) - R )/IntWidth )
        END IF
      END DO

! Neighbors array
      ni(i,j,1) = i + 1
      ni(i,j,2) = j + 1
      ni(i,j,3) = i - 1
      ni(i,j,4) = j - 1

   END DO
 END DO

!---------- Equilibrium distribution functions ---------------------------------
 DO j = yl, yu
   DO i = xl, xu

     u(i,j,:) = 0.D0
     p(i,j)   = 1.D0
     phin     = phi(i,j)

! Identify neighbors
     ie = ni(i,j,1)
     jn = ni(i,j,2)
     iw = ni(i,j,3)
     js = ni(i,j,4)

! Laplacian of the order parameter
     lapPhi = ( phi(ie,jn) + phi(ie,js) + phi(iw,jn) + phi(iw,js)   &
            + 4.D0*(phi(ie,j) + phi(iw,j) + phi(i,jn) + phi(i,js) ) &
            - 20.D0*phin )*inv6

! Chemical potential
     muPhin = alpha4*phin*(phin*phin - phistar2) - kappa*lapPhi

! Equilibrium Coeficients for f
     Af1 = 0.5D0*Gamma*muPhin
     Af0 = -2.D0*Gamma*muPhin

! Equilibrium distribution f for zero initial velocity
     f(i,j,0,now)   = Af0 + phin
     f(i,j,1:4,now) = Af1

! Equilibrium coeficients for g
     Ag1  = rhon + 3.D0*phin*muPhin
     Ag0  = 2.25D0*rhon - 1.25D0*Ag1
     Eg1A = Eg1*Ag1
     Eg2A = Eg2*Ag1

! Equilibrium distribution g for zero initial velocity
     g(i,j,0,now)   = Eg0*Ag0
     g(i,j,1:4,now) = Eg1A
     g(i,j,5:8,now) = Eg2A

! Set pressure and velocity values
     p(i,j)   = 1.D0
     u(i,j,:) = 0.D0

   END DO
 END DO


! Save initialized data
 CALL Stats
 CALL VtkPlane

 RETURN
 END SUBROUTINE Init
