!-------------------------------------------------------------------------------
! Subroutine : PostStream
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Post-Stream step for distribution functions f and g
!> @details
!! Post-Stream step for order parameter distribution function f and pressure
!! distribution function g in the D2Q9 Lee-Lin multiphase LBM.
!!
!! All differential terms are calculated using a second order centered
!! differencing scheme.

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

 SUBROUTINE PostStream

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain,    ONLY : inv12, ni, xmax, xmin, ymax, ymin
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js
 REAL(KIND = DBL) :: Ff, Fg, GX, GY, UG, v, ux, uy, U_sq
 REAL(KIND = DBL) :: rhon, tau, invTau, invTau2, invCsRho
 REAL(KIND = DBL) :: gradPsiX, gradPsiY, gradRhoSqX, gradRhoSqY
 REAL(KIND = DBL) :: gradRhoXX2, gradRhoXY2, gradRhoYX2, gradRhoYY2
 REAL(KIND = DBL) :: gradRhoU, gradPsiU, gradRhoPsi
 REAL(KIND = DBL), DIMENSION(0:8) :: feq, geq, Gam, GamW
 REAL(KIND = DBL), DIMENSION(1:8) :: gradPsiD, gradRhoD


!-------- Collision-Stream Loop ------------------------------------------------
 DO j = ymin, ymax
   DO i = xmin, xmax

! Define some local values
     rhon = rho(i,j)
     ux   = u(i,j,1)
     uy   = u(i,j,2)
     U_sq = 0.5D0*( ux*ux + uy*uy )
     tau  = tauL + (rhon - rhoL)*tauRhoStar
     invTau   = 1.D0/(1.D0 + 2.D0*tau)
     invTau2  = tau*invTau
     invCsRho = invCs_sq*rhon

! Identify neighboring cells
     ie = ni(i,j,1)
     jn = ni(i,j,2)
     iw = ni(i,j,3)
     js = ni(i,j,4)

!--------- Standard Derivatives (1/12 factor embedded in kappa_12) -------------
! Gradient of the chemical potential psi
     gradPsiX = ( 4.D0*( psi(ie,j ) - psi(iw,j ) ) + psi(ie,jn) - psi(iw,js) &
              + psi(ie,js) - psi(iw,jn) )*inv12

     gradPsiY = ( 4.D0*( psi(i ,jn) - psi(i ,js) ) + psi(ie,jn) - psi(iw,js) &
              + psi(iw,jn) - psi(ie,js) )*inv12

! Gradient of the density gradient squared
     gradRhoSqX = 4.D0*( gradRhoSq(ie,j ) - gradRhoSq(iw,j ) )           &
                + gradRhoSq(ie,jn) - gradRhoSq(iw,js) + gradRhoSq(ie,js) &
                - gradRhoSq(iw,jn)

     gradRhoSqY = 4.D0*( gradRhoSq(i ,jn) - gradRhoSq(i ,js) )           &
                + gradRhoSq(ie,jn) - gradRhoSq(iw,js) + gradRhoSq(iw,jn) &
                - gradRhoSq(ie,js)

! Second derivatives of rho
     gradRhoXX2 = 4.D0*( gradRhoXX(ie,j ) - gradRhoXX(iw,j ) )           &
                + gradRhoXX(ie,jn) - gradRhoXX(iw,js) + gradRhoXX(ie,js) &
                - gradRhoXX(iw,jn)

     gradRhoXY2 = 4.D0*( gradRhoXY(i ,jn) - gradRhoXY(i ,js) )           &
                + gradRhoXY(ie,jn) - gradRhoXY(iw,js) + gradRhoXY(iw,jn) &
                - gradRhoXY(ie,js)

     gradRhoYX2 = 4.D0*( gradRhoXY(ie,j ) - gradRhoXY(iw,j ) )           &
                + gradRhoXY(ie,jn) - gradRhoXY(iw,js) + gradRhoXY(ie,js) &
                - gradRhoXY(iw,jn)

     gradRhoYY2 = 4.D0*( gradRhoYY(i ,jn) - gradRhoYY(i ,js) )           &
                + gradRhoYY(ie,jn) - gradRhoYY(iw,js) + gradRhoYY(iw,jn) &
                - gradRhoYY(ie,js)

!--------- Directional derivatives ---------------------------------------------
! Gradient of the density rho
     gradRhoD(1) = 0.5D0*( rho(ie,j)  - rho(iw,j) )
     gradRhoD(2) = 0.5D0*( rho(i,jn)  - rho(i,js) )
     gradRhoD(3) = 0.5D0*( rho(iw,j)  - rho(ie,j) )
     gradRhoD(4) = 0.5D0*( rho(i,js)  - rho(i,jn) )
     gradRhoD(5) = 0.5D0*( rho(ie,jn) - rho(iw,js) )
     gradRhoD(6) = 0.5D0*( rho(iw,jn) - rho(ie,js) )
     gradRhoD(7) = 0.5D0*( rho(iw,js) - rho(ie,jn) )
     gradRhoD(8) = 0.5D0*( rho(ie,js) - rho(iw,jn) )

! Gradient of the chemical potential psi
     gradPsiD(1) = 0.5D0*( psi(ie,j)  - psi(iw,j) )
     gradPsiD(2) = 0.5D0*( psi(i,jn)  - psi(i,js) )
     gradPsiD(3) = 0.5D0*( psi(iw,j)  - psi(ie,j) )
     gradPsiD(4) = 0.5D0*( psi(i,js)  - psi(i,jn) )
     gradPsiD(5) = 0.5D0*( psi(ie,jn) - psi(iw,js) )
     gradPsiD(6) = 0.5D0*( psi(iw,jn) - psi(ie,js) )
     gradPsiD(7) = 0.5D0*( psi(iw,js) - psi(ie,jn) )
     gradPsiD(8) = 0.5D0*( psi(ie,js) - psi(iw,jn) )

!--------- Equilibrium values of Gamma, and the distribution functions ---------
     GamW(0) = W0C*U_sq
     Gam(0)  = W0 - GamW(0)
     feq(0)  = Gam(0)*rhon
     geq(0)  = W0C*p(i,j) - rhon*GamW(0)

     GamW(1) = W1C*(ux + 1.5D0*ux*ux - U_sq )
     Gam(1)  = W1 + GamW(1)
     feq(1)  = Gam(1)*rhon
     geq(1)  = W1C*p(i,j) + rhon*GamW(1)

     GamW(2) = W1C*( uy + 1.5D0*uy*uy - U_sq )
     Gam(2)  = W1 + GamW(2)
     feq(2)  = Gam(2)*rhon
     geq(2)  = W1C*p(i,j) + rhon*GamW(2)

     GamW(3) = W1C*( -ux + 1.5D0*ux*ux - U_sq )
     Gam(3)  = W1 + GamW(3)
     feq(3)  = Gam(3)*rhon
     geq(3)  = W1C*p(i,j) + rhon*GamW(3)

     GamW(4) = W1C*( -uy + 1.5D0*uy*uy - U_sq )
     Gam(4)  = W1 + GamW(4)
     feq(4)  = Gam(4)*rhon
     geq(4)  = W1C*p(i,j) + rhon*GamW(4)

     v       = ux + uy
     GamW(5) = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(5)  = W2 + GamW(5)
     feq(5)  = Gam(5)*rhon
     geq(5)  = W2C*p(i,j) + rhon*GamW(5)

     v      = -ux + uy
     GamW(6) = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(6)  = W2 + GamW(6)
     feq(6)  = Gam(6)*rhon
     geq(6)  = W2C*p(i,j) + rhon*GamW(6)

     v       = -ux - uy
     GamW(7) = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(7)  = W2 + GamW(7)
     feq(7)  = Gam(7)*rhon
     geq(7)  = W2C*p(i,j) + rhon*GamW(7)

     v       = ux - uy
     GamW(8) = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(8)  = W2 + GamW(8)
     feq(8)  = Gam(8)*rhon
     geq(8)  = W2C*p(i,j) + rhon*GamW(8)

! Define directional terms that are common to all (f,g) components
     gradRhoU = ux*gradRhoX(i,j) + uy*gradRhoY(i,j)
     gradPsiU = ux*gradPsiX      + uy*gradPsiY
     gradRhoPsi = gradRhoU - invCsRho*gradPsiU
     GX = kappa_12*( gradRhoSqX - gradRhoXX2 - gradRhoXY2 )
     GY = kappa_12*( gradRhoSqY - gradRhoYX2 - gradRhoYY2 )
     UG = ux*GX + uy*GY


!--------- Collision and stream for the f and g distribution functions ---------
     Ff = Gam(0)*gradRhoPsi
     Fg = GamW(0)*gradRhoU + invCs_sq*Gam(0)*UG
     f(i,j,0) = fbar(i,j,0) + invTau*( feq(0) - fbar(i,j,0) ) - invTau2*Ff
     g(i,j,0) = gbar(i,j,0) + invTau*( geq(0) - gbar(i,j,0) ) - invTau2*Fg

! Direction 1
     Ff = Gam(1)*( gradRhoD(1) - invCsRho*gradPsiD(1) - gradRhoPsi )
     Fg = GamW(1)*( gradRhoD(1) - gradRhoU ) + invCs_sq*Gam(1)*( GX - UG )
     f(i,j,1) = fbar(i,j,1) + invTau*( feq(1) - fbar(i,j,1) ) + invTau2*Ff
     g(i,j,1) = gbar(i,j,1) + invTau*( geq(1) - gbar(i,j,1) ) + invTau2*Fg

! Direction 2
     Ff = Gam(2)*( gradRhoD(2) - invCsRho*gradPsiD(2) - gradRhoPsi )
     Fg = GamW(2)*( gradRhoD(2) - gradRhoU ) + invCs_sq*Gam(2)*( GY - UG )
     f(i,j,2) = fbar(i,j,2) + invTau*( feq(2) - fbar(i,j,2) ) + invTau2*Ff
     g(i,j,2) = gbar(i,j,2) + invTau*( geq(2) - gbar(i,j,2) ) + invTau2*Fg

! Direction 3
     Ff = Gam(3)*( gradRhoD(3) - invCsRho*gradPsiD(3) - gradRhoPsi )
     Fg = GamW(3)*( gradRhoD(3) - gradRhoU ) - invCs_sq*Gam(3)*( GX + UG )
     f(i,j,3) = fbar(i,j,3) + invTau*( feq(3) - fbar(i,j,3) ) + invTau2*Ff
     g(i,j,3) = gbar(i,j,3) + invTau*( geq(3) - gbar(i,j,3) ) + invTau2*Fg

! Direction 4
     Ff = Gam(4)*( gradRhoD(4) - invCsRho*gradPsiD(4) - gradRhoPsi )
     Fg = GamW(4)*( gradRhoD(4) - gradRhoU ) - invCs_sq*Gam(4)*( GY + UG )
     f(i,j,4) = fbar(i,j,4) + invTau*( feq(4) - fbar(i,j,4) ) + invTau2*Ff
     g(i,j,4) = gbar(i,j,4) + invTau*( geq(4) - gbar(i,j,4) ) + invTau2*Fg

! Direction 5
     Ff = Gam(5)*( gradRhoD(5) - invCsRho*gradPsiD(5) - gradRhoPsi )
     Fg = GamW(5)*( gradRhoD(5) - gradRhoU ) + invCs_sq*Gam(5)*( GX + GY - UG )
     f(i,j,5) = fbar(i,j,5) + invTau*( feq(5) - fbar(i,j,5) ) + invTau2*Ff
     g(i,j,5) = gbar(i,j,5) + invTau*( geq(5) - gbar(i,j,5) ) + invTau2*Fg

! Direction 6
     Ff = Gam(6)*( gradRhoD(6) - invCsRho*gradPsiD(6) - gradRhoPsi )
     Fg = GamW(6)*( gradRhoD(6) - gradRhoU ) - invCs_sq*Gam(6)*( GX - GY + UG )
     f(i,j,6) = fbar(i,j,6) + invTau*( feq(6) - fbar(i,j,6) ) + invTau2*Ff
     g(i,j,6) = gbar(i,j,6) + invTau*( geq(6) - gbar(i,j,6) ) + invTau2*Fg

! Direction 7
     Ff = Gam(7)*( gradRhoD(7) - invCsRho*gradPsiD(7) - gradRhoPsi )
     Fg = GamW(7)*( gradRhoD(7) - gradRhoU ) - invCs_sq*Gam(7)*( GX + GY + UG )
     f(i,j,7) = fbar(i,j,7) + invTau*( feq(7) - fbar(i,j,7) ) + invTau2*Ff
     g(i,j,7) = gbar(i,j,7) + invTau*( geq(7) - gbar(i,j,7) ) + invTau2*Fg

! Direction 8
     Ff = Gam(8)*( gradRhoD(8) - invCsRho*gradPsiD(8) - gradRhoPsi )
     Fg = GamW(8)*( gradRhoD(8) - gradRhoU ) + invCs_sq*Gam(8)*( GX - GY - UG )
     f(i,j,8) = fbar(i,j,8) + invTau*( feq(8) - fbar(i,j,8) ) + invTau2*Ff
     g(i,j,8) = gbar(i,j,8) + invTau*( geq(8) - gbar(i,j,8) ) + invTau2*Fg

   END DO
 END DO

 RETURN
 END SUBROUTINE PostStream
