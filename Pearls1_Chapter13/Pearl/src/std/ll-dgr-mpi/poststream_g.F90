!-------------------------------------------------------------------------------
! Subroutine : PostStreamG
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Post-Stream step for distribution function g
!> @details
!! Post-Stream step for the pressure distribution function g in the parallel
!! dual grid D2Q9 Lee-Lin multiphase LBM.
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

 SUBROUTINE PostStreamG

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js
 REAL(KIND = DBL) :: Fs, ux, uy, UG, U_sq, v
 REAL(KIND = DBL) :: rhon, tau, invTau, invTau2
 REAL(KIND = DBL) :: gradRhoU
 REAL(KIND = DBL), DIMENSION(0:8) :: geq, Gam, GamW
 REAL(KIND = DBL), DIMENSION(1:8) :: gradRhoD


 DO j = yl, yu
   DO i = xl, xu

! Define some local values
     rhon = rho(i,j)
     ux   = u(i,j,1)
     uy   = u(i,j,2)
     U_sq = 0.5D0*( ux*ux + uy*uy )
     tau  = tauL + (rhon - rhoL)*tauRhoStar
     invTau  = 1.D0/(1.D0 + 2.D0*tau)
     invTau2 = tau*invTau

! Identify neighboring cells
     ie = ni(i,j,1)
     jn = ni(i,j,2)
     iw = ni(i,j,3)
     js = ni(i,j,4)

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

!--------- Equilibrium values of Gamma, and the distribution function g --------
! This definition of Gamma (Gam) is different from Lee-Lin because we embed the
! product by invCs_sq (needed in the collision term) in the definition itself
     GamW(0) = W0C*U_sq
     Gam(0)  = invCs_sq*( W0 - GamW(0) )
     geq(0)  = W0C*p(i,j) - rhon*GamW(0)

     GamW(1) = W1C*(ux + 1.5D0*ux*ux - U_sq )
     Gam(1)  = invCs_sq*( W1 + GamW(1) )
     geq(1)  = W1C*p(i,j) + rhon*GamW(1)

     GamW(2) = W1C*( uy + 1.5D0*uy*uy - U_sq )
     Gam(2)  = invCs_sq*( W1 + GamW(2) )
     geq(2)  = W1C*p(i,j) + rhon*GamW(2)

     GamW(3) = W1C*( -ux + 1.5D0*ux*ux - U_sq )
     Gam(3)  = invCs_sq*( W1 + GamW(3) )
     geq(3)  = W1C*p(i,j) + rhon*GamW(3)

     GamW(4) = W1C*( -uy + 1.5D0*uy*uy - U_sq )
     Gam(4)  = invCs_sq*( W1 + GamW(4) )
     geq(4)  = W1C*p(i,j) + rhon*GamW(4)

     v       = ux + uy
     GamW(5) = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(5)  = invCs_sq*( W2 + GamW(5) )
     geq(5)  = W2C*p(i,j) + rhon*GamW(5)

     v      = -ux + uy
     GamW(6) = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(6)  = invCs_sq*( W2 + GamW(6) )
     geq(6)  = W2C*p(i,j) + rhon*GamW(6)

     v       = -ux - uy
     GamW(7) = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(7)  = invCs_sq*( W2 + GamW(7) )
     geq(7)  = W2C*p(i,j) + rhon*GamW(7)

     v       = ux - uy
     GamW(8) = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(8)  = invCs_sq*( W2 + GamW(8) )
     geq(8)  = W2C*p(i,j) + rhon*GamW(8)

! Define directional terms that go with "u" and are common to all g components
     gradRhoU = ux*gradRhoX(i,j) + uy*gradRhoY(i,j)
     UG       = ux*GX(i,j)       + uy*GY(i,j)

!--------- Collision and stream for the g distribution function ----------------
     Fs = GamW(0)*gradRhoU + Gam(0)*UG
     g(i,j,0) = gbar(i,j,0) + invTau*( geq(0) - gbar(i,j,0) ) - invTau2*Fs

! Direction 1
     Fs = GamW(1)*( gradRhoD(1) - gradRhoU ) + Gam(1)*( GX(i,j) - UG )
     g(i,j,1) = gbar(i,j,1) + invTau*( geq(1) - gbar(i,j,1) ) + invTau2*Fs

! Direction 2
     Fs = GamW(2)*( gradRhoD(2) - gradRhoU ) + Gam(2)*( GY(i,j) - UG )
     g(i,j,2) = gbar(i,j,2) + invTau*( geq(2) - gbar(i,j,2) ) + invTau2*Fs

! Direction 3
     Fs = GamW(3)*( gradRhoD(3) - gradRhoU ) - Gam(3)*( GX(i,j) + UG )
     g(i,j,3) = gbar(i,j,3) + invTau*( geq(3) - gbar(i,j,3) ) + invTau2*Fs

! Direction 4
     Fs = GamW(4)*( gradRhoD(4) - gradRhoU ) - Gam(4)*( GY(i,j) + UG )
     g(i,j,4) = gbar(i,j,4) + invTau*( geq(4) - gbar(i,j,4) ) + invTau2*Fs

! Direction 5
     Fs = GamW(5)*( gradRhoD(5) - gradRhoU ) + Gam(5)*( GX(i,j) + GY(i,j) - UG )
     g(i,j,5) = gbar(i,j,5) + invTau*( geq(5) - gbar(i,j,5) ) + invTau2*Fs

! Direction 6 
     Fs = GamW(6)*( gradRhoD(6) - gradRhoU ) - Gam(6)*( GX(i,j) - GY(i,j) + UG )
     g(i,j,6) = gbar(i,j,6) + invTau*( geq(6) - gbar(i,j,6) ) + invTau2*Fs

! Direction 7
     Fs = GamW(7)*( gradRhoD(7) - gradRhoU ) - Gam(7)*( GX(i,j) + GY(i,j) + UG )
     g(i,j,7) = gbar(i,j,7) + invTau*( geq(7) - gbar(i,j,7) ) + invTau2*Fs

! Direction 8
     Fs = GamW(8)*( gradRhoD(8) - gradRhoU ) + Gam(8)*( GX(i,j) - GY(i,j) - UG )
     g(i,j,8) = gbar(i,j,8) + invTau*( geq(8) - gbar(i,j,8) ) + invTau2*Fs

   END DO
 END DO


 RETURN
 END SUBROUTINE PostStreamG
