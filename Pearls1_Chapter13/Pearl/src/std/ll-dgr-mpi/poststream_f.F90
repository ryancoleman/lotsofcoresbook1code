!-------------------------------------------------------------------------------
! Subroutine : PostStreamF
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Post-Stream step for distribution function f
!> @details
!! Post-Stream step for order parameter distribution function f in the parallel
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

 SUBROUTINE PostStreamF

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js
 REAL(KIND = DBL) :: Fs, ux, uy, U_sq, v
 REAL(KIND = DBL) :: aux, rhon, tau, invTau, invTau2, invCsRho
 REAL(KIND = DBL) :: gradPsiX_f, gradPsiY_f, gradRhoU, gradPsiU, gradRhoPsi
 REAL(KIND = DBL), DIMENSION(0:8) :: feq, Gam
 REAL(KIND = DBL), DIMENSION(1:8) :: gradPsiD, gradRhoD


 DO j = yl_f, yu_f
   DO i = xl_f, xu_f

! Define some local values
     rhon = rho_f(i,j)
     ux   = u_f(i,j,1)
     uy   = u_f(i,j,2)
     U_sq = 0.5D0*( ux*ux + uy*uy )
     tau  = tauL + (rhon - rhoL)*tauRhoStar
     invTau   = 1.D0/(1.D0 + 2.D0*tau)
     invTau2  = tau*invTau
     invCsRho = invCs_sq*rho_f(i,j)

! Identify neighboring cells
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

!--------- Standard Derivatives ------------------------------------------------
! Gradient of the chemical potential psi
     gradPsiX_f = ( 4.D0*( psi_f(ie,j ) - psi_f(iw,j ) ) + psi_f(ie,jn) &
                - psi_f(iw,js) + psi_f(ie,js) - psi_f(iw,jn) )*inv12

     gradPsiY_f = ( 4.D0*( psi_f(i ,jn) - psi_f(i ,js) ) + psi_f(ie,jn) &
                - psi_f(iw,js) + psi_f(iw,jn) - psi_f(ie,js) )*inv12

!--------- Directional derivatives ---------------------------------------------
! Gradient of the density rho
     gradRhoD(1) = 0.5D0*( rho_f(ie,j)  - rho_f(iw,j) )
     gradRhoD(2) = 0.5D0*( rho_f(i,jn)  - rho_f(i,js) )
     gradRhoD(3) = 0.5D0*( rho_f(iw,j)  - rho_f(ie,j) )
     gradRhoD(4) = 0.5D0*( rho_f(i,js)  - rho_f(i,jn) )
     gradRhoD(5) = 0.5D0*( rho_f(ie,jn) - rho_f(iw,js) )
     gradRhoD(6) = 0.5D0*( rho_f(iw,jn) - rho_f(ie,js) )
     gradRhoD(7) = 0.5D0*( rho_f(iw,js) - rho_f(ie,jn) )
     gradRhoD(8) = 0.5D0*( rho_f(ie,js) - rho_f(iw,jn) )

! Gradient of the chemical potential psi
     gradPsiD(1) = 0.5D0*( psi_f(ie,j)  - psi_f(iw,j) )
     gradPsiD(2) = 0.5D0*( psi_f(i,jn)  - psi_f(i,js) )
     gradPsiD(3) = 0.5D0*( psi_f(iw,j)  - psi_f(ie,j) )
     gradPsiD(4) = 0.5D0*( psi_f(i,js)  - psi_f(i,jn) )
     gradPsiD(5) = 0.5D0*( psi_f(ie,jn) - psi_f(iw,js) )
     gradPsiD(6) = 0.5D0*( psi_f(iw,jn) - psi_f(ie,js) )
     gradPsiD(7) = 0.5D0*( psi_f(iw,js) - psi_f(ie,jn) )
     gradPsiD(8) = 0.5D0*( psi_f(ie,js) - psi_f(iw,jn) )

!--------- Equilibrium values of Gamma, and the distribution function f --------
     aux    = W0C*U_sq
     Gam(0) = W0 - aux
     feq(0) = ( W0 - aux )*rhon

     aux    = W1C*(ux + 1.5D0*ux*ux - U_sq )
     Gam(1) = W1 + aux 
     feq(1) = ( W1 + aux )*rhon

     aux    = W1C*( uy + 1.5D0*uy*uy - U_sq )
     Gam(2) = W1 + aux
     feq(2) = ( W1 + aux )*rhon

     aux    = W1C*( -ux + 1.5D0*ux*ux - U_sq )
     Gam(3) = W1 + aux
     feq(3) = ( W1 + aux )*rhon

     aux    = W1C*( -uy + 1.5D0*uy*uy - U_sq )
     Gam(4) = W1 + aux
     feq(4) = ( W1 + aux )*rhon

     v      = ux + uy
     aux    = W2C*(v + 1.5D0*v*v - U_sq )
     Gam(5) = W2 + aux
     feq(5) = ( W2 + aux )*rhon

     v      = -ux + uy
     aux    = W2C*(v + 1.5D0*v*v - U_sq )
     Gam(6) = W2 + aux
     feq(6) = ( W2 + aux )*rhon

     v      = -ux - uy
     aux    = W2C*(v + 1.5D0*v*v - U_sq )
     Gam(7) = W2 + aux
     feq(7) = ( W2 + aux )*rhon

     v      = ux - uy
     aux    = W2C*(v + 1.5D0*v*v - U_sq )
     Gam(8) = W2 + aux
     feq(8) = ( W2 + aux )*rhon

! Directional terms common to all f components
     gradRhoU   = ux*gradRhoX_f(i,j) + uy*gradRhoY_f(i,j)
     gradPsiU   = ux*gradPsiX_f      + uy*gradPsiY_f
     gradRhoPsi = gradRhoU - invCsRho*gradPsiU

!--------- Collision and stream for the distribution function f ----------------
     Fs = Gam(0)*gradRhoPsi
     f(i,j,0) = fbar(i,j,0) + invTau*( feq(0) - fbar(i,j,0) ) - invTau2*Fs

! Direction 1
     Fs = Gam(1)*( gradRhoD(1) - invCsRho*gradPsiD(1) - gradRhoPsi )
     f(i,j,1) = fbar(i,j,1) + invTau*( feq(1) - fbar(i,j,1) ) + invTau2*Fs

! Direction 2
     Fs = Gam(2)*( gradRhoD(2) - invCsRho*gradPsiD(2) - gradRhoPsi )
     f(i,j,2) = fbar(i,j,2) + invTau*( feq(2) - fbar(i,j,2) ) + invTau2*Fs

! Direction 3
     Fs = Gam(3)*( gradRhoD(3) - invCsRho*gradPsiD(3) - gradRhoPsi )
     f(i,j,3) = fbar(i,j,3) + invTau*( feq(3) - fbar(i,j,3) ) + invTau2*Fs

! Direction 4
     Fs = Gam(4)*( gradRhoD(4) - invCsRho*gradPsiD(4) - gradRhoPsi )
     f(i,j,4) = fbar(i,j,4) + invTau*( feq(4) - fbar(i,j,4) ) + invTau2*Fs

! Direction 5
     Fs = Gam(5)*( gradRhoD(5) - invCsRho*gradPsiD(5) - gradRhoPsi )
     f(i,j,5) = fbar(i,j,5) + invTau*( feq(5) - fbar(i,j,5) ) + invTau2*Fs

! Direction 6
     Fs = Gam(6)*( gradRhoD(6) - invCsRho*gradPsiD(6) - gradRhoPsi )
     f(i,j,6) = fbar(i,j,6) + invTau*( feq(6) - fbar(i,j,6) ) + invTau2*Fs

! Direction 7
     Fs = Gam(7)*( gradRhoD(7) - invCsRho*gradPsiD(7) - gradRhoPsi )
     f(i,j,7) = fbar(i,j,7) + invTau*( feq(7) - fbar(i,j,7) ) + invTau2*Fs

! Direction 8
     Fs = Gam(8)*( gradRhoD(8) - invCsRho*gradPsiD(8) - gradRhoPsi )
     f(i,j,8) = fbar(i,j,8) + invTau*( feq(8) - fbar(i,j,8) ) + invTau2*Fs

   END DO
 END DO

 RETURN
 END SUBROUTINE PostStreamF
