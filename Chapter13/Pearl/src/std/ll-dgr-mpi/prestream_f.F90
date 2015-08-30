!-------------------------------------------------------------------------------
! Subroutine : PreStreamF
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Pre-Stream and Stream steps for distribution function f
!> @details
!! Pre-Stream and Strean steps for order parameter distribution function f in
!! the parallel dual grid D2Q9 Lee-Lin multiphase LBM.
!!
!! The pre-stream step is done in f, which is then streamed into fbar.
!!
!! Directional differentials are evaluated using a second order switching
!! differential scheme (biased/central). Other differentials are evaluated using
!! a standard second order central differencing scheme.

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

 SUBROUTINE PreStreamF

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js, ie2, iw2, jn2, js2
 REAL(KIND = DBL) :: Fs, ux, uy, U_sq, v
 REAL(KIND = DBL) :: aux, rhon, tau, invtau, AX, rho3, psi3, invCsRho
 REAL(KIND = DBL) :: gradPsiX_f, gradPsiY_f, gradRhoU, gradPsiU, gradRhoPsi
 REAL(KIND = DBL), DIMENSION(0:8) :: feq, Gam
 REAL(KIND = DBL), DIMENSION(1:8) :: gradPsiD, gradRhoD


 DO j = yl_f, yu_f
   DO i = xl_f, xu_f

! Define some local values
     rhon = rho_f(i,j)
     rho3 = 3.D0*rho_f(i,j)
     psi3 = 3.D0*psi_f(i,j)
     ux   = u_f(i,j,1)
     uy   = u_f(i,j,2)
     U_sq = 0.5D0*( ux*ux + uy*uy )
     tau  = tauL + (rhon - rhoL)*tauRhoStar
     invtau   = 0.5D0/tau
     invCsRho = invCs_sq*rho_f(i,j)

! Identify first neighbors
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

! Identify second neighbors
     ie2 = ni_f(ie,j,1)
     jn2 = ni_f(i,jn,2)
     iw2 = ni_f(iw,j,3)
     js2 = ni_f(i,js,4)

!--------- Standard Derivatives ------------------------------------------------
! Gradient of the chemical potential psi
     gradPsiX_f = ( 4.D0*( psi_f(ie,j ) - psi_f(iw,j ) ) + psi_f(ie,jn) &
                - psi_f(iw,js) + psi_f(ie,js) - psi_f(iw,jn) )*inv12

     gradPsiY_f = ( 4.D0*( psi_f(i ,jn) - psi_f(i ,js) ) + psi_f(ie,jn) &
                - psi_f(iw,js) + psi_f(iw,jn) - psi_f(ie,js) )*inv12

!--------- Directional derivatives ( Second order mixed differencing scheme ) --
! Gradient of the density rho
     AX = rho_f(ie2,j) - rho_f(i,j)
     gradRhoD(1) = 0.5D0*( - rho_f(ie2,j  ) + 4.D0*rho_f(ie,j ) - rho3 )
     IF ( AX*gradRhoD(1) < 0.D0 ) gradRhoD(1) = 0.5D0*( rho_f(ie,j) - rho_f(iw,j) )

     AX = rho_f(i,jn2) - rho_f(i,j)
     gradRhoD(2) = 0.5D0*( - rho_f(i  ,jn2) + 4.D0*rho_f(i ,jn) - rho3 )
     IF ( AX*gradRhoD(2) < 0.D0 ) gradRhoD(2) = 0.5D0*( rho_f(i,jn) - rho_f(i,js) )

     AX = rho_f(iw2,j) - rho_f(i,j)
     gradRhoD(3) = 0.5D0*( - rho_f(iw2,j  ) + 4.D0*rho_f(iw,j ) - rho3 )
     IF ( AX*gradRhoD(3) < 0.D0 ) gradRhoD(3) = 0.5D0*( rho_f(iw,j) - rho_f(ie,j) )

     AX = rho_f(i,js2) - rho_f(i,j)
     gradRhoD(4) = 0.5D0*( - rho_f(i  ,js2) + 4.D0*rho_f(i ,js) - rho3 )
     IF ( AX*gradRhoD(4) < 0.D0 ) gradRhoD(4) = 0.5D0*( rho_f(i,js) - rho_f(i,jn) )

     AX = rho_f(ie2,jn2) - rho_f(i,j)
     gradRhoD(5) = 0.5D0*( - rho_f(ie2,jn2) + 4.D0*rho_f(ie,jn) - rho3 )
     IF ( AX*gradRhoD(5) < 0.D0 ) gradRhoD(5) = 0.5D0*( rho_f(ie,jn) - rho_f(iw,js) )

     AX = rho_f(iw2,jn2) - rho_f(i,j)
     gradRhoD(6) = 0.5D0*( - rho_f(iw2,jn2) + 4.D0*rho_f(iw,jn) - rho3 )
     IF ( AX*gradRhoD(6) < 0.D0 ) gradRhoD(6) = 0.5D0*( rho_f(iw,jn) - rho_f(ie,js) )

     AX = rho_f(iw2,js2) - rho_f(i,j)
     gradRhoD(7) = 0.5D0*( - rho_f(iw2,js2) + 4.D0*rho_f(iw,js) - rho3 )
     IF ( AX*gradRhoD(7) < 0.D0 ) gradRhoD(7) = 0.5D0*( rho_f(iw,js) - rho_f(ie,jn) )

     AX = rho_f(ie2,js2) - rho_f(i,j)
     gradRhoD(8) = 0.5D0*( - rho_f(ie2,js2) + 4.D0*rho_f(ie,js) - rho3 )
     IF ( AX*gradRhoD(8) < 0.D0 ) gradRhoD(8) = 0.5D0*( rho_f(ie,js) - rho_f(iw,jn) )

! Gradient of the chemical potential psi
     AX = psi_f(ie2,j) - psi_f(i,j)
     gradPsiD(1) = 0.5D0*( - psi_f(ie2,j  ) + 4.D0*psi_f(ie,j ) - psi3 )
     IF ( AX*gradPsiD(1) < 0.D0 ) gradPsiD(1) = 0.5D0*( psi_f(ie,j) - psi_f(iw,j) )

     AX = psi_f(i,jn2) - psi_f(i,j)
     gradPsiD(2) = 0.5D0*( - psi_f(i  ,jn2) + 4.D0*psi_f(i ,jn) - psi3 )
     IF ( AX*gradPsiD(2) < 0.D0 ) gradPsiD(2) = 0.5D0*( psi_f(i,jn) - psi_f(i,js) )

     AX = psi_f(iw2,j) - psi_f(i,j)
     gradPsiD(3) = 0.5D0*( - psi_f(iw2,j  ) + 4.D0*psi_f(iw,j ) - psi3 )
     IF ( AX*gradPsiD(3) < 0.D0 ) gradPsiD(3) = 0.5D0*( psi_f(iw,j) - psi_f(ie,j) )

     AX = psi_f(i,js2) - psi_f(i,j)
     gradPsiD(4) = 0.5D0*( - psi_f(i  ,js2) + 4.D0*psi_f(i ,js) - psi3 )
     IF ( AX*gradPsiD(4) < 0.D0 ) gradPsiD(4) = 0.5D0*( psi_f(i,js) - psi_f(i,jn) )

     AX = psi_f(ie2,jn2) - psi_f(i,j)
     gradPsiD(5) = 0.5D0*( - psi_f(ie2,jn2) + 4.D0*psi_f(ie,jn) - psi3 )
     IF ( AX*gradPsiD(5) < 0.D0 ) gradPsiD(5) = 0.5D0*( psi_f(ie,jn) - psi_f(iw,js) )

     AX = psi_f(iw2,jn2) - psi_f(i,j)
     gradPsiD(6) = 0.5D0*( - psi_f(iw2,jn2) + 4.D0*psi_f(iw,jn) - psi3 )
     IF ( AX*gradPsiD(6) < 0.D0 ) gradPsiD(6) = 0.5D0*( psi_f(iw,jn) - psi_f(ie,js) )

     AX = psi_f(iw2,js2) - psi_f(i,j)
     gradPsiD(7) = 0.5D0*( - psi_f(iw2,js2) + 4.D0*psi_f(iw,js) - psi3 )
     IF ( AX*gradPsiD(7) < 0.D0 ) gradPsiD(7) = 0.5D0*( psi_f(iw,js) - psi_f(ie,jn) )

     AX = psi_f(ie2,js2) - psi_f(i,j)
     gradPsiD(8) = 0.5D0*( - psi_f(ie2,js2) + 4.D0*psi_f(ie,js) - psi3 )
     IF ( AX*gradPsiD(8) < 0.D0 ) gradPsiD(8) = 0.5D0*( psi_f(ie,js) - psi_f(iw,jn) )

!--------- Equilibrium values of Gamma, and the distribution function f --------
     aux    =  W0C*U_sq
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
     aux    = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(5) = W2 + aux
     feq(5) = ( W2 + aux )*rhon

     v      = -ux + uy
     aux    = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(6) = W2 + aux
     feq(6) = ( W2 + aux )*rhon

     v      = -ux - uy
     aux    = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(7) = W2 + aux
     feq(7) = ( W2 + aux )*rhon

     v      = ux - uy
     aux    = W2C*( v + 1.5D0*v*v - U_sq )
     Gam(8) = W2 + aux
     feq(8) = ( W2 + aux )*rhon

! Define directional terms that go with "u" and are common to all f components
     gradRhoU = ux*gradRhoX_f(i,j) + uy*gradRhoY_f(i,j)
     gradPsiU = ux*gradPsiX_f      + uy*gradPsiY_f
     gradRhoPsi = gradRhoU - invCsRho*gradPsiU

!--------- Collision and stream for the distribution function f ----------------
! Direction 0
     Fs = Gam(0)*gradRhoPsi
     fbar(i ,j ,0) = f(i,j,0) + invtau*( feq(0) - f(i,j,0) ) - 0.5D0*Fs

! Direction 1
     Fs = Gam(1)*( gradRhoD(1) - invCsRho*gradPsiD(1) - gradRhoPsi )
     fbar(ie,j ,1) = f(i,j,1) + invtau*( feq(1) - f(i,j,1) ) + 0.5D0*Fs

! Direction 2
     Fs = Gam(2)*( gradRhoD(2) - invCsRho*gradPsiD(2) - gradRhoPsi )
     fbar(i ,jn,2) = f(i,j,2) + invtau*( feq(2) - f(i,j,2) ) + 0.5D0*Fs

! Direction 3
     Fs = Gam(3)*( gradRhoD(3) - invCsRho*gradPsiD(3) - gradRhoPsi )
     fbar(iw,j ,3) = f(i,j,3) + invtau*( feq(3) - f(i,j,3) ) + 0.5D0*Fs

! Direction 4
     Fs = Gam(4)*( gradRhoD(4) - invCsRho*gradPsiD(4) - gradRhoPsi )
     fbar(i ,js,4) = f(i,j,4) + invtau*( feq(4) - f(i,j,4) ) + 0.5D0*Fs

! Direction 5
     Fs = Gam(5)*( gradRhoD(5) - invCsRho*gradPsiD(5) - gradRhoPsi )
     fbar(ie,jn,5) = f(i,j,5) + invtau*( feq(5) - f(i,j,5) ) + 0.5D0*Fs

! Direction 6
     Fs = Gam(6)*( gradRhoD(6) - invCsRho*gradPsiD(6) - gradRhoPsi )
     fbar(iw,jn,6) = f(i,j,6) + invtau*( feq(6) - f(i,j,6) ) + 0.5D0*Fs 

! Direction 7
     Fs = Gam(7)*( gradRhoD(7) - invCsRho*gradPsiD(7) - gradRhoPsi )
     fbar(iw,js,7) = f(i,j,7) + invtau*( feq(7) - f(i,j,7) ) + 0.5D0*Fs

! Direction 8
     Fs = Gam(8)*( gradRhoD(8) - invCsRho*gradPsiD(8) - gradRhoPsi )
     fbar(ie,js,8) = f(i,j,8) + invtau*( feq(8) - f(i,j,8) ) + 0.5D0*Fs

   END DO
 END DO

 RETURN
 END SUBROUTINE PreStreamF
