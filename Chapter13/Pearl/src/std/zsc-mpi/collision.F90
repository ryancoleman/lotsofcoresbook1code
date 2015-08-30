!-------------------------------------------------------------------------------
! Subroutine : Collision
! Revision   : 1.0 (2008-06-15)
! Author     : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @file
!! Collision (f,g), streaming (g), and hydrodynamics (p,phi,rho,u) calculation
!> @details
!! Collision step for order parameter distribution function f, collision and
!! streaming steps for momentum distribution function g, and calculation of
!! updated macroscopic quantities (velocity, pressure and order parameter) in
!! the parallel D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.

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

 SUBROUTINE Collision

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain,    ONLY : inv12, inv6, ni, now, nxt, xl, xu, yl, yu
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js
 REAL(KIND = DBL) :: sFx, sFy, iFs
 REAL(KIND = DBL) :: Vg, Vsq, Vg2, ux, uy
 REAL(KIND = DBL) :: rhon, invRho, phin, phin2, muPhin
 REAL(KIND = DBL) :: gradPhiX, gradPhiY, gradPhiSq, lapPhi
 REAL(KIND = DBL) :: Af0, Af1, Cfp
 REAL(KIND = DBL) :: Ag0, Ag1, Eg1A, Eg2A, Eg1R, Eg2R, geq

!--------- Hydrodynamic fields and collision in the same loop ------------------
 DO j = yl, yu
   DO i = xl, xu

     phin  = phi(i,j)
     phin2 = phin*phin 

     ie = ni(i,j,1)
     jn = ni(i,j,2)
     iw = ni(i,j,3)
     js = ni(i,j,4)

! Laplacian
     lapPhi = ( phi(ie,jn) + phi(ie,js) + phi(iw,jn) + phi(iw,js)  &
            + 4.D0*(phi(ie,j) + phi(iw,j) + phi(i,jn) + phi(i,js)) &
            - 20.D0*phin )*inv6

! Gradient of the order parameter
     gradPhiX = ( 4.D0*( phi(ie,j) - phi(iw,j) ) + phi(ie,jn) - phi(iw,js) &
              + phi(ie,js) - phi(iw,jn) )*inv12

     gradPhiY = ( 4.D0*( phi(i,jn) - phi(i,js) ) + phi(ie,jn) - phi(iw,js) &
              - phi(ie,js) + phi(iw,jn) )*inv12

     gradPhiSq = gradPhiX*gradPhiX + gradPhiY*gradPhiY

! Chemical potential
     muPhin = alpha4*phin*( phin2 - phistar2 ) - kappa*lapPhi

!--------- Hydrodynamics -------------------------------------------------------
!  Local density value
     rhon   = g(i,j,0,now) + g(i,j,1,now) + g(i,j,2,now) + g(i,j,3,now) &
            + g(i,j,4,now) + g(i,j,5,now) + g(i,j,6,now) + g(i,j,7,now) &
            + g(i,j,8,now)
     invrho = 1.D0/rhon

! Interfacial force
     sFx = muPhin*gradPhiX
     sFy = muPhin*gradPhiY

! Velocity field at each node
     u(i,j,1) = ( g(i,j,1,now) - g(i,j,3,now) + g(i,j,5,now) - g(i,j,6,now) &
              - g(i,j,7,now) + g(i,j,8,now) + 0.5D0*sFx )*invRho

     u(i,j,2) = ( g(i,j,2,now) - g(i,j,4,now) + g(i,j,5,now) + g(i,j,6,now) &
              - g(i,j,7,now) - g(i,j,8,now) + 0.5D0*sFy )*invRho

     ux = u(i,j,1)
     uy = u(i,j,2)

! Pressure at each node
     p(i,j) = alpha*( phin2*( 3.D0*phin2 - 2.D0*phiStar2 ) - phiStar4 ) &
            - kappa*( phin*lapPhi + 0.5D0*gradPhiSq ) + Cs_sq*rhon

!--------- Collision for order parameter distribution function -----------------
     Af1 = 0.5D0*Gamma*muPhin
     Af0 = -2.D0*Gamma*muPhin
     Cfp = invEta2*phin

     f(i,j,0,nxt) = f(i,j,0,now) + invTauPhi*(Af0 + phin   - f(i,j,0,now))
     f(i,j,1,now) = f(i,j,1,now) + invTauPhi*(Af1 + Cfp*ux - f(i,j,1,now))
     f(i,j,2,now) = f(i,j,2,now) + invTauPhi*(Af1 + Cfp*uy - f(i,j,2,now))
     f(i,j,3,now) = f(i,j,3,now) + invTauPhi*(Af1 - Cfp*ux - f(i,j,3,now))
     f(i,j,4,now) = f(i,j,4,now) + invTauPhi*(Af1 - Cfp*uy - f(i,j,4,now))

!--------- Collision and Stream for momentum distribution function -------------
     Ag1  = rhon + 3.D0*phin*muPhin
     Ag0  = 2.25D0*rhon - 1.25D0*Ag1
     Eg1A = Eg1*Ag1
     Eg2A = Eg2*Ag1
     Eg1R = Eg1C*rhon
     Eg2R = Eg2C*rhon
     Vsq  = 0.5D0*( ux*ux + uy*uy )

! Direction 0
     geq = Eg0*( Ag0 - rhon*3.D0*Vsq )
     iFs = Eg0T*( -ux*sFx - uy*sFy )
     g(i,j,0,nxt) = g(i,j,0,now) + invTauRho*( geq - g(i,j,0,now) ) + iFs

! Direction 1
     geq = Eg1A + Eg1R*( ux + 1.5D0*ux*ux - Vsq )
     iFs = Eg1T*( (1.D0 + 2.D0*ux)*sFx - uy*sFy )
     g(ie,j,1,nxt) = g(i,j,1,now) + invTauRho*( geq - g(i,j,1,now) ) + iFs

! Direction 2
     geq = Eg1A + Eg1R*( uy + 1.5D0*uy*uy - Vsq )
     iFs = Eg1T*( -ux*sFx + (1.D0 + 2.D0*uy)*sFy )
     g(i,jn,2,nxt) = g(i,j,2,now) + invTauRho*( geq - g(i,j,2,now) ) + iFs

! Direction 3
     geq = Eg1A + Eg1R*( -ux + 1.5D0*ux*ux - Vsq )
     iFs = Eg1T*( (2.D0*ux - 1.D0)*sFx - uy*sFy )
     g(iw,j,3,nxt) = g(i,j,3,now) + invTauRho*( geq - g(i,j,3,now) ) + iFs

! Direction 4
     geq = Eg1A + Eg1R*( -uy + 1.5D0*uy*uy - Vsq )
     iFs = Eg1T*( -ux*sFx + (2.D0*uy - 1)*sFy )
     g(i,js,4,nxt) = g(i,j,4,now) + invTauRho*( geq - g(i,j,4,now) ) + iFs

! Direction 5
     Vg  = ux + uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (1.D0 - ux) + Vg2 )*sFx + ( (1.D0 - uy) + Vg2 )*sFy )
     g(ie,jn,5,nxt) = g(i,j,5,now) + invTauRho*(geq - g(i,j,5,now)) + iFs

! Direction 6
     Vg  = -ux + uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (-1.D0 - ux) - Vg2 )*sFx + ( (1.D0 - uy) + Vg2 )*sFy )
     g(iw,jn,6,nxt) = g(i,j,6,now) + invTauRho*( geq - g(i,j,6,now) ) + iFs

! Direction 7
     Vg  = -ux - uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (-1.D0 - ux) - Vg2 )*sFx + ( (-1.D0 - uy) - Vg2 )*sFy )
     g(iw,js,7,nxt) = g(i,j,7,now) + invTauRho*( geq - g(i,j,7,now) ) + iFs

! Direction 8
     Vg  = ux - uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (1.D0 - ux) + Vg2 )*sFx + ( (-1.D0 - uy) - Vg2 )*sFy )
     g(ie,js,8,nxt) = g(i,j,8,now) + invTauRho*( geq - g(i,j,8,now) ) + iFs

   END DO
 END DO

 RETURN
 END SUBROUTINE Collision
