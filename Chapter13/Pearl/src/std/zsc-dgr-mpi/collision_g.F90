!-------------------------------------------------------------------------------
! Subroutine : CollisionG
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Collision and streaming steps for distribution function g
!> @details
!! Collision and streaming steps for the momentum distribution function g and
!! calculation of updated macroscopic quantities (velocity, pressure and order
!! parameter) in the dual grid parallel D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.

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

 SUBROUTINE CollisionG

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain,    ONLY : ni, now, nxt, xl, xu, yl, yu
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js, xs, ys
 REAL(KIND = DBL) :: sFx, sFy, iFs
 REAL(KIND = DBL) :: Vg, Vsq, Vg2, ux, uy
 REAL(KIND = DBL) :: rhon, invrho, phin, phin2, muPhi
 REAL(KIND = DBL) :: gradPhiX, gradPhiY, gradPhiSq, lapPhi
 REAL(KIND = DBL) :: Ag0, Ag1, Eg0A, Eg1A, Eg2A, Eg0R, Eg1R, Eg2R, geq


!---------- Hydrodynamics, collision and stream for g in the same loop ---------
 DO j = yl, yu
   DO i = xl, xu

! Define near neighbors for the streaming step
     ie = ni(i,j,1)
     jn = ni(i,j,2)
     iw = ni(i,j,3)
     js = ni(i,j,4)

! Define source nodes in the order parameter mesh
     xs = 2*i - 1
     ys = 2*j - 1

! Copy and re-scale the values of the order parameter and its differentials
     phin      = phi_f(xs,ys)
     phin2     = phin*phin
     lapPhi    = 4.D0*lapPhi_f(xs,ys)
     gradPhiX  = 2.D0*gradPhiX_f(xs,ys)
     gradPhiY  = 2.D0*gradPhiY_f(xs,ys)
     gradPhiSq = gradPhiX*gradPhiX + gradPhiY*gradPhiY

! Local value of the chemical potential
     muPhi = alpha4*phin*( phin2 - phistar2 ) - kappaG*lapPhi

!---------- Hydrodynamics ------------------------------------------------------
!  Local density value
     rhon   = g(i,j,0,now) + g(i,j,1,now) + g(i,j,2,now) + g(i,j,3,now) &
            + g(i,j,4,now) + g(i,j,5,now) + g(i,j,6,now) + g(i,j,7,now) &
            + g(i,j,8,now)
     invRho = 1.D0/rhon

! Interfacial force
     sFx = muPhi*gradPhiX
     sFy = muPhi*gradPhiY

! Velocity field at each node
     ux = ( g(i,j,1,now) - g(i,j,3,now) + g(i,j,5,now) - g(i,j,6,now) &
        - g(i,j,7,now) + g(i,j,8,now) + 0.5D0*sFx )*invRho

     uy = ( g(i,j,2,now) - g(i,j,4,now) + g(i,j,5,now) + g(i,j,6,now) &
        - g(i,j,7,now) - g(i,j,8,now) + 0.5D0*sFy )*invRho

     u(i,j,1) = ux
     u(i,j,2) = uy

! Pressure at each node
     p(i,j) = alpha*( phin2*( 3.D0*phin2 - 2.D0*phiStar2 ) - phiStar4 ) &
            - kappaG*( phin*lapPhi + 0.5D0*gradPhiSq ) + Cs_sq*rhon

!---------- Collision and Stream for momentum distribution function ------------
     Ag1 = rhon + 3.D0*phin*muPhi
     Ag0 = 2.25D0*rhon - 1.25D0*Ag1
     Eg0A = Eg0*Ag0
     Eg1A = Eg1*Ag1
     Eg2A = Eg2*Ag1
     Eg0R = Eg0C*rhon
     Eg1R = Eg1C*rhon
     Eg2R = Eg2C*rhon
     Vsq  = 0.5D0*( ux*ux + uy*uy )

! Direction 0
     geq = Eg0A - Eg0R*Vsq
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
     Vg = ux + uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (1.D0 - ux) + Vg2 )*sFx + ( (1.D0 - uy) + Vg2 )*sFy )
     g(ie,jn,5,nxt) = g(i,j,5,now) + invTauRho*( geq - g(i,j,5,now) ) + iFs

! Direction 6
     Vg = -ux + uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (-1.D0 - ux) - Vg2 )*sFx + ( (1.D0 - uy) + Vg2 )*sFy )
     g(iw,jn,6,nxt) = g(i,j,6,now) + invTauRho*( geq - g(i,j,6,now) ) + iFs

! Direction 7
     Vg = -ux - uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (-1.D0 - ux) - Vg2 )*sFx + ( (-1.D0 - uy) - Vg2 )*sFy )
     g(iw,js,7,nxt) = g(i,j,7,now) + invTauRho*( geq - g(i,j,7,now) ) + iFs

! Direction 8
     Vg = ux - uy
     Vg2 = invCs_sq*Vg
     geq = Eg2A + Eg2R*( Vg + 1.5D0*Vg*Vg - Vsq )
     iFs = Eg2T*( ( (1.D0 - ux) + Vg2 )*sFx + ( (-1.D0 - uy) - Vg2 )*sFy)
     g(ie,js,8,nxt) = g(i,j,8,now) + invTauRho*( geq - g(i,j,8,now) ) + iFs

   END DO
 END DO

 RETURN
 END SUBROUTINE CollisionG

