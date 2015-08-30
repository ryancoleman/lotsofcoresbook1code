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
!! the parallel D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM including gravity.

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

!  Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain,      ONLY : inv12, inv6, ni, now, nxt, xl, xu, yl, yu, zl, zu
 USE FluidParams
 USE LBMParams,   ONLY : Eg0n, Eg1n, Eg2n, EgC0n, EgC1n, EgC2n, invCs_sq, f, g
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, ie, iw, jn, js, kt, kb
 REAL(KIND = DBL) :: muPhin, phin, phin2
 REAL(KIND = DBL) :: rhon, invRhon
 REAL(KIND = DBL) :: sFx, sFy, sFz
 REAL(KIND = DBL) :: UF, ux, uy, uz, Vg, Vsq
 REAL(KIND = DBL) :: Af0, Af1, Ag0, Ag1, Cf1
 REAL(KIND = DBL) :: Eg1A, Eg2A, Eg1R, Eg2R
 REAL(KIND = DBL) :: geq1, geq2, gradPhiX, gradPhiY, gradPhiZ, lapPhi
 REAL(KIND = DBL) :: in01, in02, in03, in04, in05, in06, in07, in08, in09
 REAL(KIND = DBL) :: in10, in11, in12, in13, in14, in15, in16, in17, in18


! Order parameter update
 DO k = zl, zu
   DO j = yl, yu
     DO i = xl, xu
       phi(i,j,k) = SUM( f(i,j,k,:,now) )
     END DO
   END DO
 END DO

 DO k = zl, zu
   DO j = yl, yu
     DO i = xl, xu

! Define some local values
       phin    = phi(i,j,k)
       phin2   = phin*phin
       rhon    = SUM( g(i,j,k,:,now) )
       invRhon = 1.D0/rhon

!-------- Differential terms ---------------------------------------------------
! Identify neighbours
       ie = ni(i,j,k,1)
       iw = ni(i,j,k,2)
       jn = ni(i,j,k,3)
       js = ni(i,j,k,4)
       kt = ni(i,j,k,5)
       kb = ni(i,j,k,6)

! Nodal phi values
       in01 = phi(ie,j ,k )
       in02 = phi(iw,j ,k )
       in03 = phi(i ,jn,k )
       in04 = phi(i ,js,k )
       in05 = phi(i ,j ,kt)
       in06 = phi(i ,j ,kb)
       in07 = phi(ie,jn,k )
       in08 = phi(iw,js,k )
       in09 = phi(ie,js,k )
       in10 = phi(iw,jn,k )
       in11 = phi(ie,j ,kt)
       in12 = phi(iw,j ,kb)
       in13 = phi(ie,j ,kb)
       in14 = phi(iw,j ,kt)
       in15 = phi(i ,jn,kt)
       in16 = phi(i ,js,kb)
       in17 = phi(i ,jn,kb)
       in18 = phi(i ,js,kt)

! Laplacian of the order parameter phi
       lapPhi = ( in07 + in08 + in09 + in10 + in11 + in12 + in13 + in14 + in15 &
              + in16 + in17 + in18 + 2.0D0*( in01 + in02 + in03 + in04 + in05  &
              + in06 - 12.D0*phin ) )*inv6

! Gradient components of the order parameter phi
       gradPhiX = ( 2.0D0*( in01 - in02 ) + in07 - in08 + in09 - in10 + in11 &
                - in12 + in13 - in14 )*inv12

       gradPhiY = ( 2.0D0*( in03 - in04 ) + in07 - in08 + in10 - in09 + in15 &
                - in16 + in17 - in18 )*inv12

       gradPhiZ = ( 2.0D0*( in05 - in06 ) + in11 - in12 + in14 - in13 + in15 &
                - in16 + in18 - in17 )*inv12

!-------- Hydrodynamics (velocity, pressure, interfacial force) ----------------
! Chemical potential
       muphin = alpha4*phin*( phin2 - phiStar2 ) - kappa*lapPhi

! Interfacial force and gravity
       sFx = muPhin*gradPhiX
       sFy = muPhin*gradPhiY
       sFz = muPhin*gradPhiZ
       IF (phin > 0.D0) sFz = sFz + grav

! Velocity
       ux = ( g(i,j,k, 1,now) - g(i,j,k, 2,now) + g(i,j,k, 7,now) &
          -   g(i,j,k, 8,now) + g(i,j,k, 9,now) - g(i,j,k,10,now) &
          +   g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) &
          -   g(i,j,k,14,now) + 0.5D0*sFx )*invRhon

       uy = ( g(i,j,k, 3,now) - g(i,j,k, 4,now) + g(i,j,k, 7,now) &
          -   g(i,j,k, 8,now) - g(i,j,k, 9,now) + g(i,j,k,10,now) &
          +   g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) &
          - g(i,j,k,18,now) + 0.5D0*sFy )*invRhon

       uz = ( g(i,j,k, 5,now) - g(i,j,k, 6,now) + g(i,j,k,11,now) &
          -   g(i,j,k,12,now) - g(i,j,k,13,now) + g(i,j,k,14,now) &
          +   g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) &
          + g(i,j,k,18,now) + 0.5D0*sFz )*invRhon

!-------- Collision for order parameter distribution function f ----------------
       Af0 = -3.D0*Gamma*muPhin*invTauPhi
       Af1 = 0.5D0*Gamma*muPhin*invTauPhi
       Cf1 = invTauPhi*invEta2*phin

       f(i,j,k,0,nxt) = invTauPhi1*f(i,j,k,0,now) + Af0 + invTauPhi*phin
       f(i,j,k,1,now) = invTauPhi1*f(i,j,k,1,now) + Af1 + Cf1*ux
       f(i,j,k,2,now) = invTauPhi1*f(i,j,k,2,now) + Af1 - Cf1*ux
       f(i,j,k,3,now) = invTauPhi1*f(i,j,k,3,now) + Af1 + Cf1*uy
       f(i,j,k,4,now) = invTauPhi1*f(i,j,k,4,now) + Af1 - Cf1*uy
       f(i,j,k,5,now) = invTauPhi1*f(i,j,k,5,now) + Af1 + Cf1*uz
       f(i,j,k,6,now) = invTauPhi1*f(i,j,k,6,now) + Af1 - Cf1*uz

!-------- Collision for momentum distribution function g -----------------------
! The factor invTauRho for geq is embedded in Egn
       Ag0  = rhon - 6.D0*phin*muPhin
       Ag1  = 3.D0*phin*muPhin + rhon
       Eg1A = Eg1n*Ag1
       Eg2A = Eg2n*Ag1
       Eg1R = Eg1n*rhon
       Eg2R = Eg2n*rhon
       Vsq  = 1.5D0*( ux*ux + uy*uy + uz*uz )
       UF   = ux*sFx + uy*sFy + uz*sFz

! Rescale velocity to avoid unnecessary product operations
       ux = ux*invCs_sq
       uy = uy*invCs_sq
       uz = uz*invCs_sq

! The equilibrium g value and the force are bundled into gFs
! DIRECTION 0
       g(i,j,k,0,nxt) = invTauRhoOne*g(i,j,k,0,now) + Eg0n*( Ag0 - rhon*Vsq ) &
                      - EgC0n*UF

! DIRECTIONS 1 & 2
       geq1 = Eg1A + Eg1R*( 0.5D0*ux*ux - Vsq )
       geq2 = Eg1R*ux

       g(ie,j,k,1,nxt) = invTauRhoOne*g(i,j,k,1,now) + geq1 + geq2 &
                       + EgC1n*( ( 1.D0 + ux )*sFx - UF )

       g(iw,j,k,2,nxt) = invTauRhoOne*g(i,j,k,2,now) + geq1 - geq2 &
                       + EgC1n*( (-1.D0 + ux )*sFx - UF )

! DIRECTIONS 3 & 4
       geq1 = Eg1A + Eg1R*( 0.5D0*uy*uy - Vsq )
       geq2 = Eg1R*uy

       g(i,jn,k,3,nxt) = invTauRhoOne*g(i,j,k,3,now) + geq1 + geq2 &
                       + EgC1n*( ( 1.D0 + uy )*sFy - UF )

       g(i,js,k,4,nxt) = invTauRhoOne*g(i,j,k,4,now) + geq1 - geq2 &
                       + EgC1n*( (-1.D0 + uy )*sFy - UF )

! DIRECTIONS 5 & 6
       geq1 = Eg1A + Eg1R*( 0.5D0*uz*uz - Vsq )
       geq2 = Eg1R*uz

       g(i,j,kt,5,nxt) = invTauRhoOne*g(i,j,k,5,now) + geq1 + geq2 &
                       + EgC1n*( ( 1.D0 + uz )*sFz - UF )

       g(i,j,kb,6,nxt) = invTauRhoOne*g(i,j,k,6,now) + geq1 - geq2 &
                       + EgC1n*( (-1.D0 + uz )*sFz - UF )

! DIRECTION 7 & 8
       Vg   = ux + uy
       geq1 = Eg2A + Eg2R*( 0.5D0*Vg*Vg - Vsq )
       geq2 = Eg2R*Vg

       g(ie,jn,k,7,nxt) = invTauRhoOne*g(i,j,k,7,now) + geq1 + geq2 &
                        + EgC2n*( ( 1.D0 + Vg )*( sFx + sFy ) - UF )

       g(iw,js,k,8,nxt) = invTauRhoOne*g(i,j,k,8,now) + geq1 - geq2 &
                        + EgC2n*( (-1.D0 + Vg )*( sFx + sFy ) - UF )

! DIRECTIONS 9 & 10
       Vg   = ux - uy
       geq1 = Eg2A + Eg2R*( 0.5D0*Vg*Vg - Vsq )
       geq2 = Eg2R*Vg

       g(ie,js,k,9,nxt)  = invTauRhoOne*g(i,j,k,9,now) + geq1 + geq2  &
                         + EgC2n*( ( 1.D0 + Vg )*( sFx - sFy ) - UF )

       g(iw,jn,k,10,nxt) = invTauRhoOne*g(i,j,k,10,now) + geq1 - geq2 &
                         + EgC2n*( (-1.D0 + Vg )*( sFx - sFy ) - UF )

! DIRECTIONS 11 & 12
       Vg   = ux + uz
       geq1 = Eg2A + Eg2R*( 0.5D0*Vg*Vg - Vsq )
       geq2 = Eg2R*Vg

       g(ie,j,kt,11,nxt) = invTauRhoOne*g(i,j,k,11,now) + geq1 + geq2 &
                         + EgC2n*( ( 1.D0 + Vg )*( sFx + sFz ) - UF )

       g(iw,j,kb,12,nxt) = invTauRhoOne*g(i,j,k,12,now) + geq1 - geq2 &
                         + EgC2n*( (-1.D0 + Vg )*( sFx + sFz ) - UF )

! DIRECTIONS 13 & 14
       Vg   = ux - uz
       geq1 = Eg2A + Eg2R*( 0.5D0*Vg*Vg - Vsq )
       geq2 = Eg2R*Vg

       g(ie,j,kb,13,nxt) = invTauRhoOne*g(i,j,k,13,now) + geq1 + geq2 &
                         + EgC2n*( ( 1.D0 + Vg )*( sFx - sFz ) - UF )

       g(iw,j,kt,14,nxt) = invTauRhoOne*g(i,j,k,14,now) + geq1 - geq2 &
                         + EgC2n*( (-1.D0 + Vg )*( sFx - sFz ) - UF )
! DIRECTIONS 15 & 16
       Vg   = uy + uz
       geq1 = Eg2A + Eg2R*( 0.5D0*Vg*Vg - Vsq )
       geq2 = Eg2R*Vg

       g(i,jn,kt,15,nxt) = invTauRhoOne*g(i,j,k,15,now) + geq1 + geq2 &
                         + EgC2n*( ( 1.D0 + Vg )*( sFy + sFz ) - UF )

       g(i,js,kb,16,nxt) = invTauRhoOne*g(i,j,k,16,now) + geq1 - geq2 &
                         + EgC2n*( (-1.D0 + Vg )*( sFy + sFz ) - UF )

! DIRECTIONS 17 & 18
       Vg   = uy - uz
       geq1 = Eg2A + Eg2R*( 0.5D0*Vg*Vg - Vsq )
       geq2 = Eg2R*Vg

       g(i,jn,kb,17,nxt) = invTauRhoOne*g(i,j,k,17,now) + geq1 + geq2 &
                         + EgC2n*( ( 1.D0 + Vg )*( sFy - sFz ) - UF )

       g(i,js,kt,18,nxt) = invTauRhoOne*g(i,j,k,18,now) + geq1 - geq2 &
                         + EgC2n*( (-1.D0 + Vg )*( sFy - sFz ) - UF )

     END DO
   END DO
 END DO

 RETURN
 END SUBROUTINE Collision
