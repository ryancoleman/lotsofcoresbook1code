!-------------------------------------------------------------------------------
! Subroutine : Collision
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Collision step for order parameter distribution function f, collision and
! streaming steps for momentum distribution function g, and calculation of
! updated macroscopic quantities (velocity, pressure and order parameter) in
! the OMP D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM including gravity.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
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
 USE Domain,      ONLY : inv12, inv6, NX, NXG, NY, NYG, NZ, NZG, now, nxt
 USE FluidParams
 USE LBMParams
 USE OMP_LIB
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, itemp, j, k, m
 REAL(KIND = DBL) :: muPhin, phin, phin2
 REAL(KIND = DBL) :: rhon, invRhon
 REAL(KIND = DBL) :: sFx, sFy, sFz
 REAL(KIND = DBL) :: UF, ux, uy, uz, Vsq
 REAL(KIND = DBL) :: Af0, Af1, Ag0, Ag1, Cf1
 REAL(KIND = DBL) :: Eg1A, Eg2A, Eg1R, Eg2R
 REAL(KIND = DBL) :: geq1(1:9), geq2(1:9)
 REAL(KIND = DBL) :: in01,in02,in03,in04,in05,in06,in07,in08,in09  
 REAL(KIND = DBL) :: in10,in11,in12,in13,in14,in15,in16,in17,in18

!$OMP PARALLEL

! Order parameter update
!$OMP DO PRIVATE(i,j)
 DO k = 1, NZ
   DO j = 1, NY
     DO i = 1, NX
       phi(i,j,k) = f(i,j,k,0,now) + f(i,j,k,1,now) + f(i,j,k,2,now) + f(i,j,k,3,now) &
                  + f(i,j,k,4,now) + f(i,j,k,5,now) + f(i,j,k,6,now)
     END DO
   END DO
 END DO

!-------- Differential terms ---------------------------------------------------
!$OMP DO PRIVATE(i,j,in01,in02,in03,in04,in05,in06,in07,in08,in09,in10,in11,in12,in13,in14,in15,in16,in17,in18)
 DO k = 1, NZ
   DO j = 1, NY
     DO i = 1, NX

! Nodal phi values
   in01 = phi( i+1,j,k )
   in02 = phi( i-1,j,k )
   in03 = phi( i,j+1,k )
   in04 = phi( i,j-1,k )
   in05 = phi( i,j,k+1 )
   in06 = phi( i,j,k-1 )
   in07 = phi( i+1,j+1,k )
   in08 = phi( i-1,j-1,k ) 
   in09 = phi( i+1,j-1,k )
   in10 = phi( i-1,j+1,k )
   in11 = phi( i+1,j,k+1 )
   in12 = phi( i-1,j,k-1 )
   in13 = phi( i+1,j,k-1 )
   in14 = phi( i-1,j,k+1 )
   in15 = phi( i,j+1,k+1 )
   in16 = phi( i,j-1,k-1 )
   in17 = phi( i,j+1,k-1 )
   in18 = phi( i,j-1,k+1 )

! Laplacian of the order parameter phi
   lapPhi(i,j,k) = ( in07 + in08 + in09 + in10 + in11 + in12 + in13 + in14 + in15  &
             +   in16 + in17 + in18 + 2.0D0*( in01 + in02 + in03 + in04 + in05 &
             +   in06 - 12.D0*phi(i,j,k) ) )*inv6

! Gradient components of the order parameter phi
   gradPhiX(i,j,k) = ( 2.0D0*( in01 - in02 ) + in07 - in08 + in09 - in10 + in11 &
               - in12 + in13 - in14 )*inv12

   gradPhiY(i,j,k) = ( 2.0D0*( in03 - in04 ) + in07 - in08 + in10 - in09 + in15 &
               - in16 + in17 - in18 )*inv12

   gradPhiZ(i,j,k) = ( 2.0D0*( in05 - in06 ) + in11 - in12 + in14 - in13 + in15 &
               - in16 + in18 - in17 )*inv12

   END DO
   END DO
 END DO

!$OMP DO PRIVATE(i,j,phin,phin2,rhon,invRhon,muPhin,sFx,sFy,sFz,ux,uy,uz,Af0,Af1,Cf1,Ag0,Ag1,Eg1A, Eg2A,Eg1R,Eg2R,Vsq,UF,geq1,geq2)
 DO k = 1, NZ
   DO j = 1, NY
     DO i = 1, NX

! Define some local values
   phin    = phi(i,j,k)
   phin2   = phin*phin
   rhon    = g(i,j,k,0,now)  + g(i,j,k,1,now)  + g(i,j,k,2,now)  + g(i,j,k,3,now)  + g(i,j,k,4,now)  &
           + g(i,j,k,5,now)  + g(i,j,k,6,now)  + g(i,j,k,7,now)  + g(i,j,k,8,now)  + g(i,j,k,9,now)  &
           + g(i,j,k,10,now) + g(i,j,k,11,now) + g(i,j,k,12,now) + g(i,j,k,13,now) + g(i,j,k,14,now) &
           + g(i,j,k,15,now) + g(i,j,k,16,now) + g(i,j,k,17,now) + g(i,j,k,18,now)
   invRhon = 1.D0/rhon

!-------- Hydrodynamics (velocity, pressure, interfacial force) ----------------
! Chemical potential
   muPhin = alpha4*phin*( phin2 - phiStar2 ) - kappa*lapPhi(i,j,k)

! Interfacial force and gravity
   sFx = muPhin*gradPhiX(i,j,k)
   sFy = muPhin*gradPhiY(i,j,k)
   sFz = muPhin*gradPhiZ(i,j,k)
   IF (phin > 0.D0) sFz = sFz + grav

! Velocity
   ux = ( g(i,j,k,1,now)  - g(i,j,k,2,now)  + g(i,j,k,7,now)  - g(i,j,k,8,now)  + g(i,j,k,9,now)  &
      -   g(i,j,k,10,now) + g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) - g(i,j,k,14,now) &
      +   0.5D0*sFx )*invRhon
   uy = ( g(i,j,k,3,now)  - g(i,j,k,4,now)  + g(i,j,k,7,now)  - g(i,j,k,8,now)  - g(i,j,k,9,now)  &
      +   g(i,j,k,10,now) + g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) - g(i,j,k,18,now) &
      +   0.5D0*sFy )*invRhon
   uz = ( g(i,j,k,5,now)  - g(i,j,k,6,now)  + g(i,j,k,11,now) - g(i,j,k,12,now) - g(i,j,k,13,now) &
      +   g(i,j,k,14,now) + g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) + g(i,j,k,18,now) &
      +   0.5D0*sFz )*invRhon

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
   g(i,j,k,0,nxt) = invTauRhoOne*g(i,j,k,0,now) + Eg0n*( Ag0 - rhon*Vsq ) - EgC0n*UF

! DIRECTIONS 1 & 2
   geq1(1) = Eg1A + Eg1R*( 0.5D0*ux*ux - Vsq ) + EgC1n*( ux*sFx - UF ) 
   geq2(1) = Eg1R*ux + EgC1n*sFx

   g(i+1,j,k,1,nxt) = invTauRhoOne*g(i,j,k,1,now) + geq1(1) + geq2(1) 
   g(i-1,j,k,2,nxt) = invTauRhoOne*g(i,j,k,2,now) + geq1(1) - geq2(1) 

! DIRECTIONS 3 & 4
   geq1(2) = Eg1A + Eg1R*( 0.5D0*uy*uy - Vsq ) + EgC1n*( uy*sFy - UF )
   geq2(2) = Eg1R*uy + EgC1n*sFy

   g(i,j+1,k,3,nxt) = invTauRhoOne*g(i,j,k,3,now) + geq1(2) + geq2(2)
   g(i,j-1,k,4,nxt) = invTauRhoOne*g(i,j,k,4,now) + geq1(2) - geq2(2)

! DIRECTIONS 5 & 6
   geq1(3) = Eg1A + Eg1R*( 0.5D0*uz*uz - Vsq ) + EgC1n*( uz*sFz - UF )
   geq2(3) = Eg1R*uz + EgC1n*sFz

   g(i,j,k+1,5,nxt) = invTauRhoOne*g(i,j,k,5,now) + geq1(3) + geq2(3)
   g(i,j,k-1,6,nxt) = invTauRhoOne*g(i,j,k,6,now) + geq1(3) - geq2(3)

! DIRECTION 7 & 8
   geq1(4) = Eg2A + Eg2R*( 0.5D0*( ux + uy )*( ux + uy ) - Vsq ) &
           + EgC2n*( ( ux + uy )*( sFx + sFy ) - UF )
   geq2(4) = Eg2R*( ux + uy ) + EgC2n*( sFx + sFy )

   g(i+1,j+1,k,7,nxt) = invTauRhoOne*g(i,j,k,7,now) + geq1(4) + geq2(4)
   g(i-1,j-1,k,8,nxt) = invTauRhoOne*g(i,j,k,8,now) + geq1(4) - geq2(4)

! DIRECTIONS 9 & 10
   geq1(5) = Eg2A + Eg2R*( 0.5D0*( ux - uy )*( ux - uy ) - Vsq ) &
           + EgC2n*( ( ux - uy )*( sFx - sFy ) - UF )
   geq2(5) = Eg2R*( ux - uy )+ EgC2n*( sFx - sFy )

   g(i+1,j-1,k,9,nxt)  = invTauRhoOne*g(i,j,k,9,now)  + geq1(5) + geq2(5)
   g(i-1,j+1,k,10,nxt) = invTauRhoOne*g(i,j,k,10,now) + geq1(5) - geq2(5)

! DIRECTIONS 11 & 12
   geq1(6) = Eg2A + Eg2R*( 0.5D0*( ux + uz )*( ux + uz ) - Vsq ) &
           + EgC2n*( ( ux + uz )*( sFx + sFz ) - UF )
   geq2(6) = Eg2R*( ux + uz ) + EgC2n*( sFx + sFz )

   g(i+1,j,k+1,11,nxt) = invTauRhoOne*g(i,j,k,11,now) + geq1(6) + geq2(6)
   g(i-1,j,k-1,12,nxt) = invTauRhoOne*g(i,j,k,12,now) + geq1(6) - geq2(6)

! DIRECTIONS 13 & 14
   geq1(7) = Eg2A + Eg2R*( 0.5D0*( ux - uz )*( ux - uz ) - Vsq ) &
           + EgC2n*( ( ux - uz )*( sFx - sFz ) - UF )
   geq2(7) = Eg2R*( ux - uz ) + EgC2n*( sFx - sFz )

   g(i+1,j,k-1,13,nxt) = invTauRhoOne*g(i,j,k,13,now) + geq1(7) + geq2(7)
   g(i-1,j,k+1,14,nxt) = invTauRhoOne*g(i,j,k,14,now) + geq1(7) - geq2(7)

! DIRECTIONS 15 & 16
   geq1(8) = Eg2A + Eg2R*( 0.5D0*( uy + uz )*( uy + uz ) - Vsq ) &
           + EgC2n*( ( uy + uz )*( sFy + sFz ) - UF )
   geq2(8) = Eg2R*( uy + uz ) + EgC2n*( sFy + sFz )

   g(i,j+1,k+1,15,nxt) = invTauRhoOne*g(i,j,k,15,now) + geq1(8) + geq2(8)
   g(i,j-1,k-1,16,nxt) = invTauRhoOne*g(i,j,k,16,now) + geq1(8) - geq2(8)

! DIRECTIONS 17 & 18
   geq1(9) = Eg2A + Eg2R*( 0.5D0*( uy - uz )*( uy - uz ) - Vsq ) &
           + EgC2n*( ( uy - uz )*( sFy - sFz ) - UF )
   geq2(9) = Eg2R*( uy - uz ) + EgC2n*( sFy - sFz )

   g(i,j+1,k-1,17,nxt) = invTauRhoOne*g(i,j,k,17,now) + geq1(9) + geq2(9) 
   g(i,j-1,k+1,18,nxt) = invTauRhoOne*g(i,j,k,18,now) + geq1(9) - geq2(9)
   END DO
   END DO
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE Collision
