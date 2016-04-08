!-------------------------------------------------------------------------------
! Subroutine : Collision
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Collision step for order parameter distribution function f, collision and
! streaming steps for momentum distribution function g, and calculation of
! updated macroscopic quantities (velocity, pressure and order parameter) in
! the parallel D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM including gravity.
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
!$OMP DO PRIVATE(i,j,m)
 DO k = 1, NZ
   DO j = 1, NY
!DIR$ IVDEP
     DO i = 1, NX
       m = i + NXG*( j + NYG*k )
       phi(m) = f0(m)     + f1(m+now) + f2(m+now) + f3(m+now) &
              + f4(m+now) + f5(m+now) + f6(m+now)
     END DO
   END DO
 END DO

!-------- Differential terms ---------------------------------------------------
!$OMP DO PRIVATE(i,j,m,in01,in02,in03,in04,in05,in06,in07,in08,in09,in10,in11,in12,in13,in14,in15,in16,in17,in18)
 DO k = 1, NZ
   DO j = 1, NY
!DIR$ IVDEP
     DO i = 1, NX

       m = i + NXG*( j + NYG*k )

! Nodal phi values
   in01 = phi( m + 1                 )
   in02 = phi( m - 1                 )
   in03 = phi( m     + NXG           )
   in04 = phi( m     - NXG           )
   in05 = phi( m           + NXG*NYG )
   in06 = phi( m           - NXG*NYG )
   in07 = phi( m + 1 + NXG           )
   in08 = phi( m - 1 - NXG           ) 
   in09 = phi( m + 1 - NXG           )
   in10 = phi( m - 1 + NXG           )
   in11 = phi( m + 1       + NXG*NYG )
   in12 = phi( m - 1       - NXG*NYG )
   in13 = phi( m + 1       - NXG*NYG )
   in14 = phi( m - 1       + NXG*NYG )
   in15 = phi( m     + NXG + NXG*NYG )
   in16 = phi( m     - NXG - NXG*NYG )
   in17 = phi( m     + NXG - NXG*NYG )
   in18 = phi( m     - NXG + NXG*NYG )

! Laplacian of the order parameter phi
   lapPhi(m) = ( in07 + in08 + in09 + in10 + in11 + in12 + in13 + in14 + in15  &
             +   in16 + in17 + in18 + 2.0D0*( in01 + in02 + in03 + in04 + in05 &
             +   in06 - 12.D0*phi(m) ) )*inv6

! Gradient components of the order parameter phi
   gradPhiX(m) = ( 2.0D0*( in01 - in02 ) + in07 - in08 + in09 - in10 + in11 &
               - in12 + in13 - in14 )*inv12

   gradPhiY(m) = ( 2.0D0*( in03 - in04 ) + in07 - in08 + in10 - in09 + in15 &
               - in16 + in17 - in18 )*inv12

   gradPhiZ(m) = ( 2.0D0*( in05 - in06 ) + in11 - in12 + in14 - in13 + in15 &
               - in16 + in18 - in17 )*inv12

   END DO
   END DO
 END DO

!$OMP DO PRIVATE(i,j,m,phin,phin2,rhon,invRhon,muPhin,sFx,sFy,sFz,ux,uy,uz,Af0,Af1,Cf1,Ag0,Ag1,Eg1A, Eg2A,Eg1R,Eg2R,Vsq,UF,geq1,geq2)
 DO k = 1, NZ
   DO j = 1, NY
!DIR$ IVDEP
     DO i = 1, NX

       m = i + NXG*( j + NYG*k )

! Define some local values
   phin    = phi(m)
   phin2   = phin*phin
   rhon    = g0(m)      + g1(m+now)  + g2(m+now)  + g3(m+now)  + g4(m+now)  &
           + g5(m+now)  + g6(m+now)  + g7(m+now)  + g8(m+now)  + g9(m+now)  &
           + g10(m+now) + g11(m+now) + g12(m+now) + g13(m+now) + g14(m+now) &
           + g15(m+now) + g16(m+now) + g17(m+now) + g18(m+now)
   invRhon = 1.D0/rhon

!-------- Hydrodynamics (velocity, pressure, interfacial force) ----------------
! Chemical potential
   muPhin = alpha4*phin*( phin2 - phiStar2 ) - kappa*lapPhi(m)

! Interfacial force and gravity
   sFx = muPhin*gradPhiX(m)
   sFy = muPhin*gradPhiY(m)
   sFz = muPhin*gradPhiZ(m)
   IF (phin > 0.D0) sFz = sFz + grav

! Velocity
   ux = ( g1(m+now)  - g2(m+now)  + g7(m+now)  - g8(m+now)  + g9(m+now)  &
      -   g10(m+now) + g11(m+now) - g12(m+now) + g13(m+now) - g14(m+now) &
      +   0.5D0*sFx )*invRhon
   uy = ( g3(m+now)  - g4(m+now)  + g7(m+now)  - g8(m+now)  - g9(m+now)  &
      +   g10(m+now) + g15(m+now) - g16(m+now) + g17(m+now) - g18(m+now) &
      +   0.5D0*sFy )*invRhon
   uz = ( g5(m+now)  - g6(m+now)  + g11(m+now) - g12(m+now) - g13(m+now) &
      +   g14(m+now) + g15(m+now) - g16(m+now) - g17(m+now) + g18(m+now) &
      +   0.5D0*sFz )*invRhon

!-------- Collision for order parameter distribution function f ----------------
   Af0 = -3.D0*Gamma*muPhin*invTauPhi
   Af1 = 0.5D0*Gamma*muPhin*invTauPhi
   Cf1 = invTauPhi*invEta2*phin

   f0(m)     = invTauPhi1*f0(m)     + Af0 + invTauPhi*phin
   f1(m+now) = invTauPhi1*f1(m+now) + Af1 + Cf1*ux
   f2(m+now) = invTauPhi1*f2(m+now) + Af1 - Cf1*ux
   f3(m+now) = invTauPhi1*f3(m+now) + Af1 + Cf1*uy
   f4(m+now) = invTauPhi1*f4(m+now) + Af1 - Cf1*uy
   f5(m+now) = invTauPhi1*f5(m+now) + Af1 + Cf1*uz
   f6(m+now) = invTauPhi1*f6(m+now) + Af1 - Cf1*uz

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
   g0(m) = invTauRhoOne*g0(m) + Eg0n*( Ag0 - rhon*Vsq ) - EgC0n*UF

! DIRECTIONS 1 & 2
   geq1(1) = Eg1A + Eg1R*( 0.5D0*ux*ux - Vsq ) + EgC1n*( ux*sFx - UF ) 
   geq2(1) = Eg1R*ux + EgC1n*sFx

   g1(m+nxt + 1) = invTauRhoOne*g1(m+now) + geq1(1) + geq2(1) 
   g2(m+nxt - 1) = invTauRhoOne*g2(m+now) + geq1(1) - geq2(1) 

! DIRECTIONS 3 & 4
   geq1(2) = Eg1A + Eg1R*( 0.5D0*uy*uy - Vsq ) + EgC1n*( uy*sFy - UF )
   geq2(2) = Eg1R*uy + EgC1n*sFy

   g3(m+nxt + NXG) = invTauRhoOne*g3(m+now) + geq1(2) + geq2(2)
   g4(m+nxt - NXG) = invTauRhoOne*g4(m+now) + geq1(2) - geq2(2)

! DIRECTIONS 5 & 6
   geq1(3) = Eg1A + Eg1R*( 0.5D0*uz*uz - Vsq ) + EgC1n*( uz*sFz - UF )
   geq2(3) = Eg1R*uz + EgC1n*sFz

   g5(m+nxt + NXG*NYG) = invTauRhoOne*g5(m+now) + geq1(3) + geq2(3)
   g6(m+nxt - NXG*NYG) = invTauRhoOne*g6(m+now) + geq1(3) - geq2(3)

! DIRECTION 7 & 8
   geq1(4) = Eg2A + Eg2R*( 0.5D0*( ux + uy )*( ux + uy ) - Vsq ) &
           + EgC2n*( ( ux + uy )*( sFx + sFy ) - UF )
   geq2(4) = Eg2R*( ux + uy ) + EgC2n*( sFx + sFy )

   g7(m+nxt + 1 + NXG) = invTauRhoOne*g7(m+now) + geq1(4) + geq2(4)
   g8(m+nxt - 1 - NXG) = invTauRhoOne*g8(m+now) + geq1(4) - geq2(4)

! DIRECTIONS 9 & 10
   geq1(5) = Eg2A + Eg2R*( 0.5D0*( ux - uy )*( ux - uy ) - Vsq ) &
           + EgC2n*( ( ux - uy )*( sFx - sFy ) - UF )
   geq2(5) = Eg2R*( ux - uy )+ EgC2n*( sFx - sFy )

   g9(m+nxt + 1 - NXG)  = invTauRhoOne*g9(m+now)  + geq1(5) + geq2(5)
   g10(m+nxt - 1 + NXG) = invTauRhoOne*g10(m+now) + geq1(5) - geq2(5)

! DIRECTIONS 11 & 12
   geq1(6) = Eg2A + Eg2R*( 0.5D0*( ux + uz )*( ux + uz ) - Vsq ) &
           + EgC2n*( ( ux + uz )*( sFx + sFz ) - UF )
   geq2(6) = Eg2R*( ux + uz ) + EgC2n*( sFx + sFz )

   g11(m+nxt + 1 + NXG*NYG) = invTauRhoOne*g11(m+now) + geq1(6) + geq2(6)
   g12(m+nxt - 1 - NXG*NYG) = invTauRhoOne*g12(m+now) + geq1(6) - geq2(6)

! DIRECTIONS 13 & 14
   geq1(7) = Eg2A + Eg2R*( 0.5D0*( ux - uz )*( ux - uz ) - Vsq ) &
           + EgC2n*( ( ux - uz )*( sFx - sFz ) - UF )
   geq2(7) = Eg2R*( ux - uz ) + EgC2n*( sFx - sFz )

   g13(m+nxt + 1 - NXG*NYG) = invTauRhoOne*g13(m+now) + geq1(7) + geq2(7)
   g14(m+nxt - 1 + NXG*NYG) = invTauRhoOne*g14(m+now) + geq1(7) - geq2(7)

! DIRECTIONS 15 & 16
   geq1(8) = Eg2A + Eg2R*( 0.5D0*( uy + uz )*( uy + uz ) - Vsq ) &
           + EgC2n*( ( uy + uz )*( sFy + sFz ) - UF )
   geq2(8) = Eg2R*( uy + uz ) + EgC2n*( sFy + sFz )

   g15(m+nxt + NXG + NXG*NYG) = invTauRhoOne*g15(m+now) + geq1(8) + geq2(8)
   g16(m+nxt - NXG - NXG*NYG) = invTauRhoOne*g16(m+now) + geq1(8) - geq2(8)

! DIRECTIONS 17 & 18
   geq1(9) = Eg2A + Eg2R*( 0.5D0*( uy - uz )*( uy - uz ) - Vsq ) &
           + EgC2n*( ( uy - uz )*( sFy - sFz ) - UF )
   geq2(9) = Eg2R*( uy - uz ) + EgC2n*( sFy - sFz )

   g17(m+nxt + NXG - NXG*NYG) = invTauRhoOne*g17(m+now) + geq1(9) + geq2(9) 
   g18(m+nxt - NXG + NXG*NYG) = invTauRhoOne*g18(m+now) + geq1(9) - geq2(9)
   END DO
   END DO
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE Collision
