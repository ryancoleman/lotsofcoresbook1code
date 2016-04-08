!-------------------------------------------------------------------------------
! Subroutine : CollisionMIC
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

 SUBROUTINE CollisionMIC
!DIR$ ATTRIBUTES OFFLOAD:mic :: CollisionMIC

!  Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain
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
   DO j = 1, NY_MIC
!DIR$ IVDEP
     DO i = 1, NX
       m = i + NXG*( j + NYG_MIC*k )
       phi_mic(m) = f0_mic(m)     + f1_mic(m+now_mic) + f2_mic(m+now_mic) + f3_mic(m+now_mic) &
                  + f4_mic(m+now_mic) + f5_mic(m+now_mic) + f6_mic(m+now_mic)
     END DO
   END DO
 END DO

!-------- Differential terms ---------------------------------------------------
!$OMP DO PRIVATE(i,j,m,in01,in02,in03,in04,in05,in06,in07,in08,in09,in10,in11,in12,in13,in14,in15,in16,in17,in18)
 DO k = 1, NZ
   DO j = 1, NY_MIC
!DIR$ IVDEP
     DO i = 1, NX

       m = i + NXG*( j + NYG_MIC*k )

! Nodal phi values
   in01 = phi_mic( m + 1                     )
   in02 = phi_mic( m - 1                     )
   in03 = phi_mic( m     + NXG               )
   in04 = phi_mic( m     - NXG               )
   in05 = phi_mic( m           + NXG*NYG_MIC )
   in06 = phi_mic( m           - NXG*NYG_MIC )
   in07 = phi_mic( m + 1 + NXG               )
   in08 = phi_mic( m - 1 - NXG               ) 
   in09 = phi_mic( m + 1 - NXG               )
   in10 = phi_mic( m - 1 + NXG               )
   in11 = phi_mic( m + 1       + NXG*NYG_MIC )
   in12 = phi_mic( m - 1       - NXG*NYG_MIC )
   in13 = phi_mic( m + 1       - NXG*NYG_MIC )
   in14 = phi_mic( m - 1       + NXG*NYG_MIC )
   in15 = phi_mic( m     + NXG + NXG*NYG_MIC )
   in16 = phi_mic( m     - NXG - NXG*NYG_MIC )
   in17 = phi_mic( m     + NXG - NXG*NYG_MIC )
   in18 = phi_mic( m     - NXG + NXG*NYG_MIC )

! Laplacian of the order parameter phi
   lapPhi_mic(m) = ( in07 + in08 + in09 + in10 + in11 + in12 + in13 + in14 + in15  &
                 +   in16 + in17 + in18 + 2.0D0*( in01 + in02 + in03 + in04 + in05 &
                 +   in06 - 12.D0*phi_mic(m) ) )*inv6

! Gradient components of the order parameter phi
   gradPhiX_mic(m) = ( 2.0D0*( in01 - in02 ) + in07 - in08 + in09 - in10 + in11 &
                   - in12 + in13 - in14 )*inv12

   gradPhiY_mic(m) = ( 2.0D0*( in03 - in04 ) + in07 - in08 + in10 - in09 + in15 &
                   - in16 + in17 - in18 )*inv12

   gradPhiZ_mic(m) = ( 2.0D0*( in05 - in06 ) + in11 - in12 + in14 - in13 + in15 &
                   - in16 + in18 - in17 )*inv12

   END DO
   END DO
 END DO

!$OMP DO PRIVATE(i,j,m,phin,phin2,rhon,invRhon,muPhin,sFx,sFy,sFz,ux,uy,uz,Af0,Af1,Cf1,Ag0,Ag1,Eg1A, Eg2A,Eg1R,Eg2R,Vsq,UF,geq1,geq2)
 DO k = 1, NZ
   DO j = 1, NY_MIC
!DIR$ IVDEP
     DO i = 1, NX

       m = i + NXG*( j + NYG_MIC*k )

! Define some local values
   phin    = phi_mic(m)
   phin2   = phin*phin
   rhon    = g0_mic(m)      + g1_mic(m+now_mic)  + g2_mic(m+now_mic)  + g3_mic(m+now_mic)  + g4_mic(m+now_mic)  &
           + g5_mic(m+now_mic)  + g6_mic(m+now_mic)  + g7_mic(m+now_mic)  + g8_mic(m+now_mic)  + g9_mic(m+now_mic)  &
           + g10_mic(m+now_mic) + g11_mic(m+now_mic) + g12_mic(m+now_mic) + g13_mic(m+now_mic) + g14_mic(m+now_mic) &
           + g15_mic(m+now_mic) + g16_mic(m+now_mic) + g17_mic(m+now_mic) + g18_mic(m+now_mic)
   invRhon = 1.D0/rhon

!-------- Hydrodynamics (velocity, pressure, interfacial force) ----------------
! Chemical potential
   muPhin = alpha4*phin*( phin2 - phiStar2 ) - kappa*lapPhi_mic(m)

! Interfacial force and gravity
   sFx = muPhin*gradPhiX_mic(m)
   sFy = muPhin*gradPhiY_mic(m)
   sFz = muPhin*gradPhiZ_mic(m)
   IF (phin > 0.D0) sFz = sFz + grav

! Velocity
   ux = ( g1_mic(m+now_mic)  - g2_mic(m+now_mic)  + g7_mic(m+now_mic)  - g8_mic(m+now_mic)  + g9_mic(m+now_mic)  &
      -   g10_mic(m+now_mic) + g11_mic(m+now_mic) - g12_mic(m+now_mic) + g13_mic(m+now_mic) - g14_mic(m+now_mic) &
      +   0.5D0*sFx )*invRhon
   uy = ( g3_mic(m+now_mic)  - g4_mic(m+now_mic)  + g7_mic(m+now_mic)  - g8_mic(m+now_mic)  - g9_mic(m+now_mic)  &
      +   g10_mic(m+now_mic) + g15_mic(m+now_mic) - g16_mic(m+now_mic) + g17_mic(m+now_mic) - g18_mic(m+now_mic) &
      +   0.5D0*sFy )*invRhon
   uz = ( g5_mic(m+now_mic)  - g6_mic(m+now_mic)  + g11_mic(m+now_mic) - g12_mic(m+now_mic) - g13_mic(m+now_mic) &
      +   g14_mic(m+now_mic) + g15_mic(m+now_mic) - g16_mic(m+now_mic) - g17_mic(m+now_mic) + g18_mic(m+now_mic) &
      +   0.5D0*sFz )*invRhon

!-------- Collision for order parameter distribution function f ----------------
   Af0 = -3.D0*Gamma*muPhin*invTauPhi
   Af1 = 0.5D0*Gamma*muPhin*invTauPhi
   Cf1 = invTauPhi*invEta2*phin

   f0_mic(m)     = invTauPhi1*f0_mic(m)     + Af0 + invTauPhi*phin
   f1_mic(m+now_mic) = invTauPhi1*f1_mic(m+now_mic) + Af1 + Cf1*ux
   f2_mic(m+now_mic) = invTauPhi1*f2_mic(m+now_mic) + Af1 - Cf1*ux
   f3_mic(m+now_mic) = invTauPhi1*f3_mic(m+now_mic) + Af1 + Cf1*uy
   f4_mic(m+now_mic) = invTauPhi1*f4_mic(m+now_mic) + Af1 - Cf1*uy
   f5_mic(m+now_mic) = invTauPhi1*f5_mic(m+now_mic) + Af1 + Cf1*uz
   f6_mic(m+now_mic) = invTauPhi1*f6_mic(m+now_mic) + Af1 - Cf1*uz

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
   g0_mic(m)   = invTauRhoOne*g0_mic(m) + Eg0n*( Ag0 - rhon*Vsq ) - EgC0n*UF

! DIRECTIONS 1 & 2
   geq1(1) = Eg1A + Eg1R*( 0.5D0*ux*ux - Vsq ) + EgC1n*( ux*sFx - UF ) 
   geq2(1) = Eg1R*ux + EgC1n*sFx

   g1_mic(m+nxt_mic + 1) = invTauRhoOne*g1_mic(m+now_mic) + geq1(1) + geq2(1) 
   g2_mic(m+nxt_mic - 1) = invTauRhoOne*g2_mic(m+now_mic) + geq1(1) - geq2(1) 

! DIRECTIONS 3 & 4
   geq1(2) = Eg1A + Eg1R*( 0.5D0*uy*uy - Vsq ) + EgC1n*( uy*sFy - UF )
   geq2(2) = Eg1R*uy + EgC1n*sFy

   g3_mic(m+nxt_mic + NXG) = invTauRhoOne*g3_mic(m+now_mic) + geq1(2) + geq2(2)
   g4_mic(m+nxt_mic - NXG) = invTauRhoOne*g4_mic(m+now_mic) + geq1(2) - geq2(2)

! DIRECTIONS 5 & 6
   geq1(3) = Eg1A + Eg1R*( 0.5D0*uz*uz - Vsq ) + EgC1n*( uz*sFz - UF )
   geq2(3) = Eg1R*uz + EgC1n*sFz

   g5_mic(m+nxt_mic + NXG*NYG_MIC) = invTauRhoOne*g5_mic(m+now_mic) + geq1(3) + geq2(3)
   g6_mic(m+nxt_mic - NXG*NYG_MIC) = invTauRhoOne*g6_mic(m+now_mic) + geq1(3) - geq2(3)

! DIRECTION 7 & 8
   geq1(4) = Eg2A + Eg2R*( 0.5D0*( ux + uy )*( ux + uy ) - Vsq ) &
           + EgC2n*( ( ux + uy )*( sFx + sFy ) - UF )
   geq2(4) = Eg2R*( ux + uy ) + EgC2n*( sFx + sFy )

   g7_mic(m+nxt_mic + 1 + NXG) = invTauRhoOne*g7_mic(m+now_mic) + geq1(4) + geq2(4)
   g8_mic(m+nxt_mic - 1 - NXG) = invTauRhoOne*g8_mic(m+now_mic) + geq1(4) - geq2(4)

! DIRECTIONS 9 & 10
   geq1(5) = Eg2A + Eg2R*( 0.5D0*( ux - uy )*( ux - uy ) - Vsq ) &
           + EgC2n*( ( ux - uy )*( sFx - sFy ) - UF )
   geq2(5) = Eg2R*( ux - uy )+ EgC2n*( sFx - sFy )

   g9_mic(m+nxt_mic + 1 - NXG)  = invTauRhoOne*g9_mic(m+now_mic)  + geq1(5) + geq2(5)
   g10_mic(m+nxt_mic - 1 + NXG) = invTauRhoOne*g10_mic(m+now_mic) + geq1(5) - geq2(5)

! DIRECTIONS 11 & 12
   geq1(6) = Eg2A + Eg2R*( 0.5D0*( ux + uz )*( ux + uz ) - Vsq ) &
           + EgC2n*( ( ux + uz )*( sFx + sFz ) - UF )
   geq2(6) = Eg2R*( ux + uz ) + EgC2n*( sFx + sFz )

   g11_mic(m+nxt_mic + 1 + NXG*NYG_MIC) = invTauRhoOne*g11_mic(m+now_mic) + geq1(6) + geq2(6)
   g12_mic(m+nxt_mic - 1 - NXG*NYG_MIC) = invTauRhoOne*g12_mic(m+now_mic) + geq1(6) - geq2(6)

! DIRECTIONS 13 & 14
   geq1(7) = Eg2A + Eg2R*( 0.5D0*( ux - uz )*( ux - uz ) - Vsq ) &
           + EgC2n*( ( ux - uz )*( sFx - sFz ) - UF )
   geq2(7) = Eg2R*( ux - uz ) + EgC2n*( sFx - sFz )

   g13_mic(m+nxt_mic + 1 - NXG*NYG_MIC) = invTauRhoOne*g13_mic(m+now_mic) + geq1(7) + geq2(7)
   g14_mic(m+nxt_mic - 1 + NXG*NYG_MIC) = invTauRhoOne*g14_mic(m+now_mic) + geq1(7) - geq2(7)

! DIRECTIONS 15 & 16
   geq1(8) = Eg2A + Eg2R*( 0.5D0*( uy + uz )*( uy + uz ) - Vsq ) &
           + EgC2n*( ( uy + uz )*( sFy + sFz ) - UF )
   geq2(8) = Eg2R*( uy + uz ) + EgC2n*( sFy + sFz )

   g15_mic(m+nxt_mic + NXG + NXG*NYG_MIC) = invTauRhoOne*g15_mic(m+now_mic) + geq1(8) + geq2(8)
   g16_mic(m+nxt_mic - NXG - NXG*NYG_MIC) = invTauRhoOne*g16_mic(m+now_mic) + geq1(8) - geq2(8)

! DIRECTIONS 17 & 18
   geq1(9) = Eg2A + Eg2R*( 0.5D0*( uy - uz )*( uy - uz ) - Vsq ) &
           + EgC2n*( ( uy - uz )*( sFy - sFz ) - UF )
   geq2(9) = Eg2R*( uy - uz ) + EgC2n*( sFy - sFz )

   g17_mic(m+nxt_mic + NXG - NXG*NYG_MIC) = invTauRhoOne*g17_mic(m+now_mic) + geq1(9) + geq2(9) 
   g18_mic(m+nxt_mic - NXG + NXG*NYG_MIC) = invTauRhoOne*g18_mic(m+now_mic) + geq1(9) - geq2(9)
   END DO
   END DO
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE CollisionMIC
