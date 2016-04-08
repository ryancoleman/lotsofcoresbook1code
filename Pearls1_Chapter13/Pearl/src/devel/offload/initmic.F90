!-------------------------------------------------------------------------------
! Subroutine : InitMIC
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Initialization step for the hybrid D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM.
! Smeared interface initialized using equilibrium order parameter function for
! each drop defined in the input (in the range [R-IntWidth,R+IntWidth]). The
! distribution functions f and g are initialized to their equilibrium values
! for zero velocity.
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

 SUBROUTINE InitMIC
!DIR$ ATTRIBUTES OFFLOAD:mic :: InitMIC

!  Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, m, p
 REAL(KIND = DBL) :: r, Af0, Ag0, Af1, Ag1, Ag2
 REAL(KIND = DBL) :: muPhin, phin, rhon
 REAL(KIND = DBL) :: lapPhiInit

 now_mic = 0
 nxt_mic = NG_MIC

! Initialize counters
 iStep = 0
 tCall = 1
 eps   = 1.D0
 STAGE = 0

!$OMP PARALLEL

!--------- Set near neighbors and order parameter ------------------------------
!$OMP DO PRIVATE(i,j,m,p,R)
 DO k = 0, NZ+1
   DO j = 0, NY_MIC+1
     DO i = 0, NX+1

       m = i + NXG*( j + NYG_MIC*k )

       phi_mic(m) = -phistar
       DO p = 1, nBubbles

         R =  DSQRT( ( DBLE(i+xl) - bubbles(p,1) )**2 &
           +         ( DBLE(j+yl) - bubbles(p,2) )**2 &
           +         ( DBLE(k+zl) - bubbles(p,3) )**2 )

         IF ( R <= ( DBLE(bubbles(p,4)) + IntWidth ) ) THEN
           phi_mic(m) = phiStar*DTANH( 2.D0*( DBLE(bubbles(p,4)) - R )/IntWidth )
         END IF

       END DO

     END DO
   END DO
 END DO


!--------- Initialize distribution functions -----------------------------------
!$OMP DO PRIVATE(i,j,m,lapPhiInit,phin,rhon,muPhin,Af0,Af1,Ag0,Ag1,Ag2)
 DO k = 1, NZ
   DO j = 1, NY_MIC
!DIR$ IVDEP
     DO i = 1, NX

       m = i + NXG*( j + NYG_MIC*k )

! Laplacian of the order parameter
   lapPhiInit = ( phi_mic(m+1+NXG) + phi_mic(m-1-NXG)     + phi_mic(m+1-NXG)      &
              +   phi_mic(m-1+NXG) + phi_mic(m+1+NXG*NYG_MIC) + phi_mic(m-1-NXG*NYG_MIC)  &
              +   phi_mic(m+1-NXG*NYG_MIC) + phi_mic(m-1+NXG*NYG_MIC) + phi_mic(m+NXG+NXG*NYG_MIC) &
              +   phi_mic(m-NXG-NXG*NYG_MIC) + phi_mic(m+NXG-NXG*NYG_MIC) + phi_mic(m-NXG+NXG*NYG_MIC) &
              +   2.0D0*( phi_mic(m+1) + phi_mic(m-1) + phi_mic(m+NXG) + phi_mic(m-NXG)        &
              +   phi_mic(m+NXG*NYG_MIC) + phi_mic(m-NXG*NYG_MIC) - 12.d0*phi_mic(m) ) )*inv6

!  Chemical potential
   phin   = phi_mic(m)
   rhon   = 0.5D0*(rhoH + rhoL)
   muPhin = alpha4*phin*(phin*phin - phiStar*phiStar) - kappa*lapPhiInit

! Set distribution function f to its equilibrium value
   Af0 = -3.D0*Gamma*muPhin
   Af1 = 0.5D0*Gamma*muPhin
   f0_mic(m) = Af0 + phin
   f1_mic(m+now_mic) = Af1
   f2_mic(m+now_mic) = Af1
   f3_mic(m+now_mic) = Af1
   f4_mic(m+now_mic) = Af1
   f5_mic(m+now_mic) = Af1
   f6_mic(m+now_mic) = Af1

! Set distribution functiong to its equilibrium value
   Ag0 = Eg0*(rhon - 6.D0*phin*muPhin)
   Ag1 = Eg1*(3.D0*phin*muPhin + rhon)
   Ag2 = Eg2*(3.D0*phin*muPhin + rhon)

   g0_mic(m)  = Ag0

   g1_mic(m+now_mic)  = Ag1
   g2_mic(m+now_mic)  = Ag1
   g3_mic(m+now_mic)  = Ag1
   g4_mic(m+now_mic)  = Ag1
   g5_mic(m+now_mic)  = Ag1
   g6_mic(m+now_mic)  = Ag1

   g7_mic(m+now_mic)  = Ag2
   g8_mic(m+now_mic)  = Ag2
   g9_mic(m+now_mic)  = Ag2
   g10_mic(m+now_mic) = Ag2
   g11_mic(m+now_mic) = Ag2
   g12_mic(m+now_mic) = Ag2
   g13_mic(m+now_mic) = Ag2
   g14_mic(m+now_mic) = Ag2
   g15_mic(m+now_mic) = Ag2
   g16_mic(m+now_mic) = Ag2
   g17_mic(m+now_mic) = Ag2
   g18_mic(m+now_mic) = Ag2

 END DO
 END DO
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE InitMIC
