!-------------------------------------------------------------------------------
! Subroutine : Init
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Initialization step for the OMP D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM.
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

 SUBROUTINE Init

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
 REAL(KIND = DBL) :: in01,in02,in03,in04,in05,in06,in07,in08,in09  
 REAL(KIND = DBL) :: in10,in11,in12,in13,in14,in15,in16,in17,in18

! Initialize counters
 iStep = 0
 tCall = 1
 eps   = 1.D0
 STAGE = 0

!$OMP PARALLEL
!--------- Set near neighbors and order parameter ------------------------------
!$OMP DO PRIVATE(i,j,m,p,R)
 DO k = 0, NZ+1
   DO j = 0, NY+1
     DO i = 0, NX+1

       phi(i,j,k) = -phistar
       DO p = 1, nBubbles

         R =  DSQRT( ( DBLE(i+xl) - bubbles(p,1) )**2 &
           +         ( DBLE(j+yl) - bubbles(p,2) )**2 &
           +         ( DBLE(k+zl) - bubbles(p,3) )**2 )

         IF ( R <= ( DBLE(bubbles(p,4)) + IntWidth ) ) THEN
           phi(i,j,k) = phistar*DTANH( 2.D0*( DBLE(bubbles(p,4)) - R )/IntWidth )
         END IF

       END DO

     END DO
   END DO
 END DO


!--------- Initialize distribution functions -----------------------------------
!$OMP DO PRIVATE(i,j,m,lapPhiInit,phin,rhon,muPhin,Af0,Af1,Ag0,Ag1,Ag2)
 DO k = 1, NZ
   DO j = 1, NY
!DIR$ IVDEP
     DO i = 1, NX

! Laplacian of the order parameter
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
   lapPhiInit = ( in07 + in08 + in09 + in10 + in11 + in12 + in13 + in14 + in15  &
             +   in16 + in17 + in18 + 2.0D0*( in01 + in02 + in03 + in04 + in05 &
             +   in06 - 12.D0*phi(i,j,k) ) )*inv6

!  Chemical potential
   phin   = phi(i,j,k)
   rhon   = 0.5D0*(rhoH + rhoL)
   muPhin = alpha4*phin*(phin*phin - phiStar*phiStar) - kappa*lapPhiInit

! Set distribution function f to its equilibrium value
   Af0 = -3.D0*Gamma*muPhin
   Af1 = 0.5D0*Gamma*muPhin
   f(i,j,k,0,now) = Af0 + phin
   f(i,j,k,1,now) = Af1
   f(i,j,k,2,now) = Af1
   f(i,j,k,3,now) = Af1
   f(i,j,k,4,now) = Af1
   f(i,j,k,5,now) = Af1
   f(i,j,k,6,now) = Af1

! Set distribution functiong to its equilibrium value
   Ag0 = Eg0*(rhon - 6.D0*phin*muPhin)
   Ag1 = Eg1*(3.D0*phin*muPhin + rhon)
   Ag2 = Eg2*(3.D0*phin*muPhin + rhon)

   g(i,j,k,0,now)  = Ag0

   g(i,j,k,1,now)  = Ag1
   g(i,j,k,2,now)  = Ag1
   g(i,j,k,3,now)  = Ag1
   g(i,j,k,4,now)  = Ag1
   g(i,j,k,5,now)  = Ag1
   g(i,j,k,6,now)  = Ag1

   g(i,j,k,7,now)  = Ag2
   g(i,j,k,8,now)  = Ag2
   g(i,j,k,9,now)  = Ag2
   g(i,j,k,10,now) = Ag2
   g(i,j,k,11,now) = Ag2
   g(i,j,k,12,now) = Ag2
   g(i,j,k,13,now) = Ag2
   g(i,j,k,14,now) = Ag2
   g(i,j,k,15,now) = Ag2
   g(i,j,k,16,now) = Ag2
   g(i,j,k,17,now) = Ag2
   g(i,j,k,18,now) = Ag2

 END DO
 END DO
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE Init
