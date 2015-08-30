!-------------------------------------------------------------------------------
! Subroutine : CollisionF
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Collision step for distribution function f
!> @details
!! Collision step for the order parameter distribution function f in the dual
!! grid parallel D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM.

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

 SUBROUTINE CollisionF

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain,    ONLY : xl_f, xu_f, yl_f, yu_f
 USE FluidParams
 USE LBMParams, ONLY : f, fcol
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j
 REAL(KIND = DBL) :: phin, muPhi, ux, uy
 REAL(KIND = DBL) :: Af0, Af1, Cfp


 DO j = yl_f, yu_f
   DO i = xl_f, xu_f

! Local values of the velocity
    ux = u_f(i,j,1)
    uy = u_f(i,j,2)

! Local values of phi and the chemical potential
     phin  = phi_f(i,j)
     muPhi = alpha4*phin*( phin*phin - phiStar2 ) - kappa*lapPhi_f(i,j)

! Collision
     Af1 = 0.5D0*Gamma*muPhi
     Af0 = -2.D0*Gamma*muPhi
     Cfp = invEta2*phin
     f(i,j,0)    = f(i,j,0) + invTauPhi*( Af0 + phin   - f(i,j,0) ) 
     fcol(i,j,1) = f(i,j,1) + invTauPhi*( Af1 + Cfp*ux - f(i,j,1) )
     fcol(i,j,2) = f(i,j,2) + invTauPhi*( Af1 + Cfp*uy - f(i,j,2) )
     fcol(i,j,3) = f(i,j,3) + invTauPhi*( Af1 - Cfp*ux - f(i,j,3) )
     fcol(i,j,4) = f(i,j,4) + invTauPhi*( Af1 - Cfp*uy - f(i,j,4) )

   END DO
 END DO

 RETURN
 END SUBROUTINE CollisionF

