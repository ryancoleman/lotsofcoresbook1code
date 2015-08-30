!-------------------------------------------------------------------------------
! Subroutine : HydrodynamicsF
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Calculate chemical potential and density gradient
!> @details
!! Calculate the chemical potential and the density gradient in the order
!! parameter grid for the parallel dual grid D2Q9 Lee-Lin multiphase LBM.

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

 SUBROUTINE HydrodynamicsF

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js
 REAL(KIND = DBL) :: lapRho_f


 DO j = ylg_f, yug_f
   DO i = xlg_f, xug_f

! Identify first neighbors
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

! Gradient of the density rho
     gradRhoX_f(i,j) = ( 4.D0*( rho_f(ie,j ) - rho_f(iw,j ) ) + rho_f(ie,jn) &
                     - rho_f(iw,jn) + rho_f(ie,js) - rho_f(iw,js) )*inv12

     gradRhoY_f(i,j) = ( 4.D0*( rho_f(i ,jn) - rho_f(i ,js) ) + rho_f(ie,jn) &
                     - rho_f(ie,js) + rho_f(iw,jn) - rho_f(iw,js) )*inv12
   END DO
 END DO

!---------- Calculate differential terms in the order parameter mesh (stress formulation) -----------------------------
 DO j = ylg2_f, yug2_f
   DO i = xlg2_f, xug2_f

! Identify first neighbors
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

! Laplacian of the density rho
     lapRho_f = 4.D0*( rho_f(ie,j ) + rho_f(iw,j ) + rho_f(i ,jn)            &
              + rho_f(i ,js) ) +  rho_f(ie,jn) + rho_f(ie,js) + rho_f(iw,jn) &
              + rho_f(iw,js) - 20.D0*rho_f(i,j)

! Define Psi
     psi_f(i,j) = beta4*( rho_f(i,j) - rhoStar )*( rho_f(i,j) - rhoL )*( rho_f(i,j) - rhoH ) &
                - kappa_6*lapRho_f

   END DO
 END DO

 RETURN
 END SUBROUTINE HydrodynamicsF
