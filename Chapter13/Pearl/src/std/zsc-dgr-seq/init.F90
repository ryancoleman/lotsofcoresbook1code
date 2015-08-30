!-------------------------------------------------------------------------------
! Subroutine : Init
! Revision   : 1.0 (2008-06-15)
! Author     : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @file
!! Initialize all variables and arrays and save initialized data to file.
!> @details
!! Initialization step for the dual grid D2Q5/D2Q9 Zheng-Shu-Chew multiphase
!! LBM. Smeared interface initialized using equilibrium order parameter function
!! for each drop defined in the input (in the range [R-IntWidth,R+IntWidth]).
!! The distribution functions f and g are initialized to their equilibrium
!! values for zero velocity.

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

 SUBROUTINE Init

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams, ONLY : Eg0, Eg1, Eg2, f, g
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, ie, iw, jn, js, m, xs, ys
 REAL(KIND = DBL) :: R, rhon, phin, muPhi, lapPhi
 REAL(KIND = DBL) :: Af0, Af1
 REAL(KIND = DBL) :: Ag0, Ag1, Eg1A, Eg2A

! Initialize counters
 now   = 0
 nxt   = 1
 iStep = 0
 tCall = 1
 eps   = 1.D0

! Set vtk limits
 NX = xmax - xmin + 1
 NY = ymax - ymin + 1

 NX_f = xmax_f - xmin_f + 1
 NY_f = ymax_f - ymin_f + 1

! Set initial average density
 rhon = 0.5D0*(rhoH + rhoL)

!--------- Initialize the order parameter and near neigbor list ----------------
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

! Initialize phi assuming heavy fluid is continuous phase
      phi_f(i,j) = -phistar
      DO m = 1, nBubbles
        R =  DSQRT( (DBLE(i)-bubbles(m,1))**2 + (DBLE(j)-bubbles(m,2))**2 )
        IF ( R <= (DBLE(bubbles(m,3)) + IntWidth) ) THEN
          phi_f(i,j) = phistar*TANH( 2.D0*( DBLE(bubbles(m,3)) - R )/IntWidth )
        END IF
      END DO

! Assign near neigbors
      ni_f(i,j,1) = i + 1
      ni_f(i,j,2) = j + 1
      ni_f(i,j,3) = i - 1
      ni_f(i,j,4) = j - 1

    END DO
 END DO
! Correct near neighbours at domain edges (periodic boundary conditions)
 ni_f(xmin_f,:,3) = xmax_f
 ni_f(xmax_f,:,1) = xmin_f
 ni_f(:,ymin_f,4) = ymax_f
 ni_f(:,ymax_f,2) = ymin_f


!---------- Initialize the order parameter distribution function f -------------
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

! Local value of the order parameter
     phin  = phi_f(i,j)

! Identify near neigbors
     ie = ni_f(i,j,1)
     jn = ni_f(i,j,2)
     iw = ni_f(i,j,3)
     js = ni_f(i,j,4)

! Laplacian of the order parameter
     lapPhi_f(i,j) = ( phi_f(ie,jn) + phi_f(ie,js) + phi_f(iw,jn)    &
                   + phi_f(iw,js) + 4.D0*( phi_f(ie,j) + phi_f(iw,j) &
                   + phi_f(i,jn) + phi_f(i,js) ) - 20.D0*phin )*inv6

! Local value of the chemical potential
     muPhi = alpha4*phin*( phin*phin - phistar2 ) - kappa*lapPhi_f(i,j)

! Equilibrium coefficients
     Af1 = 0.5D0*Gamma*muPhi
     Af0 = -2.D0*Gamma*muPhi

! Equilibrium distribution for zero initial velocity
     f(i,j,0)   = Af0 + phin
     f(i,j,1:4) = Af1

   END DO
 END DO

!---------- Initialize g, neighbors in momentum mesh, pressure and velocity ----
 DO j = ymin, ymax
   DO i = xmin, xmax

! Identify near neighbors
      ni(i,j,1) = i + 1
      ni(i,j,2) = j + 1
      ni(i,j,3) = i - 1
      ni(i,j,4) = j - 1

! Define source nodes in the order parameter mesh
     xs = 2*i - 1
     ys = 2*j - 1

! Local values of the phase and the chemical potential
     phin   = phi_f(xs,ys)
     lapPhi = 4.D0*lapPhi_f(xs,ys)
     muPhi  = alpha4*phin*( phin*phin - phistar2 ) - kappaG*lapPhi

! Coefficients
      Ag1  = rhon + 3.D0*phin*muPhi
      Ag0  = 2.25D0*rhon - 1.25D0*Ag1
      Eg1A = Eg1*Ag1
      Eg2A = Eg2*Ag1

! Equilibrium distribution for zero initial velocity
      g(i,j,0,now)   = Eg0*Ag0
      g(i,j,1:4,now) = Eg1A
      g(i,j,5:8,now) = Eg2A

! Initialize pressure and velocity arrays
      p(i,j)   = 1.D0
      u(i,j,:) = 0.D0

   END DO
 END DO
! Correct near neighbours at domain edges (periodic boundary conditions)
 ni(xmin,:,3) = xmax
 ni(xmax,:,1) = xmin
 ni(:,ymin,4) = ymax
 ni(:,ymax,2) = ymin

! Calculate differential terms needed in Collision
 CALL Differentials

! Save initialized data
 CALL Stats
 CALL VtkPlane

 RETURN
 END SUBROUTINE Init
