!-------------------------------------------------------------------------------
! Subroutine : MemAlloc
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Allocate (FLAG == 1) or deallocate (FLAG != 1) memory for common arrays.
!> @details
!! If FLAG is 1 then memory is allocated for the following common arrays:
!!
!! - \b p          (pressure)
!! - \b u          (velocity)
!! - \b rho        (density)
!! - \b psi        (chemical potential)
!! - \b ni         (near neighbors)
!! - \b f, \b fbar (order parameter distribution function)
!! - \b g, \b gbar (momentum distribution function)
!! - \b gradRhoX   (x component of the gradient of the density)
!! - \b gradRhoY   (y component of the gradient of the density)
!! - \b gradRhoSq  (gradient of the density squared)
!! - \b gradRhoXX  (x component of the gradient of the density squared)
!! - \b gradRhoXY  (product of the x and y components of the density gradient)
!! - \b gradRhoYY  (y component of the gradient of the density squared)
!!
!! These are the arrays required by the D2Q9 Lee-Lin multiphase LBM,
!! and are deallocated when FLAG is different from one.

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

 SUBROUTINE MemAlloc(FLAG)

! Common Variables
 USE Domain,    ONLY : ni, xmax, xmin, ymax, ymin
 USE FluidParams
 USE LBMParams, ONLY : f, fbar, fdim, g, gbar, gdim
 IMPLICIT NONE

! Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN

! Memory allocation for hydrodynamics and neighbor arrays
 ALLOCATE(   p(xmin:xmax,ymin:ymax) )
 ALLOCATE( rho(xmin:xmax,ymin:ymax) )
 ALLOCATE( psi(xmin:xmax,ymin:ymax) )
 ALLOCATE(   u(xmin:xmax,ymin:ymax,1:2) )
 ALLOCATE(  ni(xmin:xmax,ymin:ymax,1:4) )

! Memory allocation for distribution functions
 ALLOCATE(    f(xmin:xmax,ymin:ymax,0:fdim) )
 ALLOCATE(    g(xmin:xmax,ymin:ymax,0:gdim) )
 ALLOCATE( fbar(xmin:xmax,ymin:ymax,0:fdim) )
 ALLOCATE( gbar(xmin:xmax,ymin:ymax,0:gdim) )

! Memory allocation for differential terms
 ALLOCATE(  gradRhoX(xmin:xmax,ymin:ymax) )
 ALLOCATE(  gradRhoY(xmin:xmax,ymin:ymax) )
 ALLOCATE( gradRhoSq(xmin:xmax,ymin:ymax) )
 ALLOCATE( gradRhoXX(xmin:xmax,ymin:ymax) )
 ALLOCATE( gradRhoXY(xmin:xmax,ymin:ymax) )
 ALLOCATE( gradRhoYY(xmin:xmax,ymin:ymax) )

 ELSE

! Free hydrodynamics and neighbor arrays
 DEALLOCATE( p )
 DEALLOCATE( u )
 DEALLOCATE( ni )
 DEALLOCATE( rho )
 DEALLOCATE( psi )

! Free distribution function arrays
 DEALLOCATE( f )
 DEALLOCATE( g )
 DEALLOCATE( fbar )
 DEALLOCATE( gbar )

! Free differential term arrays
 DEALLOCATE( gradRhoX )
 DEALLOCATE( gradRhoY )
 DEALLOCATE( gradRhoSq )
 DEALLOCATE( gradRhoXX )
 DEALLOCATE( gradRhoXY )
 DEALLOCATE( gradRhoYY )
 END IF

 RETURN
 END SUBROUTINE MemAlloc

