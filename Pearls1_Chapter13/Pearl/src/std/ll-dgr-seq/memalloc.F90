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
!! - \b p, \b p_f          (pressure)
!! - \b u, \b u_f          (velocity)
!! - \b rho, \b rho_f        (density)
!! - \b psi_f              (chemical potential)
!! - \b ni, \b ni_f           (near neighbors)
!! - \b f, \b fbar (order parameter distribution function)
!! - \b g, \b gbar (momentum distribution function)
!! - \b GX, \b GY  (interfacial force components in x and y directions)
!! - \b gradRhoX, \b gradRhoX_f   (x component of the gradient of the density)
!! - \b gradRhoY, \b gradRhoY_f   (y component of the gradient of the density)
!! - \b gradRhoSq_f  (gradient of the density squared)
!! - \b gradRhoXX_f  (x component of the gradient of the density squared)
!! - \b gradRhoXY_f  (product of the x and y components of the density gradient)
!! - \b gradRhoYY_f  (y component of the gradient of the density squared)
!!
!! These are the arrays required by the dual grid D2Q9 Lee-Lin multiphase LBM,
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
 USE Domain,     ONLY : ni, ni_f, xmax, xmax_f, xmin, xmin_f, ymax, ymax_f, ymin, ymin_f
 USE FluidParams
 USE LBMParams,  ONLY : f, fbar, fdim, g, gbar, gdim
 IMPLICIT NONE

! Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN

! Memory allocation for momentum/pressure grid arrays
 ALLOCATE(    p(xmin:xmax,ymin:ymax) )
 ALLOCATE(  rho(xmin:xmax,ymin:ymax) )
 ALLOCATE(    u(xmin:xmax,ymin:ymax,1:2) )
 ALLOCATE(   ni(xmin:xmax,ymin:ymax,1:4) )
 ALLOCATE(    g(xmin:xmax,ymin:ymax,0:gdim) )
 ALLOCATE( gbar(xmin:xmax,ymin:ymax,0:gdim) )

! Memory allocation for order parameter grid arrays
 ALLOCATE(   p_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE( rho_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE( psi_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE(   u_f(xmin_f:xmax_f,ymin_f:ymax_f,1:2) )
 ALLOCATE(  ni_f(xmin_f:xmax_f,ymin_f:ymax_f,1:4) )
 ALLOCATE(     f(xmin_f:xmax_f,ymin_f:ymax_f,0:fdim) )
 ALLOCATE(  fbar(xmin_f:xmax_f,ymin_f:ymax_f,0:fdim) )

! Memory allocation for differential terms in the momentum/pressure grid
 ALLOCATE(       GX(xmin:xmax,ymin:ymax) )
 ALLOCATE(       GY(xmin:xmax,ymin:ymax) )
 ALLOCATE( gradRhoX(xmin:xmax,ymin:ymax) )
 ALLOCATE( gradRhoY(xmin:xmax,ymin:ymax) )

! Memory allocation for differential terms in the order parameter grid
 ALLOCATE(  gradRhoX_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE(  gradRhoY_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE( gradRhoSq_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE( gradRhoXX_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE( gradRhoXY_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE( gradRhoYY_f(xmin_f:xmax_f,ymin_f:ymax_f) )

 ELSE

! Free memory
 DEALLOCATE( p )
 DEALLOCATE( u )
 DEALLOCATE( g )
 DEALLOCATE( ni )
 DEALLOCATE( rho )
 DEALLOCATE( gbar )
 DEALLOCATE( gradRhoX )
 DEALLOCATE( gradRhoY )
 DEALLOCATE( GX )
 DEALLOCATE( GY )


 DEALLOCATE( p_f )
 DEALLOCATE( u_f )
 DEALLOCATE( f )
 DEALLOCATE( ni_f )
 DEALLOCATE( rho_f )
 DEALLOCATE( psi_f )
 DEALLOCATE( fbar )
 DEALLOCATE( gradRhoX_f )
 DEALLOCATE( gradRhoY_f )
 DEALLOCATE( gradRhoSq_f )
 DEALLOCATE( gradRhoXX_f )
 DEALLOCATE( gradRhoXY_f )
 DEALLOCATE( gradRhoYY_f )

 END IF

 RETURN
 END SUBROUTINE MemAlloc

