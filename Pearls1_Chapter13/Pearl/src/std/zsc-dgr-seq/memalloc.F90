!-------------------------------------------------------------------------------
! Subroutine : memAlloc
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Allocate (FLAG == 1) or deallocate (FLAG != 1) memory for common arrays.
!> @details
!! If input FLAG is 1 then memory is allocated for the following common arrays:
!!
!! - \b p   (pressure)
!! - \b u   (velocity)
!! - \b phi (order parameter)
!! - \b ni  (near neighbors)
!! - \b f   (order parameter distribution function)
!! - \b g   (momentum distribution function)
!!
!! These are the arrays required by the dual grid D2Q5/D2Q9 Zheng-Shu-Chew
!! multiphase LBM, and are deallocated when FLAG is different from one.

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
 USE Domain,      ONLY : ni, ni_f, xmax, xmax_f, xmin, xmin_f, ymax, ymax_f, ymin, ymin_f
 USE FluidParams, ONLY : gradPhiX_f, gradPhiY_f, lapPhi_f, p, p_f, phi_f, u, u_f
 USE LBMParams,   ONLY : f, fcol, fdim, g, gdim
 IMPLICIT NONE

! Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN

! Hydrodynamics and neighbor arrays - Momentum/Pressure mesh
 ALLOCATE(   p(xmin:xmax,ymin:ymax) )
 ALLOCATE(   u(xmin:xmax,ymin:ymax,2) )
 ALLOCATE(  ni(xmin:xmax,ymin:ymax,1:4) )

! Hydrodynamics and neighbor arrays - Order parameter mesh
 ALLOCATE(   p_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE(   u_f(xmin_f:xmax_f,ymin_f:ymax_f,2) )
 ALLOCATE( phi_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE(  ni_f(xmin_f:xmax_f,ymin_f:ymax_f,1:4) )

! Differential term arrays - Order parameter mesh
 ALLOCATE( gradPhiX_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE( gradPhiY_f(xmin_f:xmax_f,ymin_f:ymax_f) )
 ALLOCATE(   lapPhi_f(xmin_f:xmax_f,ymin_f:ymax_f) )

! Memory allocation for distribution functions
 ALLOCATE(    f(xmin_f:xmax_f,ymin_f:ymax_f,0:fdim) )
 ALLOCATE( fcol(xmin_f:xmax_f,ymin_f:ymax_f,0:fdim) )
 ALLOCATE(    g(xmin:xmax,ymin:ymax,0:gdim,0:1) )

 ELSE

! Free momentum/pressure mesh arrays
 DEALLOCATE( p )
 DEALLOCATE( u )
 DEALLOCATE( ni )
 DEALLOCATE( g )

! Free order parameter mesh arrays
 DEALLOCATE( f )
 DEALLOCATE( fcol )
 DEALLOCATE( ni_f )
 DEALLOCATE( phi_f )
 DEALLOCATE( p_f )
 DEALLOCATE( u_f )
 DEALLOCATE( gradPhiX_f )
 DEALLOCATE( gradPhiY_f )
 DEALLOCATE( lapPhi_f )

 END IF

 RETURN
 END SUBROUTINE MemAlloc

