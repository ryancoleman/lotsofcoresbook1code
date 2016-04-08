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
!! - \b p   (pressure)
!! - \b u   (velocity)
!! - \b phi (order parameter)
!! - \b ni  (near neighbors)
!! - \b f   (order parameter distribution function)
!! - \b g   (momentum distribution function)
!!
!! These are the arrays required by the D2Q5/D2Q9 Zheng-Shu-Chew multiphase LBM,
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
 USE Domain,      ONLY : ni, xmin, xmax, ymin, ymax
 USE FluidParams, ONLY : phi, p, u
 USE LBMParams,   ONLY : f, fdim, g, gdim
 IMPLICIT NONE

! Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN

! Hydrodynamics and neighbor arrays
 ALLOCATE(   p(xmin:xmax,ymin:ymax) )
 ALLOCATE(   u(xmin:xmax,ymin:ymax,2) )
 ALLOCATE( phi(xmin:xmax,ymin:ymax) )
 ALLOCATE(  ni(xmin:xmax,ymin:ymax,1:4) )

! Memory allocation for distribution functions
 ALLOCATE(    f(xmin:xmax,ymin:ymax,0:fdim,0:1) )
 ALLOCATE(    g(xmin:xmax,ymin:ymax,0:gdim,0:1) )

 ELSE

! Free momentum/pressure mesh arrays
 DEALLOCATE( p )
 DEALLOCATE( u )
 DEALLOCATE( ni )
 DEALLOCATE( phi )
 DEALLOCATE( g )
 DEALLOCATE( f )

 END IF

 RETURN
 END SUBROUTINE MemAlloc

