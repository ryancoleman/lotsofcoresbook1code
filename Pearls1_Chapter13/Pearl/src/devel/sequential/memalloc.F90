!-------------------------------------------------------------------------------
! Subroutine : MemAlloc
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Allocate (FLAG == 1) or deallocate (FLAG != 1) memory for common arrays.
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

 SUBROUTINE MemAlloc(FLAG)

!  Common Variables
 USE Domain
 USE FluidParams, ONLY : phi, lapPhi, gradPhiX, gradPhiY, gradPhiZ
 USE LBMParams
 IMPLICIT NONE

 !   Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN
! Near neighbors and order parameter arrays
 ALLOCATE( phi( 0:NXG, 0:NYG, 0:NZG ) )
! ALLOCATE( nb(xmin:xmax,ymin:ymax,zmin:zmax,1:6) )
 ALLOCATE( lapPhi( 0:NXG, 0:NYG, 0:NZG ) )
 ALLOCATE( gradPhiX( 0:NXG, 0:NYG, 0:NZG ) )
 ALLOCATE( gradPhiY( 0:NXG, 0:NYG, 0:NZG ) )
 ALLOCATE( gradPhiZ( 0:NXG, 0:NYG, 0:NZG ) )

! Distribution function arrays
 ALLOCATE( f( 0:NXG, 0:NYG, 0:NZG, 0:6, 0:1 ) )
 ALLOCATE( g( 0:NXG, 0:NYG, 0:NZG, 0:18, 0:1 ) )


 ELSE

   DEALLOCATE( phi )
   DEALLOCATE( lapPhi )
   DEALLOCATE( gradPhiX )
   DEALLOCATE( gradPhiY )
   DEALLOCATE( gradPhiZ )
   DEALLOCATE( f )
   DEALLOCATE( g )

 END IF

 RETURN
 END SUBROUTINE MemAlloc
