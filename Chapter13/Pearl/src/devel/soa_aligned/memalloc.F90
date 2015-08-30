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
 ALLOCATE( phi( -7:NG ) )
! ALLOCATE( nb(xmin:xmax,ymin:ymax,zmin:zmax,1:6) )
 ALLOCATE( lapPhi( -7:NG ) )
 ALLOCATE( gradPhiX( -7:NG ) )
 ALLOCATE( gradPhiY( -7:NG ) )
 ALLOCATE( gradPhiZ( -7:NG ) )

! Distribution function arrays
 ALLOCATE( f0( -7:NG ) )
 ALLOCATE( f1( -7:2*NG ) )
 ALLOCATE( f2( -7:2*NG ) )
 ALLOCATE( f3( -7:2*NG ) )
 ALLOCATE( f4( -7:2*NG ) )
 ALLOCATE( f5( -7:2*NG ) )
 ALLOCATE( f6( -7:2*NG ) )

 ALLOCATE( g0( -7:NG ) )
 ALLOCATE( g1( -7:2*NG ) )
 ALLOCATE( g2( -7:2*NG ) )
 ALLOCATE( g3( -7:2*NG ) )
 ALLOCATE( g4( -7:2*NG ) )
 ALLOCATE( g5( -7:2*NG ) )
 ALLOCATE( g6( -7:2*NG ) )
 ALLOCATE( g7( -7:2*NG ) )
 ALLOCATE( g8( -7:2*NG ) )
 ALLOCATE( g9( -7:2*NG ) )
 ALLOCATE( g10( -7:2*NG ) )
 ALLOCATE( g11( -7:2*NG ) )
 ALLOCATE( g12( -7:2*NG ) )
 ALLOCATE( g13( -7:2*NG ) )
 ALLOCATE( g14( -7:2*NG ) )
 ALLOCATE( g15( -7:2*NG ) )
 ALLOCATE( g16( -7:2*NG ) )
 ALLOCATE( g17( -7:2*NG ) )
 ALLOCATE( g18( -7:2*NG ) )

 ELSE

   DEALLOCATE( phi )
   DEALLOCATE( lapPhi )
   DEALLOCATE( gradPhiX )
   DEALLOCATE( gradPhiY )
   DEALLOCATE( gradPhiZ )
   DEALLOCATE( f0 )
   DEALLOCATE( f1 )
   DEALLOCATE( f2 )
   DEALLOCATE( f3 )
   DEALLOCATE( f4 )
   DEALLOCATE( f5 )
   DEALLOCATE( f6 )
   DEALLOCATE( g0 )
   DEALLOCATE( g1 )
   DEALLOCATE( g2 )
   DEALLOCATE( g3 )
   DEALLOCATE( g4 )
   DEALLOCATE( g5 )
   DEALLOCATE( g6 )
   DEALLOCATE( g7 )
   DEALLOCATE( g8 )
   DEALLOCATE( g9 )
   DEALLOCATE( g10 )
   DEALLOCATE( g11 )
   DEALLOCATE( g12 )
   DEALLOCATE( g13 )
   DEALLOCATE( g14 )
   DEALLOCATE( g15 )
   DEALLOCATE( g16 )
   DEALLOCATE( g17 )
   DEALLOCATE( g18 )

 END IF

 RETURN
 END SUBROUTINE MemAlloc
