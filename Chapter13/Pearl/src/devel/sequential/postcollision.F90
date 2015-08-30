!-------------------------------------------------------------------------------
! Subroutine : PostCollision
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Update g in nodes at processor boundaries and outward f in ghosts.
! This is faster than using a neighbor array in collision.
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

 SUBROUTINE PostCollision

!  Common Variables
 USE Domain
 USE LBMParams
 IMPLICIT NONE

!  Local variables
 INTEGER :: i, j, k


! X direction
 DO k = 1, NZ
   DO j = 1, NY
     f(NX+1,j,k,1,now) = f(1,j,k,1,now)
   END DO
   DO j = 1, NY
     f(0,j,k,2,now) = f(NX,j,k,2,now)
   END DO
 END DO

! Y direction
 DO k = 1, NZ
   DO i = 1, NX
     f(i,NY+1,k,3,now) = f(i,1,k,3,now)
   END DO
   DO i = 1, NX
     f(i,0,k,4,now) = f(i,NY,k,4,now)
   END DO
 END DO

! Z direction
 DO j = 1, NY
   DO i = 1, NX
     f(i,j,NZ+1,5,now) = f(i,j,1,5,now)
   END DO
   DO i = 1, NX
     f(i,j,0,6,now) = f(i,j,NZ,6,now)
   END DO
 END DO

!
! At this stage we are ready to do the relaxation step in stream.f90
!


 RETURN
 END SUBROUTINE PostCollision

