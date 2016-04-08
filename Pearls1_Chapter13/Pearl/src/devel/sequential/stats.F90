!-------------------------------------------------------------------------------
! Subroutine : Stats
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Calculate average velocity, mass conservation factor and effective radius of
! the drop, and write them to file "stats.out" for the hybrid D3Q7/D3Q19
! Zheng-Shu-Chew multiphase LBM.
!
! The effective radius is calculated assuming the drop is a perfect circle with
! area given by Vol = (4/3)*Pi*R*R*R.
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

 SUBROUTINE Stats

! Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain
 USE FluidParams, ONLY : phi
 USE LBMParams
 USE OMP_LIB
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, m, IO_ERR
 REAL(KIND = DBL) :: ux, uy, uz, rhon, invRhon
 REAL(KIND = DBL) :: invVol, Ref, Vef, Vol
 REAL(KIND = DBL), DIMENSION(1:4) :: dataGlobal

! Initialize
 Vol = 0.D0
 ux  = 0.D0
 uy  = 0.D0
 uz  = 0.D0

! Loop through all nodes inside the drop for the current processor, vproc
 DO k = 1, NZ
   DO j = 1, NY
     DO i = 1, NX

! Calculate the accumulated quantities
       IF ( phi(i,j,k) >= 0.D0 ) THEN
   rhon    = g(i,j,k,0,now)  + g(i,j,k,1,now)  + g(i,j,k,2,now)  + g(i,j,k,3,now)  + g(i,j,k,4,now)  &
           + g(i,j,k,5,now)  + g(i,j,k,6,now)  + g(i,j,k,7,now)  + g(i,j,k,8,now)  + g(i,j,k,9,now)  &
           + g(i,j,k,10,now) + g(i,j,k,11,now) + g(i,j,k,12,now) + g(i,j,k,13,now) + g(i,j,k,14,now) &
           + g(i,j,k,15,now) + g(i,j,k,16,now) + g(i,j,k,17,now) + g(i,j,k,18,now)
   invRhon = 1.D0/rhon

   ux = ux + ( g(i,j,k,1,now)  - g(i,j,k,2,now)  + g(i,j,k,7,now)  - g(i,j,k,8,now)  + g(i,j,k,9,now)  &
      -   g(i,j,k,10,now) + g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) - g(i,j,k,14,now) &
       )*invRhon
   uy = uy +( g(i,j,k,3,now)  - g(i,j,k,4,now)  + g(i,j,k,7,now)  - g(i,j,k,8,now)  - g(i,j,k,9,now)  &
      +   g(i,j,k,10,now) + g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) - g(i,j,k,18,now) &
       )*invRhon
   uz = uz + ( g(i,j,k,5,now)  - g(i,j,k,6,now)  + g(i,j,k,11,now) - g(i,j,k,12,now) - g(i,j,k,13,now) &
      +   g(i,j,k,14,now) + g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) + g(i,j,k,18,now) &
       )*invRhon

         Vol = Vol + 1.D0
       END IF

     END DO
   END DO
 END DO

! Define average velocity of the drop and effective radius and volume
 invVol = 1.D0/Vol
 ux  = ux*invVol
 uy  = uy*invVol
 uz  = uz*invVol
 Ref = ( Vol*invPi*0.75D0 )**(1.D0/3.D0)
 Vef = Vol*invInitVol

! Save velocity and effective radius and volume data from the master node only
 OPEN(UNIT = 10, FILE = "stats.out", STATUS = "UNKNOWN", POSITION = "APPEND",&
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(10,'(I9,7ES19.9)')iStep,ux,uy,uz,Vef,Ref
   CLOSE(UNIT = 10)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: Unable to open output file 'stats.out'."
 END IF

 RETURN
 END SUBROUTINE Stats
