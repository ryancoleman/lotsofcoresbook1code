!-------------------------------------------------------------------------------
! Subroutine : Stats
! Revision   : 1.2 (2013/09/10)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Calculate average velocity, mass conservation factor and effective radius of
! the drop, and write them to file "stats.out" for the OMP parallel D3Q7/D3Q19
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

       m = i + NXG*( j + NYG*k ) + now

! Calculate the accumulated quantities
       IF ( phi(m-now) >= 0.D0 ) THEN
         rhon    = g0(m-now)  + g1(m)  + g2(m)  + g3(m)  + g4(m)  + g5(m)  + g6(m)  & 
                 + g7(m)  + g8(m)  + g9(m)  + g10(m) + g11(m) + g12(m) + g13(m) & 
                 + g14(m) + g15(m) + g16(m) + g17(m) + g18(m)
         invRhon = 1.D0/rhon

         ux = ux + ( g1(m)  - g2(m)  + g7(m) - g8(m) + g9(m) - g10(m) + g11(m) &
                 -   g12(m) + g13(m) - g14(m) )*invRhon
         uy = uy + ( g3(m)  - g4(m)  + g7(m) - g8(m) - g9(m) + g10(m) + g15(m) &
                 -   g16(m) + g17(m) - g18(m) )*invRhon
         uz = uz + ( g5(m)  - g6(m)  + g11(m) - g12(m) - g13(m) + g14(m) &
                 +   g15(m) - g16(m) - g17(m) + g18(m) )*invRhon

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
   IF( iStep == 0 ) THEN
      WRITE(10,'(A)')"#MP-LABS v1.2"
      WRITE(10,'(A,5(1x,A))')"#Iter","Ux","Uy","Uz","Vef","Ref"
   END IF
   WRITE(10,'(I9,5(1x,ES19.9))')iStep,ux,uy,uz,Vef,Ref
   CLOSE(UNIT = 10)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: Unable to open output file 'stats.out'."
 END IF

 RETURN
 END SUBROUTINE Stats
