!-------------------------------------------------------------------------------
! Subroutine : VtkSave
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Save order parameter, pressure, and velocity data in VTK format for the
! serial D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM. Filenames are
! 'EVOL_XXXXXXX.VTK' where XXXXXXX is the time step. These files are understood 
! by Paraview.
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

 SUBROUTINE VtkSave

! Common Variables
 USE Ntypes,      ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i, j, k, m, mg, IO_ERR
 REAL(KIND = DBL) :: ux, uy, uz, phin, phin2, pressure, rhon, invRhon
 REAL(KIND = DBL) :: gradPhiSq
 CHARACTER(20) :: filename

! Set output files formats and names
 IF( STAGE == 0 ) THEN
    WRITE( filename, '(A,I7.7,A)')"INIT_",iStep,".vtk"
 ELSE
    WRITE( filename, '(A,I7.7,A)')"EVOL_",iStep,".vtk"
 END IF  
   

 OPEN(UNIT = 12, FILE = filename, STATUS = "NEW", POSITION = "APPEND", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(12,'(A)')"# vtk DataFile Version 2.0"
   WRITE(12,'(A)')"MP-LABS v1.3"
   WRITE(12,'(A)')"ASCII"
   WRITE(12,*)
   WRITE(12,'(A)')"DATASET STRUCTURED_POINTS"
   WRITE(12,'(A,3(1x,I4))')"DIMENSIONS",NX,NY,NZ
   WRITE(12,'(A,3(1x,I4))')"ORIGIN",xmin,ymin,zmin
   WRITE(12,'(A,3(1x,I1))')"SPACING",1,1,1
   WRITE(12,*)
   WRITE(12,'(A,1x,I12)')"POINT_DATA",NX*NY*NZ
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Phi double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO k = 1, NZ
      DO j = 1, NY
         DO i = 1, NX
            WRITE(UNIT = 12,FMT = '(1(ES21.11E3))')phi( i,j,k )
         END DO
      END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Pressure double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO k = 1, NZ
      DO j = 1, NY
         DO i = 1, NX

   rhon    = g(i,j,k,0,now)  + g(i,j,k,1,now)  + g(i,j,k,2,now)  + g(i,j,k,3,now)  + g(i,j,k,4,now)  &
           + g(i,j,k,5,now)  + g(i,j,k,6,now)  + g(i,j,k,7,now)  + g(i,j,k,8,now)  + g(i,j,k,9,now)  &
           + g(i,j,k,10,now) + g(i,j,k,11,now) + g(i,j,k,12,now) + g(i,j,k,13,now) + g(i,j,k,14,now) &
           + g(i,j,k,15,now) + g(i,j,k,16,now) + g(i,j,k,17,now) + g(i,j,k,18,now)
            phin  = phi(i,j,k)
            phin2 = phin*phin

            gradPhiSq = gradPhiX(i,j,k)*gradPhiX(i,j,k) &
                      + gradPhiY(i,j,k)*gradPhiY(i,j,k) &
                      + gradPhiZ(i,j,k)*gradPhiZ(i,j,k)

! Calculate the pressure
            pressure = alpha*( phin2*( 3.D0*phin2 - 2.D0*phiStar2 ) - phiStar4 ) &
                     - kappa*( phin*lapPhi(i,j,k) - 0.5D0*gradPhiSq ) + Cs_sq*rhon

            WRITE(UNIT = 12,FMT = '(1(ES21.11E3))')pressure
         END DO
      END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"VECTORS Velocity double"
   DO k = 1, NZ
      DO j = 1, NY
         DO i = 1, NX

   rhon    = g(i,j,k,0,now)  + g(i,j,k,1,now)  + g(i,j,k,2,now)  + g(i,j,k,3,now)  + g(i,j,k,4,now)  &
           + g(i,j,k,5,now)  + g(i,j,k,6,now)  + g(i,j,k,7,now)  + g(i,j,k,8,now)  + g(i,j,k,9,now)  &
           + g(i,j,k,10,now) + g(i,j,k,11,now) + g(i,j,k,12,now) + g(i,j,k,13,now) + g(i,j,k,14,now) &
           + g(i,j,k,15,now) + g(i,j,k,16,now) + g(i,j,k,17,now) + g(i,j,k,18,now)
   invRhon = 1.D0/rhon

   ux = ( g(i,j,k,1,now)  - g(i,j,k,2,now)  + g(i,j,k,7,now)  - g(i,j,k,8,now)  + g(i,j,k,9,now)  &
      -   g(i,j,k,10,now) + g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) - g(i,j,k,14,now) &
       )*invRhon
   uy = ( g(i,j,k,3,now)  - g(i,j,k,4,now)  + g(i,j,k,7,now)  - g(i,j,k,8,now)  - g(i,j,k,9,now)  &
      +   g(i,j,k,10,now) + g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) - g(i,j,k,18,now) &
       )*invRhon
   uz = ( g(i,j,k,5,now)  - g(i,j,k,6,now)  + g(i,j,k,11,now) - g(i,j,k,12,now) - g(i,j,k,13,now) &
      +   g(i,j,k,14,now) + g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) + g(i,j,k,18,now) &
       )*invRhon

            WRITE(UNIT = 12,FMT = '(3(ES21.11E3))')ux,uy,uz
         END DO
      END DO
   END DO
   CLOSE(UNIT = 12)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: Unable to open output vtk file."
 END IF

 RETURN
 END SUBROUTINE VtkSave


