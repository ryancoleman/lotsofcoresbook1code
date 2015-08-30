!-------------------------------------------------------------------------------
! Subroutine : VtkPlane
! Revision   : 1.0 (2008-06-15)
! Author     : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @file
!! Save simulation data in VTK format
!> @details
!! Save order parameter, pressure, and velocity data in VTK format for the
!! D2Q5/D2Q9 Zheng-Shu-Chew Multiphase LBM. Filenames are 'DATA_XXXXXXX.VTK'
!! where XXXXXXX is the time step. These files are understood by Paraview.

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

 SUBROUTINE VtkPlane

!	Common Variables
 USE Domain,      ONLY : now, NX, NY, xmax, xmin, ymax, ymin
 USE FluidParams, ONLY : p, u
 USE LBMParams,   ONLY : f
 IMPLICIT NONE

!	Required functions
 CHARACTER(8) :: filer

!	Local Variables
 INTEGER :: i, j
 INTEGER :: IO_ERR
 CHARACTER(8)  :: filer1
 CHARACTER(16) :: fmtd1, fmtd3
 CHARACTER(17) :: filer2
 
!-------- Set output files formats and names -----------------------------------
 fmtd1 = '(1(ES21.11E3))'
 fmtd3 = '(3(ES21.11E3))'
 filer1 = filer()
 filer2 = "DATA"//filer1(1:8)//".vtk"
  
 OPEN(UNIT = 12, FILE = filer2, STATUS = "NEW", POSITION = "APPEND", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(12,'(A)')"# vtk DataFile Version 2.0"
   WRITE(12,'(A)')"MP-LABS v1.2"
   WRITE(12,'(A)')"ASCII"
   WRITE(12,*)
   WRITE(12,'(A)')"DATASET STRUCTURED_POINTS"
   WRITE(12,*)"DIMENSIONS",NX,NY,1
   WRITE(12,*)"ORIGIN",xmin,ymin,1
   WRITE(12,*)"SPACING",1,1,1
   WRITE(12,*)
   WRITE(12,*)"POINT_DATA",NX*NY
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Phi double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO j = ymin, ymax
       DO i = xmin, xmax
          WRITE(UNIT = 12,FMT = fmtd1)SUM(f(i,j,:,now))
       END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"SCALARS Pressure double"
   WRITE(12,'(A)')"LOOKUP_TABLE default"
   DO j = ymin, ymax
       DO i = xmin, xmax
          WRITE(UNIT = 12,FMT = fmtd1)p(i,j)
       END DO
   END DO
   WRITE(12,*)
   WRITE(12,'(A)')"VECTORS Velocity double"
   DO j = ymin, ymax
       DO i = xmin, xmax
         WRITE(UNIT = 12,FMT = fmtd3)u(i,j,1),u(i,j,2),0.D0
      END DO
   END DO
   CLOSE(UNIT = 12)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: Unable to open output vtk file."
 END IF

 RETURN
 END SUBROUTINE VtkPlane

!-------------------------------------------------------------------------------
! Function : filer
! Revision : 1.0 (2008-06-15)
! Author   : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @brief Set vtk file name as a function of processor number and time step.
 FUNCTION filer()

 USE DOMAIN, ONLY: iStep
 IMPLICIT NONE

! Function type
 CHARACTER(8) :: filer

!-------- Decide how to name the files depending on the timestep (iStep) -------
   IF (iStep < 10) THEN
     WRITE(filer,9001)iStep
   ELSE IF (iStep < 100) THEN
     WRITE(filer,9002)iStep
   ELSE IF (iStep < 1000) THEN
     WRITE(filer,9003)iStep
   ELSE IF (iStep < 10000) THEN
     WRITE(filer,9004)iStep
   ELSE IF (iStep < 100000) THEN
     WRITE(filer,9005)iStep
   ELSE IF (iStep < 1000000) THEN
     WRITE(filer,9006)iStep
   ELSE IF (iStep < 10000000) THEN
     WRITE(filer,9007)iStep
   END IF

!-------- Define the file name format depending on the timestep (iStep) --------
 9001 FORMAT('_000000',i1)
 9002 FORMAT('_00000',i2)
 9003 FORMAT('_0000',i3)
 9004 FORMAT('_000',i4)
 9005 FORMAT('_00',i5)
 9006 FORMAT('_0',i6)
 9007 FORMAT('_',i7)

 RETURN
 END FUNCTION filer
