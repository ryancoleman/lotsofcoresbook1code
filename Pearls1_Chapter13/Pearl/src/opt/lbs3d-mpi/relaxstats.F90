!-------------------------------------------------------------------------------
! Subroutine : RelaxStats
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Calculate average velocity, mass conservation factor, effective radius of the
! drop, pressure difference between the inside and the outside of the drop and
! the error with respect to the analytical value given by Laplace's equation,
! and write them to file "relax_stats.out" for the serial D3Q7/D3Q19
! Zheng-Shu-Chew multiphase LBM.
!
! The effective radius is calculated assuming the drop is a perfect circle with
! area given by Vol = (4/3)*Pi*R*R*R.
!
! The pressure inside and the pressure outside of the drop are calculated as
! the average pressures inside and outside the drop, excluding the interface
! area.
!
! The parameters calculated are defined below.
!
! data(1) : Volume
! data(2) : Ux
! data(3) : Uy
! data(4) : Uz
! data(5) : Pin
! data(6) : NodesIn
! data(7) : Pout
! data(8) : NodesOut
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

 SUBROUTINE RelaxStats

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MPIParams
 USE MPI
 USE OMP_LIB
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, m, mg, IO_ERR, MPI_ERR
 REAL(KIND = DBL) :: ux, uy, uz, phin, phin2, pressure, rhon, invRhon
 REAL(KIND = DBL) :: gradPhiSq
 REAL(KIND = DBL) :: invVol, Pdif, Perr, Pin, Pout, R, Ro, Ref, Vef, Vol
 REAL(KIND = DBL) :: NodesIn, NodesOut
 REAL(KIND = DBL), DIMENSION(1:8) :: dataLocal, dataGlobal

! Initialize
 Ro = bubbles(1,4)
 Vol        = 0.D0
 NodesIn    = 0.D0
 NodesOut   = 0.D0
 Pin        = 0.D0
 Pout       = 0.D0
 ux         = 0.D0
 uy         = 0.D0
 uz         = 0.D0
 dataLocal(:)  = 0.D0
 dataGlobal(:) = 0.D0

! Loop through all nodes inside the drop for the current processor, vproc
 DO k = 1, NZ
   DO j = 1, NY
     DO i = 1, NX

       m       = i + NXG*( j + NYG*k )
       mg      = i + NXG*( j + NYG*k ) + now
       rhon    = g0(m)   + g1(mg)  + g2(mg)  + g3(mg)  + g4(mg)  + g5(mg)  + g6(mg)  & 
               + g7(mg)  + g8(mg)  + g9(mg)  + g10(mg) + g11(mg) + g12(mg) + g13(mg) & 
               + g14(mg) + g15(mg) + g16(mg) + g17(mg) + g18(mg)
       invRhon = 1.D0/rhon
       phin    = phi(m)
       phin2   = phin*phin

       gradPhiSq = gradPhiX(m)*gradPhiX(m) &
                 + gradPhiY(m)*gradPhiY(m) &
                 + gradPhiZ(m)*gradPhiZ(m)

! Calculate the pressure
       pressure = alpha*( phin2*( 3.D0*phin2 - 2.D0*phiStar2 ) - phiStar4 ) &
                - kappa*( phin*lapPhi(m) - 0.5D0*gradPhiSq ) + Cs_sq*rhon

! Calculate the velocity
       IF ( phi(m) >= 0.D0 ) THEN
         dataLocal(1) = dataLocal(1) + 1.D0
         dataLocal(2) = dataLocal(2) + ( g1(mg)  - g2(mg)  + g7(mg)    &
                      - g8(mg)  + g9(mg) - g10(mg) + g11(mg) - g12(mg) &
                      + g13(mg) - g14(mg) )*invRhon

         dataLocal(3) = dataLocal(3) + ( g3(mg)  - g4(mg)  + g7(mg)   &
                     - g8(mg)  - g9(mg) + g10(mg) + g15(mg) - g16(mg) &
                     + g17(mg) - g18(mg) )*invRhon

         dataLocal(4) = dataLocal(4) + ( g5(mg)  - g6(mg)  + g11(mg)    &
                      - g12(mg) - g13(mg) + g14(mg) + g15(mg) - g16(mg) &
                      - g17(mg) + g18(mg) )*invRhon

       END IF

       R =  DSQRT( ( DBLE(i+xl) - bubbles(1,1) )**2 &
         +         ( DBLE(j+yl) - bubbles(1,2) )**2 &
         +         ( DBLE(k+zl) - bubbles(1,3) )**2 )

       IF ( R < (Ro - IntWidth) ) THEN
         dataLocal(5) = dataLocal(5) + pressure
         dataLocal(6) = dataLocal(6) + 1.D0
       ELSE IF ( R > (Ro + IntWidth) ) THEN
         dataLocal(7) = dataLocal(7) + pressure
         dataLocal(8) = dataLocal(8) + 1.D0
       END IF

     END DO
   END DO
 END DO

! Gather global information
 CALL MPI_ALLREDUCE( dataLocal, dataGlobal, 8, MPI_DOUBLE_PRECISION, MPI_SUM, &
                     MPI_COMM_VGRID, MPI_ERR)
 Vol    = dataGlobal(1)
 invVol = 1.D0/Vol
 IF ( iStep == 0 ) invInitVol = invVol

! Define average velocity of the drop and effective radius and volume
 ux  = dataGlobal(2)*invVol
 uy  = dataGlobal(3)*invVol
 uz  = dataGlobal(4)*invVol
 Ref = ( Vol*invPi*0.75D0 )**(1.D0/3.D0)
 Vef = Vol*invInitVol

! Global value of the average pressure inside and outside the drop
 Pin  = dataGlobal(5)/dataGlobal(6)
 Pout = dataGlobal(7)/dataGlobal(8)
 Pdif = Pin - Pout
 Perr = 0.5D0*( 2.D0*sigma/Ro - Pdif )*Ro/sigma

! Analyze convergence - Done over a certain period of time ( 10 x tStats ) to
! ensure that long range fluctuations do not give false positives
 Convergence(tCall) = Pdif
 IF ( MOD(tCall,11) == 0 ) THEN
   eps   = 0.D0
   tCall = 1
   DO i = 2, 11
     eps = eps + DABS( Convergence(i) - Convergence(i-1) )
   END DO
   eps = eps*0.1D0/Pdif
 ELSE
   tCall = tCall + 1
 END IF

! Save velocity and effective radius and volume data
#ifdef LBMIO
 IF (proc == master ) THEN
   OPEN(UNIT = 10, FILE = "relax_stats.out", STATUS = "UNKNOWN", &
        POSITION = "APPEND", IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     WRITE(10,'(I9,7ES19.9)')iStep,ux,uy,uz,Vef,Ref,Pdif,Perr
     CLOSE(UNIT = 10)
   ELSE
     CALL MemAlloc(2)
     STOP "Error: Unable to open output file 'relax_stats.out'."
   END IF
 END IF
#endif

 RETURN
 END SUBROUTINE RelaxStats
