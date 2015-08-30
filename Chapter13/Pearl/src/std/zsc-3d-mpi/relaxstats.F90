!-------------------------------------------------------------------------------
! Subroutine : RelaxStats
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Calculate intermediate simulation results and save to file 'stats.out'
!> @details
!! Calculate average velocity, mass conservation factor, effective radius of the
!! drop, pressure difference between the inside and the outside of the drop and
!! the error with respect to the analytical value given by Laplace's equation,
!! and write them to file "relax_stats.out" for the parallel D2Q5/D2Q9
!! Zheng-Shu-Chew multiphase LBM.
!!
!! The effective radius is calculated assuming the drop is a perfect circle with
!! area given by Vol = (4/3)*Pi*R*R*R.
!!
!! The pressure inside and the pressure outside of the drop are calculated as
!! the average pressures inside and outside the drop, excluding the interface
!! area.
!!
!! Parallel implementation using MPI. Variables with the "Local" ending refer to
!! quantities calculated locally in the current processor, vproc. Variables with
!! the ending "Global" refer to qunatities calculated on the complete domain.
!! The parameters stored in the MPI excahnge buffers are defined below.
!!
!! @param data(1) : Volume
!! @param data(2) : Ux
!! @param data(3) : Uy
!! @param data(4) : Uz
!! @param data(5) : Pin
!! @param data(6) : NodesIn
!! @param data(7) : Pout
!! @param data(8) : NodesOut

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

 SUBROUTINE RelaxStats

! Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams, ONLY : Cs_sq, g
 USE MPIParams, ONLY : master, MPI_COMM_VGRID, vproc
 USE MPI
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, ie, iw, jn, js, kt, kb
 INTEGER :: IO_ERR, MPI_ERR
 REAL(KIND = DBL) :: ux, uy, uz, phin, phin2, pressure, rhon, invRhon
 REAL(KIND = DBL) :: gradPhiX, gradPhiY, gradPhiZ, gradPhiSq, lapPhi
 REAL(KIND = DBL) :: in00, in01, in02, in03, in04, in05, in06, in07, in08, in09
 REAL(KIND = DBL) :: in10, in11, in12, in13, in14, in15, in16, in17, in18
 REAL(KIND = DBL) :: invVol, Pdif, Perr, Pin, Pout, R, Ro, Ref, Vef, Vol
 REAL(KIND = DBL), DIMENSION(1:8) :: dataLocal, dataGlobal

! Initialize
 dataLocal(:)  = 0.D0
 dataGlobal(:) = 0.D0
 Ro = bubbles(1,4)

! Loop through all nodes inside the drop for the current processor, vproc
 DO k = zl, zu
   DO j = yl, yu
     DO i = xl, xu

       rhon    = SUM( g(i,j,k,:,now) )
       invRhon = 1.D0/rhon
       phin    = phi(i,j,k)
       phin2   = phin*phin

! Identify neighbours
       ie = ni(i,j,k,1)
       iw = ni(i,j,k,2)
       jn = ni(i,j,k,3)
       js = ni(i,j,k,4)
       kt = ni(i,j,k,5)
       kb = ni(i,j,k,6)

! Nodal values of the order parameter
       in00 = phi(i ,j ,k )
       in01 = phi(ie,j ,k )
       in02 = phi(iw,j ,k )
       in03 = phi(i ,jn,k )
       in04 = phi(i ,js,k )
       in05 = phi(i ,j ,kt)
       in06 = phi(i ,j ,kb)
       in07 = phi(ie,jn,k )
       in08 = phi(iw,js,k )
       in09 = phi(ie,js,k )
       in10 = phi(iw,jn,k )
       in11 = phi(ie,j ,kt)
       in12 = phi(iw,j ,kb)
       in13 = phi(ie,j ,kb)
       in14 = phi(iw,j ,kt)
       in15 = phi(i ,jn,kt)
       in16 = phi(i ,js,kb)
       in17 = phi(i ,jn,kb)
       in18 = phi(i ,js,kt)

! Laplacian of the order parameter
       lapPhi = ( in07 + in08 + in09 + in10 + in11 + in12 + in13 + in14 &
              + in15 + in16 + in17 + in18 + 2.0D0*( in01 + in02 + in03  &
              + in04 + in05 + in06 - 12.D0*in00 ) )*inv6

! Components of the order parameter gradient
       gradPhiX  = ( 2.0D0*( in01 - in02 ) + in07 - in08 + in09 - in10 &
                 + in11 - in12 + in13 - in14 )*inv12

       gradPhiY  = ( 2.0D0*( in03 - in04 ) + in07 - in08 + in10 - in09 &
                 + in15 - in16 + in17 - in18 )*inv12

       gradPhiZ  = ( 2.0D0*( in05 - in06 ) + in11 - in12 + in14 - in13 &
                  + in15 - in16 + in18 - in17 )*inv12

       gradPhiSq = gradPhiX*gradPhiX + gradPhiY*gradPhiY + gradPhiZ*gradPhiZ

! Calculate the pressure
       pressure = alpha*( phin2*( 3.D0*phin2 - 2.D0*phiStar2 ) - phiStar4 ) &
                - kappa*( phin*lapPhi - 0.5D0*gradPhiSq ) + Cs_sq*rhon

! Calculate the velocity
       ux = ( g(i,j,k, 1,now) - g(i,j,k, 2,now) + g(i,j,k, 7,now) &
          -   g(i,j,k, 8,now) + g(i,j,k, 9,now) - g(i,j,k,10,now) &
          +   g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) &
          -   g(i,j,k,14,now) )*invRhon

       uy = ( g(i,j,k, 3,now) - g(i,j,k, 4,now) + g(i,j,k, 7,now) &
          -   g(i,j,k, 8,now) - g(i,j,k, 9,now) + g(i,j,k,10,now) &
          +   g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) &
          -   g(i,j,k,18,now) )*invRhon

       uz = ( g(i,j,k, 5,now) - g(i,j,k, 6,now) + g(i,j,k,11,now) &
          -   g(i,j,k,12,now) - g(i,j,k,13,now) + g(i,j,k,14,now) &
          +   g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) &
          +   g(i,j,k,18,now) )*invRhon

! Calculate the accumulated quantities
       IF ( phi(i,j,k) >= 0.D0 ) THEN
         dataLocal(1) = dataLocal(1) + 1.D0
         dataLocal(2) = dataLocal(2) + ux
         dataLocal(3) = dataLocal(3) + uy
         dataLocal(4) = dataLocal(4) + uz
       END IF

       R =  DSQRT( ( DBLE(i) - bubbles(1,1) )**2 + ( DBLE(j) - bubbles(1,2) )**2 &
         + ( DBLE(k) - bubbles(1,3) )**2 )

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
 CALL MPI_ALLREDUCE(dataLocal, dataGlobal, 8, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_VGRID, MPI_ERR)
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

! Save velocity and effective radius and volume data from the master node only
 IF( vproc == master ) THEN
   OPEN(UNIT = 10, FILE = "relax_stats.out", STATUS = "UNKNOWN", &
        POSITION = "APPEND", IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     WRITE(10,'(I9,7ES19.9)')iStep,ux,uy,uz,Vef,Ref,Pdif,Perr
     CLOSE(UNIT = 10)
   ELSE
     CALL MemAlloc(2)
     CALL MPI_FINALIZE(MPI_ERR)
     STOP "Error: Unable to open output file 'relax_stats.out'."
   END IF
 END IF

 RETURN
 END SUBROUTINE RelaxStats
