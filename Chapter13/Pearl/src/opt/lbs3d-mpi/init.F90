!-------------------------------------------------------------------------------
! Subroutine : Init
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Initialization step for the hybrid D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM.
! Smeared interface initialized using equilibrium order parameter function for
! each drop defined in the input (in the range [R-IntWidth,R+IntWidth]). The
! distribution functions f and g are initialized to their equilibrium values
! for zero velocity.
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

 SUBROUTINE Init

!  Common Variables
 USE NTypes,    ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MPIParams
 USE MPI
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, m, p
 REAL(KIND = DBL) :: r, Af0, Ag0, Af1, Ag1, Ag2
 REAL(KIND = DBL) :: muPhin, phin, rhon
 REAL(KIND = DBL) :: lapPhiInit
 INTEGER :: MPI_ERR, MPI_REQ(1:2)
 INTEGER :: MPI_STAT(2,MPI_STATUS_SIZE)

! Initialize counters
 iStep = 0
 tCall = 1
 eps   = 1.D0
 STAGE = 0
 now   = 0
 nxt   = NG

!$OMP PARALLEL

!--------- Set near neighbors and order parameter ------------------------------
!$OMP DO PRIVATE(i,j,m,p,R)
 DO k = 0, NZ+1
   DO j = 0, NY+1
     DO i = 0, NX+1

       m = i + NXG*( j + NYG*k )

       phi(m) = -phistar
       DO p = 1, nBubbles

         R =  DSQRT( ( DBLE(i+xl) - bubbles(p,1) )**2 &
           +         ( DBLE(j+yl) - bubbles(p,2) )**2 &
           +         ( DBLE(k+zl) - bubbles(p,3) )**2 )

         IF ( R <= ( DBLE(bubbles(p,3)) + IntWidth ) ) THEN
           phi(m) = phistar*TANH( 2.D0*( DBLE(bubbles(p,4)) - R )/IntWidth )
         END IF

       END DO

     END DO
   END DO
 END DO


!--------- Initialize distribution functions -----------------------------------
!$OMP DO PRIVATE(i,j,m,lapPhiInit,phin,rhon,muPhin,Af0,Af1,Ag0,Ag1,Ag2)
 DO k = 1, NZ
   DO j = 1, NY
!DIR$ IVDEP
     DO i = 1, NX

       m = i + NXG*( j + NYG*k )

! Laplacian of the order parameter
   lapPhiInit = ( phi(m+1+NXG) + phi(m-1-NXG)     + phi(m+1-NXG)      &
              +   phi(m-1+NXG) + phi(m+1+NXG*NYG) + phi(m-1-NXG*NYG)  &
              +   phi(m+1-NXG*NYG) + phi(m-1+NXG*NYG) + phi(m+NXG+NXG*NYG) &
              +   phi(m-NXG-NXG*NYG) + phi(m+NXG-NXG*NYG) + phi(m-NXG+NXG*NYG) &
              +   2.0D0*( phi(m+1) + phi(m-1) + phi(m+NXG) + phi(m-NXG)        &
              +   phi(m+NXG*NYG) + phi(m-NXG*NYG) - 12.d0*phi(m) ) )*inv6

!  Chemical potential
   phin   = phi(m)
   rhon   = 0.5D0*(rhoH + rhoL)
   muPhin = alpha4*phin*(phin*phin - phiStar*phiStar) - kappa*lapPhiInit

! Set distribution function f to its equilibrium value
   Af0 = -3.D0*Gamma*muPhin
   Af1 = 0.5D0*Gamma*muPhin
   f0(m) = Af0 + phin
   f1(m+now) = Af1
   f2(m+now) = Af1
   f3(m+now) = Af1
   f4(m+now) = Af1
   f5(m+now) = Af1
   f6(m+now) = Af1

! Set distribution functiong to its equilibrium value
   Ag0 = Eg0*(rhon - 6.D0*phin*muPhin)
   Ag1 = Eg1*(3.D0*phin*muPhin + rhon)
   Ag2 = Eg2*(3.D0*phin*muPhin + rhon)

   g0(m)  = Ag0

   g1(m+now)  = Ag1
   g2(m+now)  = Ag1
   g3(m+now)  = Ag1
   g4(m+now)  = Ag1
   g5(m+now)  = Ag1
   g6(m+now)  = Ag1

   g7(m+now)  = Ag2
   g8(m+now)  = Ag2
   g9(m+now)  = Ag2
   g10(m+now) = Ag2
   g11(m+now) = Ag2
   g12(m+now) = Ag2
   g13(m+now) = Ag2
   g14(m+now) = Ag2
   g15(m+now) = Ag2
   g16(m+now) = Ag2
   g17(m+now) = Ag2
   g18(m+now) = Ag2

 END DO
 END DO
 END DO

!----------MPI Section----------------------------------------------------------
! Exchange f values for ghost nodes (needed in "collision.f" within the main
! time loop for the calculation of phi). This takes values from the nodes where
! f is calculated, puts them into two buffers, and then sends these buffers
! to the corresponding ghost nodes in the adjacent processors, so that the ghost
! nodes are fully updated.

!-------- Copy values of outward f from real nodes into the buffers ------------
! X direction
!$OMP SINGLE
 CALL MPI_IRECV( f_east_rcv, xsize, MPI_DOUBLE_PRECISION, east, TAG1, &
                 MPI_COMM_VGRID, MPI_REQ(1), MPI_ERR )
 CALL MPI_IRECV( f_west_rcv, xsize, MPI_DOUBLE_PRECISION, west, TAG2, &
                 MPI_COMM_VGRID, MPI_REQ(2), MPI_ERR )
!$OMP END SINGLE
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     f_west_snd(j+NY*(k-1)) = f1( 1 + NXG*(j+NYG*k) + now )
     f_east_snd(j+NY*(k-1)) = f2( NX + NXG*(j+NYG*k) + now )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_SEND( f_west_snd, xsize, MPI_DOUBLE_PRECISION, west, TAG1, &
                MPI_COMM_VGRID, MPI_ERR )
 CALL MPI_SEND( f_east_snd, xsize, MPI_DOUBLE_PRECISION, east, TAG2, &
                MPI_COMM_VGRID, MPI_ERR )
 CALL MPI_WAITALL( 2, MPI_REQ, MPI_STAT, MPI_ERR )

! Y direction
 CALL MPI_IRECV( f_north_rcv, ysize, MPI_DOUBLE_PRECISION, north, TAG1, &
                 MPI_COMM_VGRID, MPI_REQ(1), MPI_ERR )
 CALL MPI_IRECV( f_south_rcv, ysize, MPI_DOUBLE_PRECISION, south, TAG2, &
                 MPI_COMM_VGRID, MPI_REQ(2), MPI_ERR )
!$OMP END SINGLE
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     f_south_snd(i+NX*(k-1)) = f3( i+NXG*( 1 + NYG*k ) + now )
     f_north_snd(i+NX*(k-1)) = f4( i + NXG*( NY + NYG*k ) + now )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_SEND( f_south_snd, ysize, MPI_DOUBLE_PRECISION, south, TAG1, &
                MPI_COMM_VGRID, MPI_ERR )
 CALL MPI_SEND( f_north_snd, ysize, MPI_DOUBLE_PRECISION, north, TAG2, &
                MPI_COMM_VGRID, MPI_ERR )
 CALL MPI_WAITALL( 2, MPI_REQ, MPI_STAT, MPI_ERR )

! Z direction
 CALL MPI_IRECV( f_top_rcv, zsize, MPI_DOUBLE_PRECISION, top, TAG1, &
                 MPI_COMM_VGRID, MPI_REQ(1), MPI_ERR )
 CALL MPI_IRECV( f_bot_rcv, zsize, MPI_DOUBLE_PRECISION, bot, TAG2, &
                 MPI_COMM_VGRID, MPI_REQ(2), MPI_ERR )
!$OMP END SINGLE
!$OMP DO PRIVATE(i)
 DO j = 1, NY
   DO i = 1, NX
     f_bot_snd(i+NX*(j-1)) = f5( i+NXG*(j+NYG) + now )
     f_top_snd(i+NX*(j-1)) = f6( i+NXG*(j+NYG*NZ) + now )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_SEND( f_bot_snd, zsize, MPI_DOUBLE_PRECISION, bot, TAG1, &
                MPI_COMM_VGRID, MPI_ERR )
 CALL MPI_SEND( f_top_snd, zsize, MPI_DOUBLE_PRECISION, top, TAG2, &
                MPI_COMM_VGRID, MPI_ERR )
 CALL MPI_WAITALL( 2, MPI_REQ, MPI_STAT, MPI_ERR )
!$OMP END SINGLE
!$OMP BARRIER


!-------- Update values of and f in ghost layers -------------------------------
! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     f2( NXG*( j + NYG*k ) + now )        = f_west_rcv(j+NY*(k-1))
     f1( NX+1 + NXG*( j + NYG*k ) + now ) = f_east_rcv(j+NY*(k-1))
   END DO
 END DO

!  Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     f4( i + NXG*NYG*k + now )           = f_south_rcv(i+NX*(k-1))
     f3( i + NXG*( NY+1 + NYG*k) + now ) = f_north_rcv(i+NX*(k-1))
   END DO
 END DO

!  Z direction
!$OMP DO PRIVATE(i)
 DO j = 1, NY
   DO i = 1, NX
     f6( i+NXG*j + now )              = f_bot_rcv(i+NX*(j-1))
     f5( i+NXG*(j+NYG*(NZ+1)) + now ) = f_top_rcv(i+NX*(j-1))
   END DO
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE Init
