!-------------------------------------------------------------------------------
! Subroutine : Init
! Revision   : 1.0 (2008-06-15)
! Author     : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!              Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Initialize all variables and arrays and save initialized data to file.
!> @details
!! Initialization step for the parallel D3Q7/D3Q19 Zheng-Shu-Chew multiphase LBM.
!! Smeared interface initialized using equilibrium order parameter function for
!! each drop defined in the input (in the range [R-IntWidth,R+IntWidth]). The
!! distribution functions f and g are initialized to their equilibrium values
!! for zero velocity.

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
 INTEGER :: count, i, j, k, m, ie, iw, jn, js, kt, kb
 REAL(KIND = DBL) :: r, Af0, Ag0, Af1, Ag1, Ag2
 REAL(KIND = DBL) :: muPhin, phin, rhon
 REAL(KIND = DBL) :: lapPhi
 INTEGER :: MPI_ERR, MPI_REQ1, MPI_REQ2
 INTEGER :: status(MPI_STATUS_SIZE)

! Initialize counters
 iStep = 0
 tCall = 1
 now   = 0
 nxt   = 1
 eps   = 1.D0

! Set vtk limits
 NX = xu - xl + 1
 NY = yu - yl + 1
 NZ = zu - zl + 1
 NPZ = 0
 NPY = 0
 NPX = 0
 DO k = zl, zu, xjump
   NPZ = NPZ + 1
 END DO
 DO j = yl, yu, xjump
   NPY = NPY + 1
 END DO
 DO i = xl, xu, xjump
   NPX = NPX + 1
 END DO

!--------- Set near neighbors and order parameter ------------------------------
 DO k = zlg, zug
   DO j = ylg, yug
     DO i = xlg, xug

       phi(i,j,k) = -phistar
       DO m = 1, nBubbles

         R =  DSQRT( ( DBLE(i) - bubbles(m,1) )**2 &
           + ( DBLE(j) - bubbles(m,2) )**2 + ( DBLE(k) - bubbles(m,3) )**2 )

         IF ( R <= ( DBLE(bubbles(m,3)) + IntWidth ) ) THEN
           phi(i,j,k) = phistar*TANH( 2.D0*( DBLE(bubbles(m,4)) - R )/IntWidth )
         END IF
       END DO

! Assign neighbor values to array
         ni(i,j,k,1) = i + 1
         ni(i,j,k,2) = i - 1
         ni(i,j,k,3) = j + 1
         ni(i,j,k,4) = j - 1
         ni(i,j,k,5) = k + 1
         ni(i,j,k,6) = k - 1

     END DO
   END DO
 END DO
!  Fix near neighbours at edges (necessary because of ghost node layers)
 IF(xl == xmin) THEN
   DO k = zlg, zug
     DO j = ylg, yug
       ni(xmin,j,k,2) = xlg
     END DO
   END DO
 END IF
 IF(xu == xmax) THEN
   DO k = zlg, zug
     DO j = ylg, yug
       ni(xmax,j,k,1) = xug
     END DO
   END DO
 END IF
 IF(yl == ymin) THEN
   DO k = zlg, zug
     DO i = xlg, xug
       ni(i,ymin,k,4) = ylg
     END DO
   END DO
 END IF
 IF(yu == ymax) THEN
   DO k = zlg, zug
     DO i = xlg, xug
       ni(i,ymax,k,3) = yug
     END DO
   END DO
 END IF
 IF(zl == zmin) THEN
   DO j = ylg, yug
     DO i = xlg, xug
       ni(i,j,zmin,6) = zlg
     END DO
   END DO
 END IF
 IF(zu == zmax) THEN
   DO j = ylg, yug
     DO i = xlg, xug
       ni(i,j,zmax,5) = zug
     END DO
   END DO
 END IF


!--------- Initialize distribution functions -----------------------------------
 DO k = zl, zu
   DO j = yl, yu
     DO i = xl, xu

! Assign the neighbours
       ie = ni(i,j,k,1)
       iw = ni(i,j,k,2)
       jn = ni(i,j,k,3)
       js = ni(i,j,k,4)
       kt = ni(i,j,k,5)
       kb = ni(i,j,k,6)

! Laplacian of the order parameter
       lapphi = ( phi(ie,jn,k) + phi(iw,js,k) + phi(ie,js,k) + phi(iw,jn,k) &
              +   phi(ie,j,kt) + phi(iw,j,kb) + phi(ie,j,kb) + phi(iw,j,kt) &
              +   phi(i,jn,kt) + phi(i,js,kb) + phi(i,jn,kb) + phi(i,js,kt) &
              +   2.0D0*( phi(ie,j,k) + phi(iw,j,k) + phi(i,jn,k)           &
              +   phi(i,js,k) + phi(i,j,kt) + phi(i,j,kb) - 12.d0*phi(i,j,k) ) )*inv6

!  Chemical potential
       phin   = phi(i,j,k)
       rhon   = 0.5D0*(rhoH + rhoL)
       muPhin = alpha4*phin*(phin*phin - phiStar*phiStar) - kappa*lapPhi

! Set distribution function f to its equilibrium value
       Af0 = -3.D0*Gamma*muPhin
       Af1 = 0.5D0*Gamma*muPhin
       f(i,j,k,0,now)   = Af0 + phin
       f(i,j,k,1:6,now) = Af1


! Set distribution functiong to its equilibrium value
       Ag0 = Eg0*(rhon - 6.D0*phin*muphin)
       Ag1 = Eg1*(3.D0*phin*muphin + rhon)
       Ag2 = Eg2*(3.D0*phin*muphin + rhon)
       g(i,j,k,0,now)    = Ag0
       g(i,j,k,1:6,now)  = Ag1
       g(i,j,k,7:18,now) = Ag2

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
 CALL MPI_IRECV(f_east_rcv, xsize, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(f_west_rcv, xsize, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 count = 0
 DO k = zl, zu
   DO j = yl, yu
     count = count + 1
     f_west_snd(count) = f(xl,j,k,1,now)
     f_east_snd(count) = f(xu,j,k,2,now)
   END DO
 END DO
 CALL MPI_SEND(f_west_snd, xsize, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(f_east_snd, xsize, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)


! Y direction
 CALL MPI_IRECV(f_north_rcv, ysize, MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(f_south_rcv, ysize, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 count = 0
 DO k = zl, zu
   DO i = xl, xu
     count = count + 1
     f_south_snd(count) = f(i,yl,k,3,now)
     f_north_snd(count) = f(i,yu,k,4,now)
   END DO
 END DO
 CALL MPI_SEND(f_south_snd, ysize, MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(f_north_snd, ysize, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)

! Z direction
 CALL MPI_IRECV(f_top_rcv, zsize, MPI_DOUBLE_PRECISION, top, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(f_bot_rcv, zsize, MPI_DOUBLE_PRECISION, bot, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 count = 0
 DO j = yl, yu
   DO i = xl, xu
     count = count + 1
     f_bot_snd(count) = f(i,j,zl,5,now)
     f_top_snd(count) = f(i,j,zu,6,now)
   END DO
 END DO
 CALL MPI_SEND(f_bot_snd, zsize, MPI_DOUBLE_PRECISION, bot, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(f_top_snd, zsize, MPI_DOUBLE_PRECISION, top, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)

!-------- Update values of and f in ghost layers -------------------------------
! X direction
 count = 0
 DO k = zl, zu
   DO j = yl, yu
     count = count + 1
     f(xlg,j,k,2,now) = f_west_rcv(count)
     f(xug,j,k,1,now) = f_east_rcv(count)
   END DO
 END DO

!  Y direction
 count = 0
 DO k = zl, zu
   DO i = xl, xu
     count = count + 1
     f(i,ylg,k,4,now) = f_south_rcv(count)
     f(i,yug,k,3,now) = f_north_rcv(count)
   END DO
 END DO

!  Z direction
 count = 0
 DO j = yl, yu
   DO i = xl, xu
     count = count + 1
     f(i,j,zlg,6,now) = f_bot_rcv(count)
     f(i,j,zug,5,now) = f_top_rcv(count)
   END DO
 END DO

!-------- Save initialized data ------------------------------------------------
 CALL Stats
 CALL VtkDump

 RETURN
 END SUBROUTINE Init
