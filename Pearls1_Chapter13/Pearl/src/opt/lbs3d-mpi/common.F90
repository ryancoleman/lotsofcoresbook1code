!-------------------------------------------------------------------------------
! File       : Common
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
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

!> @brief Definition of single and double data types
 MODULE NTypes
 IMPLICIT NONE
 SAVE

 INTEGER, PARAMETER :: SGL = KIND(1.0)
 INTEGER, PARAMETER :: DBL = KIND(1.D0)

 END MODULE NTypes

!> @brief Timing variables
 MODULE Timers
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

 INTEGER :: nIO
 REAL(KIND = DBL) :: tCurrent, tStart, tSetup, tRelax, tEvol, tIO, tRun

 END MODULE Timers

!> @brief Parameters related to the geometry and the time intervals
 MODULE Domain
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

! Time step and limits
 INTEGER :: iStep, tCall, tSave, tStat, MaxStep, RelaxStep, STAGE, &
            CompletedRelaxSteps

! Domain size
! (xmax, xmin) : computational domain size limits
! (xl  , xu  ) : local minimun and maximum values for x in a given processor
! (xlg , xug ) : local ghost node position (xlg = xl-1, xug = xu+1)
 INTEGER :: xo, yo, zo, now, nxt
 INTEGER :: NX, NY, NZ, NPX, NPy, NPZ, NXG, NYG, NZG, NG
 INTEGER :: xmin, xmax, xl, xu, xlg, xug
 INTEGER :: ymin, ymax, yl, yu, ylg, yug
 INTEGER :: zmin, zmax, zl, zu, zlg, zug
 INTEGER :: xsize, xsize2, xsize3, xsize4, xsize5, xsize6
 INTEGER :: ysize, ysize2, ysize3, ysize4, ysize5, ysize6
 INTEGER :: zsize, zsize2, zsize3, zsize4, zsize5, zsize6
 INTEGER :: xedge, xedge2, yedge, yedge2, zedge, zedge2

! ndim = spatial dimension
 INTEGER, PARAMETER :: ndim = 3

! Inverse of the initial volume of the discrete phase 
 REAL(KIND = DBL) :: invInitVol

! Define constants used in the code (inv6 = 1/6, inv12 = 1/12, invPi = 1/pi)
 REAL(KIND = DBL), PARAMETER :: inv6  = 0.16666666666666666667D0
 REAL(KIND = DBL), PARAMETER :: inv12 = 0.08333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: invPi = 0.31830988618379067154D0

! Nearest neighbor array
! REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:,:,:) :: nb

 END MODULE Domain

!> @brief Parameters related to hydrodynamics quantities
 MODULE FluidParams
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

 INTEGER :: nBubbles
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:) :: bubbles
 REAL(KIND = DBL), DIMENSION(1:11) :: Convergence
 REAL(KIND = DBL) :: sigma, IntWidth, Gamma
 REAL(KIND = DBL) :: alpha, alpha4, kappa
 REAL(KIND = DBL) :: tauPhi, invTauPhi, invTauPhi1, phiStar, phiStar2, phiStar4
 REAL(KIND = DBL) :: rhoL, rhoH
 REAL(KIND = DBL) :: tauRho, invTauRho, invTauRhoOne, invTauRhoHalf
 REAL(KIND = DBL) :: eta, eta2, invEta2
 REAL(KIND = DBL) :: gravity, grav, Eo, eps, pConv


 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi, lapPhi
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: gradPhiX, gradPhiY, gradPhiZ

 END MODULE FluidParams

!> @brief Parameters related to the LBM discretization
 MODULE LBMParams
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

! fdim = order parameter dimension - 1 (D3Q7)
! gdim = momentum dimension - 1        (D3Q19)
 INTEGER, PARAMETER :: fdim = 6
 INTEGER, PARAMETER :: gdim = 18

! Distribution functions
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f0, f1, f2, f3, f4, f5, f6

 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g0, g1, g2, g3, g4, g5, g6,    &
                                                g7, g8, g9, g10, g11, g12, g13, &
                                                g14, g15, g16, g17, g18

! D3Q7/D3Q19 Lattice speed of sound (Cs = 1/DSQRT(3), Cs_sq = 1/3, invCs_sq = 3)
 REAL(KIND = DBL), PARAMETER :: Cs       = 0.57735026918962576451D0
 REAL(KIND = DBL), PARAMETER :: Cs_sq    = 0.33333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: invCs_sq = 3.00000000000000000000D0

! Distributions weights (D3Q19: Eg0 = 1/3, Eg1 = 1/18, Eg2 = 1/36)
 REAL(KIND = DBL), PARAMETER :: Eg0 = 0.33333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: Eg1 = 0.05555555555555555556D0
 REAL(KIND = DBL), PARAMETER :: Eg2 = 0.02777777777777777778D0

! Modified distribution weight (to avoid operations: EgiC = Egi*invCs_sq)
 REAL(KIND = DBL), PARAMETER :: Eg0C = 1.00000000000000000000D0
 REAL(KIND = DBL), PARAMETER :: Eg1C = 0.16666666666666666667D0
 REAL(KIND = DBL), PARAMETER :: Eg2C = 0.08333333333333333333D0
 REAL(KIND = DBL) :: w, Eg0n, Eg1n, Eg2n, EgC0n, EgC1n, EgC2n

 END MODULE LBMParams

!> @brief MPI related parameters and information exchange buffer arrays
!> @details
!! @param nprocs  : total number of processors
!! @param proc    : processor ID
!! @param vproc   : processor ID in virtual grid
!! @param mpi_dim : virtual grid partition scheme (1->stripes, 2->boxes, 3->cubes)
 MODULE MPIParams
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

! Communication parameters
 INTEGER :: nprocs, proc, vproc, MPI_COMM_VGRID
 INTEGER :: mpi_xdim, mpi_ydim, mpi_zdim
 INTEGER :: east, west, north, south, top, bot
 INTEGER :: ne, nw, se, sw, te, tw, be, bw, tn, ts, bn, bs
 INTEGER, PARAMETER :: master  = 0
 INTEGER, PARAMETER :: mpi_dim = 3
 INTEGER, PARAMETER :: TAG1 = 1, TAG2 = 2, TAG3 = 3, TAG4 = 4, TAG5 = 5, TAG6 = 6, &
                       TAG7 = 7, TAG8 = 8, TAG9 = 9, TAG10 = 10, TAG11 = 11, &
                       TAG12 = 12, TAG13 = 13, TAG14 = 14, TAG15 = 15, TAG16 = 16, &
                       TAG17 = 17, TAG18 = 18

! Information exchange buffers (x direction)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_west_snd, f_east_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_west_rcv, f_east_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: west_snd, east_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: west_rcv, east_rcv

! Information exchange buffers (y direction)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_south_snd, f_north_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_south_rcv, f_north_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: south_snd, north_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: south_rcv, north_rcv

! Information exchange buffers (z direction)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_bot_snd, f_top_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_bot_rcv, f_top_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: bot_snd, top_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: bot_rcv, top_rcv

! Information exchange buffers (diagonals x)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: tn_snd, tn_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: bn_snd, bn_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: ts_snd, ts_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: bs_snd, bs_rcv

! Information exchange buffers (diagonals y)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: te_snd, te_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: be_snd, be_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: tw_snd, tw_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: bw_snd, bw_rcv

! Information exchange buffers (diagonals z)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: ne_snd, ne_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: se_snd, se_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: nw_snd, nw_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: sw_snd, sw_rcv


 END MODULE MPIParams

