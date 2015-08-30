C> \ingroup stpr_priv
C> @{
      SUBROUTINE stpr_clnstp(STEP,THRES)
c $Id: stpr_clnstp.F 25540 2014-05-01 22:33:27Z d3y133 $
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
      COMMON / CFACE / IWCTR,NATOM,ICALC
      COMMON / DIMS / NAT3, NAT3SQ, NAT3TR
      DIMENSION STEP(NAT3)
      DO 10 I=1,NAT3
        STEPX = STEP(I)/THRES
        ISTEPX = NINT(STEPX)
        STEPX = dble(ISTEPX)*THRES
        STEP(I) = STEPX
   10 CONTINUE
      RETURN
      END
C> @}
