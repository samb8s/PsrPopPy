      SUBROUTINE QSIMP(FUNC,A,B,S,constant,braking)
      external func
      PARAMETER (EPS=1.E-3, JMAX=10)
      OST=-1.E30
      OS= -1.E30
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,ST,J,constant,braking)
        S=(4.*ST-OST)/3.
        if(j.gt.3) then
           IF (ABS(S-OS).LT.EPS*ABS(OS).or.(s.eq.0..and.os.eq.0.))
     &          RETURN
        endif
        OS=S
        OST=ST
11    CONTINUE
      s = 1e15
      return
      END
