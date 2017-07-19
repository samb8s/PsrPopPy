      SUBROUTINE TRAPZD(FUNC,A,B,S,N,constant,braking)
      external func
      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A,constant,braking)+FUNC(B,constant,braking))
        IT=1
      ELSE
        it=2**(n-2)
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X,constant,braking)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
