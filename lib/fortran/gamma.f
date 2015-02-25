      subroutine gamma1(m,Tobs,m1input,m2input,Ps,omperinput,
     & ideg,ec,Pod, gammamavg)


       IMPLICIT NONE
       real*8 Msun,pi,G,c
       real*8 m,Tobs
       real*8 Omo,sini,ideg,m1,m2,omper,ec
       real*8 aR,a,ap
       real*8 Ps,Omp
       real*8 gammamavg,sumgamma,sumprob
       real*8 gammamarray(40)
       integer loop1
       integer lpTp1,lpTp2
       real*8 M0ar(40),E0ar(40),Tp(40)
       real*8 fcd(40),fc(40),prob(40),probn(40)
       
       integer loopPo,loopPs
       real*8 Po,Pod,hpd,hps
       integer npd,nps
 
       real*8 m1input,m2input,omperinput
CCCCCCC  constants in SI  CCCCCCCCC
       c=2.99792458D+8       
       G=6.67384D-11
       pi=3.141592653589793D0
       Msun=1.98892D+30
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC6.67300
CCC   inputs
CCC  m1input,m2input in Msun. incinput,omperinput in deg; Pod in day, Ps in sec

c      OPEN(unit=1,file="params.in",status="OLD")

       
c      read(1,*)m,Tobs,m1input,m2input,Ps,omperinput,ideg,ec,Pod
    
       

c      write(*,*) "PROGRAM IS RUNNING!"
      sini=SIN(ideg*pi/180.0D0)
       
CCCCCCCCCCCCCCCCCCCCc

       m1=m1input*Msun
       m2=m2input*Msun       
       

       omper=omperinput*pi/180.0D0
       

       Po=Pod*24.0D0*3600.0D0
       Omo=2.0D0*pi/Po
       Omp=2.0D0*pi/Ps

       aR=((G*(m1+m2))/(Omo*Omo))**0.33333333333333333333D0
       a=aR*sini
       ap=(a*m2)/(m1+m2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC -------- variation of Tp and weighted avg starts here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC ---- initialisation of the array to avoid accidental junk
       DO loop1=1,37
        gammamarray(loop1)=0.0D0
        fcd(loop1)=0.0D0
        fc(loop1)=0.0D0
        prob(loop1)=0.0D0
        E0ar(loop1)=0.0D0
        M0ar(loop1)=0.0D0
        Tp(loop1)=0.0D0
        probn(loop1)=0.0D0
        
       ENDDO
       
CC ------ probability calculation ------      
       DO lpTp1=1,37
 
        fcd(lpTp1)=0.0D0+dble(lpTp1-1)*10.0D0
        fc(lpTp1)=fcd(lpTp1)*pi/180.0D0
 
       prob(lpTp1)=(1.0D0)/
     &((1.0D0+ec*cos(fc(lpTp1)))*(1.0D0+ec*cos(fc(lpTp1))))

      E0ar(lpTp1)=
     &2.0D0*atan(SQRT((1.0D0-ec)/(1.0D0+ec))*tan(fc(lpTp1)/2.0D0))

      M0ar(lpTp1)=E0ar(lpTp1)-ec*sin(E0ar(lpTp1))

      Tp(lpTp1)=-M0ar(lpTp1)/Omo
 
       ENDDO
       

C   === normalisation of probabilities  ===
       sumgamma=0.0D0
       sumprob=0.0D0
       DO lpTp2=1,37
       probn(lpTp2)=prob(lpTp2)/prob(19)
         
        call findalpha1
     &(m,Omp,ap,omper,ec,Omo,Tp(lpTp2),Tobs,gammamarray(lpTp2))

        sumgamma=sumgamma+probn(lpTp2)*gammamarray(lpTp2)
        sumprob=sumprob+probn(lpTp2)
        ENDDO
        
        gammamavg=sumgamma/sumprob
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC -------- variation of Tp and weighted avg ends here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
       
c         write(*,*)gammamavg



c       STOP
c      RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gamma2(m,Tobs,m1input,m2input,Ps,omperinput,
     & ideg,ec,Pod, gammamavg)


       IMPLICIT NONE
       real*8 Msun,pi,G,c
       real*8 m,Tobs
       real*8 Omo,sini,ideg,m1,m2,omper,ec
       real*8 aR,a,ap
       real*8 Ps,Omp
       real*8 gammamavg,sumgamma,sumprob
       real*8 gammamarray(40)
       integer loop1
       integer lpTp1,lpTp2
       real*8 M0ar(40),E0ar(40),Tp(40)
       real*8 fcd(40),fc(40),prob(40),probn(40)
       
       integer loopPo,loopPs
       real*8 Po,Pod,hpd,hps
       integer npd,nps
     
       real*8 m1input,m2input,omperinput
CCCCCCC  constants in SI  CCCCCCCCC
       c=2.99792458D+8       
       G=6.67384D-11
       pi=3.141592653589793D0
       Msun=1.98892D+30
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC6.67300 
CCC   inputs
CCC  m1input,m2input in Msun. incinput,omperinput in deg; Pod in day, Ps in sec


      sini=SIN(ideg*pi/180.0D0)
       
    
CCCCCCCCCCC orbital and stellar parameters ccc
       m=4.0D0
       Tobs=1000.0D0

       m1=m1input*Msun
       m2=m2input*Msun       
       

       omper=omperinput*pi/180.0D0
       

       Po=Pod*24.0D0*3600.0D0
       Omo=2.0D0*pi/Po
       Omp=2.0D0*pi/Ps         
  
   
       aR=((G*(m1+m2))/(Omo*Omo))**0.33333333333333333333D0
       a=aR*sini
       ap=(a*m2)/(m1+m2) 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC -------- variation of Tp and weighted avg starts here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC ---- initialisation of the array to avoid accidental junk
       DO loop1=1,37
        gammamarray(loop1)=0.0D0
        fcd(loop1)=0.0D0
        fc(loop1)=0.0D0
        prob(loop1)=0.0D0
        E0ar(loop1)=0.0D0
        M0ar(loop1)=0.0D0
        Tp(loop1)=0.0D0
        probn(loop1)=0.0D0
        
       ENDDO
       
CC ------ probability calculation ------       
       DO lpTp1=1,37
 
        fcd(lpTp1)=0.0D0+dble(lpTp1-1)*10.0D0
        fc(lpTp1)=fcd(lpTp1)*pi/180.0D0
 
       prob(lpTp1)=(1.0D0)/
     &((1.0D0+ec*cos(fc(lpTp1)))*(1.0D0+ec*cos(fc(lpTp1))))

      E0ar(lpTp1)=
     &2.0D0*atan(SQRT((1.0D0-ec)/(1.0D0+ec))*tan(fc(lpTp1)/2.0D0))

      M0ar(lpTp1)=E0ar(lpTp1)-ec*sin(E0ar(lpTp1))

      Tp(lpTp1)=-M0ar(lpTp1)/Omo
 
       ENDDO
       

C   === normalisation of probabilities  ===
       sumgamma=0.0D0
       sumprob=0.0D0
       DO lpTp2=1,37
       probn(lpTp2)=prob(lpTp2)/prob(19)
         
        call findalpha2
     &(m,Omp,ap,omper,ec,Omo,Tp(lpTp2),Tobs,gammamarray(lpTp2))

        sumgamma=sumgamma+probn(lpTp2)*gammamarray(lpTp2)
        sumprob=sumprob+probn(lpTp2)
        ENDDO
        
        gammamavg=sumgamma/sumprob
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC -------- variation of Tp and weighted avg ends here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       

        
c         write(*,*)gammamavg


       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gamma3(m,Tobs,m1input,m2input,Ps,omperinput,
     & ideg,ec,Pod, gammamavg)
       IMPLICIT NONE
       real*8 Msun,pi,G,c
       real*8 m,Tobs
       real*8 Omo,sini,ideg,m1,m2,omper,ec
       real*8 aR,a,ap
       real*8 Ps,Omp
       real*8 gammamavg,sumgamma,sumprob
       real*8 gammamarray(40)
       integer loop1
       integer lpTp1,lpTp2
       real*8 M0ar(40),E0ar(40),Tp(40)
       real*8 fcd(40),fc(40),prob(40),probn(40)
       
       integer loopPo,loopPs
       real*8 Po,Pod,hpd,hps
       integer npd,nps
       real*8 m1input,m2input,omperinput 
CCCCCCC  constants in SI  CCCCCCCCC
       c=2.99792458D+8       
       G=6.67384D-11
       pi=3.141592653589793D0
       Msun=1.98892D+30
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC6.67300
CCC   inputs
CCC  m1input,m2input in Msun. incinput,omperinput in deg; Pod in day, Ps in sec

     
      sini=SIN(ideg*pi/180.0D0)
       
CCCCCCCCCCC orbital and stellar parameters ccc
       m=4.0D0
       Tobs=1000.0D0

       m1=m1input*Msun
       m2=m2input*Msun       
       

       omper=omperinput*pi/180.0D0
       

       Po=Pod*24.0D0*3600.0D0
       Omo=2.0D0*pi/Po
       Omp=2.0D0*pi/Ps         
  
   
       aR=((G*(m1+m2))/(Omo*Omo))**0.33333333333333333333D0
       a=aR*sini
       ap=(a*m2)/(m1+m2) 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC -------- variation of Tp and weighted avg starts here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC ---- initialisation of the array to avoid accidental junk
       DO loop1=1,37
        gammamarray(loop1)=0.0D0
        fcd(loop1)=0.0D0
        fc(loop1)=0.0D0
        prob(loop1)=0.0D0
        E0ar(loop1)=0.0D0
        M0ar(loop1)=0.0D0
        Tp(loop1)=0.0D0
        probn(loop1)=0.0D0
        
       ENDDO
       
CC ------ probability calculation ------       
       DO lpTp1=1,37
 
        fcd(lpTp1)=0.0D0+dble(lpTp1-1)*10.0D0
        fc(lpTp1)=fcd(lpTp1)*pi/180.0D0
 
       prob(lpTp1)=(1.0D0)/
     &((1.0D0+ec*cos(fc(lpTp1)))*(1.0D0+ec*cos(fc(lpTp1))))

      E0ar(lpTp1)=
     &2.0D0*atan(SQRT((1.0D0-ec)/(1.0D0+ec))*tan(fc(lpTp1)/2.0D0))

      M0ar(lpTp1)=E0ar(lpTp1)-ec*sin(E0ar(lpTp1))

      Tp(lpTp1)=-M0ar(lpTp1)/Omo
 
       ENDDO
       

C   === normalisation of probabilities  ===
       sumgamma=0.0D0
       sumprob=0.0D0
       DO lpTp2=1,37
       probn(lpTp2)=prob(lpTp2)/prob(19)
         
        call findalpha3
     &(m,Omp,ap,omper,ec,Omo,Tp(lpTp2),Tobs,gammamarray(lpTp2))

        sumgamma=sumgamma+probn(lpTp2)*gammamarray(lpTp2)
        sumprob=sumprob+probn(lpTp2)
        ENDDO
        
        gammamavg=sumgamma/sumprob
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC -------- variation of Tp and weighted avg ends here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       

       END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine 
     &findalpha1(xmy,Ompy,apy,ompery,ecy,Omoy,Tpy,Tob,gammaf)
       IMPLICIT NONE
       real*8 xmy,Ompy,apy,ompery,ecy,Omoy,Tpy,Tob       
       real*8 ta,tstep
       REAL*8 pi
       REAL*8 Ompx,apx,ecx,omperx,Omox,xm,Tpx
   
       COMMON /orbit/Ompx,apx,ecx,omperx,Omox,xm,Tpx

       integer i,ii,j,jj,k

       real*8 Marr(1001),Earr(1001),farr(1001),varrl(1001)
       real*8 Marrn(101),Earrn(101),farrn(101),varrln(101)
       real*8 alpha2rr(1001),gammaarr(1001)
       real*8 ntagarr(1001),ntagarrn(101)

       real*8 alpha2rrn(101),gammaarrn(101)
             
       real*8 tarr(1001),tarrn(101)

       integer nnn  
       real*8 tn0,tnstep

       real*8 alpha2f,gammaf 
       real*8 E0y,M0y,f0y,rl0
       real*8 fac1
C-----------------------------------------------
         pi=3.141592653589793D0
C------ This part is needed because the arguments passed by can't be put in common block
        xm=xmy
        Ompx=Ompy
        Omox=Omoy
        apx=apy
        ecx=ecy
        omperx=ompery
        Tpx=Tpy      
C----------------------------------------------------
       fac1=(Omoy)*(apy/SQRT(1.0D0-ecy*ecy))

       ta=0.0D0
       tstep=(Tob-ta)/1000.0D0
   
CCCCCC  setting array elements to zero to avoid accidental junk creeping in
        nnn=0
        DO j=1,1001
         alpha2rr(j)=0.0D0
         gammaarr(j)=0.0D0
         ntagarr(j)=0.0D0
         tarr(j)=0.0D0
         Marr(j)=0.0D0
         Earr(j)=0.0D0
         farr(j)=0.0D0
         varrl(j)=0.0D0
        ENDDO
       

         M0y=-Omoy*Tpy

        IF(M0y.GE.(2.0D0*pi))M0y=M0y-(2.0D0*pi)
      
        call keplersolve1(ecy,M0y,E0y)
        call keplersolve2(ecy,E0y,f0y)
       
      rl0=((apy*(1.0D0-ecy*ecy))/(1.0D0+ecy*cos(f0y)))*sin(f0y+ompery)
       

CCCCCCCCCCCC  initial alpha2, gamma calculation over the full observation   CCCCCC
       DO i=1,1001 
         tarr(i)=ta+dble(i-1)*tstep

         Marr(i)=Omoy*(tarr(i)-Tpy)

        IF(Marr(i).GE.(2.0D0*pi))Marr(i)=Marr(i)-(2.0D0*pi)
      
        call keplersolve1(ecy,Marr(i),Earr(i))
        call keplersolve2(ecy,Earr(i),farr(i))

       
       varrl(i)=fac1*(cos(farr(i)+ompery)+ecy*cos(ompery))

      alpha2rr(i)=varrl(i)
      call gam1(alpha2rr(i),rl0,Tob,gammaarr(i))
      ntagarr(i)=dble(i)
      
      ENDDO

       
       call sort2(1001,gammaarr,ntagarr)

C---------------------


       nnn=ntagarr(1001)
   

        IF(nnn.EQ.1)THEN
         tn0=tarr(nnn)
         tnstep=(tarr(nnn+1)-tarr(nnn))/50.0D0
        ELSEIF(nnn.EQ.1001)THEN
         tn0=tarr(nnn-1)
         tnstep=(tarr(nnn)-tarr(nnn-1))/50.0D0
        ELSE
        tn0=tarr(nnn-1)
        tnstep=(tarr(nnn+1)-tarr(nnn-1))/100.0D0
       ENDIF

 1     CONTINUE
        
CCCCCC  setting array elements to zero to avoid accidental junk creeping in
        nnn=0
        DO jj=1,101
         tarrn(jj)=0.0D0
         alpha2rrn(jj)=0.0D0
         gammaarrn(jj)=0.0D0
         ntagarrn(jj)=0.0D0
         Marrn(jj)=0.0D0
         Earrn(jj)=0.0D0
         farrn(jj)=0.0D0
         varrln(jj)=0.0D0
        ENDDO

CCC   finer search to maximize gamma
        
        DO ii=1,101     
                   
          tarrn(ii)=tn0+dble(ii-1)*tnstep

         Marrn(ii)=Omoy*(tarrn(ii)-Tpy)

        IF(Marrn(ii).GE.(2.0D0*pi))Marrn(ii)=Marrn(ii)-(2.0D0*pi)

     
        call keplersolve1(ecy,Marrn(ii),Earrn(ii))
        call keplersolve2(ecy,Earrn(ii),farrn(ii))


       varrln(ii)=fac1*(cos(farrn(ii)+ompery)+ecy*cos(ompery))

       alpha2rrn(ii)=varrln(ii)
CCCCCCCCC
          call gam1(alpha2rrn(ii),rl0,Tob,gammaarrn(ii))
          ntagarrn(ii)=dble(ii)
        ENDDO


         call sort2(101,gammaarrn,ntagarrn)

         nnn=ntagarrn(101)


        IF(nnn.EQ.1)THEN
         tn0=tarrn(nnn)
         tnstep=(tarrn(nnn+1)-tarrn(nnn))/50.0D0
        ELSEIF(nnn.EQ.101)THEN
         tn0=tarrn(nnn-1)
         tnstep=(tarrn(nnn)-tarrn(nnn-1))/50.0D0
        ELSE
        tn0=tarrn(nnn-1)
        tnstep=(tarrn(nnn+1)-tarrn(nnn-1))/100.0D0
       ENDIF


         if(tnstep.gt.1.0D-5)then
           goto 1
         else
           goto 2
        endif
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  2     CONTINUE
         gammaf=gammaarrn(101)
       
  
CCCCCCCCCCCCCCCCCCCCCC

       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
      subroutine 
     &findalpha2(xmy,Ompy,apy,ompery,ecy,Omoy,Tpy,Tob,gammaf)
       IMPLICIT NONE
       real*8 xmy,Ompy,apy,ompery,ecy,Omoy,Tpy,Tob       
       real*8 ta,tstep
       REAL*8 pi
       REAL*8 Ompx,apx,ecx,omperx,Omox,xm,Tpx
   
       COMMON /orbit/Ompx,apx,ecx,omperx,Omox,xm,Tpx

       integer i,ii,j,jj,k

       real*8 Marr(1001),Earr(1001),farr(1001),varrl(1001),acarrl(1001)
       real*8 Marrn(101),Earrn(101),farrn(101),varrln(101),acarrln(101)
       real*8 alpha1rr(1001),alpha2rr(1001),gammaarr(1001)
       real*8 ntagarr(1001),ntagarrn(101)
       real*8 alpha1rrn(101),alpha2rrn(101),gammaarrn(101)             
       real*8 tarr(1001),tarrn(101)  
       integer nnn  
       real*8 tn0,tnstep
       real*8 alpha1f,alpha2f,gammaf 
       real*8 E0y,M0y,f0y,rl0
       real*8 fac1, fac2
C--------------------------------------------------
       pi=3.141592653589793D0
C------ This part is needed because the arguments passed by can't be put in common block
       xm=xmy
       Ompx=Ompy
       Omox=Omoy
       apx=apy
       ecx=ecy
       omperx=ompery
       Tpx=Tpy      
C----------------------------------------------------
       fac1=(Omoy)*(apy/SQRT(1.0D0-ecy*ecy))
       fac2=-(Omoy*Omoy)*(apy/((1.0D0-ecy*ecy)*(1.0D0-ecy*ecy)))
 
       ta=0.0D0
       tstep=(Tob-ta)/1000.0D0
   
CCCCCC  setting array elements to zero to avoid accidental junk creeping in
        nnn=0
        DO j=1,1001
         alpha1rr(j)=0.0D0
         alpha2rr(j)=0.0D0
         gammaarr(j)=0.0D0
         ntagarr(j)=0.0D0
         tarr(j)=0.0D0
         Marr(j)=0.0D0
         Earr(j)=0.0D0
         farr(j)=0.0D0
         varrl(j)=0.0D0
         acarrl(j)=0.0D0
        ENDDO
       

         M0y=-Omoy*Tpy

        IF(M0y.GE.(2.0D0*pi))M0y=M0y-(2.0D0*pi)
      
        call keplersolve1(ecy,M0y,E0y)
        call keplersolve2(ecy,E0y,f0y)
       
      rl0=((apy*(1.0D0-ecy*ecy))/(1.0D0+ecy*cos(f0y)))*sin(f0y+ompery)
       

CCCCCCCCCCCC  initial alpha2, gamma calculation over the full observation   CCCCCC
       DO i=1,1001 
         tarr(i)=ta+dble(i-1)*tstep
        
         Marr(i)=Omoy*(tarr(i)-Tpy)

        IF(Marr(i).GE.(2.0D0*pi))Marr(i)=Marr(i)-(2.0D0*pi)
      
        call keplersolve1(ecy,Marr(i),Earr(i))
        call keplersolve2(ecy,Earr(i),farr(i))

       
       varrl(i)=fac1*(cos(farr(i)+ompery)+ecy*cos(ompery))

      acarrl(i)=fac2*sin(farr(i)+ompery)*
     &(1.0D0+ecy*cos(farr(i)))*(1.0D0+ecy*cos(farr(i)))
 

      alpha2rr(i)=varrl(i)
      alpha1rr(i)=acarrl(i)/2.0D0
      call gam2(alpha1rr(i),alpha2rr(i),rl0,Tob,gammaarr(i))
      ntagarr(i)=dble(i)
      
      ENDDO

       
       call sort2(1001,gammaarr,ntagarr)

C---------------------


       nnn=ntagarr(1001)
   

        IF(nnn.EQ.1)THEN
         tn0=tarr(nnn)
         tnstep=(tarr(nnn+1)-tarr(nnn))/50.0D0
        ELSEIF(nnn.EQ.1001)THEN
         tn0=tarr(nnn-1)
         tnstep=(tarr(nnn)-tarr(nnn-1))/50.0D0
        ELSE
        tn0=tarr(nnn-1)
        tnstep=(tarr(nnn+1)-tarr(nnn-1))/100.0D0
       ENDIF

 1     CONTINUE
        
CCCCCC  setting array elements to zero to avoid accidental junk creeping in
        nnn=0
        DO jj=1,101
         tarrn(jj)=0.0D0
         alpha1rrn(jj)=0.0D0
         alpha2rrn(jj)=0.0D0
         gammaarrn(jj)=0.0D0
         ntagarrn(jj)=0.0D0
         Marrn(jj)=0.0D0
         Earrn(jj)=0.0D0
         farrn(jj)=0.0D0
         varrln(jj)=0.0D0
         acarrln(jj)=0.0D0
        ENDDO

CCC   finer search to maximize gamma
        
        DO ii=1,101     
                   
          tarrn(ii)=tn0+dble(ii-1)*tnstep
 
         Marrn(ii)=Omoy*(tarrn(ii)-Tpy)

        IF(Marrn(ii).GE.(2.0D0*pi))Marrn(ii)=Marrn(ii)-(2.0D0*pi)

     
        call keplersolve1(ecy,Marrn(ii),Earrn(ii))
        call keplersolve2(ecy,Earrn(ii),farrn(ii))


       varrln(ii)=fac1*(cos(farrn(ii)+ompery)+ecy*cos(ompery))


       acarrln(ii)=fac2*sin(farrn(ii)+ompery)*
     &(1.0D0+ecy*cos(farrn(ii)))*(1.0D0+ecy*cos(farrn(ii)))

       alpha2rrn(ii)=varrln(ii)

       alpha1rrn(ii)=acarrln(ii)/2.0D0
     
CCCCCCCCC
          call gam2(alpha1rrn(ii),alpha2rrn(ii),rl0,Tob,gammaarrn(ii))
          ntagarrn(ii)=dble(ii)
        ENDDO


         call sort2(101,gammaarrn,ntagarrn)

         nnn=ntagarrn(101)

    
        IF(nnn.EQ.1)THEN
         tn0=tarrn(nnn)
         tnstep=(tarrn(nnn+1)-tarrn(nnn))/50.0D0
        ELSEIF(nnn.EQ.101)THEN
         tn0=tarrn(nnn-1)
         tnstep=(tarrn(nnn)-tarrn(nnn-1))/50.0D0
        ELSE
        tn0=tarrn(nnn-1)
        tnstep=(tarrn(nnn+1)-tarrn(nnn-1))/100.0D0
       ENDIF


         if(tnstep.gt.1.0D-5)then
           goto 1
         else
           goto 2
        endif
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  2     CONTINUE
         gammaf=gammaarrn(101)
  
CCCCCCCCCCCCCCCCCCCCCC

       END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine 
     &findalpha3(xmy,Ompy,apy,ompery,ecy,Omoy,Tpy,Tob,gammaf)
       IMPLICIT NONE
       real*8 xmy,Ompy,apy,ompery,ecy,Omoy,Tpy,Tob       
       real*8 ta,tstep
       REAL*8 pi
       REAL*8 Ompx,apx,ecx,omperx,Omox,xm,Tpx
   
       COMMON /orbit/Ompx,apx,ecx,omperx,Omox,xm,Tpx

       integer i,ii,j,jj,k

       real*8 Marr(1001),Earr(1001),farr(1001),varrl(1001),acarrl(1001)
       real*8 Marrn(101),Earrn(101),farrn(101),varrln(101),acarrln(101)

        real*8 jerkl(1001),jerkln(101)

       real*8 alpha0rr(1001),alpha1rr(1001),alpha2rr(1001)
       real*8 gammaarr(1001),gammaarrn(101)
       real*8 ntagarr(1001),ntagarrn(101)

       real*8 alpha0rrn(101),alpha1rrn(101),alpha2rrn(101)
             
       real*8 tarr(1001),tarrn(101)
       integer nnn  
       real*8 tn0,tnstep

       real*8 gammaf,E0y,M0y,f0y,rl0 
       real*8 fac1,fac2,fac3
C----------------------------------------------      
         pi=3.141592653589793D0
C------ This part is needed because the arguments passed by can't be put in common block
        xm=xmy
        Ompx=Ompy
        Omox=Omoy
        apx=apy
        ecx=ecy
        omperx=ompery
        Tpx=Tpy      
C----------------------------------------------------
       fac1=(Omoy)*(apy/SQRT(1.0D0-ecy*ecy))
       fac2=-(Omoy*Omoy)*(apy/((1.0D0-ecy*ecy)*(1.0D0-ecy*ecy)))
       fac3=-(Omoy*Omoy*Omoy)*(apy/((1.0D0-ecy*ecy)**3.5D0))
       ta=0.0D0
       tstep=(Tob-ta)/1000.0D0
   
CCCCCC  setting array elements to zero to avoid accidental junk creeping in
        nnn=0
        DO j=1,1001
         alpha0rr(j)=0.0D0
         alpha1rr(j)=0.0D0
         alpha2rr(j)=0.0D0
         gammaarr(j)=0.0D0
         ntagarr(j)=0.0D0
         tarr(j)=0.0D0
         Marr(j)=0.0D0
         Earr(j)=0.0D0
         farr(j)=0.0D0
         varrl(j)=0.0D0
         acarrl(j)=0.0D0         
         jerkl(j)=0.0D0
        ENDDO
       

         M0y=-Omoy*Tpy

        IF(M0y.GE.(2.0D0*pi))M0y=M0y-(2.0D0*pi)
      
        call keplersolve1(ecy,M0y,E0y)
        call keplersolve2(ecy,E0y,f0y)
       
      rl0=((apy*(1.0D0-ecy*ecy))/(1.0D0+ecy*cos(f0y)))*sin(f0y+ompery)
       

CCCCCCCCCCCC  initial alpha2, gamma calculation over the full observation   CCCCCC
       DO i=1,1001 
         tarr(i)=ta+dble(i-1)*tstep

         Marr(i)=Omoy*(tarr(i)-Tpy)

        IF(Marr(i).GE.(2.0D0*pi))Marr(i)=Marr(i)-(2.0D0*pi)
      
        call keplersolve1(ecy,Marr(i),Earr(i))
        call keplersolve2(ecy,Earr(i),farr(i))

       
       varrl(i)=fac1*(cos(farr(i)+ompery)+ecy*cos(ompery))

      acarrl(i)=fac2*sin(farr(i)+ompery)*
     &(1.0D0+ecy*cos(farr(i)))*(1.0D0+ecy*cos(farr(i))) 


       jerkl(i)=fac3*((1.0D0+ecy*cos(farr(i)))**3.0D0)*
     &(cos(farr(i)+ompery)+ecy*cos(ompery)-
     &3.0D0*ecy*sin(farr(i)+ompery)*sin(farr(i)))


      alpha2rr(i)=varrl(i)
      alpha1rr(i)=acarrl(i)/2.0D0
      alpha0rr(i)=jerkl(i)/6.0D0
      call gam3(alpha0rr(i),alpha1rr(i),alpha2rr(i),rl0,Tob,gammaarr(i))
      ntagarr(i)=dble(i)
      
      ENDDO

       
      call sort2(1001,gammaarr,ntagarr)

C---------------------


       nnn=ntagarr(1001)
   

        IF(nnn.EQ.1)THEN
         tn0=tarr(nnn)
         tnstep=(tarr(nnn+1)-tarr(nnn))/50.0D0
        ELSEIF(nnn.EQ.1001)THEN
         tn0=tarr(nnn-1)
         tnstep=(tarr(nnn)-tarr(nnn-1))/50.0D0
        ELSE
        tn0=tarr(nnn-1)
        tnstep=(tarr(nnn+1)-tarr(nnn-1))/100.0D0
       ENDIF

 1     CONTINUE
        
CCCCCC  setting array elements to zero to avoid accidental junk creeping in
        nnn=0
        DO jj=1,101
         alpha0rrn(jj)=0.0D0
         alpha1rrn(jj)=0.0D0
         alpha2rrn(jj)=0.0D0
         tarrn(jj)=0.0D0      
         gammaarrn(jj)=0.0D0
         ntagarrn(jj)=0.0D0
         Marrn(jj)=0.0D0
         Earrn(jj)=0.0D0
         farrn(jj)=0.0D0
         varrln(jj)=0.0D0
         acarrln(jj)=0.0D0        
         jerkln(jj)=0.0D0
        ENDDO

CCC   finer search to maximize gamma
        
        DO ii=1,101     
                   
          tarrn(ii)=tn0+dble(ii-1)*tnstep         
CCCCCCCC
         Marrn(ii)=Omoy*(tarrn(ii)-Tpy)

        IF(Marrn(ii).GE.(2.0D0*pi))Marrn(ii)=Marrn(ii)-(2.0D0*pi)

     
        call keplersolve1(ecy,Marrn(ii),Earrn(ii))
        call keplersolve2(ecy,Earrn(ii),farrn(ii))


       varrln(ii)=fac1*
     &(cos(farrn(ii)+ompery)+ecy*cos(ompery))


       acarrln(ii)=fac2*
     &sin(farrn(ii)+ompery)*
     &(1.0D0+ecy*cos(farrn(ii)))*(1.0D0+ecy*cos(farrn(ii)))
      

       jerkln(ii)=fac3*
     &((1.0D0+ecy*cos(farrn(ii)))**3.0D0)*
     &(cos(farrn(ii)+ompery)+ecy*cos(ompery)-
     &3.0D0*ecy*sin(farrn(ii)+ompery)*sin(farrn(ii)))

       alpha2rrn(ii)=varrln(ii)

       alpha1rrn(ii)=acarrln(ii)/2.0D0

       alpha0rrn(ii)=jerkln(ii)/6.0D0
     
CCCCCCCCC
       call gam3
     &(alpha0rrn(ii),alpha1rrn(ii),alpha2rrn(ii),rl0,Tob,gammaarrn(ii))
          ntagarrn(ii)=dble(ii)
        ENDDO


       call sort2(101,gammaarrn,ntagarrn)

         nnn=ntagarrn(101)


        IF(nnn.EQ.1)THEN
         tn0=tarrn(nnn)
         tnstep=(tarrn(nnn+1)-tarrn(nnn))/50.0D0
        ELSEIF(nnn.EQ.101)THEN
         tn0=tarrn(nnn-1)
         tnstep=(tarrn(nnn)-tarrn(nnn-1))/50.0D0
        ELSE
        tn0=tarrn(nnn-1)
        tnstep=(tarrn(nnn+1)-tarrn(nnn-1))/100.0D0
       ENDIF


         if(tnstep.gt.1.0D-5)then
           goto 1
         else
           goto 2
        endif
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  2     CONTINUE
         gammaf=gammaarrn(101)
   
CCCCCCCCCCCCCCCCCCCCCC

       RETURN
       END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine keplersolve2(eyy,Etyy,fyy)
       IMPLICIT NONE
       real*8 eyy,Etyy,fyy
       real*8 a1

       a1=SQRT((1.0D0+eyy)/(1.0D0-eyy))*tan(Etyy/2.0D0)

       fyy=2.0D0*atan(a1)

       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine keplersolve1(ECC,XMA,EA)
CCCC   modified Dunc's code to solve Kepeler's Equation  using Newton-Raphson method
      IMPLICIT NONE 
      real*8 EACC,RAD,EA,XMA,E,ECC
      integer NITS,I

      PARAMETER (EACC=1.0D-10,RAD=57.29577951308D0)
 
C
      EA=XMA + ECC*SIN(XMA)*(1.0D0+ECC*COS(XMA))
C
      DO 10 I=1,50
      E=(XMA + ECC*SIN(EA) - EA*ECC*COS(EA))/(1.0D0-ECC*COS(EA))
      IF(ABS(EA-E).LT.EACC) GO TO 20
   10 EA=E
C
   20 EA=E
      EA=(XMA + ECC*SIN(EA) - EA*ECC*COS(EA))/(1.0D0-ECC*COS(EA))
      NITS=I
 
       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE gam1(alpha2y,rl0y,t,ans)
       IMPLICIT NONE
       real*8 alpha2y,alpha2yy,rl0y,rl0yy
       real*8 t,low,ans
       real*8 f1,f2,f3
       real*8 func1,func2
       EXTERNAL func1,func2
       COMMON /alpha/alpha2yy,rl0yy

       low=0.0D0

C------ This part is needed because the arguments passed by can't be put in common block    
       alpha2yy=alpha2y
       rl0yy=rl0y
C-------------------------------------------------

       call qsimpmnj(func1,low,t,f1)
       call qsimpmnj(func2,low,t,f2)
       f3=SQRT(f1*f1+f2*f2)
       ans=f3/t

       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE gam2(alpha1y,alpha2y,rl0y,t,ans)
       IMPLICIT NONE
       real*8 alpha1y,alpha1yy
       real*8 alpha2y,alpha2yy,rl0y,rl0yy
       real*8 t,low,ans
       real*8 f1,f2,f3
       real*8 func1,func2
       EXTERNAL func1,func2
       COMMON /alpha/alpha1yy,alpha2yy,rl0yy

       low=0.0D0

C------ This part is needed because the arguments passed by can't be put in common block 
       alpha1yy=alpha1y   
       alpha2yy=alpha2y
       rl0yy=rl0y
C-------------------------------------------------

       call qsimpmnj(func1,low,t,f1)
       call qsimpmnj(func2,low,t,f2)
       f3=SQRT(f1*f1+f2*f2)
       ans=f3/t

       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE gam3(alpha0y,alpha1y,alpha2y,rl0y,t,ans)
       IMPLICIT NONE
       real*8 alpha1y,alpha1yy,alpha0y,alpha0yy
       real*8 alpha2y,alpha2yy,rl0y,rl0yy
       real*8 t,low,ans
       real*8 f1,f2,f3
       real*8 func1,func2
       EXTERNAL func1,func2
       COMMON /alpha/alpha0yy,alpha1yy,alpha2yy,rl0yy

       low=0.0D0

C------ This part is needed because the arguments passed by can't be put in common block 
       alpha0yy=alpha0y 
       alpha1yy=alpha1y   
       alpha2yy=alpha2y
       rl0yy=rl0y
C-------------------------------------------------

       call qsimpmnj(func1,low,t,f1)
       call qsimpmnj(func2,low,t,f2)
       f3=SQRT(f1*f1+f2*f2)
       ans=f3/t

       RETURN
       END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       REAL*8 function func1(z) 
       IMPLICIT NONE     
       REAL*8 c,pi       
       REAL*8 z
       REAL*8 g1,g2
       REAL*8 alpha2x,rl0x
       REAL*8 M,E,f
       REAL*8 Omp,ap,ec,omper,Omo,xxm,Tp
       COMMON /orbit/Omp,ap,ec,omper,Omo,xxm,Tp
       COMMON /alpha/alpha2x,rl0x

       c=2.99792458D+8
       pi=3.141592653589793D0 

       M=Omo*(z-Tp)

       IF(M.GE.(2.0D0*pi))M=M-(2.0D0*pi)
      
       call keplersolve1(ec,M,E)
       call keplersolve2(ec,E,f)
 
       g1=(((ap*(1.0D0-ec*ec))/(1.0D0+ec*cos(f)))*sin(f+omper))-rl0x-
     &alpha2x*z

       g2=(xxm*Omp*g1)/c       

       func1=cos(g2)
 

       RETURN
       END
         
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
       REAL*8 function func2(z)  
       IMPLICIT NONE     
       REAL*8 c,pi 
       REAL*8 z
       REAL*8 g1,g2       
       REAL*8 alpha2x,rl0x
       REAL*8 M,E,f
       REAL*8 Omp,ap,ec,omper,Omo,xxm,Tp
       COMMON /orbit/Omp,ap,ec,omper,Omo,xxm,Tp
       COMMON /alpha/alpha2x,rl0x
       
       c=2.99792458D+8
       pi=3.141592653589793D0 

       M=Omo*(z-Tp)

       IF(M.GE.(2.0D0*pi))M=M-(2.0D0*pi)
      
       call keplersolve1(ec,M,E)
       call keplersolve2(ec,E,f)
 
       g1=(((ap*(1.0D0-ec*ec))/(1.0D0+ec*cos(f)))*sin(f+omper))-rl0x-
     &alpha2x*z
       
       g2=(xxm*Omp*g1)/c

       func2=sin(g2)


       RETURN
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE qsimpmnj(func,a,b,s)
      IMPLICIT NONE
      INTEGER JMAX
      REAL*8 a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.D-6, JMAX=60)
CU    USES trapzd
      INTEGER j
      REAL*8 os,ost,st
      ost=-1.0D30
      os= -1.0D30
      DO j=1,JMAX        
        call trapzd(func,a,b,st,j)
        s=(4.0D0*st-ost)/3.0D0
        if(abs(s-os).lt.EPS*abs(os))goto 11
        os=s
        ost=st
       ENDDO
11    continue
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE trapzd(func,a,b,s,n)
       IMPLICIT NONE
       INTEGER n
       REAL*8 a,b,s,func
       EXTERNAL func
       INTEGER it,j
       REAL*8 del,sum,tnm,x       

       if (n.eq.1) then
         s=0.5D0*(b-a)*(func(a)+func(b))
       else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5D0*del
        sum=0.0D0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5D0*(s+(b-a)*sum/tnm)
       endif
       return
       END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      REAL*8 arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          temp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
