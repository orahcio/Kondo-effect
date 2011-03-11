      PROGRAM ORAHCIO
***********************
*     VERSAO 16/07/09
      
*#############PROGRAMA PRINCIPAL  #######################
*###########################################################
      
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      INTEGER DT1,DT11,DT2,DT22,DT37,DT38,DT54,DTUP,DTDN,OUT
      PARAMETER (DT1=1,DT11=11,DT2=2,DT37=37,DT38=38,DT54=54,OUT=3)
      PARAMETER (DTUP=41,DTDN=42,DT22=22)
      DIMENSION U(13),RF(8),EN(12),RC(13)
      DIMENSION IV(2),CONDV(1001,2),ADOS(10001,2),FP(2),DMI(2),DMF(2)
      DIMENSION VIC(2),DK(2),DD(2),AKK(2),ANUV(1001)
      DIMENSION AENERG(10001),BENERG(10001),ANORM(2)
      DIMENSION AMAX(2),ENG(2),BDOS(10001,2)
      DIMENSION ADEN(10001,2)
      DIMENSION SFLDOS(2),SXXF(2),SXXDS(2)
      EXTERNAL FMU
      COMMON /RESI/U,EN,RF,RC
      COMMON/DADOS/AEF,AEQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/RETAN/E00,EFF,E1,E2,DEL,DEL0,ZUI,EPSIL
      COMMON/ET/ETA,OI,OF,NT
      COMMON/BRETAN/GNT,GNF,GNQ,GDS,GNR,COMPLETEZA
      COMMON/OCC/TXF,TDS,TXFC,TXC,TOT,TSOMA,FRIED,R0F,R0C
      COMMON/CORREC/DEK
      COMMON/FANO/AP,AQ
      COMMON/FATOR/GAMA
      COMMON/FATOR_TOTAL/GAMAUP,GAMADN
      COMMON/BESSEL/AX
      COMMON/AMI/EMIN,KK
      COMMON/COMPLETEZA/ANXF,ANDS,ANSOMA
      COMMON/CONDUTANCIA/COND
      COMMON/DELTA/UDELTA
      COMMON/VOLTEXT/ANU
      COMMON/NNN/ANORMA
      COMMON/GLF/TXFLDOS
      COMMON/DOS_TOT/ADOS
      OPEN (UNIT=DT1,FILE='asaida1.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DT11,FILE='asaida2.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DT2,FILE='adosfull.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DT22,FILE='adenqd.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DT37,FILE='aconduti.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DT38,FILE='aocupacao_total.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DT54,FILE='adif.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DTUP,FILE='aaup.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DTDN,FILE='aadown.dat',STATUS='UNKNOWN')
      OPEN (UNIT=OUT,FILE='adados.out',STATUS='UNKNOWN')
      DATA PI /3.14159265358979323D0/
      EPS=1.D-12
********************
      CALL PARAMETROS
*********************
      AEF=-0.03D0
      T=1.D0/BETA
      VIC(1)=0.0
      VIC(2)=0.0
      DK(1)=-0.00290043D0
      DK(2)=-0.00031228D0
C      DMI(1)=-0.010D0
C      DMF(1)=0.01D0
C      DMI(2)=-0.00130D0
C      DMF(2)=0.00130D0
      DMI(1)=-0.012D0
      DMF(1)=0.012D0
      DMI(2)=-0.0028D0
      DMF(2)=0.0028D0
C      CALL EXATA
C     CALL OCUP3
C      DEK=Fmin(DMI,DMF,FMU,EPS)
C      AEQ=DEK
C      DIF=DABS(FRIED-R0F)
C      WRITE(6,36)DEK,TXF,TSOMA,R0F,FRIED,DIF
36      FORMAT('I=',I4,1X,'DEK=',F12.10,1X,
     *'S=',F6.4,1X,'R=',F10.4,1X,'F=',F10.4,1X,'D=',F12.6)
127   FORMAT(30F15.8)
      DMT=(OF-OI)/(NT-1)
      ALDOS0=0.D0
      AKFDN=1.D-6
      AKFDN=10.0D0
C      AKFDN=8.66065000D0
      AQ=0.01D0*DEXP(-AKFDN)
      DO IC=1,NT
      ADOS(IC,1)=0.D0
      ADOS(IC,2)=0.D0
      ENDDO
        API=0.50D0
        APF=0.20D0
        NP=1
        WRITE(6,4001)
4001    FORMAT('API=0.05',5X,'APF=-0.15',5X,'NP=2')
        READ(*,*)API,APF,NP
        DEF=0.D0
	IF(NP.EQ.1)GO TO 4000
        DP=(APF-API)/(NP-1)
4000    CONTINUE
      DO 5002 IB=1,NP
      AP=API+(IB-1)*DP
      P=AP
      DUP=D/(1.D0+P)
      DDN=D/(1.D0-P)
      RP=DSQRT((1.D0-P)/(1.D0+P))
      AKFUP=RP*AKFDN
      GAMA=0.5D0*PI*V2/D
      DD(1)=DUP
      DD(2)=DDN
      AKK(1)=AKFUP
      AKK(2)=AKFDN
      DFAT=(DUP-DDN)/DUP
        COND=0.D0
        GM0=0.D0
        AI=-0.001D0
        BI=0.0D0
        ADEKI=-0.0026D0
	ADEKF=0.0D0
        NDEK=1
        IF(NDEK.EQ.1)GO TO 5000
        DIDEK=(ADEKF-ADEKI)/(NDEK-1)
5000    CONTINUE
      AMAX(1)=0.D0
      AMAX(2)=0.D0
      OOI=-0.01D0
      OOF=0.01D0
      DT=(OOF-OOI)/(NT-1)
      ANORMA=0.D0
           DO  20 J=1,2
           D=DD(J)
           AX=AKK(J)
*****************************
C        DO 7002 II=1,NDEK
***********************
C      DEK=ADEKI+(II-1)*DIDEK
***********************************
         DEK=Fmin(DMI(J),DMF(J),FMU,EPS)
         VIC(J)=DEK
         AEQ=AMU+DEK
         CALL NORMALIZA(ANOR,SNF)
*****************************
           ANORMA=ANORMA+ANOR
           ANORM(J)=ANOR
           SFLDOS(J)=SNF
         CALL EXATA
         CALL OCUP3
*******************
           SXXF(J)=TXF
           SXXDS(J)=TDS
        DIF=DABS(FRIED-R0F)
         WRITE(6,36)II,DEK,TSOMA,R0F,FRIED,DIF
7002     CONTINUE
      TXF2=2.D0*TXF
      TXC2=2.D0*TXC
      IF(J.EQ.1)WRITE(DTUP,127)P,D,AX,DEK,DFAT,AEF,T,AMU,TSOMA,
     *TXC2,TXF2,TOT,R0F,R0C,FRIED,DIF
      IF(J.EQ.2)WRITE(DTDN,127)P,D,AX,DEK,DFAT,AEF,T,AMU,TSOMA,
     *TXC2,TXF2,TOT,R0F,R0C,FRIED,DIF
             DO 10 I=1,NT
             EW=OI+(I-1)*DMT
             AENERG(I)=EW
             CALL EXATA
             CALL GKONDO(EW,ZGF,ZGC,ZGFC,GCOND)
********************************
             ADEN(I,J)=-CC*DIMAG(ZGF)
             CALL LDOS(EW,ALDOS)
********************************
             ADOS(I,J)=ALDOS
10      CONTINUE
      DO K=1,NT
      EW=OOI+(K-1)*DT
      BENERG(K)=EW
             CALL EXATA
             CALL GKONDO(EW,ZGF,ZGC,ZGFC,GCOND)
             CALL LDOS(EW,BLDOS)
********************************
             BDOS(K,J)=BLDOS
      IF (AMAX(J).LT.DABS(BDOS(K,J))) AMAX(J)=BDOS(K,J)
      WRITE(DT22,128)EW,ADEN(K,1),ADEN(K,2),ADEN(K,1)+ADEN(K,2)
      ENDDO
      DO IC=1,NT
      IF (AMAX(J).EQ.DABS(BDOS(IC,J)))ENG(J)=BENERG(IC)
      ENDDO
      DE=ENG(1)-ENG(2)
20    CONTINUE
          DO KL=1,NT
          CALL EXATA
          CALL GKONDO(EW,ZGF,ZGC,ZGFC,GCOND)
           CALL LDOS(EW,BLDOS)
      WRITE(DT2,128)AENERG(KL),ADOS(KL,1)/ANORMA,ADOS(KL,2)/ANORMA,
     *(ADOS(KL,1)+ADOS(KL,2))/ANORMA
          ENDDO
          DO KL=1,NT
       WRITE(DT2,128)BENERG(KL),BDOS(KL,1)/ANORMA,BDOS(KL,2)/ANORMA,
     *(BDOS(KL,1)+BDOS(KL,2))/ANORMA
          ENDDO
          SUP=SXXF(1)+SXXDS(2)
          SDOWN=SXXF(2)+SXXDS(1)
      WRITE(DT54,128)AQ,P,DE,ENG(1),AMAX(1)/ANORMA,ENG(2),
     *AMAX(2)/ANORMA,VIC(1),VIC(2),SFLDOS(1)/ANORMA,SFLDOS(2)/ANORMA,
     *SXXF(1),SXXDS(1),SXXF(2),SXXDS(2),SUP,SDOWN,SUP+SDOWN,
     *AKFDN,AQ,V,AEF
C      WRITE(6,36)J,DEK,TSOMA,R0F,FRIED,DE
          ANUI=-0.010D0
	  ANUF=0.010D0
          NU=1
          IF(NU.EQ.1)GO TO 5001
          DINU=(ANUF-ANUI)/(NU-1)
5001      CONTINUE
         DO 3002 IM=1,NU
***********************
          ANU=ANUI+(IM-1)*DINU
          ANUV(IM)=ANU
         CALL CONDUT(CND)
*************************
         CONDV(IM,J)=CND
3002     CONTINUE
      DO 77 IL=1,NU
      RFILT=0.D0
      IF((CONDV(IL,1)+CONDV(IL,2)).LT.1.D-10)GOTO 75
      RFILT=(CONDV(IL,1)-CONDV(IL,2))/(CONDV(IL,1)+CONDV(IL,2))
75    CONTINUE
      COND=CONDV(IL,1)+CONDV(IL,2)
      WRITE(DT37,128)ANUV(IL),RFILT,CONDV(IL,1),CONDV(IL,2),COND
77    CONTINUE
C7002     CONTINUE
C5002  CONTINUE
C      DO K=1,NEK
5002  CONTINUE
128     FORMAT(30F20.8)
      CALL DENSI
****************
      STOP
      END

      function bissec(ax,bx,f,tol)
      real*8 ax, bx, tol
      real*8 x, xx, xa, xb, f0, f1, f2, f
      external f
      write(6,*) tol, ax, bx
      f0=f(ax)
      f2=f(bx)
      xa=ax
      xb=bx
      x=100.D0
      xx=1.D0
      do
      xx=x
      x=(xb+xa)/2.0D0
      f1=f(x)
      write(6,*) ':', xb, x, xa, f0, f1, f2
      if(f0*f1>0) then ! descarta o primeiro intervalo
        xa=x
        f0=f1
        if(f1*f2>0) then ! descarta o segundo intervalo
           write(6,*) 'Erro: raiz fora do intervalo'
           exit
        end if
      else
        xb=x
        f2=f1
      end if

      if(abs(x-xx)<tol) exit
      end do

      bissec=x
      return
      end


        SUBROUTINE BISSECAO(A,B,R,TOT,ER)

*****************************************

C       RAIZES DE UMA EQ. PELO METODO DA BISSECAO

        IMPLICIT COMPLEX*16 (X-Z)

        IMPLICIT REAL*8 (A-H,O-W)

        INTEGER DT50

        PARAMETER (DT50=50)

        COMMON/TOTAL/SNT

        OPEN (UNIT=DT50,FILE='raiz.dat',STATUS='UNKNOWN')

        DATA N /10000/

        EPS=1.E-10

        DO 10 I=1,N

        AMED=0.5D0*(A+B)

        FA=FUNCTMU(AMED)

        F1=FUNCTMU(A)

        AP1=FA*F1

        IF(AP1.GT.0)GO TO 20

        B=AMED

        GO TO 30

20      CONTINUE

        A=AMED

30      CONTINUE

        ADAB=ABS(A-B)

        IF(ADAB.LT.EPS)GO TO 40

        WRITE(DT50,50)I,A,B,AMED,ADAB

10      CONTINUE

40      CONTINUE

        R=AMED

        ER=ADAB

        TOT=SNT-F1

C        WRITE(6,50)I,A,B,AMED,ER

        WRITE(DT50,50)I,A,B,AMED,ADAB,TOT

50      FORMAT('I=',I5,5X,' A=',F10.5,5X,' B=',F15.8,5X

     *,' AMED=',F10.8,' ER=',F15.8,' TOT=',F15.8)

        RETURN

        END


      FUNCTION FUNCTMU(DEQ)
************************************
      IMPLICIT REAL*8 (A-H,O-W)
      IMPLICIT COMPLEX*16 (X-Z)
      REAL*8 MU(9),MUC(9),UU(9),UE(9)
      REAL*8 AMU
      INTEGER DT76
      PARAMETER (DT76=76)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/OCC/TXF,TDS,TXFC,TXC,TOT,TSOMA,FRIED,R0F,R0C
      OPEN (UNIT=DT76,FILE='afunc.dat',STATUS='UNKNOWN')
      EQ=AMU+DEQ
      CALL EXATA
      CALL OCUP3
********************
C      FUNCTMU=DABS(TSOMA-1.D0)
      FUNCTMU=FRIED-R0F
C      FUNCTMU=1.D0-TSOMA
C      WRITE(6,*)DEQ,FRIED,R0F
      RETURN
      END

          SUBROUTINE CONDUT(CND)
***************************
          IMPLICIT REAL*8 (A-H,O-W)
          IMPLICIT COMPLEX*16 (X-Z)
          REAL*8 MU
          EXTERNAL ACOND
          INTEGER DT50
          PARAMETER (DT50=50)
          COMMON/DADOS/EF,EQ,MU,V,V2,BETA,DS,CC,D
          COMMON/RETAN/E00,EFF,E1,E2,DEL,DEL0,ZUI,EPSIL
          COMMON/CONDUTANCIA/COND
          COMMON/VOLTEXT/ANU
          OPEN (UNIT=DT50,FILE='acondutancia.dat',STATUS='UNKNOWN')
          T=1.D0/BETA
          A=-30*T
          B=30*T
          RES=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,ACOND)
          CND=BETA*RES
          WRITE(DT50,*)ANU,CND,RES
          RETURN
          END
      FUNCTION ACOND(EW)
*************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      INTEGER DT78
      PARAMETER (DT78=78)
      DIMENSION U(13),RF(8),EN(12),RC(13)
      COMMON /RESI/U,EN,RF,RC
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/NNN/RENOR
      COMMON/FANO/AP,AQ
      COMMON/FATOR/GAMA
      COMMON/FATOR_TOTAL/GAMAUP,GAMADN
      COMMON/ET/ETA,OI,OF,NT
      COMMON/BESSEL/AX
      COMMON/VOLTEXT/ANU
      OPEN (UNIT=DT78,FILE='acondos.dat',STATUS='UNKNOWN')
      DATA PI /3.14159265358979323D0/
      EW2=EW*EW
      DELTA0=0.5D0*PI*V2/D
      DELTA=0.5D0*PI*V2/(DELTA0*DATAN(D/DELTA0))
      DELTA2=DELTA*DELTA
      RCFREE=0.5D0*DELTA/(ATAN(D/DELTA)*(DELTA2+EW2))
      ZW=DCMPLX(EW,ETA)
      ZW2=ZW*ZW
      ARG_BESSEL=AX*(1.D0+EW/D)
      ZHFAT=DELTA/(ZW-DCMPLX(0.D0,DELTA))
      CALL JY01A(ARG_BESSEL,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
************************************************************
      AJ0=BJ0
      AJ02=AJ0*AJ0
      ZZ=AX*DCMPLX(1.0D0,DELTA/D)
      CALL CJY01(ZZ,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
************************************************************
          ZHANK=DCMPLX(CBJ0,-CBY0)
          CALL EXATA
          CALL GKONDO(ZW,ZGF,ZGC,ZGFC,GCOND)
********************************************
          R0F=-CC*DIMAG(ZGF)
          R0C=CC*DIMAG(ZGC)
          ASIGMA=-1.D0/(RCFREE*AJ0)
          DELES=DATAN(DIMAG(ZGF)/DREAL(ZGF))
          AQ1=AQ
          AQ2=DREAL(ZHANK*ZHFAT)
          AQ3=-RCFREE*BY0
          AQTIL=AQ1
          DELQES=DATAN(ASIGMA*AQTIL)
          COS1=COS(DELES-DELQES)*COS(DELES-DELQES)
          COS2=COS(DELQES)*COS(DELQES)
****************************************
      ALDOS=-RCFREE*(1.D0-AJ02*(1.D0-COS1/COS2))/RENOR
      FE=FERM(EW)
      GCOND=FE*(1.D0-FE)*ALDOS
      WRITE(DT78,*)EW,DELES
      ACOND=GCOND
      RETURN
      END
      SUBROUTINE LDOS(EW,ALDOS)
**************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      INTEGER DT78
      PARAMETER (DT78=78)
      DIMENSION U(13),RF(8),EN(12),RC(13)
      COMMON /RESI/U,EN,RF,RC
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/FANO/AP,AQ
      COMMON/ET/ETA,OI,OF,NT
      COMMON/BESSEL/AX
      COMMON/NNN/ANOR
      OPEN (UNIT=DT78,FILE='aldos.dat',STATUS='UNKNOWN')
      DATA PI /3.14159265358979323D0/
      EW2=EW*EW
      DELTA0=0.5D0*PI*V2/D
      DELTA=0.5D0*PI*V2/(DELTA0*DATAN(D/DELTA0))
      DELTA2=DELTA*DELTA
C      RCFREE=0.5D0*DELTA/(ATAN(D/DELTA)*(DELTA2+EW2))
      RCFREE=0.5D0/D
      ZW=DCMPLX(EW,ETA)
      ZW2=ZW*ZW
      ZARG=(ZW-D+AMU)/(ZW+D+AMU)
      ZM=0.5D0*CDLOG(ZARG)/D
      RCFREE=DIMAG(ZM)/PI
      ZW=DCMPLX(EW,ETA)
      ZW2=ZW*ZW
C      ARG_BESSEL=AX*(1.D0+EW/D)
      ARG_BESSEL=AX
      ZHFAT=DELTA/(ZW-DCMPLX(0.D0,DELTA))
      CALL JY01A(ARG_BESSEL,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
************************************************************
      AJ0=BJ0
      AJ02=AJ0*AJ0
      ZZ=AX*DCMPLX(1.0D0,DELTA/D)
      CALL CJY01(ZZ,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
************************************************************
          ZHANK=DCMPLX(CBJ0,-CBY0)
          CALL EXATA
          CALL GKONDO(ZW,ZGF,ZGC,ZGFC,GCOND)
********************************************
          R0F=-CC*DIMAG(ZGF)
          R0C=CC*DIMAG(ZGC)
          ASIGMA=2.D0*RCFREE*AJ0
          DELES=DATAN(DIMAG(ZGF)/DREAL(ZGF))
          AQ1=AQ
          AQ2=DREAL(ZHANK*ZHFAT)
          AQ3=-RCFREE*BY0
C          AQTIL=AQ1+AQ2+AQ3
          AQTIL=AQ1
          DELQES=DATAN(-AQTIL/ASIGMA)
          COS1=COS(DELES-DELQES)*COS(DELES-DELQES)
          COS2=COS(DELQES)*COS(DELQES)
****************************************
      ALDOS=RCFREE*(1.D0-AJ02*(1.D0-COS1/COS2))
C      IF(EW.EQ.0.D0)WRITE(DT78,10)EW,ALDOS/ANOR
 10   FORMAT(10F25.8)
      RETURN
      END
      SUBROUTINE NORMALIZA(ANOR,SNF)
************************************
      IMPLICIT REAL*8 (A-H,O-W)
      IMPLICIT COMPLEX*16 (X-Z)
      REAL*8 MU
      EXTERNAL ANORM
      EXTERNAL OCUP_LDOS
      DIMENSION VA(5),VB(5)
      COMMON/DADOS/EF,EQ,MU,V,V2,BETA,DS,CC,D
      EPSIL=1.D-4
      SN=0.D0
      SF=0.D0
      DDD=0.005D0
      DD=0.01D0
C      VA(1)=-0.5D0
      VA(1)=-D-0.5D0
      VB(1)=EF-DDD
      VA(2)=EF-DDD
      VB(2)=EF+DDD
      VA(3)=EF+DDD
      VB(3)=MU-DD
      VA(4)=MU-DD
      VB(4)=MU+DD
      VA(5)=MU+DD
      VB(5)=D+0.5D0
C      VB(5)=0.5D0
      DO 10 I=1,5
         A=VA(I)
         B=VB(I)
         RES1=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,ANORM)
         RES2=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,OCUP_LDOS)
         SN=SN+RES1
         SF=SF+RES2
C         WRITE(6,*)'XXXXXXXXXXXX=',A,B,RES2
10        CONTINUE
          ANOR=SN
          SNF=SF
          RETURN
          END
      FUNCTION ANORM(EW)
*****************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/FANO/AP,AQ
      COMMON/ET/ETA,OI,OF,NT
      COMMON/BESSEL/AX
      CALL LDOS(EW,ALDOS)
*************************
      ANORM=ALDOS
      RETURN
      END
      FUNCTION OCUP_LDOS(EW)
*****************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/FANO/AP,AQ
      COMMON/ET/ETA,OI,OF,NT
      COMMON/BESSEL/AX
      CALL LDOS(EW,ALDOS)
*************************
      OCUP_LDOS=ALDOS*FERM(EW)
C       OCUP_LDOS=ALDOS
      RETURN
      END
      SUBROUTINE CONDUTIVIDADE(COND)
********************************
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        INTEGER DT56
        PARAMETER (DT56=56)
        COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
        COMMON/FATOR/GAMA
        COMMON/FANO/AP,AQ
        OPEN (UNIT=DT56,FILE='integrando.dat',STATUS='UNKNOWN')
        DATA PI /3.14159265358979323D0/
      R0=0.5D0*D
      ZW=DCMPLX(AMU,ETA)
      CALL GKONDO(ZW,ZGF,ZGC,ZGFC,GCOND)
****************************************
      GAM_TIP=2.D0*V2*DIMAG(ZGC)
      R0=0.5D0*D
      REG=DREAL(ZGF)
      AIMG=DIMAG(ZGF)
      AQ2=AQ*AQ
      ALDOS=R0*(1.D0+((1.D0-AQ2)*GAMA*AIMG+2.D0*GAMA*AQ*REG))
      ALDOS=ALDOS/(1.D0+AQ2)
      COND=GAM_TIP*ALDOS
        WRITE(DT56,*)EW,COND
        ACOND=GCOND
        RETURN
        END
      SUBROUTINE EXATA
**********************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      DIMENSION U(13),RF(8),EN(12),RC(13),DXP(12)
      COMMON /RESI/U,EN,RF,RC
      COMMON/DADOS/AEF,AEQ,AMU,VV,VV2,BETA,DS,CC,D
      COMMON/RETAN/E00,EFF,E1,E2,DEL,DEL0,ZUI,EPSIL
      COMMON/ET/ETA,OI,OF,NT
      COMMON/AMI/EMIN,KK
      DATA C1 /2.8284271247D0/
      DATA PI /3.14159265358979323D0/
      AVG4=PI*VV/(2.D0*D)
      V=VV*AVG4
      V2=V*V
      EF=AEF-AMU
      EQ=AEQ-AMU
      DELTA=DSQRT((EQ-EF)*(EQ-EF)+4.0D0*V2)
      DELTALI=DSQRT((EQ-EF)*(EQ-EF)+8.0D0*V2)
      TANGFI=2.D0*V/(EQ-EF+DELTA)
      DEN1=DSQRT(1.D0+TANGFI*TANGFI)
      SENFI=TANGFI/DEN1
      SENFI2=SENFI*SENFI
      COSFI=1.D0/DEN1
      COSFI2=COSFI*COSFI
      TANGLAM=C1*V/(EQ-EF+DELTALI)
      DEN2=DSQRT(1.D0+TANGLAM*TANGLAM)
      SENLAM=TANGLAM/DEN2
      SENLAM2=SENLAM*SENLAM
      COSLAM=1.D0/DEN2
      COSLAM2=COSLAM*COSLAM

C................AUTOVALORES DE ENERGIA..........

      EN(1)=0.D0
      EN(2)=0.5D0*(EF+EQ-DELTA)
      EN(3)=EN(2)
      EN(4)=0.5D0*(EF+EQ+DELTA)
      EN(5)=EN(4)
      EN(6)=EF+EQ
      EN(7)=EN(6)
      EN(8)=EN(6)
      EN(9)=0.5D0*(EF+3.0D0*EQ-DELTALI)
      EN(10)=0.5D0*(EF+3.0D0*EQ+DELTALI)
      EN(11)=EF+2.0D0*EQ
      EN(12)=EN(11)
 70   FORMAT(10F15.8)

      EMIN=EN(1)
      DO 50 I=2,12
         IF(EMIN.GT.EN(I)) EMIN=EN(I)
50      CONTINUE
        KK=0
        DO 51 I=1,12
           KK=KK+1
           IF(EMIN.EQ.EN(I))GO TO 52
51      CONTINUE
52      CONTINUE
C..........FUNCAO DE PARTICAO............
        
        AZ=0.D0
        DO 100 I=1,12
           EN(I)=EN(I)-EMIN
           AZ=AZ+DEXP(-BETA*EN(I))
100   CONTINUE
C........ DIFERENCAS DE ENERGIAS............
      
      U(1)= 0.5D0*(EQ+EF-DELTA)
      U(2)= 0.5D0*(EQ+EF+DELTA)
      U(3)= 0.5D0*(EQ+EF-DELTALI)
      U(4)= 0.5D0*(EQ+EF+DELTALI)
      U(5)= EQ-0.5D0*(DELTALI-DELTA)
      U(6)= EQ+0.5D0*(DELTALI+DELTA)
      U(7)= EQ-0.5D0*(DELTALI+DELTA)
      U(8)= EQ+0.5D0*(DELTALI-DELTA)
      U(9)= EQ
      U(10)= 0.5D0*(3*EQ+EF-DELTALI)
      U(11)= 0.5D0*(3*EQ+EF+DELTALI)
      U(12)= 0.5D0*(3*EQ+EF+DELTA)
      U(13)= 0.5D0*(3*EQ+EF-DELTA)
      
C............EXPONENCIAIS .........................
      DO 61 I=1,12
         DXP(I) = DEXP(-BETA*EN(I))
61     CONTINUE
       
       
       
C...............RESIDUOS .....................
       
       RF(1)=COSFI2*(DXP(1)+DXP(2)+1.5D0*(DXP(7)+DXP(5)))/AZ
       RF(2)=SENFI2*(DXP(1)+DXP(4)+1.5D0*(DXP(3)+DXP(7)))/AZ
       RF(3)=COSLAM2*(DXP(10)+DXP(11))/AZ
       RF(4)=SENLAM2*(DXP(9)+DXP(11))/AZ
       RF(5)=0.5D0*SENFI2*COSLAM2*(DXP(9)+DXP(2))/AZ
       RF(6)=0.5D0*SENFI2*SENLAM2*(DXP(2)+DXP(10))/AZ
       RF(7)=0.5D0*COSFI2*COSLAM2*(DXP(4)+DXP(9))/AZ
       RF(8)=0.5D0*COSFI2*SENLAM2*(DXP(4)+DXP(10))/AZ
       
       
       
       RC(1)=SENFI2*(DXP(1)+DXP(3)+1.5D0*(DXP(4)+DXP(8)))/AZ
       RC(2)=COSFI2*(DXP(1)+DXP(5)+1.5D0*(DXP(2)+DXP(8)))/AZ
       RC(3)=0.5D0*SENLAM2*(DXP(2)+DXP(10))/AZ
       RC(4)=0.5D0*COSLAM2*(DXP(9)+DXP(12))/AZ
       RC(5)=0.5D0*(0.5D0*C1*SENFI*SENLAM-COSFI*COSLAM)*
     *      (0.5D0*C1*SENFI*SENLAM-COSFI*COSLAM)*(DXP(2)+DXP(9))/AZ
       RC(6)=0.5D0*(0.5D0*C1*SENFI*COSLAM+COSFI*SENLAM)*
     *      (0.5D0*C1*SENFI*COSLAM+COSFI*SENLAM)*(DXP(2)+DXP(10))/AZ
       RC(7)=0.5D0*(0.5D0*C1*COSFI*SENLAM+SENFI*COSLAM)*
     *      (0.5D0*C1*COSFI*SENLAM+SENFI*COSLAM)*(DXP(4)+DXP(9))/AZ
       RC(8)=0.5D0*(0.5D0*C1*COSFI*COSLAM-SENFI*SENLAM)*
     *      (0.5D0*C1*COSFI*COSLAM-SENFI*SENLAM)*(DXP(4)+DXP(10))/AZ
       RC(9)=1.5D0*(DXP(6)+DXP(12))/AZ
       RC(10)=SENLAM2*(DXP(1)+DXP(9))/AZ
       RC(11)=COSLAM2*(DXP(1)+DXP(10))/AZ
       RC(12)=2.D0*COSFI2*(DXP(2)+DXP(12))/AZ
       RC(13)=2.D0*SENFI2*(DXP(4)+DXP(12))/AZ
C     ENERGIA LIVRE DE HELMHOLTZ
       RETURN
       END

***** FIM DA SUB EXATA ********************
      SUBROUTINE GKONDO(ZW,ZGF,ZGC,ZGFC,GCOND)
********************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      REAL*8 MU(8),MUC(13),UU(13),UE(12)
      COMMON/RESI/UU,UE,MU,MUC
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      DATA PI /3.14159265358979323D0/
      ZW2=ZW*ZW
      GAMMA=2.D0*D/PI
      DELTA=0.5D0*PI*V2/D
      DELTA2=DELTA*DELTA
      GAMMA2=GAMMA*GAMMA
      D2=D*D
      ZW2=ZW*ZW
      ETTA=0.0031D0
      ETTA=1.D-8
      AVG4=PI*V/(2.D0*D)
      VG2=V*V*AVG4*AVG4
      BW=DREAL(ZW)
      ZWD=DCMPLX(BW,ETTA)
      ZGFXAT=MU(1)/(ZWD-UU(1))+MU(2)/(ZWD-UU(2))+MU(3)/(ZWD-UU(3))+
     *MU(4)/(ZWD-UU(4))+MU(5)/(ZWD-UU(5))+MU(6)/(ZWD-UU(6))+
     *MU(7)/(ZWD-UU(7))+MU(8)/(ZWD-UU(8))
      ZGFXAT=-ZGFXAT
      ZATO=(ZWD-EQ+AMU)*(ZGFXAT)/((ZWD-EQ+AMU)-VG2*(ZGFXAT))
      ZARG=(ZW-D+AMU)/(ZW+D+AMU)
      ZM=0.5D0*V2*CDLOG(ZARG)/D
      Z1NUM=ZATO
      Z1DEN=1.D0-ZATO*ZM
      ZGF=-Z1NUM/Z1DEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ZGC=ZM/V2
      ZGC=ZM/V2+V2*ZGF/(D2-ZW2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ZGFC=ZM*ZATO/(V*Z1DEN)
      R0F=-CC*DIMAG(ZGF)
      R0C=-CC*DIMAG(ZGC)
      R0FC=CC*DIMAG(ZGFC)
C     CALCULO DA CONDUTANCIA
      Z1=(ZW-D+AMU)/(ZW+D+AMU)
      ZG0C=0.5D0*CDLOG(Z1)/D
      Z2=ZG0C
      Z3=1.D0-Z2*V2*ZGF
C     ZG00=Z2/Z3
      ZG00=ZG0C*(1.D0-ZG0C*V2*ZGF)
C     ZG00=ZG0C
      RZ1=DREAL(ZG00)
      RZ12=RZ1*RZ1
      RZ2=DIMAG(ZG00)
      RZ22=RZ2*RZ2
      AS=RZ12+RZ22
C     AS=ZG00*DCONJG(ZG00)
      FE=FERM(BW)
      GCOND=GAMMA2*BETA*FE*(1.D0-FE)*AS
      RETURN
      END

      SUBROUTINE PARAMETROS
****************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      REAL*8 MU
      COMMON/DADOS/EF,EQ,MU,V,V2,BETA,DS,CC,D
      COMMON/RETAN/E00,EFF,E1,E2,DEL,DEL0,ZUI,EPSIL
      COMMON/ET/ETA,OI,OF,NT
      COMMON/DELTA/DELTA
      DATA PI /3.14159265358979323D0/
      CC=1.D0/PI
      ZUI=(0.D0,1.D0)
      OI=-0.20D0
      OF=0.10D0
      W=PI
      D=1.D0
      DELTA=0.01D0
      NT=5001
C     WRITE(6,432)
C     432     FORMAT(' OI=-0.15D0',5X,' OF=0.15D0',5X,' D=1.D0',
C     *5X, ' NT=1001')
C     READ(*,*)OI,OF,D,NT
      EQ=0.D0
      MU=0.D0
      V=DSQRT(2.D0*D*DELTA/PI)
      V=0.12D0
      EF=-0.1D0
      T=1.D-5
      ETA=1.D-3
C     WRITE(6,433)
C     433     FORMAT(' EF=-5.D0',5X,' EQ=0.D0',/,
C     *' MU=0.D0',5X,' T=1.D-4',5X,' ETA=1.D-4')
C     READ(*,*)EF,EQ,MU,T,ETA
      BETA=1.D0/T
      V2=V*V
      E00=-10.D0
      EFF=10.D0
      DEL=5.D0
      E1=-0.1
      E2=0.1
      DEL0=PI/(1.04D0*BETA)
      EPSIL=1.D-8
C     WRITE(6,434)
C     434     FORMAT(' E00=-10.D0',5X,' EFF=10.D0',5X,' E1=-0.1D0',5X,' E2
C     *0.1D0',/,' DEL=5.D0',5X,' EPSIL=1.D-8')
C     READ(*,*)E00,EFF,E1,E2,DEL,EPSIL
C     WRITE(3,14)
C     WRITE(6,14)
C     WRITE(3,16)ETA,E00,EFF,EPSIL,E1,E2,DEL,DEL0
C     WRITE(6,16)ETA,E00,EFF,EPSIL,E1,E2,DEL,DEL0
C     WRITE(3,15)OI,OF,NT,D,EQ,ENA
C     WRITE(6,15)OI,OF,NT,D,EQ,ENA
C     WRITE(3,25)MU,EF,V,T
C     WRITE(6,25)MU,EF,V,T
C     16      FORMAT(' ETA=',F15.8,5X,' E00=',F15.8,/,
C     *' EFF=',F15.8,5X,' EPSIL=',F15.8,5X,' E1=',F15.8,/,
C     *' E2=',F15.8,5X,' DEL=',F15.8,5X,' DEL0=',F15.8,/)
C     14      FORMAT(' OS PARAMETROS DO PROBLEMA SAO:')
C     15      FORMAT(' OI=',F10.5,5X,' OF=',F10.5,' NT=',I5,5X,/,
C     *' D=',F10.5,5X,' EQ=',F10.5,5X,' ENA=',F10.5)
C     25      FORMAT(' MU=',F10.5,5X,5X,' EF=',F10.5,/,
C     *' V=',F10.5,5X,' T=',E10.5)
      RETURN
      END
      SUBROUTINE OCUP3
********************************
      IMPLICIT REAL*8 (A-H,O-W)
      IMPLICIT COMPLEX*16 (X-Z)
      REAL*8 MU
      INTEGER DT30,DT31
      PARAMETER (DT30=30,DT31=31)
      EXTERNAL GFF
      EXTERNAL GCC
      EXTERNAL GIDS
      EXTERNAL GFFC
      EXTERNAL GLDOS
      DIMENSION VA(7),VB(7)
      COMMON/DADOS/EF,EQ,MU,V,V2,BETA,DS,CC,D
      COMMON/RETAN/E00,EFF,E1,E2,DEL,DEL0,ZUI,EPSIL
      COMMON/ET/ETA,OI,OF,NT
      COMMON/OCC/TXF,TDS,TXFC,TXC,TOT,TSOMA,FRIED,R0F,R0C
      COMMON/CORREC/DEK
      COMMON/GLF/TXFLDOS
      OPEN (UNIT=DT30,FILE='ocupacao_1.dat',STATUS='UNKNOWN')
      OPEN (UNIT=DT31,FILE='aint_ocup.dat',STATUS='UNKNOWN')
      DATA PI /3.14159265358979323D0/
      SC=0.D0
      SF=0.D0
      SFC=0.D0
      SDS=0.D0
      GDOS=0.D0
      DDD=0.005D0
      DD=0.001D0
      VA(1)=-D-0.1D0
      VB(1)=EF-DDD
      VA(2)=EF-DDD
      VB(2)=EF+DDD
      VA(3)=EF+DDD
      VB(3)=MU-DD
      VA(4)=MU-DD
      VB(4)=MU+DD
      VA(5)=MU+DD
      VB(5)=D+0.1D0
      DO 10 I=1,5
         A=VA(I)
         B=VB(I)
         RES1=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,GCC)
         RES2=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,GFF)
         RES3=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,GIDS)
         RES4=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,GFFC)
C         RES5=QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,GLDOS)
         SC=SC+RES1
         SF=SF+RES2
         SDS=SDS+RES3
         SFC=SFC+RES4
         GDOS=GDOS+RES5
C         WRITE(6,*)'XXXXXX=',A,B,RES2
10        CONTINUE
	  ZZ=DCMPLX(0.0D0,ETA)
          CALL GKONDO(ZZ,ZGF,ZGC,ZGFC,GCOND)
          R0F=-CC*DIMAG(ZGF)
          R0C=CC*DIMAG(ZGC)
          TXC=CC*SC
          TXF=CC*SF
          TXFLDOS=GDOS
          TXFC=CC*SFC
          TDS=CC*SDS
          TSOMA=TXF+TDS
          TR2=TDS-TXF
          TOT=2.D0*(TXF+TXC)
          DELL=0.5D0*PI*V2/D
          DELLI=1.D0/(PI*DELL)
C          DELOR=1.D0/(PI*PI*V2*R0C)
C          DELTA=0.5D0*PI*V2/D
C          DELTA=0.5D0*PI*V2/(DELTA*DATAN(D/DELTA))
C          DELOR=1.D0/(PI*DELTA)
          FRIED=DELLI*SIN(PI*TXF)*SIN(PI*TXF)
C          FRIED=DELOR*SIN(3.1415*TXF)*SIN(3.1415*TXF)
          WRITE(DT30,31)T,MU,DEK,TSOMA,TR2,TXC,TXF,TDS,TXFC,TOT,
     *         R0F,R0C,FRIED
C     WRITE(6,30)T,DEK
C     MU,TXC,TXF,TDS,TXFC,TSOMA,TOT,DEK,R0F,R0C,FRIED
 30       FORMAT('T=',F15.8,5X,'DEK=',F15.8)
C     ,5X,*'TXC=',F15.8,/,'TXF=',F15.8,
C     *'TDS=',F15.8,5X,' NFC=',F15.8,/,'TSOMA=',F15.8,5X,'NTOT=',
C     *F15.8,5X,'DEK=',F15.8,
C     */,'R0F=',F15.8,5X,'R0C=',F15.8,5X,'FRIEDEL=',F15.8)
 31       FORMAT(20F15.8)
          RETURN
          END
      FUNCTION GLDOS(EW)
*****************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/NNN/ANORMA
      CALL LDOS(EW,ALDOS)
*************************
      GLDOS=FERM(EW)*ALDOS/ANORMA
      RETURN
      END

      FUNCTION GFF(EW)
*****************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      TT=1.D0/BETA
      ETTA=3.d0*TT
C     ETTA=1.D-5
      ZW=DCMPLX(EW,ETTA)
      CALL GKONDO(ZW,ZGF_IMP,ZGC_IMP,ZGFC_IMP,GCOND)
************************************************
      ZFE=ZFER(ZW)
      GFF=-DIMAG(ZFE*ZGF_IMP)
      RETURN
      END
      FUNCTION GCC(EW)
*****************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      ETTA=0.05D0
      ZW=DCMPLX(EW,ETTA)
      CALL GKONDO(ZW,ZGF_IMP,ZGC_IMP,ZGFC_IMP,GCOND)
************************************************
      ZFE=ZFER(ZW)
      GCC=DIMAG(ZFE*ZGC_IMP)
      RETURN
      END
      FUNCTION GIDS(EW)
*****************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      ETTA=0.05D0
      ZW=DCMPLX(EW,ETTA)
      CALL GKONDO(ZW,ZGF_IMP,ZGC_IMP,ZGFC_IMP,GCOND)
************************************************
      GIDS=-DIMAG(ZGF_IMP)
      RETURN
      END
      FUNCTION GFFC(EW)
*****************************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      ETTA=0.05D0
      ZW=DCMPLX(EW,ETTA)
      CALL GKONDO(ZW,ZGF_IMP,ZGC_IMP,ZGFC_IMP,GCOND)
************************************************
      ZFE=ZFER(ZW)
      GFFC=-DIMAG(ZFE*ZGFC_IMP)
      RETURN
      END
      SUBROUTINE DENSI
********************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      INTEGER DT10
      PARAMETER (DT10=10)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/ET/ETA,OI,OF,NT
      OPEN (UNIT=DT10,FILE='adensi.dat',STATUS='UNKNOWN')
      DMT=(OF-OI)/(NT-1)
      DO 241 I=1,NT
      AW=OI+(I-1)*DMT
      ZW=DCMPLX(AW,ETA)
      
************************************
C     CALL GREEN(ZW,ZGFEXATA,ZGCEXATA)
      CALL GKONDO(ZW,ZGF_IMP,ZGC_IMP,ZGFC_IMP,GCOND)
**********************************************
      RF_IMP=-CC*DIMAG(ZGF_IMP)
      RC_IMP=CC*DIMAG(ZGC_IMP)
      RFC_IMP=-CC*DIMAG(ZGFC_IMP)
      WRITE(DT10,121)AW,RF_IMP,RC_IMP
****************************************************
      
C........FREQUENCIA X DENSIDADE DE ESTADOS...............
C     RHO=parte imaginaria da FG exata(DENSIDADE DE ESTADOS)
C     AW=parte real da frequencia omega
      
      RHOF=CC*DIMAG(ZGFEXATA)
      RHOC=CC*DIMAG(ZGCEXATA)
C     WRITE(DT10,121)AW,RHOF,RHOC
241   CONTINUE
 121  FORMAT(10F15.8)
      RETURN
      END
      SUBROUTINE BISSECAO_OLD(A,B,R,ER)
***********************************
C       RAIZES DE UMA EQ. PELO METODO DA BISSECAO
        IMPLICIT COMPLEX*16 (X-Z)
        IMPLICIT REAL*8 (A-H,O-W)
        INTEGER DT50
        PARAMETER (DT50=50)
        COMMON/TOTAL/SNT
        OPEN (UNIT=DT50,FILE='raiz.dat',STATUS='UNKNOWN')
        DATA N /1000/
        EPS=1.E-8
        DO 10 I=1,N
        AMED=0.5D0*(A+B)
        FA=FUNCTMU(AMED)
        F1=FUNCTMU(A)
        AP1=FA*F1
        IF(AP1.GT.0)GO TO 20
        B=AMED
        GO TO 30
20      CONTINUE
        A=AMED
30      CONTINUE
        ADAB=ABS(A-B)
        IF(ADAB.LT.EPS)GO TO 40
10      CONTINUE
40      CONTINUE
        R=AMED
        WRITE(6,*)'ZZZZZZZZZZZZZZZ=',R
        ER=ADAB
        RETURN
        END
      FUNCTION FUNCTMU_OLD(DEQ)
************************************
      IMPLICIT REAL*8 (A-H,O-W)
      IMPLICIT COMPLEX*16 (X-Z)
      REAL*8 MU(9),MUC(9),UU(9),UE(9)
      REAL*8 AMU
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/OCC/TXF,TDS,TXFC,TXC,TOT,TSOMA,FRIED,R0F,R0C
      EQ=AMU+DEQ
      CALL EXATA
      CALL OCUP3
********************
      FUNCTMU=FRIED-R0F
C      FUNCTMU=TSOMA-1.D0
      WRITE(6,*)'vvvvvvvvvvvvvvvvvvvvvvv=',FUNCTMU
      RETURN
      END
        FUNCTION FMU(DEQ)
************************************
      IMPLICIT REAL*8 (A-H,O-W)
      IMPLICIT COMPLEX*16 (X-Z)
      REAL*8 MU(9),MUC(9),UU(9),UE(9)
      REAL*8 AMU
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      COMMON/OCC/TXF,TDS,TXFC,TXC,TOT,TSOMA,FRIED,R0F,R0C
      EQ=AMU+DEQ
      CALL EXATA
      CALL OCUP3
********************
      FMU=DABS(FRIED-R0F)
C      FMU=DABS(TSOMA-1.D0)
      RETURN
      END
      FUNCTION FERM(AX)
C     *****************************
      IMPLICIT COMPLEX*16 (X-Z)
      IMPLICIT REAL*8 (A-H,O-W)
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      BEXP=350.D0
      T=1.D0/BETA
      ABSX=DABS(AX)
      IF(ABSX - T*BEXP)250,200,200
 200  IF(AX)210,220,230
 210  FERM=1.D0
      RETURN
 220  FERM=0.5D0
      RETURN
 230  FERM=0.D0
      RETURN
 250  ARRG=AX*BETA
      FERM=1.D0/(1.D0+DEXP(ARRG))
 300  RETURN
      END
      FUNCTION ZFER(ZW)
      
***************************
      
      IMPLICIT REAL*8 (A-H,O-W)
      
      IMPLICIT COMPLEX*16 (X-Z)
      
      REAL*8 MU
      
      COMMON/DADOS/EF,EQ,AMU,V,V2,BETA,DS,CC,D
      
      DATA BEXP/35.D0/
      
      DATA ZUM/(1.D0,0.D0)/
      
      Z=ZW
      
      ZEX=Z*BETA
      
      T=1.D0/BETA
      
      RX=DREAL(Z)
      
      RY=DIMAG(Z)
      
      IF(RY.NE.0.D0)GO TO 100
      
      ZFER=DCMPLX(FERM(RX),0.D0)
      
      RETURN
      
 100  ABSX=DABS(RX)
      
      IF(ABSX - T*BEXP) 250,200,200
      
 200  IF(RX)210,220,230
      
 210  ZFER= (1.0D0,0.D0)
      
      RETURN
      
 220  ZFER=(0.5D0,0.D0)
      
      RETURN
      
 230  ZFER=(0.D0,0.D0)
      
      RETURN
      
250       CONTINUE
          
          ZFER=ZUM/(ZUM+CDEXP(ZEX))
          
          RETURN
          
          END

      
      Function Fmin(ax,bx,f,Tol)
      Implicit double precision (a-h,o-z)
      External f
      C=0.5d0*(3.0d0-dsqrt(5.0d0))
      Eps=1.0d0
 10   Eps=Eps/2.0d0
      Tol1=1.0d0+Eps
      If(tol1.gt.1.0d0) go to 10
      Eps=dsqrt(Eps)
      A=ax
      B=bx
      V=A+C*(B-A)
      W=V
      X=V
      E=0.0d0
      FX=f(X)
      FV=FX
      FW=FX
 20   XM=0.5d0*(A+B)
      TOL1=Eps*Dabs(X)+tol/3.0d0
      TOL2=2.0d0*TOL1
      
      
      If(Dabs(X-XM).le.(TOL2-0.5d0*(B-A))) go to 90
      
      If(Abs(E).le.Tol1) go to 40
      
      R=(X-W)*(FX-FV)
      Q=(X-V)*(FX-FW)
      P=(X-V)*Q-(X-W)*R
      Q=2.0d0*(Q-R)
      If(Q.gt.0.0d0) P=-P
      Q=Dabs(Q)
      R=E
      E=D
c     
c     ¨ Parabola es aceptable ?
c     
 30   If(dabs(P).ge.dabs((0.5d0*Q*R))) go to 40
      If(P.ge.Q*(A-X)) go to 40
      If(P.ge.Q*(B-X)) go to 40
c     
c     El paso de la interpolacion parabolica
c     
      D=P/Q
      U=X+D
c     
c     F no se permite calcular demasiado cerca al AX o BX
c     
      If((U-A).lt.Tol2) D=Dsign(Tol1,XM-X)
      If((B-U).lt.Tol2) D=Dsign(Tol1,XM-X)
      Go to 50
c     
c     El paso de seccion de oro
c     
 40   If(X.ge.XM) E=A-X
      If(X.lt.XM) E=B-X
      D=C*E
c     
c     F no se permite calcular demasiado cerca al X
c     
 50   If(Dabs(D).ge.Tol1) U=X+D
      If(Dabs(D).lt.Tol1) U=X+Dsign(Tol1,D)
      FU=f(u)
c     
c     Asignar valores nuevos a los parametros A,B,V,W y X
c     
      If(FU.gt.FX) go to 60
      If(U.ge.X) A=X
      If(U.lt.X) B=X
      V=W
      FV=FW
      W=X
      FW=FX
      X=U
      FX=FU
      go to 20
 60   if(U.lt.X) A=U
      If(U.ge.X) B=U
      If(FU.le.FW) go to 70
      If(W.eq.X) go to 70
      If(FU.le.FV) go to 80
      If(V.eq.X) go to 80
      If(V.eq.W) go to 80
      go to 20
 70   V=W
      FV=FW
      W=U
      FW=FU
      go to 20
 80   V=U
      FV=FU
      go to 20
c     
c     Fin del ciclo principal
c     
 90   Fmin=X
      Return
      End
      SUBROUTINE QUAD(A,B,RESULT,K,EPSIL,NPTS,ICHECK,F)
********************************************************
      IMPLICIT REAL*8 (A-H,P-Z)
      DIMENSION FUNCT(127),P(381),RESULT(8)
C     C      THIS SUBROUTINE ATTEMPTS TO CALCULATE THE INTEGRAL OF F(X)
C     C      OVER THE INTERVAL  (A,B) WITH RELATIVE ERROR NOT
C     C      EXCEEDING  EPSIL.
C     C      THE RESULT IS OBTAINED USING A SEQUENCE OF 1, 3, 7, 15, 31, 63,
C     C      127 AND 255 POINT INTERLACING FORMALAE (NO INTEGRAND
C     C      EVALUATIONS ARE WASTED) OF RESPECTIVE DEGREE 1,5,11,23,
C     C      47,95,191 AND 383. THE FORMULAE ARE BASED ON THE OPTIMAL
C     C      EXTENSION OF THE 3-POINT GAUSS FORMULA. DETAILS  OF
C     C      THE FORMULAE ARE GIVEN IN THE OPTIMUM ADDITION OF POINTS
C     C      TO QUADRATURE FORMULAE  BY T. N. L. PATRERSON,MATHS.COMP.
C     VOL. 22,847-856,1968.
C     ******INPUT *****
C     C      A          LOWER LIMIT OF INTEGRATION.
C     B          UPPER LIMIT OF INTEGRATION.
C     EPSIL      RELATIVE ACCURANCY REQUIRED. WHEN THE RELATIVE
C     DIFFERENCE OF TWO SUCCESSIVE FORMULAE DOES NOT
C     EXCEED  *EPSL* THE LAST FORMULA COMPUTED IS TAKEN
C     AS THE RESULT.
C     F          F(X) IS THE INTEGRAND.
C     ******OUTPUT *****
C     RESULT     THIS ARRAY,WHICH SHOULD BE DECLARED TO HAVE AT
C     LEAST 8 ELEMENTS,HOLDS THE RESULTS OBTAINED BY
C     THE 1,3,7, ETC., POINT FORMULAE.THE NUMBER OF
C     FORMULAE COMPUTED DEPENDS ON *EPSIL*
C     K           RESULT(K) HOLDS THE VALUE OF THE INTERVAL TO THE
C     SPECIFIED RELATIVE ACCURANCY.
C     NPTS       NUMBER INTEGRAND EVALUATIONS.
C     ICHECK     ON  EXIT NORMALLY ICHECK =0. HOWEVER IF CONVERGENCE
C     TO THE ACCURANCY REQUESTED IS NOT ACHIEVED ICHECK=1
C     ON EXIT
C     C      ABSCISSAE AND WEIGHTS OF QUADRATURE RULES ARE STACKED IN
C     ARRAY *P* IN THE ORDER IN WHICH THEY ARE NEEDED.
      DATA
     *     P( 1),P( 2),P( 3),P( 4),P( 5),P( 6),P( 7),
     *     P( 8),P( 9),P(10),P(11),P(12),P(13),P(14),
     *     P(15),P(16),P(17),P(18),P(19),P(20),P(21),
     *     P(22),P(23),P(24),P(25),P(26),P(27),P(28)/
     *     0.77459666924148337704D+00,0.55555555555555555556D+00,
     *     0.88888888888888888889D+00,0.26848808986833344073D+00,
     *     0.96049126870802028342D+00,0.10465622602646726519D+00,
     *     0.43424374934680255800D+00,0.40139741477596222291D+00,
     *     0.45091653865847414235D+00,0.13441525524378422036D+00,
     *     0.51603282997079739697D-01,0.20062852937698902103D+00,
     *     0.99383196321275502221D+00,0.17001719629940260339D-01,
     *     0.88845923287225699889D+00,0.92927195315124537686D-01,
     *     0.62110294673722640294D+00,0.17151190913639138079D+00,
     *     0.22338668642896688163D+00,0.21915685840158749640D+00,
     *     0.22551049979820668739D+00,0.67207754295990703540D-01,
     *     0.25807598096176653565D-01,0.10031427861179557877D+00,
     *     0.84345657393211062463D-02,0.46462893261757986541D-01,
     *     0.85755920049990351154D-01,0.10957842105592463824D+00/
      DATA
     *     P(29),P(30),P(31),P(32),P(33),P(34),P(35),
     *     P(36),P(37),P(38),P(39),P(40),P(41),P(42),
     *     P(43),P(44),P(45),P(46),P(47),P(48),P(49),
     *     P(50),P(51),P(52),P(53),P(54),P(55),P(56)/
     *     0.99909812496766759766D+00,0.25447807915618744154D-02,
     *     0.98153114955374010687D+00,0.16446049854387810934D-01,
     *     0.92965485742974005667D+00,0.35957103307129322097D-01,
     *     0.83672593816886873550D+00,0.56979509494123357412D-01,
     *     0.70249620649152707861D+00,0.76879620499003531043D-01,
     *     0.53131974364437562397D+00,0.93627109981264473617D-01,
     *     0.33113539325797683309D+00,0.10566989358023480974D+00,
     *     0.11248894313318662575D+00,0.11195687302095345688D+00,
     *     0.11275525672076869161D+00,0.33603877148207730542D-01,
     *     0.12903800100351265626D-01,0.50157139305899537414D-01,
     *     0.42176304415588548391D-02,0.23231446639910269443D-01,
     *     0.42877960025007734493D-01,0.54789210527962865032D-01,
     *     0.12651565562300680114D-02,0.82230079572359296693D-02,
     *     0.17978551568128270333D-01,0.28489754745833548613D-01/
      DATA
     *     P(57),P(58),P(59),P(60),P(61),P(62),P(63),
     *     P(64),P(65),P(66),P(67),P(68),P(69),P(70),
     *     P(71),P(72),P(73),P(74),P(75),P(76),P(77),
     *     P(78),P(79),P(80),P(81),P(82),P(83),P(84)/
     *     0.38439810249455532039D-01,0.46813554990628012403D-01,
     *     0.52834946790116519862D-01,0.55978436510476319408D-01,
     *     0.99987288812035761194D+00,0.36322148184553065969D-03,
     *     0.99720625937222195908D+00,0.25790497946856882724D-02,
     *     0.98868475754742947994D+00,0.61155068221172463397D-02,
     *     0.97218287474858179658D+00,0.10498246909621321898D-01,
     *     0.94634285837340290515D+00,0.15406750466559497802D-01,
     *     0.91037115695700429250D+00,0.20594233915912711149D-01,
     *     0.86390793819369047715D+00,0.25869679327214746911D-01,
     *     0.80694053195021761186D+00,0.31073551111687964880D-01,
     *     0.73975604435269475868D+00,0.36064432780782572640D-01,
     *     0.66290966002478059546D+00,0.40715510116944318934D-01,
     *     0.57719571005204581484D+00,0.44914531653632197414D-01,
     *     0.48361802694584102756D+00,0.48564330406673198716D-01/
      DATA
     *     P(85),P(86),P(87),P(88),P(89),P(90),P(91),
     *     P(92),P(93),P(94),P(95),P(96),P(97),P(98),
     *     P(99),P(100),P(101),P(102),P(103),P(104),P(105),
     *     P(106),P(107),P(108),P(109),P(110),P(111),P(112)/
     *     0.38335932419873034692D 00,0.51583253952048458777D-01,
     *     0.27774982202182431507D 00,0.53905499335266063927D-01,
     *     0.16823525155220746498D 00,0.55481404356559363988D-01,
     *     0.56344313046592789972D-01,0.56277699831254301273D-01,
     *     0.56377628360384717388D-01,0.16801938574103865271D-01,
     *     0.64519000501757369228D-02,0.25078569652949768707D-01,
     *     0.21088152457266328793D-02,0.11615723319955134727D-01,
     *     0.21438980012503867246D-01,0.27394605263981432516D-01,
     *     0.63260731936263354422D-03,0.41115039786546930472D-02,
     *     0.89892757840641357233D-02,0.14244877372916774306D-01,
     *     0.19219905124727766019D-01,0.23406777495314006201D-01,
     *     0.26417473395058259931D-01,0.27989218255238159704D-01,
     *     0.18073956444538835782D-03,0.12895240826104173921D-02,
     *     0.30577534101755311361D-02,0.52491234548088591251D-02/
      DATA
     *     P(113),P(114),P(115),P(116),P(117),P(118),P(119),
     *     P(120),P(121),P(122),P(123),P(124),P(125),P(126),
     *     P(127),P(128),P(129),P(130),P(131),P(132),P(133),
     *     P(134),P(135),P(136),P(137),P(138),P(139),P(140)/
     *     0.77033752332797418482D-02,0.10297116957956355524D-01,
     *     0.12934839663607373455D-01,0.15536775555843982440D-01,
     *     0.18032216390391286320D-01,0.20357755058472159467D-01,
     *     0.22457265826816098707D-01,0.24282165203336599358D-01,
     *     0.25791626976024229388D-01,0.26952749667633031963D-01,
     *     0.27740702178279681994D-01,0.28138849915627150636D-01,
     *     0.99998243035489159858D 00,0.50536095207862517625D-04,
     *     0.99959879967191068325D 00,0.37774664632698466027D-03,
     *     0.99831663531840739253D 00,0.93836984854238150079D-03,
     *     0.99572410469840718851D 00,0.16811428654214699063D-02,
     *     0.99149572117810613240D 00,0.25687649437940203731D-02,
     *     0.98537149959852037111D 00,0.35728927835172996494D-02,
     *     0.97714151463970571416D 00,0.46710503721143217474D-02,
     *     0.96663785155841656709D 00,0.58434498758356395076D-02/
      DATA
     *     P(141),P(142),P(143),P(144),P(145),P(146),P(147),
     *     P(148),P(149),P(150),P(151),P(152),P(153),P(154),
     *     P(155),P(156),P(157),P(158),P(159),P(160),P(161),
     *     P(162),P(163),P(164),P(165),P(166),P(167),P(168)/
     *     0.95373000642576113641D 00,0.70724899954335554680D-02,
     *     0.93832039777959288365D 00,0.83428387539681577056D-02,
     *     0.92034002547001242073D 00,0.96411777297025366953D-02,
     *     0.89974489977694003664D 00,0.10955733387837901648D-01,
     *     0.87651341448470526974D 00,0.12275830560082770087D-01,
     *     0.85064449476835027976D 00,0.13591571009765546790D-01,
     *     0.82215625436498040737D 00,0.14893641664815182035D-01,
     *     0.79108493379984836143D 00,0.16173218729577719942D-01,
     *     0.75748396638051363793D 00,0.17421930159464173747D-01,
     *     0.72142308537009891548D 00,0.18631848256138790186D-01,
     *     0.68298743109107922809D 00,0.19795495048097499488D-01,
     *     0.64227664250975951377D 00,0.20905851445812023852D-01,
     *     0.59940393024224289297D 00,0.21956366305317824939D-01,
     *     0.55449513263193254887D 00,0.22940964229387748761D-01/
      DATA
     *P(169),P(170),P(171),P(172),P(173),P(174),P(175),
     *     P(176),P(177),P(178),P(179),P(180),P(181),P(182),
     *     P(183),P(184),P(185),P(186),P(187),P(188),P(189),
     *     P(190),P(191),P(192),P(193),P(194),P(195),P(196)/
     *     0.50768775753371660215D 00,0.23854052106038540080D-01,
     *     0.45913001198983233287D 00,0.24690524744487676909D-01,
     *     0.40897982122988867241D 00,0.25445769965464765813D-01,
     *     0.35740383783153215238D 00,0.26115673376706097680D-01,
     *     0.30457644155671404334D 00,0.26696622927450359906D-01,
     *     0.25067873030348317661D 00,0.27185513229624791819D-01,
     *     0.19589750271110015392D 00,0.27579749566481873035D-01,
     *     0.14042423315256017459D 00,0.27877251476613701609D-01,
     *     0.84454040083710883710D-01,0.28076455793817246607D-01,
     *     0.28184648949745694339D-01,0.28176319033016602131D-01,
     *     0.28188814180192358694D-01,0.84009692870519326354D-02,
     *     0.32259500250878684614D-02,0.12539284826474884353D-01,
     *     0.10544076228633167722D-02,0.58078616599775673635D-02,
     *     0.10719490006251933623D-01,0.13697302631990716258D-01/
      DATA
     *     P(197),P(198),P(199),P(200),P(201),P(202),P(203),
     *     P(204),P(205),P(206),P(207),P(208),P(209),P(210),
     *     P(211),P(212),P(213),P(214),P(215),P(216),P(217),
     *     P(218),P(219),P(220),P(221),P(222),P(223),P(224)/
     *     0.31630366082226447689D-03,0.20557519893273465236D-02,
     *     0.44946378920320678616D-02,0.71224386864583871532D-02,
     *     0.96099525623638830097D-02,0.11703388747657003101D-01,
     *     0.13208736697529129966D-01,0.13994609127619079852D-01,
     *     0.90372734658751149261D-04,0.64476204130572477933D-03,
     *     0.15288767050877655684D-02,0.26245617274044295626D-02,
     *     0.38516876166398709241D-02,0.51485584789781777618D-02,
     *     0.64674198318036867274D-02,0.77683877779219912200D-02,
     *     0.90161081951956431600D-02,0.10178877529236079733D-01,
     *     0.11228632913408049354D-01,0.12141082601668299679D-01,
     *     0.12895813488012114694D-01,0.13476374833816515982D-01,
     *     0.13870351089139840997D-01,0.14069424957813575318D-01,
     *     0.25157870384280661489D-04,0.18887326450650491366D-03,
     *     0.46918492424785040975D-03,0.84057143271072246365D-03/
      DATA
     *     P(225),P(226),P(227),P(228),P(229),P(230),P(231),
     *     P(232),P(233),P(234),P(235),P(236),P(237),P(238),
     *     P(239),P(240),P(241),P(242),P(243),P(244),P(245),
     *     P(246),P(247),P(248),P(249),P(250),P(251),P(252)/
     *     0.12843824718970101768D-02,0.17864463917586498247D-02,
     *     0.23355251860571608737D-02,0.29217249379178197538D-02,
     *     0.35362449977167777340D-02,0.41714193769840788528D-02,
     *     0.48205888648512683476D-02,0.54778666939189508240D-02,
     *     0.61379152800413850435D-02,0.67957855048827733948D-02,
     *     0.74468208324075910174D-02,0.80866093647888599710D-02,
     *     0.87109650797320868736D-02,0.93159241280693950932D-02,
     *     0.98977475240487497440D-02,0.10452925722906011926D-01,
     *     0.10978183152658912470D-01,0.11470482114693874380D-01,
     *     0.11927026053019270040D-01,0.12345262372243838455D-01,
     *     0.12722884982732382906D-01,0.13057836688353048840D-01,
     *     0.13348311463725179953D-01,0.13592756614812395910D-01,
     *     0.13789874783240936517D-01,0.13938625738306850804D-01,
     *     0.14038227896908623303D-01,0.14088159516508301065D-01/
      DATA
     *     P(253),P(254),P(255),P(256),P(257),P(258),P(259),
     *     P(260),P(261),P(262),P(263),P(264),P(265),P(266),
     *     P(267),P(268),P(269),P(270),P(271),P(272),P(273),
     *     P(274),P(275),P(276),P(277),P(278),P(279),P(280)/
     *     0.99999759637974846462D 00,0.69379364324108267170D-05,
     *     0.99994399620705437576D 00,0.53275293669780613125D-04,
     *     0.99976049092443204733D 00,0.13575491094922871973D-03,
     *     0.99938033802502358193D 00,0.24921240048299729402D-03,
     *     0.99874561446809511470D 00,0.38974528447328229322D-03,
     *     0.99780535449595727456D 00,0.55429531493037471492D-03,
     *     0.99651414591489027385D 00,0.74028280424450333046D-03,
     *     0.99483150280062100052D 00,0.94536151685852538246D-03,
     *     0.99272134428278861533D 00,0.11674841174299594077D-02,
     *     0.99015137040077015918D 00,0.14049079956551446427D-02,
     *     0.98709252795403406719D 00,0.16561127281544526052D-02,
     *     0.98351865757863272876D 00,0.19197129710138724125D-02,
     *     0.97940628167086268381D 00,0.21944069253638388388D-02,
     *     0.97473445975240266776D 00,0.24789582266575679307D-02/
      DATA
     *     P(281),P(282),P(283),P(284),P(285),P(286),P(287),
     *     P(288),P(289),P(290),P(291),P(292),P(293),P(294),
     *     P(295),P(296),P(297),P(298),P(299),P(300),P(301),
     *     P(302),P(303),P(304),P(305),P(306),P(307),P(308)/
     *     0.96948465950245923177D 00,0.27721957645934509940D-02,
     *     0.96364062156981213252D 00,0.30730184347025783234D-02,
     *     0.95718821610986096274D 00,0.33803979910869203823D-02,
     *     0.95011529752129487656D 00,0.36933779170256508183D-02,
     *     0.94241156519108305981D 00,0.40110687240750233989D-02,
     *     0.93406843615772578800D 00,0.43326409680929828545D-02,
     *     0.92507893290707565236D 00,0.46573172997568547773D-02,
     *     0.91543758715576504064D 00,0.49843645647655386012D-02,
     *     0.90514035881326159519D 00,0.53130866051870565663D-02,
     *     0.89418456833555902286D 00,0.56428181013844441585D-02,
     *     0.88256884024734190684D 00,0.59729195655081658049D-02,
     *     0.87029305554811390585D 00,0.63027734490857587172D-02,
     *     0.85735831088623215653D 00,0.66317812429018878941D-02,
     *     0.84376688267270860104D 00,0.69593614093904229394D-02/
      DATA
     *     P(309),P(310),P(311),P(312),P(313),P(314),P(315),
     *     P(316),P(317),P(318),P(319),P(320),P(321),P(322),
     *     P(323),P(324),P(325),P(326),P(327),P(328),P(329),
     *     P(330),P(331),P(332),P(333),P(334),P(335),P(336)/
     *     0.82952219463740140018D 00,0.72849479805538070639D-02,
     *     0.81462878765513741344D 00,0.76079896657190565832D-02,
     *     0.79909229096084140180D 00,0.79279493342948491103D-02,
     *     0.78291939411828301639D 00,0.82443037630328680306D-02,
     *     0.76611781930376009072D 00,0.85565435613076896192D-02,
     *     0.74869629361693660282D 00,0.88641732094824942641D-02,
     *     0.73066452124218126133D 00,0.91667111635607884067D-02,
     *     0.71203315536225203459D 00,0.94636899938300652943D-02,
     *     0.69281376977911470289D 00,0.97546565363174114611D-02,
     *     0.67301883023041847920D 00,0.10039172044056840798D-01,
     *     0.65266166541001749610D 00,0.10316812330947621682D-01,
     *     0.63175643771119423041D 00,0.10587167904885197931D-01,
     *     0.61031811371518640016D 00,0.10849844089337314099D-01,
     *     0.58836243444766254143D 00,0.11104461134006926537D-01/
      DATA
     *     P(337),P(338),P(339),P(340),P(341),P(342),P(343),
     *     P(344),P(345),P(346),P(347),P(348),P(349),P(350),
     *     P(351),P(352),P(353),P(354),P(355),P(356),P(357),
     *     P(358),P(359),P(360),P(361),P(362),P(363),P(364)/
     *     0.56590588542365442262D 00,0.11350654315980596602D-01,
     *     0.54296566649831149049D 00,0.11588074033043952568D-01,
     *     0.51955966153745702199D 00,0.11816385890830235763D-01,
     *     0.49570640791876146017D 00,0.12035270785279562630D-01,
     *     0.47142506587165887693D 00,0.12244424981611985899D-01,
     *     0.44673538766202847374D 00,0.12443560190714035263D-01,
     *     0.42165768662616330006D 00,0.12632403643542078765D-01,
     *     0.39621280605761593918D 00,0.12810698163877361967D-01,
     *     0.37042208795007823014D 00,0.12978202239537399286D-01,
     *     0.34430734159943802278D 00,0.13134690091960152836D-01,
     *     0.31789081206847668318D 00,0.13279951743930530650D-01,
     *     0.29119514851824668196D 00,0.13413793085110098513D-01,
     *     0.26424337241092676194D 00,0.13536035934956213614D-01,
     *     0.23705884558982972721D 00,0.13646518102571291428D-01/
      DATA
     *     P(365),P(366),P(367),P(368),P(369),P(370),P(371),
     *     P(372),P(373),P(374),P(375),P(376),P(377),P(378),
     *     P(379),P(380),P(381)/
     *     0.20966523824318119477D 00,0.13745093443001896632D-01,
     *     0.18208649675925219825D 00,0.13831631909506428676D-01,
     *     0.15434681148137810869D 00,0.13906019601325461264D-01,
     *     0.12647058437230196685D 00,0.13968158806516938516D-01,
     *     0.98482396598119202090D-01,0.14017968039456608810D-01,
     *     0.70406976042855179063D-01,0.14055382072649964277D-01,
     *     0.42269164765363603212D-01,0.14080351962553661325D-01,
     *     0.14093886410782462614D-01,0.14092845069160408355D-01,
     *     0.14094407090096179347D-01/
      ZERO=0.D0
      DOIS=2.D0
      ICHECK=0
C     CHECK FOR TRIVIAL CASE
      IF(A.EQ.B)GO TO 70
C       SCALE FACTORS
      SUM=(B+A)/DOIS
      DIFF=(B-A)/DOIS
C     1-POIN GAUSS
      FZERO=F(SUM)
      RESULT(1)=DOIS*FZERO*DIFF
      I=0
      IOLD=0
      INEW=1
      K=2
      ACUM=ZERO
      GO TO 30
 10   IF(K.EQ.8)GO TO 50
      K=K+1
      ACUM=ZERO
C     CONTRIBUTION FROM FUNCTION VALUES ALREADY COMPUTED
      DO 20 J=1,IOLD
         I=I+1
         ACUM=ACUM+P(I)*FUNCT(J)
20      CONTINUE
C     CONTRIBUTION FROM NEW FUNCTION VALUES
 30     IOLD=IOLD+INEW
        DO 40 J=INEW,IOLD
           I=I+1
           X=P(I)*DIFF
           FUNCT(J)=F(SUM+X)+F(SUM-X)
           I=I+1
           ACUM=ACUM+P(I)*FUNCT(J)
40      CONTINUE
        INEW=IOLD+1
        I=I+1
        RESULT(K)=(ACUM+P(I)*FZERO)*DIFF
C     CHECK FOR CONVERGENCE
        IF(DABS(RESULT(K)-RESULT(K-1))-EPSIL*DABS(RESULT(K)))60,
     *       60,10
C     CONVERGENCE NOT ACHIEVED
 50     ICHECK=1
C     NORMAL TERMINATION
 60     NPTS=INEW+IOLD
        RETURN
C     TRIVIAL CASE
 70     K=2
        RESULT(1)=ZERO
        RESULT(2)=ZERO
        NPTS=0
        RETURN
        END
      FUNCTION QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,F)
C     C
C     C
C     C        THIS FUNCTION ROUTINE PERFORMS AUTOMATIC INTEGRATION
C       OVER A FINITE INTERVAL USING THE BASIC INTEGRATION
C     ALGORITHM QUAD TOGETHER WIHT, IF NECESSARY AN ADAPTIVE
C     SUBDIVISION PROCESS. IT IS GENERALLY MORE EFFICIENT THAN
C     THE NON-ADAPTIVE ALGORITHM QSUB BUT IS LIKILY TO BE LESS
C     RELIABLE(SEE COMP. J., 14,189,1971).
C     THE CALL TAKES THE FORM
C     QSUBA(A,B,EPSIL,NPTS,ICHECK,RELERR,F)
C     AND CAUSES F(X) TO BE  INTEGRATED OER (A,B) WIHT RELATIVE
C     ERROR HOPEFULLY NOT  EXCCEEDING EPSIL. SHOULD QUAD CONVERGE
C     (ICHECK=0) THEN QSUBA WILL RETURN THE VALUE OBTAINED BY IT
C     OTHERWISE SUBDIVISION  WILL BE INVOKED AS A RESCUE
C     OPERATION IN AN ADAPTIVE MANNER.THE ARGUMENT RELERR GIVES
C     A CRUDE ESTIMATE OF THE ACTUAL  RELATIVE ERROR OBATAINED.
C     THE SUBDIVISION STRATEGY IS AS FOLLOWS
C     AT EACH STAGE OF THE PROCESS AN INTERVAL IS PRESENTED FOR
C     SUBDIVISION (INITIALLY THIS WILL BE THE WHOLE INTERVAL
C     (A,B)) THE INTERVAL IS HALVED AND QUAD APPLIED TO EACH
C     SUBINTERVAL. SHOULD QUAD FAIL ON THE FIRST SUBINTERVAL
C     THE SUBINTERVAL IS STACKED FOR FUTURE SUBDIVISION AND  THE
C     SECOND SUBINTERVAL IMMEDIATILY EXAMINED. SHOULD QUAD FAIL
C     ON THE SECOND SUBINTERVAL THE INTERVAL IS
C     IMMEDIATILY SUBDIVIDED AND THE WHOLE PROCESS REPEATED
C     EACH TIME A CONVERGED RESULT IS OBTAINED IT IS
C     ACCUMULATED AS THE PARTIAL VALUE OF THE INTEGRAL. WHEN
C     QUAD CONVERGES ON BOTH SUBINTERVALS THE INTERVAL LAST
C     STACKED IS CHOSEN NEXT FOR SUBDIVISION AND THE PROCESS
C     REPEATED. A SUBINTERVAL IS NOT EXAMINED AGAIN ONCE A
C     CONVERGED RESULT IS OBTAINED FOR IT SO THAT A  SPURIOUS
C     CONVERGENGE IS MORE LIKELY TO SLIP THROUGH THAN FOR THE
C     NON-ADAPTIVE ALGORITHM QSUB.
C     THE CONVERGE CRITERION OF QUAD IS SLIGHTLY RELAXED
C     IN THAT A PANEL IS DEEMED TO HAVE BEEN SUCCESSFULLY
C     INTEGRATED IF EITHER QUAD CONVERGES OR THE ESTIMATED
C     ABSOLUTE ERRROR  COMMITTED ON THIS PANEL DOES NOT EXCCEED
C     EPSIL TIMES THE ESTIMATED ABSOLUTE VALUE OF THE INTEGRAL
C     OVER (A,B). THIS RELAXATION IS TO TRY TO TAKES ACCOUNT OF
C     A COMMON SITUATION WHERE ONE PARTICULAR PANEL CAUSES
C     SPECIAL DIFFICULTY ,PERHAPS DUE TO A SINGULARITY OF
C     SOME TYPE. IN THIS CASE QUAD COULD OBTAIN NEARLY EXACT
C     ANSWERS ON ALL OTHER PANELS AND SO THE RELATIVE ERROR FOR
C     THE TOTAL INTEGRATION WOULD BE ALMOST ENTIRELY DUE TO THE
C     DELINQUENT PANEL. WITHOUT THIS CONDITION THE COMPUTATION
C     MIGHT CONTINUE DESPITE THE REQUESTED RELATIVE ERROR BEING
C     ACHIEVED. IF THIS RELAXED CONVERGENCE CRITERION IS APPLIED TOO
C     MANY TIMES,I.E. ICHECK = 2 IN MANY PANELS, THE RELATIVE
C     ERROR WOULD BE LARGER THAN EPSIL. IN THIS CASE CHECK
C     EPSIL  VS  RELERR .
C     THE OUTCOME OF INTEGRATION IS INDICATED BY ICHECK.
C     ICHECK=0  CONVERGENCE OBTAINED WITHOUT INVOKING SUB-
C     DIVISION. THIS WOULD CORRESPOND TO THE
C     DIRECT USE OF QUAD.
C     ICHECK=1   RESULT OBTAINED AFTER INVOKING SUBDIVISION.
C     ICHECK=2   AS FOR ICHECK=1,BUT AT SOME POINT THE
C     RELAXED CONVERGENCE CRITERION WAS USED.
C     THE RISK OF UNDERESTIMATING THE RELATIVE
C     ERROR WILL BE INCREASED. IF NECESSARY,
C     CONFIDENCE MAY BE RESTORED BY CHEKING
C     EPSIL AND RELERR FOR A SERIOUS DISCREPANCY.
C     ICHECK NEGATIVE
C     IF DURING THE SUBDIVISION PROCESS THE STACK
C     OF DELINQUENT INTERVALS BECOMES FULL(IT IS
C     PRESENTLY SET TO HOLD AT MOST 100 NUMBERS)
C     A RESULT IS OBTAINED BY CONTINUING THE
C     INTEGRATION IGNORING CONVERGENCE FAILURES
C     WHICH CANNOT BE ACCOMODATED ON THE STACK.
C     THIS OCCURRENCE IS FLAGGED BY RETURNING
C     ICHECK WITH NEGATIVE SIGN.
C     THE RELIABILITY OF THE  ALGORITHM WILL DECREASE FOR LARGE
C     VALUES OF EPSIL.IT IS RECOMMENDED THAT EPSIL SHOULD
C     GENERALY BE LESS THAN ABOUT 0.001.
        IMPLICIT REAL*8 (A-H,P-Z)
        DIMENSION RESULT(8),STACK(100)
        EXTERNAL F
        ABCTR = 0.D0
        VALUE = QSUBC(A,B,EPSIL,NPTS,ICHECK,RELERR,F,ABCTR)
        QSUBA = VALUE
        RETURN
        END
      FUNCTION QSUBC(A,B,EPSIL,NPTS,ICHECK,RELERR,F,ABCTR)
C     
C     
C
C     ESTA SUBRUTINA TIENE UNA VARIABLE MAS:  ABCTR
C     SI  ABCTR = 0.D0  ES IGUAL A  QSUBA.
C     TIENE COMO OBJETO COMPARAR  EL RESULTADO DE
C     LA PRIMERA PASADA POR  QUAD  CON  EPSIL*ABCTR
C     CUANDO  DABS(QSUBC) .LT. EPSIL*DABS(ABCTR)  NO SUBDIVIDE.
C     DABS(QSUBC).LT.DABS(ABCTR)  ->  ESTIM = DABS(ABCTR)*EPSIL.
C     ESTA VARIANTE ES INTERESANTE QUANDO SE HA SUBDIVIDIDO LA
C     INTEGRAL EXTERNAMENTE A  QSUB2 , Y PUEDE PRETENDERSE CALCULAR
C     UNA DE LAS PARTES QUE ES RELATIVAMENTE MUY PEQUENA CON
C     DEMASIADA APROXIMACION.
C     HAY QUE TENER MUCHO CUIDADO CUANDO DIFERENTES PARTES TIENEN
C     SIGNO DIFERENTE, PUES EN ESE CASO SE COMPENSAN Y EL ERROR
C     RELATIVO PODRIA SER INADEQUADO. EN ESE CASO HAY QUE SEPARAR
C     LAS PARTES CON CUIDADO, PUES EL ORDEN DE SEPARACION AFECTARA
C     EL RESULTADO.
C     
C     
C     
C     
      IMPLICIT REAL*8 (A-H,P-Z)
      DIMENSION RESULT(8),STACK(100)
      EXTERNAL F
      DATA ISMAX/100/
      ABCTRQ = DABS(ABCTR)
      QSBMIN = EPSIL*ABCTRQ
      ZERO=0.D0
      CALL QUAD(A,B,RESULT,K,EPSIL,NPTS,ICHECK,F)
      QSUBC=RESULT(K)
      IF(ICHECK.NE.0)  GO TO 5
C     HAY CONVERGENCIA SIN SUBDIVISIONES ( ICHECK = 0 )  Y CALCULA  RELERR
C     PARA DEVOLVERLO A SUB. LLAMADA.
      RELERR = 0.D0
      IF(QSUBC.NE.ZERO)
     1     RELERR=DABS((RESULT(K)-RESULT(K-1))/QSUBC)
      RETURN
C     
C*****CHECK IF SUBDIVISION IS NEEDED
C     SI ABS(QSUBC) < ABCTR*EPSIL  SUSPENDE SUBDIVISIONES, PUES
C     EL VALOR ESTIMADO DE LA INTEGRAL ES DEMASIADO PEQUENO.
    5 ABQSUBC = DABS(QSUBC)
      IF(ABQSUBC.LT.QSBMIN) RETURN
C       SUBDIVIDED
      ESTIM=DABS(QSUBC*EPSIL)
C     CUANDO  ABS(QSUBC) < ABCTR  USA  ESTIM = EPSIL*ABCTR EN LUGAR
C     DE  ESTIM = EPSIL*ABC(QSUBC) PARA APLICAR CRITERIO RELAJADO
C     DE CONVERGENCIA. ESTO PUEDE EVITAR CALCULOS INNECESARIOS.
      IF(ABQSUBC.LT.ABCTRQ) ESTIM = QSBMIN
C     USA  RELERR PARA ACUMULAR EL ERROR ABSOLUTO. SOLO AL FINAL
C     LO DIVIDE POR  QSUBC  PARA OBTENER EL VERDADERO  RELERR.
      RELERR=0.0D0
      QSUBC=0.0D0
      IS=1
      IC=1
      SUB1=A
      SUB3=B
 10   SUB2=(SUB1+SUB3)*0.5D0
      CALL QUAD(SUB1,SUB2,RESULT,K,EPSIL,NF,ICHECK,F)
      NPTS=NPTS+NF
C     COMP ES EL VALOR ABSOLUTO DE LA DIFERENCIA ENTRE LAS DOS ULTIMAS
C     SUBDIVISIONES DE QUAD, QUE USA COMO VALOR ABSOLUTO Y ACUMULA
C     EN  RELERR
      COMP=DABS(RESULT(K)-RESULT(K-1))
      IF(ICHECK.EQ.0)GO TO 30
      IF(COMP.LE.ESTIM)GO TO 70
      IF(IS.GE.ISMAX)GO TO 20
C     STACK SUBINTERVAL (SUB1,SUB2) FOR FURTURE EXAMINATION
      STACK(IS)=SUB1
      IS=IS+1
      STACK(IS)=SUB2
      IS=IS+1
      GO TO 40
 20   IC=-IABS(IC)
 30   QSUBC=QSUBC+RESULT(K)
C     ACUMULA ERRORES ABSOLUTOS EN  RELERR
      RELERR=RELERR+COMP
 40   CALL QUAD(SUB2,SUB3,RESULT,K,EPSIL,NF,ICHECK,F)
      NPTS=NPTS+NF
      COMP=DABS(RESULT(K)-RESULT(K-1))
      IF(ICHECK.EQ.0)GO TO 50
      IF(COMP.LE.ESTIM)GO TO 80
C     SUBDIVIDE INTEVAL (SUB2,SUB3)
      SUB1=SUB2
      GO TO 10
 50   QSUBC=QSUBC+RESULT(K)
      RELERR=RELERR+COMP
      IF(IS.EQ.1) GO TO 60
C     SUBDIVIDE THE DELINQUENT INTERVAL LAST STACKED
      IS=IS-1
      SUB3=STACK(IS)
      IS=IS-1
      SUB1=STACK(IS)
      GO TO 10
C     SUBDIVISION RESULT
 60   ICHECK=IC
      IF(QSUBC.NE.0.D0)  GO TO 65
C     SI EL ERROR ABSOLUTO ESTIMADO NO ES NULO, Y  QSUBC = 0.D0,
C     DEVUELVE ERROR ABSOLUTO EN  RELERR CON SIGNO CAMBIADO, Y
C     AVISA EN LA TELA.
      IF(RELERR.EQ.0.D0)  RETURN
      WRITE(*,*) ' ATENCION !  QSUBC = ',QSUBC,'  RELERR = ',RELERR
      RELERR = - RELERR
      WRITE(*,*) ' CAMBIA SIGNO DE  RELERR EN LA SALIDA! '
 65   RELERR=RELERR/DABS(QSUBC)
      RETURN
C     RELAXED CONVERGENCE
 70   IC=ISIGN(2,IC)
      GO TO 30
 80   IC=ISIGN(2,IC)
      GO TO 50
      END
C-----------------------------------------------------------------------
C     
C     COMPUTER            - VAX/DOUBLE
C     
C     LATEST REVISION     - JANUARY 1, 1978
C     
C     PURPOSE             - ZERO OF A FUNCTION WHICH CHANGES SIGN IN A
C     GIVEN INTERVAL (BRENT ALGORITHM)
C     
C     USAGE               - CALL ZBRENT (F,EPS,NSIG,A,B,MAXFN,IER)
C     
C     ARGUMENTS    F      - AN EXTERNAL FUNCTION SUBPROGRAM F(X)
C     PROVIDED BY THE USER WHICH COMPUTES F FOR
C     ANY X IN THE INTERVAL (A,B). (INPUT)
C     F MUST APPEAR IN AN EXTERNAL STATEMENT IN
C     THE CALLING PROGRAM
C     EPS    - FIRST CONVERGENCE CRITERION (INPUT).  A ROOT,
C     B, IS ACCEPTED IF ABS(F(B)) IS LESS THAN OR
C     EQUAL TO EPS.  EPS MAY BE SET TO ZERO.
C     NSIG   - SECOND CONVERGENCE CRITERION (INPUT).  A ROOT,
C     B, IS ACCEPTED IF THE CURRENT APPROXIMATION
C     AGREES WITH THE TRUE SOLUTION TO NSIG
C     SIGNIFICANT DIGITS.
C     A,B    - ON INPUT, THE USER MUST SUPPLY TWO POINTS, A
C     AND B, SUCH THAT F(A) AND F(B) ARE OPPOSITE
C     IN SIGN.
C     ON OUTPUT, BOTH A AND B ARE ALTERED.  B
C     WILL CONTAIN THE BEST APPROXIMATION TO THE
C     ROOT OF F. SEE REMARK 1.
C     MAXFN  - ON INPUT, MAXFN SHOULD CONTAIN AN UPPER BOUND
C     ON THE NUMBER OF FUNCTION EVALUATIONS
C     REQUIRED FOR CONVERGENCE.  ON OUTPUT, MAXFN
C     WILL CONTAIN THE ACTUAL NUMBER OF FUNCTION
C     EVALUATIONS USED.
C     IER    - ERROR PARAMETER. (OUTPUT)
C     TERMINAL ERROR
C     IER = 129 INDICATES THE ALGORITHM FAILED TO
C     CONVERGE IN MAXFN EVALUATIONS.
C     IER = 130 INDICATES F(A) AND F(B) HAVE THE
C     SAME SIGN.
C     
C     PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C     - SINGLE/H36,H48,H60
C     
C     REQD. IMSL ROUTINES - UERTST,UGETIO
C     
C     NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C     CONVENTIONS IS AVAILABLE IN THE MANUAL
C     INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C     
C     REMARKS  1.  ON EXIT FROM ZBRENT, WHEN IER=0, A AND B SATISFY THE
C     FOLLOWING,
C     F(A)*F(B) .LE.0,
C     ABS(F(B)) .LE. ABS(F(A)), AND
C     EITHER ABS(F(B)) .LE. EPS OR
C     ABS(A-B) .LE. MAX(ABS(B),0.1)*10.0**(-NSIG).
C     THE PRESENCE OF 0.1 IN THIS ERROR CRITERION CAUSES
C     LEADING ZEROES TO THE RIGHT OF THE DECIMAL POINT TO BE
C     COUNTED AS SIGNIFICANT DIGITS. SCALING MAY BE REQUIRED
C     IN ORDER TO ACCURATELY DETERMINE A ZERO OF SMALL
C     MAGNITUDE.
C            2.  ZBRENT IS GUARANTEED TO REACH CONVERGENCE WITHIN
C     K = (ALOG((B-A)/D)+1.0)**2 FUNCTION EVALUATIONS WHERE
C     D=MIN(OVER X IN (A,B) OF
C     MAX(ABS(X),0.1)*10.0**(-NSIG)).
C                THIS IS AN UPPER BOUND ON THE NUMBER OF EVALUATIONS.
C     RARELY DOES THE ACTUAL NUMBER OF EVALUATIONS USED BY
C     ZBRENT EXCEED SQRT(K). D CAN BE COMPUTED AS FOLLOWS,
C     P = AMIN1(ABS(A),ABS(B))
C     P = AMAX1(0.1,P)
C     IF ((A-0.1)*(B-0.1).LT.0.0) P = 0.1
C     D = P*10.0**(-NSIG)
C     
C     COPYRIGHT           - 1977 BY IMSL, INC. ALL RIGHTS RESERVED.
C     
C     WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C     APPLIED TO THIS CODE. NO OTHER WARRANTY,
C     EXPRESSED OR IMPLIED, IS APPLICABLE.
C     
C-----------------------------------------------------------------------
C     
      SUBROUTINE ZBRENT (F,EPS,NSIG,A,B,MAXFN,IER)
C     SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NSIG,MAXFN,IER
      DOUBLE PRECISION   F,EPS,A,B
C     SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IC
      DOUBLE PRECISION   ZERO,HALF,ONE,THREE,TEN,
     1     T,FA,FB,C,FC,D,E,TOL,RM,S,P,Q,R,RONE,TEMP
      DATA               ZERO/0.0D0/,HALF/.5D0/,ONE/1.0D0/,THREE/3.0D0/,
     1     TEN/10.0D0/
C     FIRST EXECUTABLE STATEMENT
      IER = 0
      T = TEN**(-NSIG)
      IC = 2
      S = A
      FA = F(S)
      S = B
      FB = F(S)
C     TEST FOR SAME SIGN
      IF (FA*FB.GT.ZERO) GO TO 50
    5 C = A
      FC = FA
      D = B-C
      E = D
 10   IF (DABS(FC).GE.DABS(FB)) GO TO 15
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
   15 CONTINUE
      TOL = T*DMAX1(DABS(B),0.1D0)
      RM = (C-B)*HALF
C     TEST FOR FIRST CONVERGENCE CRITERIA
      IF (DABS(FB).LE.EPS) GO TO 40
C     TEST FOR SECOND CONVERGENCE CRITERIA
      IF (DABS(C-B).LE.TOL) GO TO 40
C     CHECK EVALUATION COUNTER
      IF (IC.GE.MAXFN) GO TO 45
C     IS BISECTION FORCED
      IF (DABS(E).LT.TOL) GO TO 30
      IF (DABS(FA).LE.DABS(FB)) GO TO 30
      S = FB/FA
      IF (A.NE.C) GO TO 20
C     LINEAR INTERPOLATION
      P = (C-B)*S
      Q = ONE-S
      GO TO 25
C     INVERSE QUADRATIC INTERPOLATION
 20   Q = FA/FC
      R = FB/FC
      RONE = R-ONE
      P = S*((C-B)*Q*(Q-R)-(B-A)*RONE)
      Q = (Q-ONE)*RONE*(S-ONE)
 25   IF (P.GT.ZERO) Q = -Q
      IF (P.LT.ZERO) P = -P
      S = E
      E = D
C     IF ABS(P/Q).GE.75*ABS(C-B) THEN
C     FORCE BISECTION
      IF (P+P.GE.THREE*RM*Q) GO TO 30
C     IF ABS(P/Q).GE..5*ABS(S) THEN FORCE
C     BISECTION. S = THE VALUE OF P/Q
C     ON THE STEP BEFORE THE LAST ONE
      IF (P+P.GE.DABS(S*Q)) GO TO 30
      D = P/Q
      GO TO 35
C     BISECTION
 30   E = RM
      D = E
C     INCREMENT B
 35   A = B
      FA = FB
      TEMP = D
      IF (DABS(TEMP).LE.HALF*TOL) TEMP = DSIGN(HALF*TOL,RM)
      B = B+TEMP
      S = B
      FB = F(S)
      IC = IC+1
      IF (FB*FC.LE.ZERO) GO TO 10
      GO TO 5
C     CONVERGENCE OF B
 40   A = C
      MAXFN = IC
      GO TO 9005
C     MAXFN EVALUATIONS
 45   IER = 129
      A = C
      MAXFN = IC
      GO TO 9000
C     TERMINAL ERROR - F(A) AND F(B) HAVE
C     THE SAME SIGN
 50   IER = 130
      MAXFN = IC
 9000 CONTINUE
      CALL UERTST (IER,6HZBRENT)
 9005 RETURN
      END
C     IMSL ROUTINE NAME   - UERTST
C     
C-----------------------------------------------------------------------
C     
C     COMPUTER            - VAX/SINGLE
C     
C     LATEST REVISION     - JUNE 1, 1982
C     
C     PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION
C     
C     USAGE               - CALL UERTST (IER,NAME)
C     
C     ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
C     IER = I+J WHERE
C     I = 128 IMPLIES TERMINAL ERROR MESSAGE,
C     I =  64 IMPLIES WARNING WITH FIX MESSAGE,
C     I =  32 IMPLIES WARNING MESSAGE.
C     J = ERROR CODE RELEVANT TO CALLING
C     ROUTINE.
C     NAME   - A CHARACTER STRING OF LENGTH SIX PROVIDING
C     THE NAME OF THE CALLING ROUTINE. (INPUT)
C     
C     PRECISION/HARDWARE  - SINGLE/ALL
C     
C     REQD. IMSL ROUTINES - UGETIO,USPKD
C     
C     NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C     CONVENTIONS IS AVAILABLE IN THE MANUAL
C     INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C     
C     REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
C     TO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
C     NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
C     FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
C     THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
C     UGETIO AS FOLLOWS..
C     NIN = 0
C     NOUT = NEW OUTPUT UNIT NUMBER
C     CALL UGETIO(3,NIN,NOUT)
C     SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
C     
C     COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C     
C     WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C     APPLIED TO THIS CODE. NO OTHER WARRANTY,
C     EXPRESSED OR IMPLIED, IS APPLICABLE.
C     
C-----------------------------------------------------------------------
C     
      SUBROUTINE UERTST (IER,NAME)
C     SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IER
      INTEGER            NAME(1)
C     SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IEQ,IEQDF,IOUNIT,LEVEL,LEVOLD,NAMEQ(6),
     *     NAMSET(6),NAMUPK(6),NIN,NMTB
      DATA               NAMSET/1HU,1HE,1HR,1HS,1HE,1HT/
      DATA               NAMEQ/6*1H /
      DATA               LEVEL/4/,IEQDF/0/,IEQ/1H=/
C     UNPACK NAME INTO NAMUPK
C     FIRST EXECUTABLE STATEMENT
      CALL USPKD (NAME,6,NAMUPK,NMTB)
C     GET OUTPUT UNIT NUMBER
      CALL UGETIO(1,NIN,IOUNIT)
C     CHECK IER
      IF (IER.GT.999) GO TO 25
      IF (IER.LT.-32) GO TO 55
      IF (IER.LE.128) GO TO 5
      IF (LEVEL.LT.1) GO TO 30
C     PRINT TERMINAL MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,35) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,35) IER,NAMUPK
      GO TO 30
    5 IF (IER.LE.64) GO TO 10
      IF (LEVEL.LT.2) GO TO 30
C     PRINT WARNING WITH FIX MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,40) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,40) IER,NAMUPK
      GO TO 30
 10   IF (IER.LE.32) GO TO 15
C     PRINT WARNING MESSAGE
      IF (LEVEL.LT.3) GO TO 30
      IF (IEQDF.EQ.1) WRITE(IOUNIT,45) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,45) IER,NAMUPK
      GO TO 30
   15 CONTINUE
C     CHECK FOR UERSET CALL
      DO 20 I=1,6
         IF (NAMUPK(I).NE.NAMSET(I)) GO TO 25
   20 CONTINUE
      LEVOLD = LEVEL
      LEVEL = IER
      IER = LEVOLD
      IF (LEVEL.LT.0) LEVEL = 4
      IF (LEVEL.GT.4) LEVEL = 4
      GO TO 30
   25 CONTINUE
      IF (LEVEL.LT.4) GO TO 30
C     PRINT NON-DEFINED MESSAGE
      IF (IEQDF.EQ.1) WRITE(IOUNIT,50) IER,NAMEQ,IEQ,NAMUPK
      IF (IEQDF.EQ.0) WRITE(IOUNIT,50) IER,NAMUPK
 30   IEQDF = 0
      RETURN
 35   FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,
     1     20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
 40   FORMAT(27H *** WARNING WITH FIX ERROR,2X,7H(IER = ,I3,
     1     20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
 45   FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,
     1     20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
 50   FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,
     1     20H) FROM IMSL ROUTINE ,6A1,A1,6A1)
C     
C     SAVE P FOR P = R CASE
C     P IS THE PAGE NAMUPK
C     R IS THE ROUTINE NAMUPK
 55   IEQDF = 1
      DO 60 I=1,6
 60      NAMEQ(I) = NAMUPK(I)
 65      RETURN
         END
C     IMSL ROUTINE NAME   - UGETIO
C     
C-----------------------------------------------------------------------
C     
C     COMPUTER            - VAX/SINGLE
C     
C     LATEST REVISION     - JUNE 1, 1981
C     
C     PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
C     VALUES FOR INPUT AND OUTPUT UNIT
C     IDENTIFIERS.
C     
C     USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
C     
C     ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
C     IF IOPT=1, THE CURRENT INPUT AND OUTPUT
C     UNIT IDENTIFIER VALUES ARE RETURNED IN NIN
C     AND NOUT, RESPECTIVELY.
C     IF IOPT=2, THE INTERNAL VALUE OF NIN IS
C     RESET FOR SUBSEQUENT USE.
C     IF IOPT=3, THE INTERNAL VALUE OF NOUT IS
C     RESET FOR SUBSEQUENT USE.
C     NIN    - INPUT UNIT IDENTIFIER.
C     OUTPUT IF IOPT=1, INPUT IF IOPT=2.
C     NOUT   - OUTPUT UNIT IDENTIFIER.
C     OUTPUT IF IOPT=1, INPUT IF IOPT=3.
C     
C     PRECISION/HARDWARE  - SINGLE/ALL
C     
C     REQD. IMSL ROUTINES - NONE REQUIRED
C     
C     NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C     CONVENTIONS IS AVAILABLE IN THE MANUAL
C     INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C     
C     REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
C     OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT
C     IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR
C     IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.
C     SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
C     
C     COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C     
C     WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C     APPLIED TO THIS CODE. NO OTHER WARRANTY,
C     EXPRESSED OR IMPLIED, IS APPLICABLE.
C     
C-----------------------------------------------------------------------
C     
      SUBROUTINE UGETIO(IOPT,NIN,NOUT)
C     SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IOPT,NIN,NOUT
C     SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            NIND,NOUTD
      DATA               NIND/1/,NOUTD/2/
C     FIRST EXECUTABLE STATEMENT
      IF (IOPT.EQ.3) GO TO 10
      IF (IOPT.EQ.2) GO TO 5
      IF (IOPT.NE.1) GO TO 9005
      NIN = NIND
      NOUT = NOUTD
      GO TO 9005
    5 NIND = NIN
      GO TO 9005
 10   NOUTD = NOUT
 9005 RETURN
      END
C     IMSL ROUTINE NAME   - USPKD
C     
C-----------------------------------------------------------------------
C     
C     COMPUTER            - VAX/SINGLE
C     
C     LATEST REVISION     - JUNE 1, 1982
C     
C     PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES THAT HAVE
C     CHARACTER STRING ARGUMENTS
C     
C     USAGE               - CALL USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C     
C     ARGUMENTS    PACKED - CHARACTER STRING TO BE UNPACKED.(INPUT)
C     NCHARS - LENGTH OF PACKED. (INPUT)  SEE REMARKS.
C     UNPAKD - INTEGER ARRAY TO RECEIVE THE UNPACKED
C     REPRESENTATION OF THE STRING. (OUTPUT)
C     NCHMTB - NCHARS MINUS TRAILING BLANKS. (OUTPUT)
C     
C     PRECISION/HARDWARE  - SINGLE/ALL
C     
C     REQD. IMSL ROUTINES - NONE
C     
C     REMARKS  1.  USPKD UNPACKS A CHARACTER STRING INTO AN INTEGER ARRAY
C     IN (A1) FORMAT.
C     2.  UP TO 129 CHARACTERS MAY BE USED.  ANY IN EXCESS OF
C     THAT ARE IGNORED.
C     
C     COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C     
C     WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C     APPLIED TO THIS CODE. NO OTHER WARRANTY,
C     EXPRESSED OR IMPLIED, IS APPLICABLE.
C     
C-----------------------------------------------------------------------
      SUBROUTINE USPKD  (PACKED,NCHARS,UNPAKD,NCHMTB)
C     SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NC,NCHARS,NCHMTB
C     
      INTEGER            UNPAKD(1),IBLANK
      INTEGER            PACKED(1)
      DATA               IBLANK /1H /
C     INITIALIZE NCHMTB
      NCHMTB = 0
C     RETURN IF NCHARS IS LE ZERO
      IF(NCHARS.LE.0) RETURN
C     SET NC=NUMBER OF CHARS TO BE DECODED
      NC = MIN0 (129,NCHARS)
C     DECODE (NC,150,PACKED) (UNPAKD(I),I=1,NC)
      WRITE(PACKED(1),150) (UNPAKD(I),I=1,NC)
 150  FORMAT (129A1)
C     CHECK UNPAKD ARRAY AND SET NCHMTB
C     BASED ON TRAILING BLANKS FOUND
      DO 200 N = 1,NC
         NN = NC - N + 1
         IF(UNPAKD(NN) .NE. IBLANK) GO TO 210
  200 CONTINUE
 210  NCHMTB = NN
      RETURN
      END
      
C              PROGRAM MJY01A
C
C       =========================================================
C       Purpose: This program computes the Bessel functions
C                Jn(x) and Yn(x) ( n=0,1 ) and their derivatives
C                using subroutine JY01A
C       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
C       Output:  BJ0 --- J0(x)
C                DJ0 --- J0'(x)
C                BJ1 --- J1(x)
C                DJ1 --- J1'(x)
C                BY0 --- Y0(x)
C                DY0 --- Y0'(x)
C                BY1 --- Y1(x)
C                DY1 --- Y1'(x)
C       Example:
C
C        x       J0(x)        J0'(x)       J1(x)        J1'(x)
C       ---------------------------------------------------------
C        1     .76519769   -.44005059    .44005059    .32514710
C        5    -.17759677    .32757914   -.32757914   -.11208094
C       10    -.24593576   -.04347275    .04347275   -.25028304
C       20     .16702466   -.06683312    .06683312    .16368301
C       30    -.08636798    .11875106   -.11875106   -.08240961
C       40     .00736689   -.12603832    .12603832    .00421593
C       50     .05581233    .09751183   -.09751183    .05776256
C
C        x       Y0(x)        Y0'(x)       Y1(x)        Y1'(x)
C      ---------------------------------------------------------
C        1     .08825696    .78121282   -.78121282    .86946979
C        5    -.30851763   -.14786314    .14786314   -.33809025
C       10     .05567117   -.24901542    .24901542    .03076962
C       20     .06264060    .16551161   -.16551161    .07091618
C       30    -.11729573   -.08442557    .08442557   -.12010992
C       40     .12593642    .00579351   -.00579351    .12608125
C       50    -.09806500    .05679567   -.05679567   -.09692908
C       =========================================================
C
C        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        WRITE(*,*)'Please enter x'
C        READ(*,*)X
C        WRITE(*,20)X
C        WRITE(*,*)'  x          J0(x)          J0''(x)         J1(x)',
C     &            '          J1''(x)'
C        WRITE(*,*)'------------------------------------------',
C     &            '----------------------------'
C        CALL JY01A(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
C        WRITE(*,10)X,BJ0,DJ0,BJ1,DJ1
C        WRITE(*,*)
C        WRITE(*,*)'  x          Y0(x)          Y0''(x)         Y1(x)',
C     &            '          Y1''(x)'
C        WRITE(*,*)'------------------------------------------',
C     &            '----------------------------'
C        WRITE(*,10)X,BY0,DY0,BY1,DY1
C10      FORMAT(1X,F5.1,4E16.8)
C20      FORMAT(3X,'x =',F5.1)
C        END


        SUBROUTINE JY01A(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
C
C       =======================================================
C       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
C                Y1(x), and their derivatives
C       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ò 0 )
C       Output:  BJ0 --- J0(x)
C                DJ0 --- J0'(x)
C                BJ1 --- J1(x)
C                DJ1 --- J1'(x)
C                BY0 --- Y0(x)
C                DY0 --- Y0'(x)
C                BY1 --- Y1(x)
C                DY1 --- Y1'(x)
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(12),B(12),A1(12),B1(12)
        PI=3.141592653589793D0
        RP2=0.63661977236758D0
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.0D0
           DJ0=0.0D0
           DJ1=0.5D0
           BY0=-1.0D+300
           BY1=-1.0D+300
           DY0=1.0D+300
           DY1=1.0D+300
           RETURN
        ENDIF
        IF (X.LE.12.0D0) THEN
           BJ0=1.0D0
           R=1.0D0
           DO 5 K=1,30
              R=-0.25D0*R*X2/(K*K)
              BJ0=BJ0+R
              IF (DABS(R).LT.DABS(BJ0)*1.0D-15) GO TO 10
5          CONTINUE
10         BJ1=1.0D0
           R=1.0D0
           DO 15 K=1,30
              R=-0.25D0*R*X2/(K*(K+1.0D0))
              BJ1=BJ1+R
              IF (DABS(R).LT.DABS(BJ1)*1.0D-15) GO TO 20
15         CONTINUE
20         BJ1=0.5D0*X*BJ1
           EC=DLOG(X/2.0D0)+0.5772156649015329D0
           CS0=0.0D0
           W0=0.0D0
           R0=1.0D0
           DO 25 K=1,30
              W0=W0+1.0D0/K
              R0=-0.25D0*R0/(K*K)*X2
              R=R0*W0
              CS0=CS0+R
              IF (DABS(R).LT.DABS(CS0)*1.0D-15) GO TO 30
25         CONTINUE
30         BY0=RP2*(EC*BJ0-CS0)
           CS1=1.0D0
           W1=0.0D0
           R1=1.0D0
           DO 35 K=1,30
              W1=W1+1.0D0/K
              R1=-0.25D0*R1/(K*(K+1))*X2
              R=R1*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS1=CS1+R
              IF (DABS(R).LT.DABS(CS1)*1.0D-15) GO TO 40
35         CONTINUE
40         BY1=RP2*(EC*BJ1-1.0D0/X-0.25D0*X*CS1)
        ELSE
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01,
     &            -.1100171402692467D+03,.3038090510922384D+04,
     &            -.1188384262567832D+06,.6252951493434797D+07,
     &            -.4259392165047669D+09,.3646840080706556D+11,
     &            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02,
     &             .5513358961220206D+03,-.1825775547429318D+05,
     &             .8328593040162893D+06,-.5006958953198893D+08,
     &             .3836255180230433D+10,-.3649010818849833D+12,
     &             .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01,
     &             .1215978918765359D+03,-.3302272294480852D+04,
     &             .1276412726461746D+06,-.6656367718817688D+07,
     &             .4502786003050393D+09,-.3833857520742790D+11,
     &             .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02,
     &             -.6038440767050702D+03,.1971837591223663D+05,
     &             -.8902978767070678D+06,.5310411010968522D+08,
     &             -.4043620325107754D+10,.3827011346598605D+12,
     &             -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (X.GE.35.0) K0=10
           IF (X.GE.50.0) K0=8
           T1=X-0.25D0*PI
           P0=1.0D0
           Q0=-0.125D0/X
           DO 45 K=1,K0
              P0=P0+A(K)*X**(-2*K)
45            Q0=Q0+B(K)*X**(-2*K-1)
           CU=DSQRT(RP2/X)
           BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
           BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
           T2=X-0.75D0*PI
           P1=1.0D0
           Q1=0.375D0/X
           DO 50 K=1,K0
              P1=P1+A1(K)*X**(-2*K)
50            Q1=Q1+B1(K)*X**(-2*K-1)
           CU=DSQRT(RP2/X)
           BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
           BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
        ENDIF
        DJ0=-BJ1
        DJ1=BJ0-BJ1/X
        DY0=-BY1
        DY1=BY0-BY1/X
        RETURN
        END
        SUBROUTINE CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
C
C       =======================================================
C       Purpose: Compute Bessel functions J0(z), J1(z), Y0(z),
C                Y1(z), and their derivatives for a complex
C                argument
C       Input :  z --- Complex argument
C       Output:  CBJ0 --- J0(z)
C                CDJ0 --- J0'(z)
C                CBJ1 --- J1(z)
C                CDJ1 --- J1'(z)
C                CBY0 --- Y0(z)
C                CDY0 --- Y0'(z)
C                CBY1 --- Y1(z)
C                CDY1 --- Y1'(z)
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION A(12),B(12),A1(12),B1(12)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        RP2=2.0D0/PI
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBJ0=(1.0D0,0.0D0)
           CBJ1=(0.0D0,0.0D0)
           CDJ0=(0.0D0,0.0D0)
           CDJ1=(0.5D0,0.0D0)
           CBY0=-(1.0D300,0.0D0)
           CBY1=-(1.0D300,0.0D0)
           CDY0=(1.0D300,0.0D0)
           CDY1=(1.0D300,0.0D0)
           RETURN
        ENDIF
        IF (DREAL(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CBJ0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,40
              CR=-0.25D0*CR*Z2/(K*K)
              CBJ0=CBJ0+CR
              IF (CDABS(CR).LT.CDABS(CBJ0)*1.0D-15) GO TO 15
10         CONTINUE
15         CBJ1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+1.0D0))
              CBJ1=CBJ1+CR
              IF (CDABS(CR).LT.CDABS(CBJ1)*1.0D-15) GO TO 25
20         CONTINUE
25         CBJ1=0.5D0*Z1*CBJ1
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(0.0D0,0.0D0)
           DO 30 K=1,40
              W0=W0+1.0D0/K
              CR=-0.25D0*CR/(K*K)*Z2
              CP=CR*W0
              CS=CS+CP
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 35
30         CONTINUE
35         CBY0=RP2*(CDLOG(Z1/2.0D0)+EL)*CBJ0-RP2*CS
           W1=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(1.0D0,0.0D0)
           DO 40 K=1,40
              W1=W1+1.0D0/K
              CR=-0.25D0*CR/(K*(K+1))*Z2
              CP=CR*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS=CS+CP
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 45
40         CONTINUE
45         CBY1=RP2*((CDLOG(Z1/2.0D0)+EL)*CBJ1-1.0D0/Z1-.25D0*Z1*CS)
        ELSE
           DATA A/-.703125D-01,.112152099609375D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01,
     &            -.1100171402692467D+03,.3038090510922384D+04,
     &            -.1188384262567832D+06,.6252951493434797D+07,
     &            -.4259392165047669D+09,.3646840080706556D+11,
     &            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .732421875D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02,
     &             .5513358961220206D+03,-.1825775547429318D+05,
     &             .8328593040162893D+06,-.5006958953198893D+08,
     &             .3836255180230433D+10,-.3649010818849833D+12,
     &             .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875D+00,-.144195556640625D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01,
     &             .1215978918765359D+03,-.3302272294480852D+04,
     &             .1276412726461746D+06,-.6656367718817688D+07,
     &             .4502786003050393D+09,-.3833857520742790D+11,
     &             .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02,
     &             -.6038440767050702D+03,.1971837591223663D+05,
     &             -.8902978767070678D+06,.5310411010968522D+08,
     &             -.4043620325107754D+10,.3827011346598605D+12,
     &             -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           CT1=Z1-.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 50 K=1,K0
50            CP0=CP0+A(K)*Z1**(-2*K)
           CQ0=-0.125D0/Z1
           DO 55 K=1,K0
55            CQ0=CQ0+B(K)*Z1**(-2*K-1)
           CU=CDSQRT(RP2/Z1)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CT2=Z1-.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 60 K=1,K0
60            CP1=CP1+A1(K)*Z1**(-2*K)
           CQ1=0.375D0/Z1
           DO 65 K=1,K0
65            CQ1=CQ1+B1(K)*Z1**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
        ENDIF
        IF (DREAL(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBY0=CBY0-2.0D0*CI*CBJ0
           IF (DIMAG(Z).GT.0.0) CBY0=CBY0+2.0D0*CI*CBJ0
           IF (DIMAG(Z).LT.0.0) CBY1=-(CBY1-2.0D0*CI*CBJ1)
           IF (DIMAG(Z).GT.0.0) CBY1=-(CBY1+2.0D0*CI*CBJ1)
           CBJ1=-CBJ1
        ENDIF
        CDJ0=-CBJ1
        CDJ1=CBJ0-1.0D0/Z*CBJ1
        CDY0=-CBY1
        CDY1=CBY0-1.0D0/Z*CBY1
        RETURN
        END

