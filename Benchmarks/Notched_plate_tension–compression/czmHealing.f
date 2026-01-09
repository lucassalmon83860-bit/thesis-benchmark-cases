cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Constitutive Behavior for an Interface model for crack propagation and healing (3D version)
c	 With heat power evolution
c
c    Author: Stephane Lejeunes
c    Institute: LMA-CNRS UMR7031
c    Date: 19/01/24     
c      
c    PROPS={Kplus,Kmoins,Kt,Kf,eta1,eta2,Gc,Am,Ea,n,mu,Eps_crit} 
c
c    Size(PROPS)=13 (en 3D) 9 (en 2D)  Size(DEPVAR)=10+2
c    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     1 JLTYP,TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION FLUX(2), TIME(2), COORDS(3)
      CHARACTER*80 SNAME

	  power_max = 382.e6	  
	  time_change = 72000.
	  IF (KSTEP.LE.1) Then
		FLUX(1) = (power_max*(TIME(2)))/time_change
	  ELSE
		FLUX(1) = power_max
	  ENDIF
c	  WRITE(*,*)"Flux : ",FLUX(1)

      RETURN
      END
	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine uexpan(expan,dexpandt,temp,time,dtime,predef,dpred,
     $     statev,cmname,nstatv,noel)
c                    
      include 'aba_param.inc'
c                    
      character*80 cmname
c                    
      dimension expan(*),dexpandt(*),temp(2),time(2),predef(*),
     $     dpred(*),statev(nstatv)
c                    
         expanmat = 7.5e-06
		 expangonfl = 4.e-8
		 tempcrit = 2000.
		 time_change = 72000.
		 time_stop = 2000000.
		 timemax = 720000
		 IF (time(2).LE.time_change) Then
			expan(1) = expanmat*temp(2)
		 ELSE
			IF (time(2).LE.(time_change+time_stop)) Then
				expan(1) = expanmat*temp(2)+expangonfl*dtime
			ELSE
				expan(1) = expanmat*temp(2)
			ENDIF
		 ENDIF

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE GAPCON(AK,D,FLOWM,TEMP,PREDEF,TIME,CINAME,SLNAME,
     1 MSNAME,COORDS,NOEL,NODE,NPRED,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CINAME,SLNAME,MSNAME
C
      DIMENSION AK(5),D(2), FLOWM(2), TEMP(2), PREDEF(2,*),
     1 TIME(2), COORDS(3)
		 temperature = (TEMP(1)+TEMP(2))/2.
c 		 -------Conductivite------
c		 conductivite  = 0.274
		 conductivite = 1.
c		 conductivite = 1. * PP(1.-(D(1)/5.e-3))
		 IF (D(1).GE.1.e-7) Then
		 	h_cond = conductivite/D(1)
		 ELSE 
			h_cond  = conductivite /1.e-7
		 ENDIF
c		 h_cond = 100 * PP(1-(D(1)/5.e-3))
c 		 -------Radiation------
		 emissivite_MOX = 0.9
		 emissivite_clad = 0.85
		 const_StefanBoltzman = 5.67032e-8
		 IF ((ABS(TEMP(1)-TEMP(2))).GE.1.e-3) Then
		 	h_rad = (const_StefanBoltzman/((1/emissivite_MOX) + (1/emissivite_MOX) - 1))*((ABS(TEMP(1)**4-TEMP(2)**4))/(ABS(TEMP(1)-TEMP(2))))
		 ELSE 
			h_rad  = 0
		 ENDIF
c		 h_rad = (const_StefanBoltzman/((1/emissivite_MOX) + (1/emissivite_MOX) - 1))*((TEMP(1)-297)+(TEMP(2)-297))*((TEMP(1)-297)**2+(TEMP(2)-297)**2)
		 
		
c		 WRITE(*,*)"Time : ",TIME(2) 
c		 WRITE(*,*)"h_cond : ",h_cond 
c		 WRITE(*,*)"h_rad : ",h_rad 
c		 WRITE(*,*)"dT : ",ABS(TEMP(1)-TEMP(2))
		 
		 AK(1) = h_cond + h_rad
		 
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),H1(NTENS,NTENS),
     4 DFGRD0(3,3),DFGRD1(3,3),H3(NTENS,NTENS),H4(NTENS,NTENS),
     5 H2(NTENS,NTENS),H5(NTENS,NTENS),H6(NTENS,NTENS),
     6 H7(NTENS,NTENS),EPSTi(NTENS-1),EPSTi0(NTENS-1),
     7 EPSTe(NTENS-1),EPSTeOld(NTENS-1),SIGMAT(NTENS-1),
     8 EPSTe0(NTENS-1),DN(NTENS-1),Fid(NTENS-1,NTENS-1)
	  
      DKplus=PROPS(1)
      DKminus=PROPS(2)
      DKt=PROPS(3)
      DKf=PROPS(4)
      Dnu=PROPS(11)
      SMALL_K=1.D-8

C     INITIALIZATION OF INTERNAL VARIABLES
      IF(KINC.EQ.1.AND.STATEV(NSTATV).NE.KINC.AND.KSTEP.EQ.1) THEN
        DO I=1,NSTATV
          STATEV(I)=0.0
        ENDDO 
        STATEV(NSTATV)=KINC
      ENDIF
C     GET INITIAL AND CURRENT VALUES
      Heal0=STATEV(3)
      Da0=STATEV(4) 
      Da=STATEV(2)
      Heal=STATEV(1)
      DO I=1,NTENS-1
        DO J=1,NTENS-1
          Fid(I,J)=0.
         ENDDO
      ENDDO    
      DO I=1,NTENS-1      
        Fid(I,I)=1.
        EPSTi(I)=STATEV(5+3*(I-1))
        EPSTi0(I)=STATEV(6+3*(I-1))
        EPSTeOld(I)=STATEV(7+3*(I-1))
      ENDDO
      DO IJ=1,NTENS
       DO KL=1,NTENS  
        DDSDDE(IJ,KL)=0.D0
       ENDDO
      ENDDO
      
C     UPDATE  
      IF(KINC.NE.STATEV(NSTATV)) THEN
           Heal0=Heal
           Da0=Da
           DO I=1,NTENS-1
             EPSTi0(I)=EPSTi(I)
           ENDDO
      END IF
      EPSNplus=0.5*(ABS(STRAN(1)+DSTRAN(1))+STRAN(1)+DSTRAN(1))
      EPSNplus0=0.5*(ABS(STRAN(1))+STRAN(1))    
      EPSNminus=0.5*(STRAN(1)+DSTRAN(1)-ABS(STRAN(1)+DSTRAN(1)))
      EPSNminus0=0.5*(STRAN(1)-ABS(STRAN(1)))
      DO I=1,NTENS-1
        EPSTe(I)=STRAN(I+1)+DSTRAN(I+1)
        EPSTe0(I)=STRAN(I+1)        
      ENDDO  
      IF(Da>=0.9999) then
        DO I=1,NTENS-1
          EPSTeOld(I)=EPSTe(I)
         enddo
      endif     
     
   
      
C     INTEGRATION OF THE EVOLUTIONS (DAMAGE AND CHEMISTRY AND FRICTION) 
      IF((KINC>0 )) THEN
		CALL EULER_INTEGRATOR(EPSNplus,EPSNplus0,EPSNminus,EPSNminus0,
     1   EPSTe,EPSTe0,EPSTi,EPSTi0,Heal0,Heal,Da0,Da,
     2   DTIME,TIME(2),TEMP+DTEMP,TEMP,PROPS,1.D-4,NPROPS,NTENS,EPSTeOld)
      ENDIF
	  
	  IF (Da.NE.Da) Then
		Da = Da0
		Heal = Heal0
		PNEWDT = 0.5
	  ENDIF
      damageFunc=(1.-Da)**2
	  
      dNN=0
      DO I=1,NTENS-1
        dNN=dNN+DKf*(EPSTe(I)-EPSTi(I))*DKf*(EPSTe(I)-EPSTi(I))
      ENDDO  
      dNN=SQRT(dNN)
      Af=Da*dNN+Dnu*DKminus*EPSNminus
      IF(dNN.EQ.0) then
         dNN=1
      ENDIF    
      DO I=1,NTENS-1
        DN(I)=DKf*(EPSTe(I)-EPSTi(I))/dNN          
      ENDDO
      DHeavisideF=0
      if(Af>0) then        
        DHeavisideF=1
      endif
      DHeavisideComp=1
      if(EPSNplus>0) then
        DHeavisideComp=0
      endif	     
      SIGMAN=damageFunc*DKplus*EPSNplus+DKminus*EPSNminus
	  
      DO I=1,NTENS-1
        SIGMAT(I)=damageFunc*DKt*EPSTe(I)+DHeavisideComp*DKf*(EPSTe(I)-EPSTi(I))+SMALL_K*DKt*EPSTe(I)
      ENDDO  
      
      STRESS(1)=SIGMAN
      DO I=1,NTENS-1
        STRESS(I+1)=SIGMAT(I)
      ENDDO  
      IF(EPSNplus.GT.0) THEN
         DDSDDE(1,1)=damageFunc*DKplus+SMALL_K*DKplus
         DO I=1,NTENS-1
           DDSDDE(I+1,I+1)=SMALL_K*DKt
         ENDDO            
         DO I=1,NTENS-1
           DO J=1,NTENS-1
             DDSDDE(I+1,J+1)=DDSDDE(I+1,J+1)+(damageFunc*DKt+DHeavisideComp*DKf)*Fid(I,J)-
     1        DHeavisideComp*Dkf*DTIME/(PROPS(5)*PROPS(3))*(DHeavisideF*Da*(DN(I)*DN(J))+
     2        0.5*(Af+ABS(Af))*((Fid(I,J)/dNN)+DN(I)*DN(J)))
           ENDDO
         ENDDO   
      ELSE 
         DDSDDE(1,1)=DKminus
         DO I=1,NTENS-1
           DDSDDE(I+1,I+1)=SMALL_K*DKt
         ENDDO  
         DO I=1,NTENS-1
           DO J=1,NTENS-1
             DDSDDE(I+1,J+1)=DDSDDE(I+1,J+1)+(damageFunc*DKt+DHeavisideComp*DKf)*Fid(I,J)-
     1        DHeavisideComp*Dkf*DTIME/(PROPS(5)*PROPS(3))*(DHeavisideF*Da*(DN(I)*DN(J))+
     2        0.5*(Af+ABS(Af))*((Fid(I,J)/dNN)+DN(I)*DN(J)))
           ENDDO
         ENDDO 
      ENDIF
C     INTERNAL VARIABLE STORAGE
      STATEV(1)=Heal
      STATEV(2)=Da
      STATEV(3)=Heal0
      STATEV(4)=Da0
      DO I=1,NTENS-1
        STATEV(5+3*(I-1))=EPSTi(I)
        STATEV(6+3*(I-1))=EPSTi0(I)
        STATEV(7+3*(I-1))=EPSTeOld(I)
      ENDDO  
      STATEV(8)=EPSNplus
      STATEV(9)=EPSTe(1)
C	   STATEV(9)=EPSNminus

c	  STATEV(10)= damageFunc*DKplus*EPSNplus+DKminus*EPSNminus
C	  STATEV(11)= PROPS(7)*Da 
C	  STATEV(12)= STATEV(10)+STATEV(11)
C	  STATEV(13)= SIGMAN
      STATEV(NSTATV)=KINC

c     Calcul de l'évolution de l'énergie dissipée
	  Eta1=PROPS(5)
	  Eta2=PROPS(6)
      Gc=PROPS(7)
	  at=PROPS(8)*Exp(-PROPS(9)/(8.34*TEMP))
	  m=PROPS(10)
	  Eps_crit = -1e-15
	  coeff = 1.
c 	  a remplacer
	  
c     Forces thermodynamiques 
	  Ad1=-2*(1.-Da)*((EPSNplus**2.*DKplus/2.) + (coeff*(EPSTe(1)**2.*DKt/2.)))
	  Ad2=2*Da*at*((1-Heal)**(m+1))/(m+1)
	  Dpsi_Dda = Ad1 + Ad2
	  Dpsi_Dh= -(Da**2)*at*(1-Heal)**m
	  
c     lois de comportement	  
	  dd = (1/Eta1)*(PP(-Ad1-(Gc*Da))*HEAV(-Dpsi_Dda)-(Ad2*HEAV(Eps_crit-EPSNminus)*HEAV(Dpsi_Dda)))
	  dh = (1/Eta2)*(-Dpsi_Dh*HEAV(Eps_crit-EPSNminus))
	  
	  dd_2 = (Da-Da0)/DTIME
	  dh_2 = (Heal-Heal0)/DTIME
	  
	  
	  IF (Da0.LE.0.99999) Then
c		Energie dissipée par fissuration
		STATEV(11) = - Dpsi_Dda*dd_2
		STATEV(13) = STATEV(13) + STATEV(11)
C	  	SPD = SPD - Dpsi_Dda*dd
c     	Energie dissipée par guérison
		STATEV(12) =  - Dpsi_Dh*dh_2
		STATEV(10) = STATEV(10) + STATEV(12)
C	  	SCD = SCD - Dpsi_Dh*dh
c     	Energie totale dissipée
c		STATEV(13) =  - Dpsi_Dda*dd - Dpsi_Dh*dh
c	  	SSE = SSE + damageFunc*(EPSNplus**2.*DKplus/2. + coeff*(EPSTe(1)**2.*DKt/2.))+ EPSNminus**2.*DKminus/2.
c     	Energie élastique
c		STATEV(NSTATV)= damageFunc*((EPSNplus**2.*DKplus/2.) + (coeff*(EPSTe(1)**2.*DKt/2.))) + (EPSNminus**2.*DKminus/2.)
	  ELSE
c		Energie dissipée par fissuration
		STATEV(11) = - 0
C	  	SPD = SPD - Dpsi_Dda*dd
c     	Energie dissipée par guérison
		STATEV(12) =  - 0
C	  	SCD = SCD - Dpsi_Dh*dh
c     	Energie totale dissipée
c		STATEV(13) =  - 0 - 0
c	  	SSE = SSE + damageFunc*(EPSNplus**2.*DKplus/2. + coeff*(EPSTe(1)**2.*DKt/2.))+ EPSNminus**2.*DKminus/2.
c     	Energie élastique
c		STATEV(NSTATV)= damageFunc*((EPSNplus**2.*DKplus/2.) + (coeff*(EPSTe(1)**2.*DKt/2.))) + (EPSNminus**2.*DKminus/2.)
	  ENDIF
c     
	  
	  
      RETURN
      CLOSE(99)      
      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C     INTEGRATOR BACKWARD EULER FOR DAMAGE and CHEMISTRY 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE EULER_INTEGRATOR(EPSNplus,EPSNplus0,EPSNminus,EPSNminus0,
     1  EPSTe,EPSTe0,EPSTi,EPSTi0,Heal0,Heal,Da0,Da,DTIME,T0,
     2  Theta,Theta0,PROPS,Dprec,NPROPS,NTENS,EPSTeOld)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION PROPS(NPROPS),EPSTe(NTENS-1),EPSTe0(NTENS-1),EPSTi(NTENS-1),
     1 EPSTi0(NTENS-1),EPSTii(NTENS-1),EPSTet(NTENS-1),EPSTedemi(NTENS-1),
     2 EPSTiS(NTENS-1),EPSTidemi(NTENS-1),EPSTeOld(NTENS-1)

      Tcur=T0
      Ti=T0
      Tfin=T0+DTIME
      Dai=Da0
      Heali=Heal0
      DO I=1,NTENS-1
       EPSTii(I)=EPSTi0(I)
      ENDDO 
      
      hMax=0.5*DTIME
      hMin=DTIME/1000.
      h0=(hMin+hMax)/2.
      h=h0
	  compteur = 0
      DO WHILE(Tfin>Tcur)
10      Tcur=Ti+h
        Tdemi=Ti+h/2.
        IF(Tcur>Tfin) Then
          Tcur=Tfin+1.D-9
          h=Tcur-Ti
          Tdemi=Ti+h/2. 
        ENDIF

        EPSNplust=EPSNplus0+(EPSNplus-EPSNplus0)*(Tcur-T0)/DTIME
        EPSNplusdemi=EPSNplus0+(EPSNplus-EPSNplus0)*(Tdemi-T0)/DTIME
        DO I=1,NTENS-1
          EPSTet(I)=EPSTe0(I)+(EPSTe(I)-EPSTe0(I))*(Tcur-T0)/DTIME-EPSTeOld(I)
          EPSTedemi(I)=EPSTe0(I)+(EPSTe(I)-EPSTe0(I))*(Tdemi-T0)/DTIME-EPSTeOld(I)
        ENDDO  
        EPSNminust=EPSNminus0+(EPSNminus-EPSNminus0)*(Tcur-T0)/DTIME
        EPSNminusdemi=EPSNminus0+(EPSNminus-EPSNminus0)*(Tdemi-T0)/DTIME
        Thetat=Theta0+(Theta-Theta0)*(Tcur-T0)/DTIME
        Thetatdemi=Theta0+(Theta-Theta0)*(Tdemi-T0)/DTIME
        SigmaMinust=PROPS(2)*EPSNminust
        SigmaMinusdemi=PROPS(2)*EPSNminusdemi
        
		
        CALL LOCAL_NEWTON(EPSNplust,EPSNminust,EPSTet,Heali,Heal,Dai,Da,h,Thetat,PROPS,NPROPS,NTENS)
        CALL LOCAL_NEWTON_FRICTION(EPSNplust,EPSTet,EPSTii,EPSTi,Heal,Da,h,Thetat,
     1       SigmaMinust,PROPS,NPROPS,NTENS)
        CALL LOCAL_NEWTON(EPSNplusdemi,EPSNminusdemi,EPSTedemi,Heali,Healdemi,Dai,Dademi,
     1   0.5*h,Thetatdemi,PROPS,NPROPS,NTENS)
        CALL LOCAL_NEWTON_FRICTION(EPSNplusdemi,EPSTedemi,EPSTii,EPSTidemi,Heal,Da,0.5*h,
     1     Thetatdemi,SigmaMinusdemi,PROPS,NPROPS,NTENS)
        CALL LOCAL_NEWTON(EPSNplust,EPSNminust,EPSTet,Healdemi,HealS,Dademi,DaS,0.5*h,
     1   Thetat,PROPS,NPROPS,NTENS)     
        CALL LOCAL_NEWTON_FRICTION(EPSNplust,EPSTet,EPSTidemi,EPSTiS,Heal,Da,0.5*h,
     1     Thetat,SigmaMinust,PROPS,NPROPS,NTENS)
     
        if(NTENS==3) then
           dResi=Sqrt((HealS-Heal)**2.+(DaS-Da)**2.+(EPSTiS(1)-EPSTi(1))**2.+(EPSTiS(2)-EPSTi(2))**2.) 
        else
           dResi=Sqrt((HealS-Heal)**2.+(DaS-Da)**2.+(EPSTiS(1)-EPSTi(1))**2.) 
        endif
        dErr=dResi/((Dprec*MAX(Da,Heal))+Dprec)
        If((dErr>1.0 .and. h>hMin) .and. Tcur<Tfin)  then
           h=Max(h*0.5,hMin)
           h0=h
           goto 10
        else if(dResi.NE.dResi) then
c           WRITE(*,*)"Error in Euler Integration",Da,Heal,Da0,Heal0,dErr,dResi  
c           Da=Da0
           Heal=Heal0
           goto 20  
        else

           Dai=MAX(MIN(Da,1.),0.)
           Heali=MAX(MIN(Heal,1.),0.)
           EPSTii=EPSTi
           Ti=Tcur
           if(dErr<1) then
	     h=(0.8)*(Dprec/(dResi+1.D-9))**(0.3)*h0
             h=MIN(MAX(h,hMin),hMax)
           endif  
           h0=h
        end if   
      ENDDO
20    IF((Da.LT.0).OR.(Da.GT.1).OR.(Heal.LT.0).OR.(Heal.GT.1)) THEN
c         WRITE(*,*)"Error in Euler Integration",Da,Heal
      ENDIF
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccc
C     positive Part
ccccccccccccccccccccccccccccccccccccccccccc
      REAL*8 FUNCTION PP(a)
      REAL*8 a
      PP=(a+ABS(A))*0.5
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccc
C     Heaviside 
ccccccccccccccccccccccccccccccccccccccccccc      
      REAL*8 FUNCTION HEAV(a)
      REAL*8 a
      HEAV=0.
      if(a>=0) HEAV=1.
      RETURN
      END   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     SOLVING THE NON LINEAR LOCAL SYSTEM OF EQUATIONS CHEMISTRY AND DAMAGE
c    PROPS={Kplus,Kmoins,Kt,Kf,eta1,eta2,Gc,Am,Ea,n,mu} 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LOCAL_NEWTON(EPSNplus,EPSNminus,EPSTe,Heal0,Heal,Da0,Da,h,Theta,PROPS,NPROPS,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION PROPS(NPROPS),R(2),Dk(2,2),Dkinv(2,2),DeltaY(2),EPSTe(NTENS-1)
      Eta1=PROPS(5)
      Eta2=PROPS(6)
      Gc=PROPS(7)
      DKplus=PROPS(1)
      DKt=PROPS(3)
      at=PROPS(8)*Exp(-PROPS(9)/(8.34*Theta))
      Dn=PROPS(10)
      Da=Da0
      Heal=Heal0
C	  Eps_crit = PROPS(12)
      Eps_crit = -1e-15
	  coeff = 1.
C     ITERATIONS DE NEWTON    
      dR0=1.0D10
      dNEPSTe=0
      DO I=1,NTENS-1
        dNEPSTe=dNEPSTe+EPSTe(I)*EPSTe(I)
      ENDDO
      DO IT=1,100
      
        Ad1=2*(1.-Da)*(EPSNplus**2.*DKplus/2. + coeff*(DKt/2.*dNEPSTe))-Gc*Da
        Ad2=2*Da*at*((1-Heal)**(Dn+1))/(Dn+1)
C        Ah=(Da**2)*at*(1-Heal)**Dn-(1-Da)**2*(EPSNplus**2.*DKplus + DKt*dNEPSTe)
        Ah=(Da**2)*at*(1-Heal)**Dn
        Dpsi_Dda=-2*(1.-Da)*(EPSNplus**2.*DKplus/2. + coeff*(DKt/2.*dNEPSTe)) + 2*Da*at*((1-Heal)**Dn)/(Dn+1)
        
        DAd1_Dda=-2*(EPSNplus**2.*DKplus/2. + coeff*(DKt/2.*dNEPSTe))-Gc
        DAd2_Dda=2*at*((1-Heal)**(Dn+1))/(Dn+1)
        DAh_Dda=2*Da*at*(1-Heal)**Dn
        
        DAd1_Dh=0
        DAd2_Dh=2*Da*at*(1-Heal)**Dn
        DAh_Dh=-(Da**2)*at*Dn*(1-Heal)**(Dn-1)
        
        R(1)=Da-Da0-(h/Eta1)*(PP(Ad1)*HEAV(-Dpsi_Dda)-Ad2*HEAV(Eps_crit-EPSNminus)*HEAV(Dpsi_Dda))
        R(2)=Heal-Heal0-(h/Eta2)*HEAV(Eps_crit-EPSNminus)*Ah
C        if(EPSNplus==0 .and. Da>0) write(*,*)IT,"Da",Da,"Heal",Heal,"Ah",Ah,"-Dpsi_Dda",-Dpsi_Dda
        dR=SQRT(R(1)**2+R(2)**2)
        IF(dR.LE.1.E-5) GOTO 10  
   
        IF(dR.GE.1.D6*dR0 .OR. dR.NE.dR) THEN
c          WRITE(*,*)"DIVERGENCE NEWTON...",dR,dR0," ",IT," ",h," ",R(1),R(2),Heal,Da," ",NOEL
C          call XIT()
           Da=sqrt(-Da)
          GOTO 10
        ENDIF
        dR0=dR 
C       Correction ...
        Dk(1,1)=1.-(h/Eta1)*(HEAV(Ad1)*HEAV(-Dpsi_Dda)*DAd1_Dda-HEAV(Eps_crit-EPSNminus)*HEAV(Dpsi_Dda)*DAd2_Dda)
        Dk(1,2)=-(h/Eta1)*(-DAd2_Dh*HEAV(Eps_crit-EPSNminus)*HEAV(Dpsi_Dda))
        Dk(2,1)=-(h/Eta2)*(DAh_Dda*HEAV(Eps_crit-EPSNminus))
        Dk(2,2)=1.-(h/Eta2)*(DAh_Dh*HEAV(Eps_crit-EPSNminus))

        CALL INV2(Dk,Dkinv)     
        DO IJ=1,2
         S=0
         DO KL=1,2
           S=S+Dkinv(IJ,KL)*R(KL)
         ENDDO
         DeltaY(IJ)=-S
        ENDDO    
        Da=MIN(1.,MAX(0.,Da+DeltaY(1)))
        Heal=MIN(1.,MAX(0.,Heal+DeltaY(2)))        
        if(Heal+1.e-10<Heal0)then
c          write(*,*)"Pb Convergence Newton:",Da,Da0,Heal,Heal0
          Da=sqrt(-Da)
          goto 10
        endif 
      ENDDO
10    CONTINUE
C      IF(IT>10) WRITE(*,*)"SLOW CONV NEWTON...",dR," IT",IT
      RETURN
      END
      
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     SOLVING THE NON LINEAR LOCAL SYSTEM FOR FRICTION
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LOCAL_NEWTON_FRICTION(EPSNplus,EPSTe,EPSTi0,EPSTi,Heal,Da,h,Theta,
     1         SigmaMinust,PROPS,NPROPS,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION PROPS(NPROPS),EPSTe(NTENS-1),EPSTi0(NTENS-1),EPSTi(NTENS-1),
     1 Af(NTENS-1),R(NTENS-1),Dk(NTENS-1,NTENS-1),FId(NTENS-1,NTENS-1),
     2 Dkinv(NTENS-1,NTENS-1),DeltaY(NTENS-1)
      Eta=PROPS(5)*PROPS(3)
      DKf=PROPS(4)
      Dnu=PROPS(11)
         
      DO I=1,NTENS-1    
       EPSTi(I)=EPSTi0(I)
       DO J=1,NTENS-1
        FId(I,J)=0.
       ENDDO
       FId(I,I)=1.               
      ENDDO 
      
      If(EPSNplus>0) then
        GOTO 10
      endif
C     ITERATIONS DE NEWTON    
      dR0=1.0D10
      DO IT=1,100
        dNAf=0      
        DO I=1,NTENS-1
          Af(I)=Da*Dkf*(EPSTe(I)-EPSTi(I))
          dNAf=dNAf+Af(I)*Af(I)
        ENDDO
        dNAf=SQRT(dNAf)
         
        DO I=1,NTENS-1
         R(I)=EPSTi(I)-EPSTi0(I)-(h/Eta)*(ABS(dNAf+Dnu*SigmaMinust)+dNAf+Dnu*SigmaMinust)*0.5*(EPSTe(I)-EPSTi(I))
         IF(dNAf>0) THEN
           R(I)=R(I)/dNAf
         ENDIF  
        ENDDO 
         
        dR=0
        DO I=1,NTENS-1
          dR=dR+R(I)*R(I)
        ENDDO
        dR=SQRT(dR)
        
        IF(dR.LE.1.E-6) THEN
		  GOTO 10
		ENDIF
        IF(dR.GE.dR0 .OR. dR.NE.dR) THEN
c          WRITE(*,*)"DIVERGENCE NEWTON FRICTION...",dR," ",IT," ",h," ",R(1),R(2)
C          call XIT()
          GOTO 10
        ENDIF
        dR0=dR 
        DHeavisideF=0
        dFcrit=(ABS(dNAf+Dnu*SigmaMinust)+dNAf+Dnu*SigmaMinust)*0.5
        
        if(dNAf+Dnu*SigmaMinust>0) then        
	   DHeavisideF=1
        endif	     
C       Correction ...
        If(dNAf<1D-8) THEN
          dNAf=1
        endif  
        
        DO I=1,NTENS-1
          DO J=1,NTENS-1
           Dk(I,J)=FId(I,J)+DHeavisideF*(h/Eta)*DKf*Da*Af(I)*Af(J)/(dNaf*dNaf)+
     1          DKf*dFcrit*(h/Eta)*(FId(I,J)/dNaf+Af(I)*Af(J)/dNaf**2)
          ENDDO
        ENDDO    
        CALL INV2(Dk,Dkinv) 
        DO IJ=1,NTENS-1
         S=0
         DO KL=1,NTENS-1
           S=S+Dkinv(IJ,KL)*R(KL)
         ENDDO
         DeltaY(IJ)=-S
         EPSTi(IJ)=EPSTi(IJ)+DeltaY(IJ)
        ENDDO    
      ENDDO
10    CONTINUE
C      IF(IT>10) WRITE(*,*)"SLOW CONV NEWTON...",dR," IT",IT
      RETURN
      END      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     INVERSE OF A 2x2 MATRIX
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE INV2(A,Ainv)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(2,2),Ainv(2,2)
      DETA=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF(DETA.EQ.0) WRITE(*,*) 'NON-INVERTIBLE TENSOR'
      Ainv(1,1)=(A(2,2))/DETA  
      Ainv(2,1)=-(A(2,1))/DETA  

      Ainv(1,2)=-(A(1,2))/DETA  
      Ainv(2,2)=(A(1,1))/DETA  
      RETURN
      END
      

      
