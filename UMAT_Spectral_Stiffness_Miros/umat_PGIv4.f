     
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1           RPL,DDSDDT,DRPLDE,DRPLDT,
     2           STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3           NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4           CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      integer getpid, test_id, myid, ierr, numcpu
      
       parameter NCRYS=1000,SUP_IND=94,PI=3.14159265358979323846264,
     # C11=168.0e3,C12=121.4e3,C44=75.4e3 
      integer TensorInd(21,4),indMap(6,2),a,b,c,d,r,counter
      real*8 Phi1(NCRYS),PHI(NCRYS),Phi2(NCRYS),Super_set_ODF(SUP_IND,3)
     #  ,Fo_C_real(21,SUP_IND),Fo_C_imag(21,SUP_IND),stiffnessSE(6,6),
     # Q_c_s(NCRYS,3,3),stiffness(6,6),ELAS(6,6),Q(3,3),Ctensor(3,3,3,3)
     # ,temp,Id(3,3),C_eq(6,6),C_diff(6,6),C_diff_norm,stiffness1(6,6),
     # CtensorSE(3,3,3,3),temp66(6,6),Ctensor_avg(3,3,3,3)
      complex Fo_C(21,SUP_IND)
      integer index_El_r(732),index_El_r_aux(22),index_El_i(711),
     # index_El_i_aux(22),index_ss_12_r(104),index_ss_12_i(76),
     # index_ws_12_r(183),index_ws_12_i(134),index_GD_r(123),
     # index_GD_i(56)
      character*256 outdir
      double precision ROT(3,3),aux33(3,3)
      YM=128000.0D0
      PV=0.33D0     
     
c      call GETOUTDIR(outdir,lenoutdir)


C
C     Elastic tangent modulus
C

 

c      Phi1=[12.0,5.0,33.0]
c      PHI=[333.0,120.0,22.0]
c      Phi2=[145.0,60.0,0.0]

      C44x2=C44*2.0
      
c     Phi1=[4.284783313646079]
c     PHI=[0.898844564777080]
c     Phi2=[0.092502450355699]
      
       open(unit=20,file='txfft16_1000gr_rad.txt'
     #,status='old')
       do i=1,NCRYS
           read(20,*) Phi1(i),PHI(i),Phi2(i)
       enddo
      close(20)

c ##### Reorient Crystals using DROT from ABAQUS

      Do i=1,NCRYS
         ph =Phi1(i)
         th =PHI(i)
         om =Phi2(i)
          call euler(2,ph,th,om,aux33)
          call orient(aux33,DROT)
          call euler(1,ph,th,om,aux33)
         Phi1(i)= ph
         PHI(i)= th
         Phi2(i)= om
      Enddo


c     Phi1=Phi1*180.0/PI
c     PHI=PHI*180.0/PI
c     Phi2=Phi2*180.0/PI
      
      ! loading of super set, coefs and indices
      !----------------
      !load super set
      open(unit=3,file='Super_set_ODF.txt',
     # status='old')
      do i=1,SUP_IND
          read(3,*) Super_set_ODF(i,:)
      enddo
      close(3)
      
      !load Fourier coefs
      !real part
      open(unit=3,file='Fo_C_real.txt',
     # status='old')
      do i=1,21
          read(3,*) Fo_C_real(i,:)
      enddo
      close(3)
      
      !imag part
      open(unit=3,file='Fo_C_imag.txt',
     # status='old')
      do i=1,21
          read(3,*) Fo_C_imag(i,:)
      enddo
      close(3)
      
      !total
       Fo_C=cmplx(Fo_C_real,Fo_C_imag)
c      Fo_C=dconjg(Fo_C)!if commented out, conjugated coefs are used
      
      index_El_r=0
      index_El_r_aux=0
      k=0
      do i=1,21
        index_El_r_aux(i)=k !determines which portion of index_El_r is belonging to which row of Fo_C_real
        do j=1,SUP_IND
          if(abs(Fo_C_real(i,j)).gt.0.0)then
              k=k+1
              index_El_r(k)=j !contains indices of all non-zero elements of Fo_C_real for each row
          endif
        enddo
      enddo
      index_El_r_aux(i)=k
      
      index_El_i=0
      index_El_i_aux=0
      k=0
      do i=1,21
        index_El_i_aux(i)=k !determines which portion of index_El_i is belonging to which row of Fo_C_imag
        do j=1,SUP_IND
          if(abs(Fo_C_imag(i,j)).gt.0.0)then
              k=k+1
              index_El_i(k)=j !contains indices of all non-zero elements of Fo_C_imag for each row
          endif
        enddo
      enddo
      index_El_i_aux(i)=k
      
      !load tensor indices
      open(unit=3,file='TensorInd.txt',
     # status='old')
      do i=1,21
          read(3,*) TensorInd(i,:)
      enddo
      close(3)

      call spectral_el(Phi1,PHI,Phi2,NCRYS,C11,C12,C44x2,
     # Super_set_ODF,TensorInd,index_El_r,index_El_r_aux,
     # index_El_i,index_El_i_aux,Fo_C_real,Fo_C_imag,
     # stiffness)

      
 
       
     
      DDSDDE(1,1)=stiffness(1,1)
      DDSDDE(1,2)=stiffness(1,2)
      DDSDDE(1,3)=stiffness(1,3)
      DDSDDE(1,4)=stiffness(1,4)
      DDSDDE(1,5)=stiffness(1,5)
      DDSDDE(1,6)=stiffness(1,6)
      DDSDDE(2,1)=stiffness(2,1)
      DDSDDE(2,2)=stiffness(2,2)
      DDSDDE(2,3)=stiffness(2,3)
      DDSDDE(2,4)=stiffness(2,4)
      DDSDDE(2,5)=stiffness(2,5)
      DDSDDE(2,6)=stiffness(2,6)
      DDSDDE(3,1)=stiffness(3,1)
      DDSDDE(3,2)=stiffness(3,2)
      DDSDDE(3,3)=stiffness(3,3)
      DDSDDE(3,4)=stiffness(3,4)
      DDSDDE(3,5)=stiffness(3,5)
      DDSDDE(3,6)=stiffness(3,6)
      DDSDDE(4,1)=stiffness(4,1)
      DDSDDE(4,2)=stiffness(4,2)
      DDSDDE(4,3)=stiffness(4,3)
      DDSDDE(4,4)=stiffness(4,4)
      DDSDDE(4,5)=stiffness(4,5)
      DDSDDE(4,6)=stiffness(4,6)
      DDSDDE(5,1)=stiffness(5,1)
      DDSDDE(5,2)=stiffness(5,2)
      DDSDDE(5,3)=stiffness(5,3)
      DDSDDE(5,4)=stiffness(5,4)
      DDSDDE(5,5)=stiffness(5,5)
      DDSDDE(5,6)=stiffness(5,6)
      DDSDDE(6,1)=stiffness(6,1)
      DDSDDE(6,2)=stiffness(6,2)
      DDSDDE(6,3)=stiffness(6,3)
      DDSDDE(6,4)=stiffness(6,4)
      DDSDDE(6,5)=stiffness(6,5)
      DDSDDE(6,6)=stiffness(6,6)
   

C

 

      Do i=1,6
      Do j=1,6
      STRESS(i)=STRESS(i)+DDSDDE(i,j)*DSTRAN(j)
      Enddo
      Enddo 
      f=2
 

C
C     new stress at the end of increment
C 

 
     
      epx=0.0D0
c
      return 
      End SUBROUTINE UMAT
    
C
C**********************************************************************
C********************************************************************** 
C
      
      
      subroutine spectral_el(Phi1_in,PHI_in,Phi2_in,NCRYS,C11,C12,C44x2,
     # Super_set_ODF,TensorInd,index_El_r,index_El_r_aux,
     # index_El_i,index_El_i_aux,Fo_C_real,Fo_C_imag,C)
      !--------------
      !Input:
      !   Phi1,PHI,Phi2 - euler angles in degrees
      !   NCRYS - number of crystals
      !   Super_set_ODF - super set for ODF transform
      !   Fo_C - coefs for stiffness
      !   TensorInd - indices of tensor
      !Output:
      !   C - average stiffness in matrix representation
      !-----------------
      parameter SUP_IND=94,CI=(0.0,1.0),PI=3.141592653589793238462643
      integer NCRYS,b1(NCRYS),b2(NCRYS),b3(NCRYS),frequency(NCRYS),
     # TensorInd(21,4),indMap(6,2)
     # ,signInd(8,4)!test
      real*8 Phi1(NCRYS),PHI(NCRYS),Phi2(NCRYS),tempSuper(SUP_IND)
     # ,Super_set_ODF(SUP_IND,3),Ctensor(3,3,3,3),C(6,6),Id(3,3)
     # ,CtensorS(3,3,3,3),Phi1_in(NCRYS),PHI_in(NCRYS),Phi2_in(NCRYS),
     # C11,C12,C44,C44x2,temp66(6,6)
      complex FTerm(SUP_IND)!,Fo_C(21,SUP_IND)!(eff3)
      !(eff3)
      integer index_El_r(732),index_El_r_aux(22),index_El_i(711),
     # index_El_i_aux(22)
      real*8 Fo_C_real(21,SUP_IND),Fo_C_imag(21,SUP_IND)
      real*8 sum,FTerm_G(NCRYS,SUP_IND ),tempSuper_G(NCRYS,SUP_IND ),
     #Super_set_ODF_G(NCRYS,SUP_IND,3 ) 
       integer ii
      complex
     #Fterm1,
     #Fterm2,
     #Fterm3,
     #Fterm4,
     #Fterm5,
     #Fterm6,
     #Fterm7,
     #Fterm8,
     #Fterm9,
     #Fterm10,
     #Fterm11,
     #Fterm12,
     #Fterm13,
     #Fterm14,
     #Fterm15,
     #Fterm16,
     #Fterm17,
     #Fterm18,
     #Fterm19,
     #Fterm20,
     #Fterm21,
     #Fterm22,
     #Fterm23,
     #Fterm24,
     #Fterm25,
     #Fterm26,
     #Fterm27,
     #Fterm28,
     #Fterm29,
     #Fterm30,
     #Fterm31,
     #Fterm32,
     #Fterm33,
     #Fterm34,
     #Fterm35,
     #Fterm36,
     #Fterm37,
     #Fterm38,
     #Fterm39,
     #Fterm40,
     #Fterm41,
     #Fterm42,
     #Fterm43,
     #Fterm44,
     #Fterm45,
     #Fterm46,
     #Fterm47,
     #Fterm48,
     #Fterm49,
     #Fterm50,
     #Fterm51,
     #Fterm52,
     #Fterm53,
     #Fterm54,
     #Fterm55,
     #Fterm56,
     #Fterm57,
     #Fterm58,
     #Fterm59,
     #Fterm60,
     #Fterm61,
     #Fterm62,
     #Fterm63,
     #Fterm64,
     #Fterm65,
     #Fterm66,
     #Fterm67,
     #Fterm68,
     #Fterm69,
     #Fterm70,
     #Fterm71,
     #Fterm72,
     #Fterm73,
     #Fterm74,
     #Fterm75,
     #Fterm76,
     #Fterm77,
     #Fterm78,
     #Fterm79,
     #Fterm80,
     #Fterm81,
     #Fterm82,
     #Fterm83,
     #Fterm84,
     #Fterm85,
     #Fterm86,
     #Fterm87,
     #Fterm88,
     #Fterm89,
     #Fterm90,
     #Fterm91,
     #Fterm92,
     #Fterm93,
     #Fterm94 
      
      real*8 
     #tempSuper1,
     #tempSuper2,
     #tempSuper3,
     #tempSuper4,
     #tempSuper5,
     #tempSuper6,
     #tempSuper7,
     #tempSuper8,
     #tempSuper9,
     #tempSuper10,
     #tempSuper11,
     #tempSuper12,
     #tempSuper13,
     #tempSuper14,
     #tempSuper15,
     #tempSuper16,
     #tempSuper17,
     #tempSuper18,
     #tempSuper19,
     #tempSuper20,
     #tempSuper21,
     #tempSuper22,
     #tempSuper23,
     #tempSuper24,
     #tempSuper25,
     #tempSuper26,
     #tempSuper27,
     #tempSuper28,
     #tempSuper29,
     #tempSuper30,
     #tempSuper31,
     #tempSuper32,
     #tempSuper33,
     #tempSuper34,
     #tempSuper35,
     #tempSuper36,
     #tempSuper37,
     #tempSuper38,
     #tempSuper39,
     #tempSuper40,
     #tempSuper41,
     #tempSuper42,
     #tempSuper43,
     #tempSuper44,
     #tempSuper45,
     #tempSuper46,
     #tempSuper47,
     #tempSuper48,
     #tempSuper49,
     #tempSuper50,
     #tempSuper51,
     #tempSuper52,
     #tempSuper53,
     #tempSuper54,
     #tempSuper55,
     #tempSuper56,
     #tempSuper57,
     #tempSuper58,
     #tempSuper59,
     #tempSuper60,
     #tempSuper61,
     #tempSuper62,
     #tempSuper63,
     #tempSuper64,
     #tempSuper65,
     #tempSuper66,
     #tempSuper67,
     #tempSuper68,
     #tempSuper69,
     #tempSuper70,
     #tempSuper71,
     #tempSuper72,
     #tempSuper73,
     #tempSuper74,
     #tempSuper75,
     #tempSuper76,
     #tempSuper77,
     #tempSuper78,
     #tempSuper79,
     #tempSuper80,
     #tempSuper81,
     #tempSuper82,
     #tempSuper83,
     #tempSuper84,
     #tempSuper85,
     #tempSuper86,
     #tempSuper87,
     #tempSuper88,
     #tempSuper89,
     #tempSuper90,
     #tempSuper91,
     #tempSuper92,
     #tempSuper93,
     #tempSuper94 
     
      !----------------
      C44=C44x2*0.5
      
      !identity
      Id=0.0
      do i=1,3
          do j=1,3
              if(i.eq.j)Id(i,j)=1.0
          enddo
      enddo

      !euler angles to degrees
c      Phi1=Phi1_in*180.0/PI
c      PHI=PHI_in*180.0/PI
c      Phi2=Phi2_in*180.0/PI   
      Phi1=Phi1_in
      PHI=PHI_in
      Phi2=Phi2_in

      !to grid points
      b1=NINT(Phi1)
      b2=NINT(PHI)
      b3=NINT(Phi2)
      
      !frequencies of orientations
      call genFrequencyCount(b1,b2,b3,NCRYS,frequency) 
      
      !calculate ODF Fourier coefs
      FTerm=0.0
      FTerm_G=0.0
      sum=0.0
      Fterm1=FTerm(1)
      Fterm2=FTerm(2)
      Fterm3=FTerm(3)
      Fterm4=FTerm(4)
      Fterm5=FTerm(5)
      Fterm6=FTerm(6)
      Fterm7=FTerm(7)
      Fterm8=FTerm(8)
      Fterm9=FTerm(9)
      Fterm10=FTerm(10)
      Fterm11=FTerm(11)
      Fterm12=FTerm(12)
      Fterm13=FTerm(13)
      Fterm14=FTerm(14)
      Fterm15=FTerm(15)
      Fterm16=FTerm(16)
      Fterm17=FTerm(17)
      Fterm18=FTerm(18)
      Fterm19=FTerm(19)
      Fterm20=FTerm(20)
      Fterm21=FTerm(21)
      Fterm22=FTerm(22)
      Fterm23=FTerm(23)
      Fterm24=FTerm(24)
      Fterm25=FTerm(25)
      Fterm26=FTerm(26)
      Fterm27=FTerm(27)
      Fterm28=FTerm(28)
      Fterm29=FTerm(29)
      Fterm30=FTerm(30)
      Fterm31=FTerm(31)
      Fterm32=FTerm(32)
      Fterm33=FTerm(33)
      Fterm34=FTerm(34)
      Fterm35=FTerm(35)
      Fterm36=FTerm(36)
      Fterm37=FTerm(37)
      Fterm38=FTerm(38)
      Fterm39=FTerm(39)
      Fterm40=FTerm(40)
      Fterm41=FTerm(41)
      Fterm42=FTerm(42)
      Fterm43=FTerm(43)
      Fterm44=FTerm(44)
      Fterm45=FTerm(45)
      Fterm46=FTerm(46)
      Fterm47=FTerm(47)
      Fterm48=FTerm(48)
      Fterm49=FTerm(49)
      Fterm50=FTerm(50)
      Fterm51=FTerm(51)
      Fterm52=FTerm(52)
      Fterm53=FTerm(53)
      Fterm54=FTerm(54)
      Fterm55=FTerm(55)
      Fterm56=FTerm(56)
      Fterm57=FTerm(57)
      Fterm58=FTerm(58)
      Fterm59=FTerm(59)
      Fterm60=FTerm(60)
      Fterm61=FTerm(61)
      Fterm62=FTerm(62)
      Fterm63=FTerm(63)
      Fterm64=FTerm(64)
      Fterm65=FTerm(65)
      Fterm66=FTerm(66)
      Fterm67=FTerm(67)
      Fterm68=FTerm(68)
      Fterm69=FTerm(69)
      Fterm70=FTerm(70)
      Fterm71=FTerm(71)
      Fterm72=FTerm(72)
      Fterm73=FTerm(73)
      Fterm74=FTerm(74)
      Fterm75=FTerm(75)
      Fterm76=FTerm(76)
      Fterm77=FTerm(77)
      Fterm78=FTerm(78)
      Fterm79=FTerm(79)
      Fterm80=FTerm(80)
      Fterm81=FTerm(81)
      Fterm82=FTerm(82)
      Fterm83=FTerm(83)
      Fterm84=FTerm(84)
      Fterm85=FTerm(85)
      Fterm86=FTerm(86)
      Fterm87=FTerm(87)
      Fterm88=FTerm(88)
      Fterm89=FTerm(89)
      Fterm90=FTerm(90)
      Fterm91=FTerm(91)
      Fterm92=FTerm(92)
      Fterm93=FTerm(93)
      Fterm94=FTerm(94)



             
#ifdef ScalarVersion
 




#ifdef OMP 
!$omp parallel do  
!$omp& private(tempSuper1  )
!$omp& private(tempSuper2  )
!$omp& private(tempSuper3  )
!$omp& private(tempSuper4  )
!$omp& private(tempSuper5  )
!$omp& private(tempSuper6  )
!$omp& private(tempSuper7  )
!$omp& private(tempSu per8  )
!$omp& private(tempSuper9  )
!$omp& private(tempSuper10 )
!$omp& private(tempSuper11 )
!$omp& private(tempSuper12 )
!$omp& private(tempSuper13 )
!$omp& private(tempSuper14 )
!$omp& private(tempSuper15 )
!$omp& private(tempSuper16 )
!$omp& private(tempSuper17 )
!$omp& private(tempSuper18 )
!$omp& private(tempSuper19 )
!$omp& private(tempSuper20 )
!$omp& private(tempSuper21 )
!$omp& private(tempSuper22 )
!$omp& private(tempSuper23 )
!$omp& private(tempSuper24 )
!$omp& private(tempSuper25 )
!$omp& private(tempSuper26 )
!$omp& private(tempSuper27 )
!$omp& private(tempSuper28 )
!$omp& private(tempSuper29 )
!$omp& private(tempSuper30 )
!$omp& private(tempSuper31 )
!$omp& private(tempSuper32 )
!$omp& private(tempSuper33 )
!$omp& private(tempSuper34 )
!$omp& private(tempSuper35 )
!$omp& private(tempSuper36 )
!$omp& private(tempSuper37 )
!$omp& private(tempSuper38 )
!$omp& private(tempSuper39 )
!$omp& private(tempSuper40 )
!$omp& private(tempSuper41 )
!$omp& private(tempSuper42 )
!$omp& private(tempSuper43 )
!$omp& private(tempSuper44 )
!$omp& private(tempSuper45 )
!$omp& private(tempSuper46 )
!$omp& private(tempSuper47 )
!$omp& private(tempSuper48 )
!$omp& private(tempSuper49 )
!$omp& private(tempSuper50 )
!$omp& private(tempSuper51 )
!$omp& private(tempSuper52 )
!$omp& private(tempSuper53 )
!$omp& private(tempSuper54 )
!$omp& private(tempSuper55 )
!$omp& private(tempSuper56 )
!$omp& private(tempSuper57 )
!$omp& private(tempSuper58 )
!$omp& private(tempSuper59 )
!$omp& private(tempSuper60 )
!$omp& private(tempSuper61 )
!$omp& private(tempSuper62 )
!$omp& private(tempSuper63 )
!$omp& private(tempSuper64 )
!$omp& private(tempSuper65 )
!$omp& private(tempSuper66 )
!$omp& private(tempSuper67 )
!$omp& private(tempSuper68 )
!$omp& private(tempSuper69 )
!$omp& private(tempSuper70 )
!$omp& private(tempSuper71 )
!$omp& private(tempSuper72 )
!$omp& private(tempSuper73 )
!$omp& private(tempSuper74 )
!$omp& private(tempSuper75 )
!$omp& private(tempSuper76 )
!$omp& private(tempSuper77 )
!$omp& private(tempSuper78 )
!$omp& private(tempSuper79 )
!$omp& private(tempSuper80 )
!$omp& private(tempSuper81 )
!$omp& private(tempSuper82 )
!$omp& private(tempSuper83 )
!$omp& private(tempSuper84 )
!$omp& private(tempSuper85 )
!$omp& private(tempSuper86 )
!$omp& private(tempSuper87 )
!$omp& private(tempSuper88 )
!$omp& private(tempSuper89 )
!$omp& private(tempSuper90 )
!$omp& private(tempSuper91 )
!$omp& private(tempSuper92 )
!$omp& private(tempSuper93 )
!$omp& private(tempSuper94 )


!$omp& reduction(+:Fterm1 )
!$omp& reduction(+:Fterm2 )
!$omp& reduction(+:Fterm3 )
!$omp& reduction(+:Fterm4 )
!$omp& reduction(+:Fterm5 )
!$omp& reduction(+:Fterm6 )
!$omp& reduction(+:Fterm7 )
!$omp& reduction(+:Fterm8 )
!$omp& reduction(+:Fterm9 )
!$omp& reduction(+:Fterm10 )
!$omp& reduction(+:Fterm11 )
!$omp& reduction(+:Fterm12 )
!$omp& reduction(+:Fterm13 )
!$omp& reduction(+:Fterm14 )
!$omp& reduction(+:Fterm15 )
!$omp& reduction(+:Fterm16 )
!$omp& reduction(+:Fterm17 )
!$omp& reduction(+:Fterm18 )
!$omp& reduction(+:Fterm19 )
!$omp& reduction(+:Fterm20 )
!$omp& reduction(+:Fterm21 )
!$omp& reduction(+:Fterm22 )
!$omp& reduction(+:Fterm23 )
!$omp& reduction(+:Fterm24 )
!$omp& reduction(+:Fterm25 )
!$omp& reduction(+:Fterm26 )
!$omp& reduction(+:Fterm27 )
!$omp& reduction(+:Fterm28 )
!$omp& reduction(+:Fterm29 )
!$omp& reduction(+:Fterm30 )
!$omp& reduction(+:Fterm31 )
!$omp& reduction(+:Fterm32 )
!$omp& reduction(+:Fterm33 )
!$omp& reduction(+:Fterm34 )
!$omp& reduction(+:Fterm35 )
!$omp& reduction(+:Fterm36 )
!$omp& reduction(+:Fterm37 )
!$omp& reduction(+:Fterm38 )
!$omp& reduction(+:Fterm39 )
!$omp& reduction(+:Fterm40 )
!$omp& reduction(+:Fterm41 )
!$omp& reduction(+:Fterm42 )
!$omp& reduction(+:Fterm43 )
!$omp& reduction(+:Fterm44 )
!$omp& reduction(+:Fterm45 )
!$omp& reduction(+:Fterm46 )
!$omp& reduction(+:Fterm47 )
!$omp& reduction(+:Fterm48 )
!$omp& reduction(+:Fterm49 )
!$omp& reduction(+:Fterm50 )
!$omp& reduction(+:Fterm51 )
!$omp& reduction(+:Fterm52 )
!$omp& reduction(+:Fterm53 )
!$omp& reduction(+:Fterm54 )
!$omp& reduction(+:Fterm55 )
!$omp& reduction(+:Fterm56 )
!$omp& reduction(+:Fterm57 )
!$omp& reduction(+:Fterm58 )
!$omp& reduction(+:Fterm59 )
!$omp& reduction(+:Fterm60 )
!$omp& reduction(+:Fterm61 )
!$omp& reduction(+:Fterm62 )
!$omp& reduction(+:Fterm63 )
!$omp& reduction(+:Fterm64 )
!$omp& reduction(+:Fterm65 )
!$omp& reduction(+:Fterm66 )
!$omp& reduction(+:Fterm67 )
!$omp& reduction(+:Fterm68 )
!$omp& reduction(+:Fterm69 )
!$omp& reduction(+:Fterm70 )
!$omp& reduction(+:Fterm71 )
!$omp& reduction(+:Fterm72 )
!$omp& reduction(+:Fterm73 )
!$omp& reduction(+:Fterm74 )
!$omp& reduction(+:Fterm75 )
!$omp& reduction(+:Fterm76 )
!$omp& reduction(+:Fterm77 )
!$omp& reduction(+:Fterm78 )
!$omp& reduction(+:Fterm79 )
!$omp& reduction(+:Fterm80 )
!$omp& reduction(+:Fterm81 )
!$omp& reduction(+:Fterm82 )
!$omp& reduction(+:Fterm83 )
!$omp& reduction(+:Fterm84 )
!$omp& reduction(+:Fterm85 )
!$omp& reduction(+:Fterm86 )
!$omp& reduction(+:Fterm87 )
!$omp& reduction(+:Fterm88 )
!$omp& reduction(+:Fterm89 )
!$omp& reduction(+:Fterm90 )
!$omp& reduction(+:Fterm91 )
!$omp& reduction(+:Fterm92 )
!$omp& reduction(+:Fterm93 )
!$omp& reduction(+:Fterm94 )
                     
#elif GPU
!$ACC kernels

!$ACC loop 

!$ACC& private(tempSuper1 )
!$ACC& private(tempSuper2 )
!$ACC& private(tempSuper3 )
!$ACC& private(tempSuper4 )
!$ACC& private(tempSuper5 )
!$ACC& private(tempSuper6 )
!$ACC& private(tempSuper7 )
!$ACC& private(tempSuper8 )
!$ACC& private(tempSuper9 )
!$ACC& private(tempSuper10)
!$ACC& private(tempSuper11)
!$ACC& private(tempSuper12)
!$ACC& private(tempSuper13)
!$ACC& private(tempSuper14)
!$ACC& private(tempSuper15)
!$ACC& private(tempSuper16)
!$ACC& private(tempSuper17)
!$ACC& private(tempSuper18)
!$ACC& private(tempSuper19)
!$ACC& private(tempSuper20)
!$ACC& private(tempSuper21)
!$ACC& private(tempSuper22)
!$ACC& private(tempSuper23)
!$ACC& private(tempSuper24)
!$ACC& private(tempSuper25)
!$ACC& private(tempSuper26)
!$ACC& private(tempSuper27)
!$ACC& private(tempSuper28)
!$ACC& private(tempSuper29)
!$ACC& private(tempSuper30)
!$ACC& private(tempSuper31)
!$ACC& private(tempSuper32)
!$ACC& private(tempSuper33)
!$ACC& private(tempSuper34)
!$ACC& private(tempSuper35)
!$ACC& private(tempSuper36)
!$ACC& private(tempSuper37)
!$ACC& private(tempSuper38)
!$ACC& private(tempSuper39)
!$ACC& private(tempSuper40)
!$ACC& private(tempSuper41)
!$ACC& private(tempSuper42)
!$ACC& private(tempSuper43)
!$ACC& private(tempSuper44)
!$ACC& private(tempSuper45)
!$ACC& private(tempSuper46)
!$ACC& private(tempSuper47)
!$ACC& private(tempSuper48)
!$ACC& private(tempSuper49)
!$ACC& private(tempSuper50)
!$ACC& private(tempSuper51)
!$ACC& private(tempSuper52)
!$ACC& private(tempSuper53)
!$ACC& private(tempSuper54)
!$ACC& private(tempSuper55)
!$ACC& private(tempSuper56)
!$ACC& private(tempSuper57)
!$ACC& private(tempSuper58)
!$ACC& private(tempSuper59)
!$ACC& private(tempSuper60)
!$ACC& private(tempSuper61)
!$ACC& private(tempSuper62)
!$ACC& private(tempSuper63)
!$ACC& private(tempSuper64)
!$ACC& private(tempSuper65)
!$ACC& private(tempSuper66)
!$ACC& private(tempSuper67)
!$ACC& private(tempSuper68)
!$ACC& private(tempSuper69)
!$ACC& private(tempSuper70)
!$ACC& private(tempSuper71)
!$ACC& private(tempSuper72)
!$ACC& private(tempSuper73)
!$ACC& private(tempSuper74)
!$ACC& private(tempSuper75)
!$ACC& private(tempSuper76)
!$ACC& private(tempSuper77)
!$ACC& private(tempSuper78)
!$ACC& private(tempSuper79)
!$ACC& private(tempSuper80)
!$ACC& private(tempSuper81)
!$ACC& private(tempSuper82)
!$ACC& private(tempSuper83)
!$ACC& private(tempSuper84)
!$ACC& private(tempSuper85)
!$ACC& private(tempSuper86)
!$ACC& private(tempSuper87)
!$ACC& private(tempSuper88)
!$ACC& private(tempSuper89)
!$ACC& private(tempSuper90)
!$ACC& private(tempSuper91)
!$ACC& private(tempSuper92)
!$ACC& private(tempSuper93)
!$ACC& private(tempSuper94)




!$ACC& reduction(+:Fterm1 )
!$ACC& reduction(+:Fterm2 )
!$ACC& reduction(+:Fterm3 )
!$ACC& reduction(+:Fterm4 )
!$ACC& reduction(+:Fterm5 )
!$ACC& reduction(+:Fterm6 )
!$ACC& reduction(+:Fterm7 )
!$ACC& reduction(+:Fterm8 )
!$ACC& reduction(+:Fterm9 )
!$ACC& reduction(+:Fterm10 )
!$ACC& reduction(+:Fterm11 )
!$ACC& reduction(+:Fterm12 )
!$ACC& reduction(+:Fterm13 )
!$ACC& reduction(+:Fterm14 )
!$ACC& reduction(+:Fterm15 )
!$ACC& reduction(+:Fterm16 )
!$ACC& reduction(+:Fterm17 )
!$ACC& reduction(+:Fterm18 )
!$ACC& reduction(+:Fterm19 )
!$ACC& reduction(+:Fterm20 )
!$ACC& reduction(+:Fterm21 )
!$ACC& reduction(+:Fterm22 )
!$ACC& reduction(+:Fterm23 )
!$ACC& reduction(+:Fterm24 )
!$ACC& reduction(+:Fterm25 )
!$ACC& reduction(+:Fterm26 )
!$ACC& reduction(+:Fterm27 )
!$ACC& reduction(+:Fterm28 )
!$ACC& reduction(+:Fterm29 )
!$ACC& reduction(+:Fterm30 )
!$ACC& reduction(+:Fterm31 )
!$ACC& reduction(+:Fterm32 )
!$ACC& reduction(+:Fterm33 )
!$ACC& reduction(+:Fterm34 )
!$ACC& reduction(+:Fterm35 )
!$ACC& reduction(+:Fterm36 )
!$ACC& reduction(+:Fterm37 )
!$ACC& reduction(+:Fterm38 )
!$ACC& reduction(+:Fterm39 )
!$ACC& reduction(+:Fterm40 )
!$ACC& reduction(+:Fterm41 )
!$ACC& reduction(+:Fterm42 )
!$ACC& reduction(+:Fterm43 )
!$ACC& reduction(+:Fterm44 )
!$ACC& reduction(+:Fterm45 )
!$ACC& reduction(+:Fterm46 )
!$ACC& reduction(+:Fterm47 )
!$ACC& reduction(+:Fterm48 )
!$ACC& reduction(+:Fterm49 )
!$ACC& reduction(+:Fterm50 )
!$ACC& reduction(+:Fterm51 )
!$ACC& reduction(+:Fterm52 )
!$ACC& reduction(+:Fterm53 )
!$ACC& reduction(+:Fterm54 )
!$ACC& reduction(+:Fterm55 )
!$ACC& reduction(+:Fterm56 )
!$ACC& reduction(+:Fterm57 )
!$ACC& reduction(+:Fterm58 )
!$ACC& reduction(+:Fterm59 )
!$ACC& reduction(+:Fterm60 )
!$ACC& reduction(+:Fterm61 )
!$ACC& reduction(+:Fterm62 )
!$ACC& reduction(+:Fterm63 )
!$ACC& reduction(+:Fterm64 )
!$ACC& reduction(+:Fterm65 )
!$ACC& reduction(+:Fterm66 )
!$ACC& reduction(+:Fterm67 )
!$ACC& reduction(+:Fterm68 )
!$ACC& reduction(+:Fterm69 )
!$ACC& reduction(+:Fterm70 )
!$ACC& reduction(+:Fterm71 )
!$ACC& reduction(+:Fterm72 )
!$ACC& reduction(+:Fterm73 )
!$ACC& reduction(+:Fterm74 )
!$ACC& reduction(+:Fterm75 )
!$ACC& reduction(+:Fterm76 )
!$ACC& reduction(+:Fterm77 )
!$ACC& reduction(+:Fterm78 )
!$ACC& reduction(+:Fterm79 )
!$ACC& reduction(+:Fterm80 )
!$ACC& reduction(+:Fterm81 )
!$ACC& reduction(+:Fterm82 )
!$ACC& reduction(+:Fterm83 )
!$ACC& reduction(+:Fterm84 )
!$ACC& reduction(+:Fterm85 )
!$ACC& reduction(+:Fterm86 )
!$ACC& reduction(+:Fterm87 )
!$ACC& reduction(+:Fterm88 )
!$ACC& reduction(+:Fterm89 )
!$ACC& reduction(+:Fterm90 )
!$ACC& reduction(+:Fterm91 )
!$ACC& reduction(+:Fterm92 )
!$ACC& reduction(+:Fterm93 )
!$ACC& reduction(+:Fterm94 )
#endif
 

 
      do i=1,NCRYS
    
      tempSuper1=
     #(Super_set_ODF(1,3)-1.0)*b1(i)+
     #(Super_set_ODF(1,2)-1.0)*b2(i)+
     #(Super_set_ODF(1,1)-1.0)*b3(i)
      tempSuper2=
     #(Super_set_ODF(2,3)-1.0)*b1(i)+
     #(Super_set_ODF(2,2)-1.0)*b2(i)+
     #(Super_set_ODF(2,1)-1.0)*b3(i)
      tempSuper3=
     #(Super_set_ODF(3,3)-1.0)*b1(i)+
     #(Super_set_ODF(3,2)-1.0)*b2(i)+
     #(Super_set_ODF(3,1)-1.0)*b3(i)
      tempSuper4=
     #(Super_set_ODF(4,3)-1.0)*b1(i)+
     #(Super_set_ODF(4,2)-1.0)*b2(i)+
     #(Super_set_ODF(4,1)-1.0)*b3(i)
      tempSuper5=
     #(Super_set_ODF(5,3)-1.0)*b1(i)+
     #(Super_set_ODF(5,2)-1.0)*b2(i)+
     #(Super_set_ODF(5,1)-1.0)*b3(i)
      tempSuper6=
     #(Super_set_ODF(6,3)-1.0)*b1(i)+
     #(Super_set_ODF(6,2)-1.0)*b2(i)+
     #(Super_set_ODF(6,1)-1.0)*b3(i)
      tempSuper7=
     #(Super_set_ODF(7,3)-1.0)*b1(i)+
     #(Super_set_ODF(7,2)-1.0)*b2(i)+
     #(Super_set_ODF(7,1)-1.0)*b3(i)
      tempSuper8=
     #(Super_set_ODF(8,3)-1.0)*b1(i)+
     #(Super_set_ODF(8,2)-1.0)*b2(i)+
     #(Super_set_ODF(8,1)-1.0)*b3(i)
      tempSuper9=
     #(Super_set_ODF(9,3)-1.0)*b1(i)+
     #(Super_set_ODF(9,2)-1.0)*b2(i)+
     #(Super_set_ODF(9,1)-1.0)*b3(i)
      tempSuper10=
     #(Super_set_ODF(10,3)-1.0)*b1(i)+
     #(Super_set_ODF(10,2)-1.0)*b2(i)+
     #(Super_set_ODF(10,1)-1.0)*b3(i)
      tempSuper11=
     #(Super_set_ODF(11,3)-1.0)*b1(i)+
     #(Super_set_ODF(11,2)-1.0)*b2(i)+
     #(Super_set_ODF(11,1)-1.0)*b3(i)
      tempSuper12=
     #(Super_set_ODF(12,3)-1.0)*b1(i)+
     #(Super_set_ODF(12,2)-1.0)*b2(i)+
     #(Super_set_ODF(12,1)-1.0)*b3(i)
      tempSuper13=
     #(Super_set_ODF(13,3)-1.0)*b1(i)+
     #(Super_set_ODF(13,2)-1.0)*b2(i)+
     #(Super_set_ODF(13,1)-1.0)*b3(i)
      tempSuper14=
     #(Super_set_ODF(14,3)-1.0)*b1(i)+
     #(Super_set_ODF(14,2)-1.0)*b2(i)+
     #(Super_set_ODF(14,1)-1.0)*b3(i)
      tempSuper15=
     #(Super_set_ODF(15,3)-1.0)*b1(i)+
     #(Super_set_ODF(15,2)-1.0)*b2(i)+
     #(Super_set_ODF(15,1)-1.0)*b3(i)
      tempSuper16=
     #(Super_set_ODF(16,3)-1.0)*b1(i)+
     #(Super_set_ODF(16,2)-1.0)*b2(i)+
     #(Super_set_ODF(16,1)-1.0)*b3(i)
      tempSuper17=
     #(Super_set_ODF(17,3)-1.0)*b1(i)+
     #(Super_set_ODF(17,2)-1.0)*b2(i)+
     #(Super_set_ODF(17,1)-1.0)*b3(i)
      tempSuper18=
     #(Super_set_ODF(18,3)-1.0)*b1(i)+
     #(Super_set_ODF(18,2)-1.0)*b2(i)+
     #(Super_set_ODF(18,1)-1.0)*b3(i)
      tempSuper19=
     #(Super_set_ODF(19,3)-1.0)*b1(i)+
     #(Super_set_ODF(19,2)-1.0)*b2(i)+
     #(Super_set_ODF(19,1)-1.0)*b3(i)
      tempSuper20=
     #(Super_set_ODF(20,3)-1.0)*b1(i)+
     #(Super_set_ODF(20,2)-1.0)*b2(i)+
     #(Super_set_ODF(20,1)-1.0)*b3(i)
      tempSuper21=
     #(Super_set_ODF(21,3)-1.0)*b1(i)+
     #(Super_set_ODF(21,2)-1.0)*b2(i)+
     #(Super_set_ODF(21,1)-1.0)*b3(i)
      tempSuper22=
     #(Super_set_ODF(22,3)-1.0)*b1(i)+
     #(Super_set_ODF(22,2)-1.0)*b2(i)+
     #(Super_set_ODF(22,1)-1.0)*b3(i)
      tempSuper23=
     #(Super_set_ODF(23,3)-1.0)*b1(i)+
     #(Super_set_ODF(23,2)-1.0)*b2(i)+
     #(Super_set_ODF(23,1)-1.0)*b3(i)
      tempSuper24=
     #(Super_set_ODF(24,3)-1.0)*b1(i)+
     #(Super_set_ODF(24,2)-1.0)*b2(i)+
     #(Super_set_ODF(24,1)-1.0)*b3(i)
      tempSuper25=
     #(Super_set_ODF(25,3)-1.0)*b1(i)+
     #(Super_set_ODF(25,2)-1.0)*b2(i)+
     #(Super_set_ODF(25,1)-1.0)*b3(i)
      tempSuper26=
     #(Super_set_ODF(26,3)-1.0)*b1(i)+
     #(Super_set_ODF(26,2)-1.0)*b2(i)+
     #(Super_set_ODF(26,1)-1.0)*b3(i)
      tempSuper27=
     #(Super_set_ODF(27,3)-1.0)*b1(i)+
     #(Super_set_ODF(27,2)-1.0)*b2(i)+
     #(Super_set_ODF(27,1)-1.0)*b3(i)
      tempSuper28=
     #(Super_set_ODF(28,3)-1.0)*b1(i)+
     #(Super_set_ODF(28,2)-1.0)*b2(i)+
     #(Super_set_ODF(28,1)-1.0)*b3(i)
      tempSuper29=
     #(Super_set_ODF(29,3)-1.0)*b1(i)+
     #(Super_set_ODF(29,2)-1.0)*b2(i)+
     #(Super_set_ODF(29,1)-1.0)*b3(i)
      tempSuper30=
     #(Super_set_ODF(30,3)-1.0)*b1(i)+
     #(Super_set_ODF(30,2)-1.0)*b2(i)+
     #(Super_set_ODF(30,1)-1.0)*b3(i)
      tempSuper31=
     #(Super_set_ODF(31,3)-1.0)*b1(i)+
     #(Super_set_ODF(31,2)-1.0)*b2(i)+
     #(Super_set_ODF(31,1)-1.0)*b3(i)
      tempSuper32=
     #(Super_set_ODF(32,3)-1.0)*b1(i)+
     #(Super_set_ODF(32,2)-1.0)*b2(i)+
     #(Super_set_ODF(32,1)-1.0)*b3(i)
      tempSuper33=
     #(Super_set_ODF(33,3)-1.0)*b1(i)+
     #(Super_set_ODF(33,2)-1.0)*b2(i)+
     #(Super_set_ODF(33,1)-1.0)*b3(i)
      tempSuper34=
     #(Super_set_ODF(34,3)-1.0)*b1(i)+
     #(Super_set_ODF(34,2)-1.0)*b2(i)+
     #(Super_set_ODF(34,1)-1.0)*b3(i)
      tempSuper35=
     #(Super_set_ODF(35,3)-1.0)*b1(i)+
     #(Super_set_ODF(35,2)-1.0)*b2(i)+
     #(Super_set_ODF(35,1)-1.0)*b3(i)
      tempSuper36=
     #(Super_set_ODF(36,3)-1.0)*b1(i)+
     #(Super_set_ODF(36,2)-1.0)*b2(i)+
     #(Super_set_ODF(36,1)-1.0)*b3(i)
      tempSuper37=
     #(Super_set_ODF(37,3)-1.0)*b1(i)+
     #(Super_set_ODF(37,2)-1.0)*b2(i)+
     #(Super_set_ODF(37,1)-1.0)*b3(i)
      tempSuper38=
     #(Super_set_ODF(38,3)-1.0)*b1(i)+
     #(Super_set_ODF(38,2)-1.0)*b2(i)+
     #(Super_set_ODF(38,1)-1.0)*b3(i)
      tempSuper39=
     #(Super_set_ODF(39,3)-1.0)*b1(i)+
     #(Super_set_ODF(39,2)-1.0)*b2(i)+
     #(Super_set_ODF(39,1)-1.0)*b3(i)
      tempSuper40=
     #(Super_set_ODF(40,3)-1.0)*b1(i)+
     #(Super_set_ODF(40,2)-1.0)*b2(i)+
     #(Super_set_ODF(40,1)-1.0)*b3(i)
      tempSuper41=
     #(Super_set_ODF(41,3)-1.0)*b1(i)+
     #(Super_set_ODF(41,2)-1.0)*b2(i)+
     #(Super_set_ODF(41,1)-1.0)*b3(i)
      tempSuper42=
     #(Super_set_ODF(42,3)-1.0)*b1(i)+
     #(Super_set_ODF(42,2)-1.0)*b2(i)+
     #(Super_set_ODF(42,1)-1.0)*b3(i)
      tempSuper43=
     #(Super_set_ODF(43,3)-1.0)*b1(i)+
     #(Super_set_ODF(43,2)-1.0)*b2(i)+
     #(Super_set_ODF(43,1)-1.0)*b3(i)
      tempSuper44=
     #(Super_set_ODF(44,3)-1.0)*b1(i)+
     #(Super_set_ODF(44,2)-1.0)*b2(i)+
     #(Super_set_ODF(44,1)-1.0)*b3(i)
      tempSuper45=
     #(Super_set_ODF(45,3)-1.0)*b1(i)+
     #(Super_set_ODF(45,2)-1.0)*b2(i)+
     #(Super_set_ODF(45,1)-1.0)*b3(i)
      tempSuper46=
     #(Super_set_ODF(46,3)-1.0)*b1(i)+
     #(Super_set_ODF(46,2)-1.0)*b2(i)+
     #(Super_set_ODF(46,1)-1.0)*b3(i)
      tempSuper47=
     #(Super_set_ODF(47,3)-1.0)*b1(i)+
     #(Super_set_ODF(47,2)-1.0)*b2(i)+
     #(Super_set_ODF(47,1)-1.0)*b3(i)
      tempSuper48=
     #(Super_set_ODF(48,3)-1.0)*b1(i)+
     #(Super_set_ODF(48,2)-1.0)*b2(i)+
     #(Super_set_ODF(48,1)-1.0)*b3(i)
      tempSuper49=
     #(Super_set_ODF(49,3)-1.0)*b1(i)+
     #(Super_set_ODF(49,2)-1.0)*b2(i)+
     #(Super_set_ODF(49,1)-1.0)*b3(i)
      tempSuper50=
     #(Super_set_ODF(50,3)-1.0)*b1(i)+
     #(Super_set_ODF(50,2)-1.0)*b2(i)+
     #(Super_set_ODF(50,1)-1.0)*b3(i)
      tempSuper51=
     #(Super_set_ODF(51,3)-1.0)*b1(i)+
     #(Super_set_ODF(51,2)-1.0)*b2(i)+
     #(Super_set_ODF(51,1)-1.0)*b3(i)
      tempSuper52=
     #(Super_set_ODF(52,3)-1.0)*b1(i)+
     #(Super_set_ODF(52,2)-1.0)*b2(i)+
     #(Super_set_ODF(52,1)-1.0)*b3(i)
      tempSuper53=
     #(Super_set_ODF(53,3)-1.0)*b1(i)+
     #(Super_set_ODF(53,2)-1.0)*b2(i)+
     #(Super_set_ODF(53,1)-1.0)*b3(i)
      tempSuper54=
     #(Super_set_ODF(54,3)-1.0)*b1(i)+
     #(Super_set_ODF(54,2)-1.0)*b2(i)+
     #(Super_set_ODF(54,1)-1.0)*b3(i)
      tempSuper55=
     #(Super_set_ODF(55,3)-1.0)*b1(i)+
     #(Super_set_ODF(55,2)-1.0)*b2(i)+
     #(Super_set_ODF(55,1)-1.0)*b3(i)
      tempSuper56=
     #(Super_set_ODF(56,3)-1.0)*b1(i)+
     #(Super_set_ODF(56,2)-1.0)*b2(i)+
     #(Super_set_ODF(56,1)-1.0)*b3(i)
      tempSuper57=
     #(Super_set_ODF(57,3)-1.0)*b1(i)+
     #(Super_set_ODF(57,2)-1.0)*b2(i)+
     #(Super_set_ODF(57,1)-1.0)*b3(i)
      tempSuper58=
     #(Super_set_ODF(58,3)-1.0)*b1(i)+
     #(Super_set_ODF(58,2)-1.0)*b2(i)+
     #(Super_set_ODF(58,1)-1.0)*b3(i)
      tempSuper59=
     #(Super_set_ODF(59,3)-1.0)*b1(i)+
     #(Super_set_ODF(59,2)-1.0)*b2(i)+
     #(Super_set_ODF(59,1)-1.0)*b3(i)
      tempSuper60=
     #(Super_set_ODF(60,3)-1.0)*b1(i)+
     #(Super_set_ODF(60,2)-1.0)*b2(i)+
     #(Super_set_ODF(60,1)-1.0)*b3(i)
      tempSuper61=
     #(Super_set_ODF(61,3)-1.0)*b1(i)+
     #(Super_set_ODF(61,2)-1.0)*b2(i)+
     #(Super_set_ODF(61,1)-1.0)*b3(i)
      tempSuper62=
     #(Super_set_ODF(62,3)-1.0)*b1(i)+
     #(Super_set_ODF(62,2)-1.0)*b2(i)+
     #(Super_set_ODF(62,1)-1.0)*b3(i)
      tempSuper63=
     #(Super_set_ODF(63,3)-1.0)*b1(i)+
     #(Super_set_ODF(63,2)-1.0)*b2(i)+
     #(Super_set_ODF(63,1)-1.0)*b3(i)
      tempSuper64=
     #(Super_set_ODF(64,3)-1.0)*b1(i)+
     #(Super_set_ODF(64,2)-1.0)*b2(i)+
     #(Super_set_ODF(64,1)-1.0)*b3(i)
      tempSuper65=
     #(Super_set_ODF(65,3)-1.0)*b1(i)+
     #(Super_set_ODF(65,2)-1.0)*b2(i)+
     #(Super_set_ODF(65,1)-1.0)*b3(i)
      tempSuper66=
     #(Super_set_ODF(66,3)-1.0)*b1(i)+
     #(Super_set_ODF(66,2)-1.0)*b2(i)+
     #(Super_set_ODF(66,1)-1.0)*b3(i)
      tempSuper67=
     #(Super_set_ODF(67,3)-1.0)*b1(i)+
     #(Super_set_ODF(67,2)-1.0)*b2(i)+
     #(Super_set_ODF(67,1)-1.0)*b3(i)
      tempSuper68=
     #(Super_set_ODF(68,3)-1.0)*b1(i)+
     #(Super_set_ODF(68,2)-1.0)*b2(i)+
     #(Super_set_ODF(68,1)-1.0)*b3(i)
      tempSuper69=
     #(Super_set_ODF(69,3)-1.0)*b1(i)+
     #(Super_set_ODF(69,2)-1.0)*b2(i)+
     #(Super_set_ODF(69,1)-1.0)*b3(i)
      tempSuper70=
     #(Super_set_ODF(70,3)-1.0)*b1(i)+
     #(Super_set_ODF(70,2)-1.0)*b2(i)+
     #(Super_set_ODF(70,1)-1.0)*b3(i)
      tempSuper71=
     #(Super_set_ODF(71,3)-1.0)*b1(i)+
     #(Super_set_ODF(71,2)-1.0)*b2(i)+
     #(Super_set_ODF(71,1)-1.0)*b3(i)
      tempSuper72=
     #(Super_set_ODF(72,3)-1.0)*b1(i)+
     #(Super_set_ODF(72,2)-1.0)*b2(i)+
     #(Super_set_ODF(72,1)-1.0)*b3(i)
      tempSuper73=
     #(Super_set_ODF(73,3)-1.0)*b1(i)+
     #(Super_set_ODF(73,2)-1.0)*b2(i)+
     #(Super_set_ODF(73,1)-1.0)*b3(i)
      tempSuper74=
     #(Super_set_ODF(74,3)-1.0)*b1(i)+
     #(Super_set_ODF(74,2)-1.0)*b2(i)+
     #(Super_set_ODF(74,1)-1.0)*b3(i)
      tempSuper75=
     #(Super_set_ODF(75,3)-1.0)*b1(i)+
     #(Super_set_ODF(75,2)-1.0)*b2(i)+
     #(Super_set_ODF(75,1)-1.0)*b3(i)
      tempSuper76=
     #(Super_set_ODF(76,3)-1.0)*b1(i)+
     #(Super_set_ODF(76,2)-1.0)*b2(i)+
     #(Super_set_ODF(76,1)-1.0)*b3(i)
      tempSuper77=
     #(Super_set_ODF(77,3)-1.0)*b1(i)+
     #(Super_set_ODF(77,2)-1.0)*b2(i)+
     #(Super_set_ODF(77,1)-1.0)*b3(i)
      tempSuper78=
     #(Super_set_ODF(78,3)-1.0)*b1(i)+
     #(Super_set_ODF(78,2)-1.0)*b2(i)+
     #(Super_set_ODF(78,1)-1.0)*b3(i)
      tempSuper79=
     #(Super_set_ODF(79,3)-1.0)*b1(i)+
     #(Super_set_ODF(79,2)-1.0)*b2(i)+
     #(Super_set_ODF(79,1)-1.0)*b3(i)
      tempSuper80=
     #(Super_set_ODF(80,3)-1.0)*b1(i)+
     #(Super_set_ODF(80,2)-1.0)*b2(i)+
     #(Super_set_ODF(80,1)-1.0)*b3(i)
      tempSuper81=
     #(Super_set_ODF(81,3)-1.0)*b1(i)+
     #(Super_set_ODF(81,2)-1.0)*b2(i)+
     #(Super_set_ODF(81,1)-1.0)*b3(i)
      tempSuper82=
     #(Super_set_ODF(82,3)-1.0)*b1(i)+
     #(Super_set_ODF(82,2)-1.0)*b2(i)+
     #(Super_set_ODF(82,1)-1.0)*b3(i)
      tempSuper83=
     #(Super_set_ODF(83,3)-1.0)*b1(i)+
     #(Super_set_ODF(83,2)-1.0)*b2(i)+
     #(Super_set_ODF(83,1)-1.0)*b3(i)
      tempSuper84=
     #(Super_set_ODF(84,3)-1.0)*b1(i)+
     #(Super_set_ODF(84,2)-1.0)*b2(i)+
     #(Super_set_ODF(84,1)-1.0)*b3(i)
      tempSuper85=
     #(Super_set_ODF(85,3)-1.0)*b1(i)+
     #(Super_set_ODF(85,2)-1.0)*b2(i)+
     #(Super_set_ODF(85,1)-1.0)*b3(i)
      tempSuper86=
     #(Super_set_ODF(86,3)-1.0)*b1(i)+
     #(Super_set_ODF(86,2)-1.0)*b2(i)+
     #(Super_set_ODF(86,1)-1.0)*b3(i)
      tempSuper87=
     #(Super_set_ODF(87,3)-1.0)*b1(i)+
     #(Super_set_ODF(87,2)-1.0)*b2(i)+
     #(Super_set_ODF(87,1)-1.0)*b3(i)
      tempSuper88=
     #(Super_set_ODF(88,3)-1.0)*b1(i)+
     #(Super_set_ODF(88,2)-1.0)*b2(i)+
     #(Super_set_ODF(88,1)-1.0)*b3(i)
      tempSuper89=
     #(Super_set_ODF(89,3)-1.0)*b1(i)+
     #(Super_set_ODF(89,2)-1.0)*b2(i)+
     #(Super_set_ODF(89,1)-1.0)*b3(i)
      tempSuper90=
     #(Super_set_ODF(90,3)-1.0)*b1(i)+
     #(Super_set_ODF(90,2)-1.0)*b2(i)+
     #(Super_set_ODF(90,1)-1.0)*b3(i)
      tempSuper91=
     #(Super_set_ODF(91,3)-1.0)*b1(i)+
     #(Super_set_ODF(91,2)-1.0)*b2(i)+
     #(Super_set_ODF(91,1)-1.0)*b3(i)
      tempSuper92=
     #(Super_set_ODF(92,3)-1.0)*b1(i)+
     #(Super_set_ODF(92,2)-1.0)*b2(i)+
     #(Super_set_ODF(92,1)-1.0)*b3(i)
      tempSuper93=
     #(Super_set_ODF(93,3)-1.0)*b1(i)+
     #(Super_set_ODF(93,2)-1.0)*b2(i)+
     #(Super_set_ODF(93,1)-1.0)*b3(i)
      tempSuper94=
     #(Super_set_ODF(94,3)-1.0)*b1(i)+
     #(Super_set_ODF(94,2)-1.0)*b2(i)+
     #(Super_set_ODF(94,1)-1.0)*b3(i)


    
    
    
    
        Fterm1=Fterm1+dfloat(frequency(1))
     #*exp(-2.0*PI*CI/360.0*tempSuper1)
      Fterm2=Fterm2+dfloat(frequency(2))
     #*exp(-2.0*PI*CI/360.0*tempSuper2)
      Fterm3=Fterm3+dfloat(frequency(3))
     #*exp(-2.0*PI*CI/360.0*tempSuper3)
      Fterm4=Fterm4+dfloat(frequency(4))
     #*exp(-2.0*PI*CI/360.0*tempSuper4)
      Fterm5=Fterm5+dfloat(frequency(5))
     #*exp(-2.0*PI*CI/360.0*tempSuper5)
      Fterm6=Fterm6+dfloat(frequency(6))
     #*exp(-2.0*PI*CI/360.0*tempSuper6)
      Fterm7=Fterm7+dfloat(frequency(7))
     #*exp(-2.0*PI*CI/360.0*tempSuper7)
      Fterm8=Fterm8+dfloat(frequency(8))
     #*exp(-2.0*PI*CI/360.0*tempSuper8)
      Fterm9=Fterm9+dfloat(frequency(9))
     #*exp(-2.0*PI*CI/360.0*tempSuper9)
      Fterm10=Fterm10+dfloat(frequency(10))
     #*exp(-2.0*PI*CI/360.0*tempSuper10)
      Fterm11=Fterm11+dfloat(frequency(11))
     #*exp(-2.0*PI*CI/360.0*tempSuper11)
      Fterm12=Fterm12+dfloat(frequency(12))
     #*exp(-2.0*PI*CI/360.0*tempSuper12)
      Fterm13=Fterm13+dfloat(frequency(13))
     #*exp(-2.0*PI*CI/360.0*tempSuper13)
      Fterm14=Fterm14+dfloat(frequency(14))
     #*exp(-2.0*PI*CI/360.0*tempSuper14)
      Fterm15=Fterm15+dfloat(frequency(15))
     #*exp(-2.0*PI*CI/360.0*tempSuper15)
      Fterm16=Fterm16+dfloat(frequency(16))
     #*exp(-2.0*PI*CI/360.0*tempSuper16)
      Fterm17=Fterm17+dfloat(frequency(17))
     #*exp(-2.0*PI*CI/360.0*tempSuper17)
      Fterm18=Fterm18+dfloat(frequency(18))
     #*exp(-2.0*PI*CI/360.0*tempSuper18)
      Fterm19=Fterm19+dfloat(frequency(19))
     #*exp(-2.0*PI*CI/360.0*tempSuper19)
      Fterm20=Fterm20+dfloat(frequency(20))
     #*exp(-2.0*PI*CI/360.0*tempSuper20)
      Fterm21=Fterm21+dfloat(frequency(21))
     #*exp(-2.0*PI*CI/360.0*tempSuper21)
      Fterm22=Fterm22+dfloat(frequency(22))
     #*exp(-2.0*PI*CI/360.0*tempSuper22)
      Fterm23=Fterm23+dfloat(frequency(23))
     #*exp(-2.0*PI*CI/360.0*tempSuper23)
      Fterm24=Fterm24+dfloat(frequency(24))
     #*exp(-2.0*PI*CI/360.0*tempSuper24)
      Fterm25=Fterm25+dfloat(frequency(25))
     #*exp(-2.0*PI*CI/360.0*tempSuper25)
      Fterm26=Fterm26+dfloat(frequency(26))
     #*exp(-2.0*PI*CI/360.0*tempSuper26)
      Fterm27=Fterm27+dfloat(frequency(27))
     #*exp(-2.0*PI*CI/360.0*tempSuper27)
      Fterm28=Fterm28+dfloat(frequency(28))
     #*exp(-2.0*PI*CI/360.0*tempSuper28)
      Fterm29=Fterm29+dfloat(frequency(29))
     #*exp(-2.0*PI*CI/360.0*tempSuper29)
      Fterm30=Fterm30+dfloat(frequency(30))
     #*exp(-2.0*PI*CI/360.0*tempSuper30)
      Fterm31=Fterm31+dfloat(frequency(31))
     #*exp(-2.0*PI*CI/360.0*tempSuper31)
      Fterm32=Fterm32+dfloat(frequency(32))
     #*exp(-2.0*PI*CI/360.0*tempSuper32)
      Fterm33=Fterm33+dfloat(frequency(33))
     #*exp(-2.0*PI*CI/360.0*tempSuper33)
      Fterm34=Fterm34+dfloat(frequency(34))
     #*exp(-2.0*PI*CI/360.0*tempSuper34)
      Fterm35=Fterm35+dfloat(frequency(35))
     #*exp(-2.0*PI*CI/360.0*tempSuper35)
      Fterm36=Fterm36+dfloat(frequency(36))
     #*exp(-2.0*PI*CI/360.0*tempSuper36)
      Fterm37=Fterm37+dfloat(frequency(37))
     #*exp(-2.0*PI*CI/360.0*tempSuper37)
      Fterm38=Fterm38+dfloat(frequency(38))
     #*exp(-2.0*PI*CI/360.0*tempSuper38)
      Fterm39=Fterm39+dfloat(frequency(39))
     #*exp(-2.0*PI*CI/360.0*tempSuper39)
      Fterm40=Fterm40+dfloat(frequency(40))
     #*exp(-2.0*PI*CI/360.0*tempSuper40)
      Fterm41=Fterm41+dfloat(frequency(41))
     #*exp(-2.0*PI*CI/360.0*tempSuper41)
      Fterm42=Fterm42+dfloat(frequency(42))
     #*exp(-2.0*PI*CI/360.0*tempSuper42)
      Fterm43=Fterm43+dfloat(frequency(43))
     #*exp(-2.0*PI*CI/360.0*tempSuper43)
      Fterm44=Fterm44+dfloat(frequency(44))
     #*exp(-2.0*PI*CI/360.0*tempSuper44)
      Fterm45=Fterm45+dfloat(frequency(45))
     #*exp(-2.0*PI*CI/360.0*tempSuper45)
      Fterm46=Fterm46+dfloat(frequency(46))
     #*exp(-2.0*PI*CI/360.0*tempSuper46)
      Fterm47=Fterm47+dfloat(frequency(47))
     #*exp(-2.0*PI*CI/360.0*tempSuper47)
      Fterm48=Fterm48+dfloat(frequency(48))
     #*exp(-2.0*PI*CI/360.0*tempSuper48)
      Fterm49=Fterm49+dfloat(frequency(49))
     #*exp(-2.0*PI*CI/360.0*tempSuper49)
      Fterm50=Fterm50+dfloat(frequency(50))
     #*exp(-2.0*PI*CI/360.0*tempSuper50)
      Fterm51=Fterm51+dfloat(frequency(51))
     #*exp(-2.0*PI*CI/360.0*tempSuper51)
      Fterm52=Fterm52+dfloat(frequency(52))
     #*exp(-2.0*PI*CI/360.0*tempSuper52)
      Fterm53=Fterm53+dfloat(frequency(53))
     #*exp(-2.0*PI*CI/360.0*tempSuper53)
      Fterm54=Fterm54+dfloat(frequency(54))
     #*exp(-2.0*PI*CI/360.0*tempSuper54)
      Fterm55=Fterm55+dfloat(frequency(55))
     #*exp(-2.0*PI*CI/360.0*tempSuper55)
      Fterm56=Fterm56+dfloat(frequency(56))
     #*exp(-2.0*PI*CI/360.0*tempSuper56)
      Fterm57=Fterm57+dfloat(frequency(57))
     #*exp(-2.0*PI*CI/360.0*tempSuper57)
      Fterm58=Fterm58+dfloat(frequency(58))
     #*exp(-2.0*PI*CI/360.0*tempSuper58)
      Fterm59=Fterm59+dfloat(frequency(59))
     #*exp(-2.0*PI*CI/360.0*tempSuper59)
      Fterm60=Fterm60+dfloat(frequency(60))
     #*exp(-2.0*PI*CI/360.0*tempSuper60)
      Fterm61=Fterm61+dfloat(frequency(61))
     #*exp(-2.0*PI*CI/360.0*tempSuper61)
      Fterm62=Fterm62+dfloat(frequency(62))
     #*exp(-2.0*PI*CI/360.0*tempSuper62)
      Fterm63=Fterm63+dfloat(frequency(63))
     #*exp(-2.0*PI*CI/360.0*tempSuper63)
      Fterm64=Fterm64+dfloat(frequency(64))
     #*exp(-2.0*PI*CI/360.0*tempSuper64)
      Fterm65=Fterm65+dfloat(frequency(65))
     #*exp(-2.0*PI*CI/360.0*tempSuper65)
      Fterm66=Fterm66+dfloat(frequency(66))
     #*exp(-2.0*PI*CI/360.0*tempSuper66)
      Fterm67=Fterm67+dfloat(frequency(67))
     #*exp(-2.0*PI*CI/360.0*tempSuper67)
      Fterm68=Fterm68+dfloat(frequency(68))
     #*exp(-2.0*PI*CI/360.0*tempSuper68)
      Fterm69=Fterm69+dfloat(frequency(69))
     #*exp(-2.0*PI*CI/360.0*tempSuper69)
      Fterm70=Fterm70+dfloat(frequency(70))
     #*exp(-2.0*PI*CI/360.0*tempSuper70)
      Fterm71=Fterm71+dfloat(frequency(71))
     #*exp(-2.0*PI*CI/360.0*tempSuper71)
      Fterm72=Fterm72+dfloat(frequency(72))
     #*exp(-2.0*PI*CI/360.0*tempSuper72)
      Fterm73=Fterm73+dfloat(frequency(73))
     #*exp(-2.0*PI*CI/360.0*tempSuper73)
      Fterm74=Fterm74+dfloat(frequency(74))
     #*exp(-2.0*PI*CI/360.0*tempSuper74)
      Fterm75=Fterm75+dfloat(frequency(75))
     #*exp(-2.0*PI*CI/360.0*tempSuper75)
      Fterm76=Fterm76+dfloat(frequency(76))
     #*exp(-2.0*PI*CI/360.0*tempSuper76)
      Fterm77=Fterm77+dfloat(frequency(77))
     #*exp(-2.0*PI*CI/360.0*tempSuper77)
      Fterm78=Fterm78+dfloat(frequency(78))
     #*exp(-2.0*PI*CI/360.0*tempSuper78)
      Fterm79=Fterm79+dfloat(frequency(79))
     #*exp(-2.0*PI*CI/360.0*tempSuper79)
      Fterm80=Fterm80+dfloat(frequency(80))
     #*exp(-2.0*PI*CI/360.0*tempSuper80)
      Fterm81=Fterm81+dfloat(frequency(81))
     #*exp(-2.0*PI*CI/360.0*tempSuper81)
      Fterm82=Fterm82+dfloat(frequency(82))
     #*exp(-2.0*PI*CI/360.0*tempSuper82)
      Fterm83=Fterm83+dfloat(frequency(83))
     #*exp(-2.0*PI*CI/360.0*tempSuper83)
      Fterm84=Fterm84+dfloat(frequency(84))
     #*exp(-2.0*PI*CI/360.0*tempSuper84)
      Fterm85=Fterm85+dfloat(frequency(85))
     #*exp(-2.0*PI*CI/360.0*tempSuper85)
      Fterm86=Fterm86+dfloat(frequency(86))
     #*exp(-2.0*PI*CI/360.0*tempSuper86)
      Fterm87=Fterm87+dfloat(frequency(87))
     #*exp(-2.0*PI*CI/360.0*tempSuper87)
      Fterm88=Fterm88+dfloat(frequency(88))
     #*exp(-2.0*PI*CI/360.0*tempSuper88)
      Fterm89=Fterm89+dfloat(frequency(89))
     #*exp(-2.0*PI*CI/360.0*tempSuper89)
      Fterm90=Fterm90+dfloat(frequency(90))
     #*exp(-2.0*PI*CI/360.0*tempSuper90)
      Fterm91=Fterm91+dfloat(frequency(91))
     #*exp(-2.0*PI*CI/360.0*tempSuper91)
      Fterm92=Fterm92+dfloat(frequency(92))
     #*exp(-2.0*PI*CI/360.0*tempSuper92)
      Fterm93=Fterm93+dfloat(frequency(93))
     #*exp(-2.0*PI*CI/360.0*tempSuper93)
      Fterm94=Fterm94+dfloat(frequency(94))
     #*exp(-2.0*PI*CI/360.0*tempSuper94)
         
         
         
         
    
    
    
    
    
    
    
   
      
      enddo
       
 
#ifdef OMP
!$omp end parallel do
#elif GPU
!$acc end kernels
#endif


      FTerm(1)=Fterm1
      FTerm(2)=Fterm2
      FTerm(3)=Fterm3
      FTerm(4)=Fterm4
      FTerm(5)=Fterm5
      FTerm(6)=Fterm6
      FTerm(7)=Fterm7
      FTerm(8)=Fterm8
      FTerm(9)=Fterm9
      FTerm(10)=Fterm10
      FTerm(11)=Fterm11
      FTerm(12)=Fterm12
      FTerm(13)=Fterm13
      FTerm(14)=Fterm14
      FTerm(15)=Fterm15
      FTerm(16)=Fterm16
      FTerm(17)=Fterm17
      FTerm(18)=Fterm18
      FTerm(19)=Fterm19
      FTerm(20)=Fterm20
      FTerm(21)=Fterm21
      FTerm(22)=Fterm22
      FTerm(23)=Fterm23
      FTerm(24)=Fterm24
      FTerm(25)=Fterm25
      FTerm(26)=Fterm26
      FTerm(27)=Fterm27
      FTerm(28)=Fterm28
      FTerm(29)=Fterm29
      FTerm(30)=Fterm30
      FTerm(31)=Fterm31
      FTerm(32)=Fterm32
      FTerm(33)=Fterm33
      FTerm(34)=Fterm34
      FTerm(35)=Fterm35
      FTerm(36)=Fterm36
      FTerm(37)=Fterm37
      FTerm(38)=Fterm38
      FTerm(39)=Fterm39
      FTerm(40)=Fterm40
      FTerm(41)=Fterm41
      FTerm(42)=Fterm42
      FTerm(43)=Fterm43
      FTerm(44)=Fterm44
      FTerm(45)=Fterm45
      FTerm(46)=Fterm46
      FTerm(47)=Fterm47
      FTerm(48)=Fterm48
      FTerm(49)=Fterm49
      FTerm(50)=Fterm50
      FTerm(51)=Fterm51
      FTerm(52)=Fterm52
      FTerm(53)=Fterm53
      FTerm(54)=Fterm54
      FTerm(55)=Fterm55
      FTerm(56)=Fterm56
      FTerm(57)=Fterm57
      FTerm(58)=Fterm58
      FTerm(59)=Fterm59
      FTerm(60)=Fterm60
      FTerm(61)=Fterm61
      FTerm(62)=Fterm62
      FTerm(63)=Fterm63
      FTerm(64)=Fterm64
      FTerm(65)=Fterm65
      FTerm(66)=Fterm66
      FTerm(67)=Fterm67
      FTerm(68)=Fterm68
      FTerm(69)=Fterm69
      FTerm(70)=Fterm70
      FTerm(71)=Fterm71
      FTerm(72)=Fterm72
      FTerm(73)=Fterm73
      FTerm(74)=Fterm74
      FTerm(75)=Fterm75
      FTerm(76)=Fterm76
      FTerm(77)=Fterm77
      FTerm(78)=Fterm78
      FTerm(79)=Fterm79
      FTerm(80)=Fterm80
      FTerm(81)=Fterm81
      FTerm(82)=Fterm82
      FTerm(83)=Fterm83
      FTerm(84)=Fterm84
      FTerm(85)=Fterm85
      FTerm(86)=Fterm86
      FTerm(87)=Fterm87
      FTerm(88)=Fterm88
      FTerm(89)=Fterm89
      FTerm(90)=Fterm90
      FTerm(91)=Fterm91
      FTerm(92)=Fterm92
      FTerm(93)=Fterm93
      FTerm(94)=Fterm94




#elif ArrayVersion






#ifdef OMP 
!$omp parallel do private(tempSuper) reduction(+:FTerm)
#elif GPU
#if 1
!$ACC kernels
           
#else
!$ACC parallel
#endif
!$ACC loop private(tempSuper,FTerm) 
#endif
     
      do i=1,NCRYS
      do ii=1,SUP_IND 
     
          !Twiddle-Factors
          tempSuper(ii)=(Super_set_ODF(ii,3)-1.0)*b1(i)+
     #     (Super_set_ODF(ii,2)-1.0)*b2(i)+
     # (Super_set_ODF(ii,1)-1.0)*b3(i)
          
          !sumation of ODF Fourier coefs of single orientations
c          if(frequency(i).ne.0)then
           
#if 0
!$ACC atomic
#endif    

          

           FTerm(ii)=FTerm(ii)+dfloat(frequency(i))*
     #exp(-2.0*PI*CI/360.0*tempSuper(ii))
      
c             FTerm(ii)=2.0/NCRYS*FTerm(ii)


#if 0
!$ACC end atomic
#endif
c          endif
      enddo
      enddo
 
#ifdef OMP
!$omp end parallel do
#elif GPU
#if 1
!$acc end kernels
#else
!$acc end parallel
#endif
#endif
  


#endif
  

 
       do ii=1,SUP_IND
       FTerm(ii)=2.0/NCRYS*FTerm(ii) !averaging, ODF Fourier coefs for polycrystal 
       Enddo    
   

    
      !calculate volume average
      do i=1,21
       
        
          
          !(eff3)
c          CtensorS(TensorInd(i,1),TensorInd(i,2),TensorInd(i,3),
c     #     TensorInd(i,4))=1.0/360.0**3*
c     #     real(dot_product(Fo_C(i,:),FTerm))

          !(eff3) sum only non-zero elements, only real parts
          CtensorS(TensorInd(i,1),TensorInd(i,2),TensorInd(i,3),
     #     TensorInd(i,4))=1.0/360.0**3*(
     #     dot_product(
     #     Fo_C_real( !first term dot product
     #     i,index_El_r(index_El_r_aux(i)+1:index_El_r_aux(i+1))
     #     ), 
     #     real(FTerm( !second term dot product
     #     index_El_r(index_El_r_aux(i)+1:index_El_r_aux(i+1)))
     #     )
     #     )
     #     + dot_product(
     #     Fo_C_imag(  !first term dot product
     #     i,index_El_i(index_El_i_aux(i)+1:index_El_i_aux(i+1))
     #     ),
     #     aimag(FTerm( !second term dot product
     #     index_El_i(index_El_i_aux(i)+1:index_El_i_aux(i+1)))
     #     )
     #     )
     #     )

       
      enddo

 

      do i=1,3
          do j=i,3
              do k=1,3
                  do l=k,3
                      Ctensor(i,j,k,l)=C12*Id(i,j)*Id(k,l)
     #                +C44*(Id(i,k)*Id(j,l)+Id(i,l)*Id(j,k))+
     #                (C11-C12-2*C44)*CtensorS(i,j,k,l)
                      
                  enddo
              enddo
          enddo
      enddo
      
      !matrix representation of stiffness
      indMap(:,1)=[1,2,3,2,1,1]
      indMap(:,2)=[1,2,3,3,3,2]
      do i=1,6
          do j=i,6
              C(i,j)=Ctensor(indMap(i,1),indMap(i,2),indMap(j,1),
     #         indMap(j,2))
              C(j,i)=Ctensor(indMap(i,1),indMap(i,2),indMap(j,1),
     #         indMap(j,2))
c              write(*,*)i,j
c              write(*,*)indMap(i,1),indMap(i,2),indMap(j,1),
c     #         indMap(j,2)
c              write(*,*)'*************'
          enddo
      enddo
      
      !to abaqus notation
      temp66=C
      C(4,:)=temp66(6,:)
      C(6,:)=temp66(4,:)
      temp66=C
      C(:,4)=temp66(:,6)
      C(:,6)=temp66(:,4)

      end subroutine spectral_el
C
C**********************************************************************
C********************************************************************** 
C
      subroutine genFrequencyCount(b1,b2,b3,NCRYS,frequency)
      integer NCRYS,b1(NCRYS),b2(NCRYS),b3(NCRYS),frequency(NCRYS)
      
      frequency=1
c#ifdef UGPU
c!$ACC kernels
c!$ACC loop independent
c!$ACC& Private(i,j,frequency)
c#endif      
      do i=1,NCRYS
          do j=i+1,NCRYS
              if(b1(i).eq.b1(j).and.b2(i).eq.b2(j)
     #         .and.b3(i).eq.b3(j))then
                  frequency(i)=frequency(i)+1
                  frequency(j)=0
              endif
          enddo
      enddo
c#ifdef UGPU
c!$ACC end kernels
c#endif
      end subroutine genFrequencyCount
C
C**********************************************************************
C********************************************************************** 
C
      
!*******************************************
!*******************************************      
      
      		SUBROUTINE euler2Q (NCRYS, Phi1, PHI, Phi2, Q_c_s)

!		Input: Euler angles in degrees.
!		Output: 3d array of transformation matrices crystal to sample
		REAL*8 PI
		INTEGER i,NCRYS
		PARAMETER (PI = 3.1415926535897932384626433832795)
		REAL*8 Phi1(NCRYS), PHI(NCRYS), Phi2(NCRYS), Q_c_s(NCRYS,3,3)  
		REAL*8 SPhi1(NCRYS), CPhi1(NCRYS), SPHI(NCRYS), CPHI(NCRYS)      
		REAL*8 SPhi2(NCRYS), CPhi2(NCRYS), to_rad, to_deg

          to_rad = PI / 180.d0!here was error
		to_deg = 1.d0 / to_rad
    
		Phi1 = to_rad * Phi1 
		PHI  = to_rad * PHI
		Phi2 = to_rad * Phi2

		SPhi1 = DSIN(Phi1)
		CPhi1 = DCOS(Phi1)
		SPHI  = DSIN(PHI)
		CPHI  = DCOS(PHI)
		SPhi2 = DSIN(Phi2)
		CPhi2 = DCOS(Phi2)
    
		DO i = 1,NCRYS
        
			Q_c_s(i,1,1) = CPhi1(i)*CPhi2(i)-SPhi1(i)*CPHI(i)*SPhi2(i)
			Q_c_s(i,1,2) = -CPhi1(i)*SPhi2(i)-SPhi1(i)*CPHI(i)*CPhi2(i)
			Q_c_s(i,1,3) = SPhi1(i)*SPHI(i) 
			Q_c_s(i,2,1) = SPhi1(i)*CPhi2(i)+CPhi1(i)*CPHI(i)*SPhi2(i)
			Q_c_s(i,2,2) = -SPhi1(i)*SPhi2(i)+CPhi1(i)*CPHI(i)*CPhi2(i)
			Q_c_s(i,2,3) = -CPhi1(i)*SPHI(i)
			Q_c_s(i,3,1) = SPHI(i)*SPhi2(i)
			Q_c_s(i,3,2) = SPHI(i)*CPhi2(i)
			Q_c_s(i,3,3) = CPHI(i)
    
		ENDDO
    
			Phi1 = to_deg * Phi1 
			PHI  = to_deg * PHI
			Phi2 = to_deg * Phi2

		RETURN
    
		END SUBROUTINE euler2Q



    
      SUBROUTINE ELAST_BCC(Q, C11, C12, C44, ELAS)
	  implicit none
C
C      INCLUDE 'C:\Milan\PhD\Research\1000_model\COMM.f'

        REAL*8 C11, C12, C44, Q(3,3), ELAS(6,6), HC44, HALF


C      INTEGER INTP
C      
      HALF = 0.5D0      
      HC44 = HALF*C44 ! IN INPUT FILE C44 = 2*C44
C
C
      ELAS(1,1) = C12 + C44 + (C11 - C12 - C44)*
     &   (Q(1,1)*Q(1,1)*Q(1,1)*Q(1,1) + Q(1,2)*Q(1,2)*Q(1,2)*Q(1,2) +
     &    Q(1,3)*Q(1,3)*Q(1,3)*Q(1,3))
C
C
      ELAS(2,2) = C12 + C44 + (C11 - C12 - C44)*
     &  (Q(2,1)*Q(2,1)*Q(2,1)*Q(2,1) + Q(2,2)*Q(2,2)*Q(2,2)*Q(2,2) +
     &  Q(2,3)*Q(2,3)*Q(2,3)*Q(2,3))
C
C
      ELAS(3,3) = C12 + C44 + (C11 - C12 - C44)* 
     &  (Q(3,1)*Q(3,1)*Q(3,1)*Q(3,1) + Q(3,2)*Q(3,2)*Q(3,2)*Q(3,2) +
     &   Q(3,3)*Q(3,3)*Q(3,3)*Q(3,3))
C
C
      ELAS(1,2) = C12 + (C11 - C12 - C44)* 
     & (Q(1,1)*Q(1,1)*Q(2,1)*Q(2,1) + Q(1,2)*Q(1,2)*Q(2,2)*Q(2,2) + 
     &  Q(1,3)*Q(1,3)*Q(2,3)*Q(2,3))
C
C
      ELAS(1,3) = C12 + (C11 - C12 - C44)*
     & (Q(1,1)*Q(1,1)*Q(3,1)*Q(3,1) + Q(1,2)*Q(1,2)*Q(3,2)*Q(3,2) + 
     &  Q(1,3)*Q(1,3)*Q(3,3)*Q(3,3))
C
C
      ELAS(2,3) = C12 + (C11 - C12 - C44)* 
     & (Q(2,1)*Q(2,1)*Q(3,1)*Q(3,1) + Q(2,2)*Q(2,2)*Q(3,2)*Q(3,2) + 
     &  Q(2,3)*Q(2,3)*Q(3,3)*Q(3,3))
C
C
C
        ELAS(2,1) = ELAS(1,2)
        ELAS(3,1) = ELAS(1,3)
        ELAS(3,2) = ELAS(2,3)
C
C
      ELAS(4,4) = HC44 + (C11 - C12 - C44)* 
     & (Q(1,1)*Q(2,1)*Q(1,1)*Q(2,1) + Q(1,2)*Q(2,2)*Q(1,2)*Q(2,2) + 
     &  Q(1,3)*Q(2,3)*Q(1,3)*Q(2,3))
C
C
      ELAS(5,5) = HC44 + (C11 - C12 - C44)* 
     & (Q(1,1)*Q(3,1)*Q(1,1)*Q(3,1) + Q(1,2)*Q(3,2)*Q(1,2)*Q(3,2) + 
     &  Q(1,3)*Q(3,3)*Q(1,3)*Q(3,3))
C
C
      ELAS(6,6) = HC44 + (C11 - C12 - C44)* 
     & (Q(2,1)*Q(3,1)*Q(2,1)*Q(3,1) + Q(2,2)*Q(3,2)*Q(2,2)*Q(3,2) + 
     &  Q(2,3)*Q(3,3)*Q(2,3)*Q(3,3))
C
C
      ELAS(4,5) = (C11 - C12 - C44)* 
     & (Q(1,1)*Q(2,1)*Q(1,1)*Q(3,1) + Q(1,2)*Q(2,2)*Q(1,2)*Q(3,2) + 
     &  Q(1,3)*Q(2,3)*Q(1,3)*Q(3,3))   
C
C
      ELAS(4,6) = (C11 - C12 - C44)* 
     & (Q(1,1)*Q(2,1)*Q(2,1)*Q(3,1) + Q(1,2)*Q(2,2)*Q(2,2)*Q(3,2) + 
     &  Q(1,3)*Q(2,3)*Q(2,3)*Q(3,3))  
C
C
      ELAS(5,6) = (C11 - C12 - C44)* 
     & (Q(1,1)*Q(3,1)*Q(2,1)*Q(3,1) + Q(1,2)*Q(3,2)*Q(2,2)*Q(3,2) + 
     &  Q(1,3)*Q(3,3)*Q(2,3)*Q(3,3))
C
C
C
        ELAS(5,4) = ELAS(4,5)
        ELAS(6,4) = ELAS(4,6)
        ELAS(6,5) = ELAS(5,6)!correction
C
C
      ELAS(1,4) = (C11 - C12 - C44)* 
     & (Q(1,1)*Q(1,1)*Q(1,1)*Q(2,1) + Q(1,2)*Q(1,2)*Q(1,2)*Q(2,2) + 
     &  Q(1,3)*Q(1,3)*Q(1,3)*Q(2,3))
C
C
      ELAS(1,5) = (C11 - C12 - C44)* 
     & (Q(1,1)*Q(1,1)*Q(1,1)*Q(3,1) + Q(1,2)*Q(1,2)*Q(1,2)*Q(3,2) + 
     &  Q(1,3)*Q(1,3)*Q(1,3)*Q(3,3))
C
C
      ELAS(1,6) = (C11 - C12 - C44)* 
     & (Q(1,1)*Q(1,1)*Q(2,1)*Q(3,1) + Q(1,2)*Q(1,2)*Q(2,2)*Q(3,2) + 
     &  Q(1,3)*Q(1,3)*Q(2,3)*Q(3,3))
C
C
      ELAS(2,4) = (C11 - C12 - C44)* 
     & (Q(2,1)*Q(2,1)*Q(1,1)*Q(2,1) + Q(2,2)*Q(2,2)*Q(1,2)*Q(2,2) + 
     &  Q(2,3)*Q(2,3)*Q(1,3)*Q(2,3))
C
C
      ELAS(2,5) = (C11 - C12 - C44)* 
     & (Q(2,1)*Q(2,1)*Q(1,1)*Q(3,1) + Q(2,2)*Q(2,2)*Q(1,2)*Q(3,2) +  
     &  Q(2,3)*Q(2,3)*Q(1,3)*Q(3,3))
C
C
      ELAS(2,6) = (C11 - C12 - C44)* 
     & (Q(2,1)*Q(2,1)*Q(2,1)*Q(3,1) + Q(2,2)*Q(2,2)*Q(2,2)*Q(3,2) + 
     &  Q(2,3)*Q(2,3)*Q(2,3)*Q(3,3))
C
C
      ELAS(3,4) = (C11 - C12 - C44)* 
     & (Q(3,1)*Q(3,1)*Q(1,1)*Q(2,1) + Q(3,2)*Q(3,2)*Q(1,2)*Q(2,2) + 
     &  Q(3,3)*Q(3,3)*Q(1,3)*Q(2,3))
C
C
      ELAS(3,5) = (C11 - C12 - C44)* 
     & (Q(3,1)*Q(3,1)*Q(1,1)*Q(3,1) + Q(3,2)*Q(3,2)*Q(1,2)*Q(3,2) + 
     &  Q(3,3)*Q(3,3)*Q(1,3)*Q(3,3))
C
C
      ELAS(3,6) = (C11 - C12 - C44)* 
     & (Q(3,1)*Q(3,1)*Q(2,1)*Q(3,1) + Q(3,2)*Q(3,2)*Q(2,2)*Q(3,2) + 
     &  Q(3,3)*Q(3,3)*Q(2,3)*Q(3,3))
C
C
C        
         ELAS(4,1) = ELAS(1,4)
         ELAS(4,2) = ELAS(2,4)
         ELAS(4,3) = ELAS(3,4)
         ELAS(5,1) = ELAS(1,5)
         ELAS(5,2) = ELAS(2,5)
         ELAS(5,3) = ELAS(3,5)
         ELAS(6,1) = ELAS(1,6)
         ELAS(6,2) = ELAS(2,6)
         ELAS(6,3) = ELAS(3,6)        
C      
C      DO K1 = 1,6
C       DO K2 = 1,6
C        IF (K2.GT.3) THEN
C         ELAS(K1,K2) = ELAS(K1,K2,ICRYS,KINTK)*2.0D0 ! comment: what is KINTK?
C        ENDIF
C       ENDDO
C      ENDDO
C        
      RETURN
      END SUBROUTINE ELAST_BCC
c#####################################################      
c#####################################################
      
       subroutine euler(iopt,ph,th,tm,a)
      dimension a(3,3)
      pi=4.*atan(1.d0)
c
c     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
c     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
c     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
c     ph,th,om ARE THE EULER ANGLES OF ca REFERRED TO sa.
c
      if(iopt.eq.1) then
        th=acos(a(3,3))
        if(abs(a(3,3)).ge.0.9999) then
          tm=0.
          ph=atan2(a(1,2),a(1,1))
        else
          sth=sin(th)
          tm=atan2(a(1,3)/sth,a(2,3)/sth)
          ph=atan2(a(3,1)/sth,-a(3,2)/sth)
        endif
        th=th*180./pi
        ph=ph*180./pi
        tm=tm*180./pi
         
      else if(iopt.eq.2) then
        sph=sin(ph)
        cph=cos(ph)
        sth=sin(th)
        cth=cos(th)
        stm=sin(tm)
        ctm=cos(tm)
        a(1,1)=ctm*cph-sph*stm*cth
        a(2,1)=-stm*cph-sph*ctm*cth
        a(3,1)=sph*sth
        a(1,2)=ctm*sph+cph*stm*cth
        a(2,2)=-sph*stm+cph*ctm*cth
        a(3,2)=-sth*cph
        a(1,3)=sth*stm
        a(2,3)=ctm*sth
        a(3,3)=cth
      endif

      return
      end
c
      
      
c#####################################################      
c#####################################################     
      
      
      
      
      
      subroutine orient(a,c)
      dimension a(3,3),c(3,3),th2(3,3),v(3),vbar(3)
      dimension th(3,3)
      dimension rot(3,3),anew(3,3)

c     BUILD ROTATION TENSOR BASED ON RODRIGUES FORMULA

c#if 0
c      v(1)=c(3,2)
c      v(2)=c(1,3)
c      v(3)=c(2,1)
cc      v(4)=c(1,1)
cc      v(5)=c(2,2)
cc      v(6)=c(3,3)
      
c      snorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3) 
cc     #      +v(4)*v(4)+v(5)*v(5)+v(6)*v(6) )
c      snorm1=tan(snorm/2.)
c      if(snorm.gt.1.e-06) go to 97
c      snorm=1.
c97    do 20 i=1,3
c      vbar(i)=snorm1*v(i)/snorm
c20    continue
c      snorm=vbar(1)*vbar(1)+vbar(2)*vbar(2)+vbar(3)*vbar(3)
cc     # +vbar(4)*vbar(4)+vbar(5)*vbar(5)+vbar(6)*vbar(6)
c      th(3,2)=vbar(1)
c      th(1,3)=vbar(2)
c      th(2,1)=vbar(3)
c      th(2,3)=-vbar(1)
c      th(3,1)=-vbar(2)
c      th(1,2)=-vbar(3)
c      do 40 i=1,3
c40    th(i,i)=0.
c      do 30 i=1,3
c      do 30 j=1,3
c      th2(i,j)=0.
c      do 50 k=1,3
c50    th2(i,j)=th2(i,j)+th(i,k)*th(k,j)
c30    continue
c      do 60 i=1,3
c      do 60 j=1,3
c60    rot(i,j)=(i/j)*(j/i)+2.*(th(i,j)+th2(i,j))/(1.+snorm)
      
c#endif     
      do 70 i=1,3
      do 70 j=1,3
      anew(i,j)=0.
      do 80 k=1,3
80    anew(i,j)=anew(i,j)+c(i,k)*a(k,j)
70    continue
      do 90 i=1,3
      do 90 j=1,3
90    a(i,j)=anew(i,j)
      f=2
      return
      end
c
C *****************************************************************************       
  
    
      
       
