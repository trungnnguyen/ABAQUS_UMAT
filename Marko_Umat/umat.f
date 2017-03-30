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
      real*8 temp
      integer iii,jjj

#if 1
!$ACC parallel do  
       Do iii=1,100000
       temp=0.0
       Enddo
!$ACC end parallel do
#endif 
 
      YM=128000.0D0
      PV=0.33D0
 

C
C     Elastic tangent modulus
C
      
      s1=(1.D0+pv)*(1.0D0-2.0D0*pv)
      s2=(1.D0-pv)/s1
      s3=ym/(1.D0+pv)
      s4=s3*0.5d0
      s5=s3/(1.D0-2.D0*pv)
      s6=s5*pv
      s7=s5-s6
     

     
      DDSDDE(1,1)=s7
      DDSDDE(1,2)=s6
      DDSDDE(1,3)=s6
      DDSDDE(1,4)=0.D0
      DDSDDE(1,5)=0.D0
      DDSDDE(1,6)=0.D0
      DDSDDE(2,1)=s6
      DDSDDE(2,2)=s7
      DDSDDE(2,3)=s6
      DDSDDE(2,4)=0.D0
      DDSDDE(2,5)=0.D0
      DDSDDE(2,6)=0.D0
      DDSDDE(3,1)=s6
      DDSDDE(3,2)=s6
      DDSDDE(3,3)=s7
      DDSDDE(3,4)=0.D0
      DDSDDE(3,5)=0.D0
      DDSDDE(3,6)=0.D0
      DDSDDE(4,1)=0.D0
      DDSDDE(4,2)=0.D0
      DDSDDE(4,3)=0.D0
      DDSDDE(4,4)=s4
      DDSDDE(4,5)=0.D0
      DDSDDE(4,6)=0.D0
      DDSDDE(5,1)=0.D0
      DDSDDE(5,2)=0.D0
      DDSDDE(5,3)=0.D0
      DDSDDE(5,4)=0.D0
      DDSDDE(5,5)=s4
      DDSDDE(5,6)=0.D0
      DDSDDE(6,1)=0.D0
      DDSDDE(6,2)=0.D0
      DDSDDE(6,3)=0.D0
      DDSDDE(6,4)=0.D0
      DDSDDE(6,5)=0.D0
      DDSDDE(6,6)=s4

      h1=DSTRAN(1)
      h2=DSTRAN(2)
      h3=DSTRAN(3)
      h4=DSTRAN(4)
      h5=DSTRAN(5)
      h6=DSTRAN(6)

C
      t1=STRESS(1)+s7*h1+s6*h2+s6*h3
      t2=STRESS(2)+s6*h1+s7*h2+s6*h3
      t3=STRESS(3)+s6*h1+s6*h2+s7*h3
      t4=STRESS(4)+s4*h4
      t5=STRESS(5)+s4*h5
      t6=STRESS(6)+s4*h6   
C
C     new stress at the end of increment
C 
      STRESS(1)=t1
      STRESS(2)=t2
      STRESS(3)=t3
      STRESS(4)=t4
      STRESS(5)=t5
      STRESS(6)=t6

     
      epx=0.0D0
c
      return
      end
