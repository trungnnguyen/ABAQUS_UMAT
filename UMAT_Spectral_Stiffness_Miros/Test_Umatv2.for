      Program Test_Umat 
       
      parameter(ninc=40) !number of increments
      parameter(nstatv=100000) !number of state variables, needs to be calculated
      parameter(NDI=3,NSHR=3,NTENS=6,NOEL=1,NPT=1,NPROPS=1)
      DIMENSION STRESS(NTENS),STATEV(NSTATV),STATEV_old(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      
      double precision  dtime,strain_rate 
      double  precision t1,t2
       DROT(1,2)=0.
       DROT(1,3)=0.
       DROT(2,3)=0.
       DROT(1,1)=1.
       DROT(2,2)=1.
       DROT(3,3)=1. 
      strain_rate=0.001
      dtime=.01
      STRESS=1.
      DSTRAN=strain_rate*dtime
      print*,"Applied Strain Rate: "
      print*
      print*, strain_rate
       
      print*,"=========================="
      print*,"Time Increment: "
      print*
      
      print*, dtime
      print*,"=========================="
       
      print*,"Stress before calling UMAT: "
      print*
      print*, STRESS
      print*,"=========================="
      print*
      print*,"Calling UMAT ..." 
      print*
      print*

       call CPU_TIME(t1)
       
      Do i=1,20
      Do j=1,100 

      call UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1           RPL,DDSDDT,DRPLDE,DRPLDT,
     2           STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3           NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4           CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      Enddo
      Enddo       

       call CPU_TIME(t2)
         
     

c      f=2 
      
            
      print*,"UMAT Job finished! "
      print* 
      print*,"=========================="
      print*,"Updated Stress From UMAT "
      print*
      print*, STRESS
      print*,"=========================="
      print*,"Jacobian Obtained From UMAT "
      print*
      print*, DDSDDE
c      pause
      
      print*,"Elapsed Time:"
      print*
      print*,(t2-t1)
      
      End program Test_Umat
      
 
   