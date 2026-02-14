USE FULIOUMULTI
USE GENERATE_FULIOU_LEVELS ,only : gflq, generate_level_scheme
USE EXTRAS       ,only : getatmosphere, aer_scale_hgt
USE CALIPSO_OUTPUT, only : pack_sky,print_pack_sky,skyp,SKYP_TYPE

USE ICEDIRSFC,only: tau_uc !! Debug Diagnostic
implicit none
 real  aot_by_band,aotf,wlf
 common /aotbyband/ aot_by_band(18) 
 common /tau_spline_aot/ aotf(18),wlf(18)


TYPE (SKYP_TYPE) ut,tu

real psfc,deltau,delcf,delReff,DelGCR,a,b,c,d,e,f,g,h,dmin,PI
integer kk,i,is,j,na,nb,imodis,IL(12),Nfiles,X,mid_day_month(12)
real psel(6),DATA(6,36),temp(20),lev(20),r3(20),r4(20),aod_exp(360,180,12)

character(LEN=140)  foo

PI=4.D0*DATAN(1.D0)

!DATA USED IN MAKING TEMPERATURE PROFILE BELOW
 OPEN(10,FILE='./jhsv.dat') !read top of file
  DO i=1,20 
     read(10,*) a,b,c
      lev(I)=a
      r3(I)=b
      r4(I)=c
  END DO  
! write(*,*) lev
! write(*,*) r3
! write(*,*) r4
! STOP

OPEN(99,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/src/simple/Exp_aod.dat') !read file
read(99,*) aod_exp

!READNUMBER OF FILE LINES IN DATA FILES (BELOW)
OPEN(10,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PETM_M_lines.dat') !read file
 DO i=1,12 
   read(10,*) na,nb
   IL(I)=nb
   END DO  
write(*,*) IL
!write(*,*) r3
!write(*,*) r4

DO Nfiles=1,12!LOOP OVER DATAFILES
X = Nfiles
!READ IN CLOUD MICROPHYSICS VARIABLES FOUND FROM MODIS DATA (HSV)I
     IF (X  .EQ. 1) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_1.dat') !Read DATA files 11
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT1e2aod1H.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 2) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_2.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT2e2aod1H.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 3) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_3.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT3e2aod1H.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 4) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_4.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT4e2aod1H.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 5) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_5.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT5e2aod1H.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 6) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_6.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT6e2aod1H.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 7) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_7.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT7e2aod1H.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 8) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_8.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT8e2aod1H.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 9) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_9.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT9e2aod1H.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 10) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_10.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT10e2aodH.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 11) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_11.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT11e2aodH.dat') !OUTPUT FILE WITH RESULTS
ELSE IF (X .EQ. 12) THEN
OPEN(11,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/Input_MONTHLYtheory_12.dat') !read top of file
OPEN(12,FILE='/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/PT_MODIS_RESULT12e2aodH.dat') !OUTPUT FILE WITH RESULTS
END IF
    DO i=1,10  
       read(11,*)  foo
       WRITE(*,*)  foo
    END DO
   Mid_day_month =[ 14,44,73,104,134,165,195,226,257,287,318,348] 

!STOP
!********************************************************  
!MAIN LOOP OVER GRIT CELLS WITH MODIS CLOUD MICROPHYSICS*
!********************************************************
DO IMODIS=1,IL(Nfiles) !(LON,LAT)


READ(11,'(2(I5,2x),8(F12.6,2x))') na,nb,a,b,c,d,e,f,g,h
     WRITE(*,*) na,nb,a,b,c,d,e,f,g,h
     READ(11,'(20(F12.6,2x))') temp
!MAKE TEMPERATURE PROFILE FILE and write to file=jnew.lev
     open(13,FILE='../../testatms/jnew.lay')
      write(13,*) temp(20) , 19
     DO i=1,20
      write(13,*) lev(I),temp(I),r3(I),r4(I), 0
     END DO
    CLOSE(13) 
!    STOP
!     WRITE(*,*) temp
!STOP
!CHange in cloud microphysics parameterized from Forbush decreases
deltau  = -2.87 !+1.08
delcf   = -5.5  !+ 1.47
delReff =  0.7  !- 0.41
!::::::::::::::::::::::::::::::
 call set_default_options_fu ! Sets some of the more obsure inputs to reasonable values.
fi%lscm   = .false. 
fi%lscm(1:2)   = .true. 
fi%lscm   = .true. 
!InPut Profile Assignment
 call getatmosphere('../../testatms/jnew.lay ', &
! call getatmosphere('./testatms/jmls.lay ',&
 FI%VI%nlev,&
 FI%VI%pp,&
 FI%VI%pt,&
 FI%VI%ph,&
 FI%VI%po,&
 FI%pts)
 FI%VI%nlev = FI%VI%nlev+1  ! LAYER(getatm) to LEVEL

 FI%VI%hsfc = 0.00 !! SURFACE GEOPOTENTIAL OF FI%VI profile
! FI%VI%hsfc = 1600 !! SURFACE GEOPOTENTIAL OF FI%VI profile

 !gflq%hsfc = 1500. !Meters Surface elev. of ACTUAL FOV... to nearest 120m Multiple
  gflq%hsfc = 2.0   ! 


 gflq%mode = 'CALIP'
 gflq%mode = 'CERE3'
 
!CERES INSTRUMENT SPECTRAL WINDOWS FOR FLIGHT MODEL 3 HSV
fi%instrument= 3
! gflq%nld =4
! gflq%internal_levels(1:4) = (/70.,200.,500.,850./)

fi%HYBRID_SW_SOLVER =.true. !200130802 SYNI Ed4 $S Clear , 2S HOMO CLD , GWTSA INHOM Cld
fi%HYBRID_SW_SOLVER =.false. !checks isksolve,fourssl

fi%isksolve= 0  ! Solver Method (0=fu 1=gwtsa) 
fi%fourssl = .false.
!withincrease in solar irradiation 1.001 ~ 0.1%
fi%ss      = 1365.*(1.000110 + 0.034221*COS(2.*PI*Mid_day_month(I)/365.)+& ! Solar Constant wm-2+
                               0.001280*SIN(2.*PI*Mid_day_month(I)/365.)+& ! Eccentricity
                               0.000719*COS(4.*PI*Mid_day_month(I)/365.)+&
                               0.000077*SIN(4.*PI*Mid_day_month(I)/365.)  ) 
fi%u0      =  a   ! Cosine Solar Zenith Angle
fi%ur      =  0.8 ! Cosine View Zenith Angle (for IR Radiance)
fi%umco2   =  400 ! CO2 concentration in ppm


!-------Cnd 2
fi%wp_hgt_flag = 0  ! Constant lwc with height
!fi%wp_hgt_flag = 1  ! Water Cloud  Top thicker than Base
!fi%wp_hgt_flag = 2  ! Ice Cloud  Bottom thicker than top
!*********************************************************************************
!*********************************************************************************

DO j =1,50

DelGCR = FLOAT(j-1)/5.


fi%fc%dpi%ldpi  = .false.

fi%fc(1)%cldfrac   = b*(1.0 + DelGCR*Delcf*aod_exp(na,nb,Nfiles)/100.)  ! Cloud Fraction Liquid (0-1) 
fi%fc(2)%cldfrac   = c*(1.0 + DelGCR*Delcf*aod_exp(na,nb,Nfiles)/100.)  !HSV ICE
!i1fi%fc(1)%novl      =   1 !layer 2 clouds
fi%fc(1)%novl      =   1 !First layer clouds used in primary test
fi%fc(2)%novl      =   1 !HSV layer 2 clouds   
!  WRITE, DATA(2,i)
!FI%VD%cldpres(1:2, 1,1) = (/200,400/)
!FI%VD%cldpres(1:2, 1,2) = (/704,725/)

! SETTING CLOUD BOTTOM PRESSURE (NOT MEASURED)
dmin = 1000.
IF (d < 800.)  THEN
   dmin = d+200.  
ELSE IF (d > 800. .and. d < 900. )  THEN
   dmin = d + 100 
ELSE IF (d > 900. .and. d < 950)   THEN
        dmin = 940.
ELSE IF (d > 950. .and. d < 980)   THEN
        dmin = 990
ENDIF

FI%VD%cldpres(1:2, 1,1) = (/d,dmin/) !cloud extension used in primary test
!write(*,*) FI%VD%cldpres(1:2, 1,1),d,dmin
!STOP
FI%VD%cldpres(1:2, 2,1) = (/200.,300./) !HSV

fi%fc(1)%rphase(1)    =  1.0    ! Cloud Phase 1=Water 2=Ice
fi%fc(2)%rphase(1)    =  2.0    ! Cloud Phase 1=Water 2=Ice
fi%fc(1)%re(1) = e*(1.0 + DelGCR*DelReff*aod_exp(na,nb,Nfiles)/100.) 
fi%fc(2)%de(1) = f
!!!    WRITE(*,*), DATA(5,i)
!fi%fc(1)%asp(1) = exp(iasp*0.1) !! Fu 20006 Ice AspectRatio !!!!! NEW FOR 20010130


fi%fc(1)%tau_vis(1)       = g*(1.0 + DelGCR*Deltau*aod_exp(na,nb,Nfiles)/100.)  ! Cloud Visible Optical Depth ( Minnis)
fi%fc(2)%tau_vis(1)       = h  ! HSV Cloud Visible Optical Depth ( Minnis)
fi%fc(1)%sc(1)%mn_lin_tau =  fi%fc(1)%tau_vis(1) *1.15
fi%fc(2)%sc(1)%mn_lin_tau =  fi%fc(2)%tau_vis(1) *1.15 !HSV
!	    WRITE(*,*), DATA(2,i)
!-----
!fi%fc(1)%rphase(2)    =  2.0    ! Cloud Phase 1=Water 2=Ice
!fi%fc(1)%re(2) = 10.

!HSV ***********************************************
!fi%fc(1)%asp(1) = exp(iasp*0.1) !! Fu 20006 Ice AspectRatio !!!!! NEW FOR 20010130

!fi%fc(1)%tau_vis(2)       = 30    ! Cloud Visible Optical Depth ( Minnis)
!fi%fc(1)%sc(2)%mn_lin_tau =  fi%fc(1)%tau_vis(2) 
!
!fi%fc(1)%tau_vis(2)       = 1E-20    ! Cloud Visible Optical Depth ( Minnis)
!fi%fc(1)%sc(2)%mn_lin_tau = 1E-20 !fi%fc(1)%tau_vis(2) 


!Surface Properties --------------------------------------------------

!Allow different albedos for Aerosol Vs. NO Aerosol cases , And for each Clear/Cloud Conditions
fi%sfcalb(1:18,1,0)  = 0.0 ! Clear sky -Spectral Surface Albedo SW
fi%sfcalb(1:18,2,0)  = 0.0 ! Pristine sky -Spectral Surface Albedo SW
fi%sfcalb(1:18,1,1:)  = 0.0  ! CLOUDY w/AOT  sky -Spectral Surface Albedo SW
fi%sfcalb(1:18,2,1:)  = 0.0  ! CLOUDY w/o AOT sky -Spectral Surface Albedo SW

fi%ee(1:12)  = 0.99 ! Spectral Surface Emissivity LW

!Aerosols ------------------------------------------------------------
fi%nac         = 2   ! 2 aerosol types 
fi%itps(1)     = 2   ! Continental see types (1-18)
fi%itps(2)     = 1   
!fi%itps(2)	     = 11	   ! Soot	  see types (1-18)

fi%n_atau      = 1   ! 1 wavelength input for aerosols
fi%a_wli(1)    = 0.641   ! AOT wavelength(microns) of a_taus
fi%a_taus(1,1) =  0.80   ! AOT for constituent 1

fi%a_taus(1,2) =  0.20   ! AOT for constituent 2

!----------------------------------------------------------------------


 call generate_level_scheme !! Define model Fixed layer structure pre-cloud by fixed DZ intervals...
! call print_vla_in 

 call prepare_model_profile_fu !! CALL After all FI%VD and FI%VI structures are defined.
 call vla_interface_fu     ! uses FI%VO !! Assign Model ATM Profile and CLD Levels
! call print_vla_out

!Aerosol Profile (after fi%pp is created )-----------------------------

 call aer_scale_hgt(fi%nv,fi%pp,3.0,fi%aprofs(1:fi%nv,1) )
 call aer_scale_hgt(fi%nv,fi%pp,3.0,fi%aprofs(1:fi%nv,2) )
! RADIATVE TRANSFER --------------------------------------------------

! if (i == 1 )
 call print_in_fu   ! PRINTS INPUTS  AS ASCII 
 
 call rad_multi_fu  ! CALL THE CODE !!!
WRITE(12,'(2(I5,2x),11(F12.6,2x))') na ,nb, delGCR, ftoa(2)%swup,ftoa(2)%olr,a,&
        b*(1.0 + DelGCR*Delcf/100.),c,d,e*(1.0 +DelGCR*DelReff/100.),f,g*(1.0 + DelGCR*Deltau/100.),h 
!:q call print_out_fu		   ! PRINTS Lots of OUTPUTS  AS ASCII
print'(3f10.2,f10.3)',&
!fi%fc(1)%de(1),&
ftoa(2)%olr,&
!ftoa(2)%swdn,&
ftoa(2)%swup!, &
!!ftoa(2)%swup/ftoa(2)%swdn
!fsfc(2)%swdn,&
!fsfc(2)%swdir,&
!fsfc(2)%swdif,&
!tau_uc(:,1)
 call pack_sky
 
! call print_pack_sky
 END DO ! GCR LOOP
 END DO ! MODIS LOOP (LON,LAT)
 END DO ! FILE MONTH MODIS LOOP
stop ' Simple.f90 normal end'

end
