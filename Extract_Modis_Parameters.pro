pro  ExtractModisParameters
  ;*********************************************************************************
  ;Parameters for the Fu Liou program
  ;One file for each month
  ;Reads MODIS files and extract parameters for Fu Liou program
  ;Data are sent to '/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS20XXMONTH/'
  ;Have 2003 used monthly data                                          
  ;*********************************************************************************
  SET_PLOT,'ps'
    close,/all
  dir = '/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/'
  dirdat = '/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS2003MONTH/'

  pos = cgLayout([1,1], OXMargin=[10, 10], OYMargin=[10, 10], XGap=5, YGap=6)
  
  Nteoriflag = 0 ;USING TEORETICAL ZENITH ANGLES:1 OR MODIS MEASURED :0 
  
  ;Find middle day of all 12 month. Used in cosz
  days_month=[31,28,31,30,31,30,31,31,30,31,30,31]
  ;month = ['jan1','feb2','mar3','apr4','may5','jun6','jul7','aug8','sep9','oct10','nov11','dec12']
  month = ['1','2','3','4','5','6','7','8','9','10','11','12']
   days_s   = fltarr(12) 
   Phi      = !pi*(findgen(180)-89.5)/180.
   delta    = 23.45/180.*!pi*sin(!pi*360.*(284. + findgen(365)+1)/365./180)
  FOR I=0,11 DO BEGIN
    ; print,days_month(I)/2
    days_s(I)=FIX(julday(I+1,days_month(I)/2,2003) - julday(1,1,2003))
  ENDFOR
  
  
  FILENAME = FILE_SEARCH(dirdat,'*.hdf'); FIND FILES

  P = FLTARR(360,180,8)
  Np = fltarr(12)
  N_parm  = [0,12,13,10,4,5,2,3] ;relevant parameters

  FOR Nm=0,11 DO BEGIN ;LOOP OVER FILES (12 MONTH)
 
  OPENW,1,STRCOMPRESS(dirdat+'Input_MONTHLY_mu2003' + month(Nm)+'.dat') ; OPEN FILE TO PUT EXTRACTED PARAMETERS 
  ; get hdf data
  sd_id   = HDF_SD_START(FILENAME(NM), /READ )
  PARAM   = ['Solar_Zenith_Mean_Mean','Cloud_Fraction_Mean_Mean','Cloud_Optical_Thickness_Liquid_Mean_Mean','Cloud_Optical_Thickness_Ice_Mean_Mean','Cloud_Effective_Radius_Liquid_Mean_Mean',$
    'Cloud_Effective_Radius_Ice_Mean_Mean','Cirrus_Fraction_Infrared','High_Cloud_Fraction_Infrared','Cloud_Water_Path_Liquid_Mean_Mean',$
    'Cloud_Water_Path_Ice_Mean_Mean','Cloud_Top_Pressure_Mean_Mean','Retrieved_Temperature_Profile_Mean_Mean','Cloud_Retrieval_Fraction_Liquid_FMean','Cloud_Retrieval_Fraction_Ice_FMean']
  PARMfac = [0.01,1.E-4,0.01,0.01,0.01,0.01,1.E-4,1.E-4,1.0,1.0,0.1,0.01,1.E-4,1.E-4]

  HDF_SD_GETDATA, HDF_SD_SELECT(sd_id, HDF_SD_NAMETOINDEX(sd_id, 'Ocean_Fraction_Day_FMean')), Dummy
  Ocean = reverse(DUMMY,2)


  FOR II=0,7 DO BEGIN ;READ HDF FILE PARAMETERS GIVEN BY PARAM(N_parm)
    HDF_SD_GETDATA, HDF_SD_SELECT(sd_id, HDF_SD_NAMETOINDEX(sd_id, PARAM(N_parm(II)))), Dummy
    P(*,*,II) = reverse(DUMMY,2)*PARMfac(N_parm(II))
  ENDFOR

  ;READ TEMPERATURE-PROFILE
  HDF_SD_GETDATA, HDF_SD_SELECT(sd_id, HDF_SD_NAMETOINDEX(sd_id, PARAM(11))), Dummy
  DUMMY = reverse(DUMMY,2)
  TEMP  = PARMfac(11)*(DUMMY+15000)

  ;WRITING THE EXTRACTED PARAMETERS TO FILE
  printf,1,'lon',' lat'
  FOR I=0,7 DO printf,1,param(N_parm(I)) ;LOOP OVER PARAMETERS names
  printf,1,'Retrieved_Temperature_Profile_Mean'
  FOR I=0,359 DO BEGIN                   ;LOOP LON
    FOR J=0,179 DO BEGIN                 ;LOOP LAT
      A = WHERE(P(I,J,*) GE 0.,COUNT)
      B = WHERE(DUMMY(I,J,*) gt 0,COUNT2)
      IF COUNT EQ 8 and OCEAN(I,J) EQ 10000 and COUNT2 EQ 20 THEN BEGIN
        IF Nteoriflag EQ 0 THEN  BEGIN 
         PRINTf,1,I,J,COS(P(I,J,0)*!pi/180.),$
          P(I,J,1),P(I,J,2),P(I,J,3),P(I,J,4),P(I,J,5),P(I,J,6),P(I,J,7),FORMAT='(2(I5,2x),8(F12.6,2x))'
        ENDIF
        IF Nteoriflag EQ 1 THEN  BEGIN
         PRINTf,1,I,J,insolation_weighted_cosz( phi(J),delta(days_s(Nm))),$
                 P(I,J,1),P(I,J,2),P(I,J,3),P(I,J,4),P(I,J,5),P(I,J,6),P(I,J,7),FORMAT='(2(I5,2x),8(F12.6,2x))'
        ENDIF           
         PRINTf,1,TEMP(I,J,*),FORMAT='(20(F12.6,2x))'
      ENDIF
    ENDFOR                               ;END LAT
  ENDFOR                                 ;END 



  close,/all
  OPENR,1,STRCOMPRESS(dirdat+'PETM_MONTHLY' + STRING(NM+1)+'.dat')
  foo =''
  Nlines= FILE_lines(STRCOMPRESS(dirdat+'PETM_MONTHLY' + STRING(NM+1)+'.dat'))
  PRINT,NM+1,(Nlines-10)/2 ;number of grit cells with parameters ;10 lines of header and 2 lines pr grid point
  Np(Nm) = (Nlines-10)/2
  
  CLOSE,1
  ENDFOR ;LOOP OVER FILES
  
  OPENW,2,STRCOMPRESS(dirdat+'PETM_M_lines.dat')
  FOR I=0,11 DO PRINTF,2,I+1,Np(I),FORMAT='(2(I6,2x))'
  CLOSE,2 
  STOP
END
FUNCTION  Hour_angle_ss , phi, delt
  ;FINDING THE HOUR ANGLE
  IF delt*phi LT 0.D0                    THEN h0 = 0.D0
  IF delt*phi GT 0.D0                    THEN h0 =!pi
  IF (abs(delt)-!pi/2.+abs(phi)) LT 0.D0 THEN h0 = ACOS(-tan(phi)*tan(delt))
  return,h0
END
FUNCTION  insolation_weighted_cosz ,phi ,delta
  h0 = Hour_angle_ss(phi, delta)
  X1 = h0*sin(delta)*sin(phi)+cos(delta)*cos(phi)*sin(h0)
  X2 = h0*(2.D0*(sin(delta)*sin(phi))^2 + (cos(delta)*cos(phi))^2) + $
    cos(delta)*cos(phi)*sin(h0)*(cos(delta)*cos(phi)*cos(h0)+4.D0*sin(delta)*sin(phi))
  cosz = X2/X1/2.D0

  return,cosz
END



