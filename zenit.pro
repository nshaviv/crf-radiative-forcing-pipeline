pro CERESmapFIG_GRL

SET_PLOT,'ps'
DEVICE,bits=8,file='//Users/hesv/Documents/Work/CERES/CERESmapFIG_GRL.ps'
DEVICE,xsize=20,ysize=24
DEVICE,xoffset=0,yoffset=3
DEVICE,/Color


!p.charthick = 1.0
;!x.thick     = 3.0
;!y.thick     = 3.0
;!p.thick     = 3.0
!P.CHARSIZE=1.2
!x.charsize=0.9
!y.charsize=0.9
!p.multi=[0,2,1]
!p.font=0

SPAWN, 'killall Skim'
 

;Make plot symbol (circle)
amp = 1.5
t = findgen(21)/20.*360.
X =  amp*sin(t/180.*!pi)
Y =  amp*cos(t/180.*!pi)
USERSYM, X, Y ,/fill
;ct = 16 ;color table used
;mkparula
;cgLoadCT, 28, /Brewer, NColors=ncolors, Bottom=1, Clip=[10, 210]
mkviridis
;fake_parula
;magma
  Z_map_all = fltarr(180,6)

;FOR Nsign =0,1 DO BEGIN 
; kp=-1
FOR mainloop=0,0 DO BEGIN ;MAIN LOOP TO END OF PROGRAM *******************************************************************
 
  dir = '/Users/hesv/Documents/Work/CERES/'


  name = ['CERES_SYN1deg-200003-200111.nc',$ ;0 #days 640
          'CERES_SYN1deg-200112-200308.nc',$ ;1 #days 639
          'CERES_SYN1deg-200309-200505.nc',$ ;2 #days 639
          'CERES_SYN1deg-200506-200601.nc',$ ;3 #days 245
          'CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20060201-20080804.nc',$ ;4  #days 916
          'CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20080805-20110206.nc',$ ;5  #days 916
          'CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20110207-20130810.nc',$ ;6  #days 916
          'CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20130811-20160212.nc',$ ;7  #days 916
          'CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20160213-20161231.nc',$ ;8  #days 323
          'CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20240401-20240531.nc','CERES_SYN1deg-201202-201205.nc']

;Open the netcdf file,  ncid is an idl string
nday    = 36   ;Total number of days
mdel =  0     ; additional days before day -15
DEL  = Nday-36L ;fixing to make longer time intervals than 36 days 
;NFdays = [640,639,639,245]
NFdays = [640,639,639,245,916,916,916,916,323]
Ndays_total =  TOTAL(NFdays, /INTEGER)
Fsum    = fltarr(Ndays_total)
Fsum_z  = fltarr(180,Ndays_total)
Fsum_lon  = fltarr(360,Ndays_total)
Fmap    = fltarr(360,180,Ndays_total)
SUM_map = fltarr(360,180,nday)
;nt      = [0,640,639,639,245]
;NFdays = [640,639,639,245]

nt      = [0,640,639,639,245,916,916,916,916,323]
NFdays  =   [640,639,639,245,916,916,916,916,323]

ntsum=FIX(total(nt,/cumulative)) 


;print,Ndays_total
;stop


Pacific_mask =land_sea_mask(360,180,/pacific) ;sea mask for centerd at the Pacific as the CERES data are


;**********************************************************************************
FOR K=0,8 DO BEGIN ;LOOP FOR OPENING CERES FILES

ncid = NCDF_OPEN(dir+name(k),/nowrite)

  filename=dir+name(k)
    nstru = NCDF_INQUIRE(ncid) ;getting information of the number of variables 
   print, nstru

    dims = lonarr(nstru.ndims)
    names = strarr(nstru.ndims)
    
    for i=0L,nstru.ndims-1 do begin  ;names of the dimensions
      NCDF_DIMINQ,ncid,i,dname,dsize
      names[i] = dname
      dims[i] = dsize
;        print,'Dimension '+strtrim(i,2)+': ',dname,dsize
    endfor 
    for i=0L,nstru.nvars-1 do begin  ;getting names of variables
    vardesc = NCDF_VARINQ(ncid, i)
    print,vardesc
    endfor

 PARM = ['toa_net_all_daily','toa_sw_all_daily','toa_lw_all_daily']   

varid3 = PARM(mainloop)
out_dat = dir + STRCOMPRESS(varid3+'_2000_2005_zonal.dat')  ;data out file

NCDF_VARGET, ncid,  varid3, var_name3
help,var_name3

;Close the netcdf file
NCDF_CLOSE, ncid

FOR ai=0,179 DO BEGIN
  FOR aj=0,359 DO BEGIN
    Fmap(aj,ai,Ntsum(k):Ntsum(k)+Nt(k+1)-1) =   var_name3(aj,ai,*)
  ENDFOR
ENDFOR


ENDFOR ; END READING CERES FILES LOOP******************************************************
;

print,max(Fmap),MIN(Fmap)


ZZ    = global_sum(Fmap)
Zonal = zonal_sum(Fmap);,mask=Pacific_mask);,/invert)
Lon   = longitude_sum(Fmap)
co = cos(!pi*(findgen(180)-90)/180.)
;Zonal_area = fltarr(180,2163)
Zonal_area = fltarr(180,Ndays_total )

For i=0,179 do Zonal_area(i,*) = Zonal(i,*)*co(i)
Fsum(*)    =    total(Zonal_area(30:105,*),1)/85.
;Fsum(Ntsum(k):Ntsum(k)+Nt(k+1)-1) =    ZZ

FOR i=0,179 do Fsum_z(i,*) =    Zonal(i,*)
FOR i=0,359 do Fsum_lon(i,*) =    Lon(i,*)

close,1
;READ IN NEUTRON DATA IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII;
dir = '/Users/hesv/Documents/Work/Isccp/MODIS/MODIS_long/'
dat_file = 'OULU2000_2005II.dat'
openr,1,dir+dat_file
Y = 0.
M = 0.
D = 0.
Fsum_ou   = fltarr(2053+50)
FOR i =0,2052 DO BEGIN
  readf,1,Y,M,D,Neu
  Fsum_ou(i)   = Neu
ENDFOR
close,1
;ØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØØ

nfig =0

FOR jk = 0,1 DO BEGIN

kp = jk

J0  =  julday(3,1,2000)
;JJd = [julday(10,31,2003),$;0   
;       julday(1, 19,2005),$;1
;       julday(9, 13,2005),$;2
;       julday(7, 16,2000),$;3
;       julday(4, 12,2001),$;4
;       julday(11,10,2004),$;5
;       julday(9, 26,2001),$;6
;       julday(7, 17,2005),$;7
;       julday(7, 27,2004),$;8       
;       julday(5, 31,2003),$;9              
;       julday(11,25,2001),$;10       
;       julday(5, 15,2005),$;11
;       julday(8, 28,2001) ];12
       
     
       
Nfd   = 5
Nsign = jk ;0 negativ declination.  1 positiv declination
      
JJd = [julday(10,31,2003),$;7 negativ dec -15
       julday(1, 19,2005),$;10             -20
       julday(11,10,2004),$;9             -18
       julday(11,25,2001),$;3             -21
       julday(12,15,2006)];14             -23


;test with small declinations
;     JJd = [julday(3,9,2012),$;             -4.8
;            julday(9,26,2001),$;         dec -2.2
;            julday(9, 13,2005)];             3.0



       
       
       
 IF Nsign EQ 1 THEN BEGIN      
   Nfd   = 6
;  positiv declination
;  JJd = [julday(9, 13,2005),$  ;positiv dec
;         julday(7, 16,2000),$
;         julday(4, 12,2001),$
;         julday(7, 17,2005),$;7
;         julday(7, 27,2004),$;8 
;         julday(5, 31,2003),$;9 
;         julday(5, 15,2005),$;11
;         julday(8, 28,2001) ];12    
         
       JJd = [$  ;positiv dec (all 6)
         julday(7, 16,2000),$;1  21
         julday(7, 17,2005),$;12 21
         julday(7, 27,2004),$;8  18
         julday(5, 31,2003),$;6  21
         julday(5, 15,2005),$;11 18
         julday(6, 25,2015) $;16 23. 
          ];12             
ENDIF

;  JJd = [julday(9, 13,2005),$;2
;  julday(9, 26,2001)];6


;JJd  = J0+50 + FIX(1800.*randomu(seed,Nfd)) ;take random FD dates
Fs   = 0    
FLAG = 0  ; 0 use lattitude.  1 use Solar zenit lattitude
Sum   = fltarr(nday)
Sum_z = fltarr(180,nday)
Sum_lon = fltarr(360,nday)
Sum_zs = fltarr(180,nday)
; Neutron monitor data
Sum_ou    = fltarr(nday)

Del_gcr = fltarr(Nfd)
;GLOBAL NEUTRON MONITOR DATA
;FOR I=0,Nfd-1 DO BEGIN
;  SUM_ou(*) = SUM_ou(*) + Fsum_ou(JJd(I)-J0-15L-mdel:JJd(I)-J0+20L+DEL-mdel)
;  foo = Fsum_ou(JJd(I)-J0-15L-mdel:JJd(I)-J0+20L+DEL-mdel)
;  
;  Del_gcr(I) = (mean(foo) -MIN(foo(13:20)))/MEAN(foo)*100
;  
;Zen_deg = Zenit(JJd(I))
;I0 =  (Zen_deg)
;print,' I_0 = ',I0, '  F#  ',I
;ENDFOR


;;;print, 'Del_gcr [%]  ',del_gcr,mean(Del_gcr(0:3)),stddev(Del_gcr(0:3))
;STOP
;Del_gcr [%] neg       19.9571      12.7167      9.62968      6.82617
;12.2824      5.65380
;Del_gcr [%]  pos      12.0170      7.33527      7.68268      7.49456      6.47751
;8.20141      2.18230

;stop

;GLOBAL CERES DATA
print,JJd(0)-J0-15L-mdel,JJd(0)-J0+20L+DEL-mdel


    FOR I=0,Nfd-1 DO BEGIN
      FOR J=0,179 DO BEGIN                                                                        
        FOR k=0,359 DO SUM_map(k,j,*) = SUM_map(k,j,*) + Fmap(k,j,JJd(I+Fs)-J0-15L-mdel:JJd(I+Fs)-J0+20L+DEL-mdel); * Del_gcr(I)/total(Del_gcr(*))*float(Nfd)     
      ENDFOR
    ENDFOR
    xtit='!7Lattitude [deg]'
    xtit2='!7Longitude [deg]'

SUM_map = SUM_map/float(Nfd)



;LEVELS and labels of CONTOUR PLOT
if varid3 eq 'toa_net_all_daily' then begin
  lev = [-3.,-1,0.,1.,3.]
  bartit = '!9D!7NET  [W/m!u2!n]'
endif
if varid3 eq 'toa_sw_all_daily'  then begin
  lev = [-4.,-2.,0.,2.,4.]
  bartit = '!9D!7SW  [W/m!u2!n]'
endif
if varid3 eq 'toa_lw_all_daily'  then begin
  lev = [-4.,-2.,0.,2.,4.]
  bartit = '!9D!7LW  [W/m!u2!n]'
endif


;ENDFOR ;JK



;REMOVE LINEAR TREND
  For i=0,179 DO BEGIN
   For j=0,359 DO BEGIN
    P=linfit(findgen(nday),SUM_map(j,i,*))
    SUM_map(j,i,*) = SUM_map(j,i,*) - POLY(findgen(nday),P)
   endfor
  endfor


;****************************************************************************************
;showing the areas where there are maximum/minimum in effect
;LEVELS and labels of CONTOUR PLOT
if varid3 eq 'toa_net_all_daily' then begin
  lev = [6,8]
  bartit = '!9D!7NET  [W/m!u2!n]'
endif
if varid3 eq 'toa_sw_all_daily'  then begin
  lev = [-12,-8]
  bartit = '!7TOA SW [W/m!u2!n] '
endif
if varid3 eq 'toa_lw_all_daily'  then begin
  lev = [6,8]
  bartit = '!7TOA LW [W/m!u2!n] '
endif


X0 = 0.
Y0 = 0.2
X1 = 1.
Y1 = 1.

lat = findgen(180)-90 +0.5
lon = findgen(360) +0.5

;*****************************************************************************************
pos = cgLayout([2,3], OXMargin=[8, 1], OYMargin=[15, 5], XGap=10, YGap=12,iymargin=[1,1])
;stop
DELx = 0.0
dely = 0.00


xtit = ['','','','','!7longitude','!7longitude']
ytit = ['!7latitude','','!7latitude','','!7latitude','']
Fig_num = ['a)','d)','b)','e)','c)','f)']
tit  = ['Average day -15 to 0','Average day 3 to 12','','','',''] 
xyou = ['NET before FD','NET before FD','NET after FD','NET after FD','LW before FD','LW after FD']
xsty = [7,7,7,7,1,1]
ysty = [1,7,1,7,1,7]



FOR jm=0,1 DO begin

  Max_map = fltarr(360,180)
  
  
  ;Average of period before and after FD's  
  days_a = 9    ;nom 9
  L0_0   = 18   ;nom =18, correspond to day 3
  dd     =  7    ;nom = 7,  
  
 
  print,'jm  = ',jm
  L0 =  L0_0*(jm mod 2)
  for l=0,days_a-1 + dd*((jm+1) mod 2) DO MAX_map(*,*) = MAX_map(*,*) + SMOOTH(Sum_map(*,*,L0+l),[10,5])
  Max_map = max_map/(days_a + dd*((jm+1) mod 2)) ;9 days
  foo_map = Max_map
  print,mean(foo_map),jm
  IF jm eq 0 then print,'day0',L0,' to ',days_a-1 + dd*((jm+1) mod 2)
  IF jm eq 1 then print,'day0',L0,' to ',L0+days_a-1 + dd*((jm+1) mod 2),days_a-1 + dd*((jm+1) mod 2)

  
  Pacific_mask(*,0)=0. ;First row in antarctica set to zero as it should be 
  
  A = WHERE(Pacific_mask EQ 0)
  Foo_map(A)     = -999
  Z_map_all(*,jm)= zonal_sum0(foo_map,360,180)
  A = WHERE(Z_map_all(*,jm) LT -900)
  Z_map_all(A,jm) = 0.
  
  ;END FD's AVERAGE

;colar tables
mkparula,0,255
;cgloadct,72,/reverse

     
    kp = jm*2 + jk 
       
       print, '  KP = ',kp

   cgplot,lon,lat,/nodata,yrange=[-90,90],xrange=[0,360],xstyle=1,ystyle=1,POSITION=pos[*,kp],xtitle='lon',ytitle='lat'$
    ,ticklen=-0.02,yticks=6,xticks=6,NoErase=jm NE 0,title=tit(jm)

   cgImage, bytscl(foo_map*Pacific_mask,min=-10,max=10) , Margin=1.0, /Scale,  /Save ,/overplot,NoErase=1 NE 0

   cgContour, foo_map, LEVELS=lev, /OnImage,C_Color=cgcolor('orange','black'),NoErase=jm NE 0
    
    
  FOR i=1,5 DO oplot,[0,360],30*[(i),(i)],li=1,col=cgcolor('charcoal')
  FOR i=1,5 DO oplot,60*[(i),(i)],[0,180],li=1,col=cgcolor('charcoal')
  cgMap_Set,0,180, /Cylindrical, /NoBorder, /NoErase, $
        Limit=[-90, 0, 90, 360], POSITION=pos[*,kp]
;    Limit=[-90, 0, 90, 360], POSITION=[0.15, 0.3, 0.95, 1.]
  cgMap_Continents, Color=cgcolor('black'),/fill
 PRINT,'jmx ',kp
 
 XYOUTS, 10,75,Fig_num(kp) ,Color=cgcolor('white')
 nfig = nfig+1
 XYOUTS, 10,-85, STRING(xyou(kp)),Color=cgcolor('orange')
 ENDFOR ;jm loop
 jm = 1

 mkviridis
 mkparula,0,255
; loadct,72




 kp = (jm+1)*2 + jk
 print,'KP test ',KP,jm,jk 

 cgplot,findgen(180)-90,SMOOTH(Z_map_all(*,1)+0.5,20,/edge_truncate),yrange=[-2,8],noerase=1, POSITION=pos[*,kp],$
     xtitle ='lat [deg]',ytitle='NET [W/m!u2!n]',xticks=6,xrange=[-90,90],thick=6
     oplot,findgen(180)-90,SMOOTH(Z_map_all(*,0)+0.5,20,/edge_truncate),li=2,thick=6
     
     tek_color    
     cgLegend, Location=[-70, 7.5], colors=[0,0,2],li=[0,2,0], Titles=['Forbush decrease response',$
      'Reference before Forbush decrease','GCR modified MODIS observations '], $
       Length=0.04,charsize=0.8,thick=6,/data,VSpace=1. 
       xytitle = ['mean(dec) = -18.5!uo!n','mean(gcr) =  -12.3%!n','mean(dec) = 19.8!uo!n','mean(gcr) =  -8.2%!n']
  ;   xyouts,-75,6,xytitle(0+Nsign*2), CHARSIZE=1.    
  ;   xyouts,0,6,xytitle(1+Nsign*2) , CHARSIZE=1.
     
     XYOUTS,-85,7,Fig_num(kp) ,Color=cgcolor('black')
      
     close,1
     openr,1,'//Users/hesv/Documents/Work/CERES/cloudresponce_neg_dec_5e2.dat'
     IF Nsign EQ 1 THEN BEGIN
      close,1
      openr,1,'//Users/hesv/Documents/Work/CERES/cloudresponce_pos_dec_3e2.dat'
    
     ENDIF
   
   kp=(kp+2)
      

     lat = findgen(180)-89.5
     Szonal = fltarr(180)
     a=0.
     b=0.
     FOR I = 0,179 DO BEGIN
       readf,1,a,b
       Szonal(I)=b
     ENDFOR
     close,1
     tek_color

     oplot,lat,SMOOTH(Szonal,1,/edge_truncate),col=cgcolor('red'),thick=6
     PRINT,'jmx  out ',kp
 
 ENDFOR ;JK
 
; 0.0250000
; 
mkparula,0,255
;cgloadct,72,/reverse
fac = 0.06
 cgCOLORBAR, RANGE=[-8,8],POSITION=[0.3, 0.00+fac, 0.8,0.025+fac],title=bartit,$
  Tlocation='bottom',Color='charcoal',tickinterval=2
 

ENDFOR ;Mainloop******************************************************************************************
;ENDFOR ;Nsign *******************************************************************************************

 device,/close
  PDFFIG ,'//Users/hesv/Documents/Work/CERES/', 'CERESmapFIG_GRL.ps'
  SPAWN, 'open /Users/hesv/Documents/Work/CERES/CERESmapFIG_GRL.ps'
 ; SPAWN, 'open /Users/hesv/Documents/Work/CERES/CERESmapFIG_GRL.pdf'
  STOP

 
; openw,1,'/Users/hesv/Documents/Work/CERES/AOD_MODIS/MAX_MAP.bin'
; writeu,1,Max_map
; close,1
 
end
function longitude_sum,x,MASK=mask,Invert=invert
  ;PROGRAM TO MAKE longitudnal SUM ON A SPHERE
  ; X  : is the field that is summed X(lon,lat,time)
  ;mask: is a field used to map out certain regions
  ;is has the values 1 or zero
  ;invert : Keyword makes it possible to take the complementary
  ;area of mask
  ;HSV 2018
  Tmask   = keyword_set(mask)
  Tinvert = keyword_set(invert)
  S     = SIZE(X)
  ntime = S(3)
  nlon  = S(1)
  nlat  = S(2)
  foo   = dblarr(nlon,nlat)
  foon  = dblarr(nlon,nlat)
  y     = dblarr(nlon,nlat,ntime)
  yn    = dblarr(nlon,nlat,ntime)
  Y     = X
  z     = fltarr(nlon,Ntime)
  phi = ((FINDGEN(nlat)+0.5D0)*180.D0/nlat -90.0D0)*!PI/180D0     ;LAT ARRAY
  FOR i=0,nlat-1   DO BEGIN  ;Area correction
    y(*,i,*)=Y(*,i,*)*cos(phi(i))*(!pi/2.D0)
    yn(*,i,*)=cos(phi(i))*(!pi/2.D0) ;normalization factor
  ENDFOR 
  IF Tmask then BEGIN
      A=where(mask eq 1,count,complement=A_C) 
      IF Tinvert then A=A_C
      if count eq 0 then BEGIN
        print,'mask error: NO 1'
        STOP
      ENDIF    
     FOR i=0,ntime-1 do begin
       FOR j=0,nlon-1 do begin
          foo  = y(j,*,i)
          foon = yn(j,*,i)  
          z(j,i) = total(foo(A),/nan)/total(foon(A),/nan)
       ENDFOR
     ENDFOR
   ENDIF ELSE BEGIN
      FOR i=0,ntime-1 do begin
       FOR j=0,nlon-1 do begin
         foo  = y(j,*,i)
         foon = yn(j,*,i) 
         z(j,i) = total(foo,/nan)/total(foon,/nan)
       ENDFOR
      ENDFOR 
    ENDELSE
  RETURN,Z
END
;function zonal_sum,x
;  ;HSV
;  ntime = N_ELEMENTS(x(0,0,*))
;  nlon = N_ELEMENTS(x(*,0,0))
;  nlat = N_ELEMENTS(x(0,*,0))
;  print,'ntime=  ',ntime
;  z = fltarr(nlat,Ntime)
;  FOR i=0,ntime-1 do begin
;    FOR j=0,nlat-1 do begin
;      z(j,i) = TOTAL(x(*,j,i))/360.
;    ENDFOR
;  ENDFOR
;  RETURN,Z
;END
FUNCTION zonal_sum0,x,nlon,nlat
  ;print,'***GLOBAL_SUM********************************************
  ;EX = 0 latband around equtor
  ;EX = 1 Complementary (extratropics)
  z = fltarr(nlat)
  
  FOR j=0,nlat-1 DO BEGIN
    foo = x(*,j)  ;Area correction
    A = WHERE(foo GT -100.0,count)
    IF COUNT GE 0 THEN  BEGIN
      Z(J) = MEAN(foo(A))
    ENDIF ELSE BEGIN
      Z(J) = 0.
    ENDELSE
  ENDFOR
 RETURN,Z
END
function zonal_sum,x,MASK=mask,Invert=invert
;PROGRAM TO MAKE ZONAL SUM ON A SPHERE
; X  : is the field that is summed X(lon,lat,time)
;mask: is a field used to map out certain regions 
;is has the values 1 or zero
;invert : Keyword makes it possible to take the complementary 
;area of mask
;HSV 2018 
  S       = SIZE(X)
  ntime   = S(3)
  nlon    = S(1)
  nlat    = S(2)
  y       = dblarr(nlon,nlat,ntime)
  foo     = dblarr(nlon)
  z       = fltarr(nlat,Ntime)
  Y       = X
  Tmask   = keyword_set(mask)
  Tinvert = keyword_set(invert)
  IF Tmask then BEGIN
;      IF Tinvert then A = A_C
;      if count eq 0 then BEGIN
;        print,'mask error: NO 1'
;        STOP
;      ENDIF    
    FOR i=0,ntime-1 do begin
      FOR j=0,nlat-1 do begin
           foo = y(*,j,i)
           A = WHERE(mask(*,j) eq 1)
           IF Tinvert THEN A = WHERE(mask(*,j) ne 1)
        z(j,i) = MEAN(foo(A),/nan)
      ENDFOR
    ENDFOR 
    ENDIF ELSE BEGIN       
      FOR i=0,ntime-1 do begin
       FOR j=0,nlat-1 do begin
        z(j,i) = MEAN(Y(*,j,i))
       ENDFOR
      ENDFOR
    ENDELSE
  RETURN,Z
END
function global_sum,x,MASK=mask,Invert=invert
  ;PROGRAM TO MAKE Global SUM ON A SPHERE
  ; X  : is the field that is summed X(lon,lat,time)
  ;mask: is a field used to map out certain regions
  ;is has the values 1 or zero
  ;invert : Keyword makes it possible to take the complementary
  ;area of mask
  ;HSV 2018
  S     = SIZE(X)
  ntime = S(3)
  nlon  = S(1)
  nlat  = S(2)
  z    = dblarr(Ntime)
  y    = dblarr(nlon,nlat,ntime)
  yn   = dblarr(nlon,nlat,ntime)
  foo  = dblarr(nlon,nlat)
  foon = dblarr(nlon,nlat)
  Y=X
  phi = ((FINDGEN(nlat)+0.5D0)*180.D0/nlat -90.0D0)*!PI/180D0     ;LAT ARRAY
  FOR i=0,nlat-1   DO BEGIN  ;Area correction
    y(*,i,*)=Y(*,i,*)*cos(phi(i))*(!pi/2.D0)  
    yn(*,i,*)=cos(phi(i))*(!pi/2.D0) ;normalization factor
  ENDFOR    
  Tmask   = keyword_set(mask)
  Tinvert = keyword_set(invert)
  IF Tmask then BEGIN
      A=where(mask eq 1,count,complement=A_C) 
      IF Tinvert then A=A_C
    if count eq 0 then BEGIN
      print,'mask error' 
      STOP
    endif
   FOR i=0,ntime-1 do begin
     foo = y(*,*,i)
     foon= yn(*,*,i)
     z(i) = total(foo(A),/nan)/total(foon(A),/nan)
   ENDFOR
  ENDIF else begin
       FOR i=0,ntime-1 do begin
        z(i) = total(y(*,*,i),/nan)/total(yn(*,*,i),/Nan)        
       ENDFOR 
    endelse
  RETURN,Z
END
FUNCTION FD_strength,nfd
  list2 = [119,87, 83, 75, 70, 70, 64, 56, 54, 54, 53, 50, 48, 47, 45, 45, 44, 44, 39, 38, 37, 36, 35, 33, 33, 28]
  IF N_elements(nfd) eq 1 THEN list2 = list2(0:nfd-1)
  IF N_elements(nfd) gt 1 THEN list2 = list2(nfd)   
  RETURN,list2
END
FUNCTION LIN_CONFIDENCE, xi,yi,xc,meanp=mea
R = LINFIT(xi,yi)
delY =yi-POLY(xi,R)
sigma = sqrt(1./(N_elements(xi)-2.)*total(delY^2))
xmean = mean(xi)
varxi2 = total((xi-xmean)^2)
SEyc  = sigma*SQRT(1.+1./N_elements(xi) + (xc-xmean)^2/varxi2)
IF Keyword_Set(mea) THEN  SEyc  = sigma*SQRT(1./N_elements(xi) + (xc-xmean)^2/varxi2)
return,SEyc
end
FUNCTION Zenit, Jday
;calculate solar zenit
month1 = 0
day1 = 0
year1 = 0
CALDAT, Jday, Month1, Day1, Year1
DoY = Jday - julday(1,1,Year1)+1
Zenit = 23.45*sin( (DoY+284.)/365.*2*!pi)
return,Zenit
end
