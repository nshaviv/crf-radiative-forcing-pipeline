pro FIG2_GRL
  ;*****************************************************************************
  ; PROGRAM THAT MAKES FIGURE 2
  ;
  ; This program:
  ;  - Reads CERES SYN1deg daily TOA radiation data (Net, SW, LW)
  ;  - Composites the data around selected Forbush Decrease (FD) events
  ;  - Separates cases of negative and positive solar declination
  ;  - Averages pre- and post-FD periods
  ;  - Produces Figure 2 as a PostScript file
  ;
  ; External data used:
  ;  - CERES netCDF files
  ;  - Fu-Liou model output files:
  ;       cloudresponce_neg_dec_5e2.dat
  ;       cloudresponce_pos_dec_3e2.dat
  ;  - MODIS-based cloud input parameters (used in Fu-Liou simulations)
  ;
  ; Output:
  ;  - FIG2_GRL.ps
  ;*****************************************************************************
  
  ;-------------------------------------------------------------
  ; Define working directory and output PostScript device
  ;-------------------------------------------------------------
DIR = '//Users/hesv/Documents/Work/CERES/GRLSvensmark/'
SET_PLOT,'ps'
DEVICE,bits=8,file=DIR+'FIG2_GRL.ps'
DEVICE,xsize=20,ysize=24
DEVICE,xoffset=0,yoffset=3
DEVICE,/Color

;-------------------------------------------------------------
; Plot appearance settings
;-------------------------------------------------------------
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
;-------------------------------------------------------------
; Define a custom circular plot symbol
;------------------------------------------------------------- 
;Make plot symbol (circle)
amp = 1.5
t = findgen(21)/20.*360.
X =  amp*sin(t/180.*!pi)
Y =  amp*cos(t/180.*!pi)
USERSYM, X, Y ,/fill

mkviridis

;-------------------------------------------------------------
; Allocate arrays
;-------------------------------------------------------------
  Z_map_all = fltarr(180,6)
  
;-------------------------------------------------------------
; CERES file names
;-------------------------------------------------------------
 ;READ CERES FILES
  dirC = '/Users/hesv/Documents/Work/CERES/'
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

;-------------------------------------------------------------
; Time and array initialization
;-------------------------------------------------------------
;Open the netcdf file,  ncid is an idl string
nday    = 36   ;Total number of days
mdel =  0     ; additional days before day -15
DEL  = Nday-36L ;fixing to make longer time intervals than 36 days 

NFdays = [640,639,639,245,916,916,916,916,323]
Ndays_total =  TOTAL(NFdays, /INTEGER)
Fsum    = fltarr(Ndays_total)
Fsum_z  = fltarr(180,Ndays_total)
Fsum_lon  = fltarr(360,Ndays_total)
Fmap    = fltarr(360,180,Ndays_total)
SUM_map = fltarr(360,180,nday)

nt      = [0,640,639,639,245,916,916,916,916,323]
NFdays  =   [640,639,639,245,916,916,916,916,323]

ntsum=FIX(total(nt,/cumulative)) 


Pacific_mask =land_sea_mask(360,180,/pacific) ;sea mask for centerd at the Pacific as the CERES data are

;*************************************************************
; READ CERES FILES
;*************************************************************
FOR K=0,8 DO BEGIN ;LOOP FOR OPENING CERES FILES

ncid = NCDF_OPEN(dirC+name(k),/nowrite)

  filename=dirC+name(k)
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

varid3 = PARM(0) ;NET
out_dat = dir + STRCOMPRESS(varid3+'_2000_2005_zonal.dat')  ;data out file

NCDF_VARGET, ncid,  varid3, var_name3
help,var_name3

;Close the netcdf file
NCDF_CLOSE, ncid

FOR ai=0,179 DO BEGIN
  FOR aj=0,359 DO BEGIN
    Fmap(aj,ai,Ntsum(k):Ntsum(k)+Nt(k+1)-1) =  var_name3(aj,ai,*)
  ENDFOR
ENDFOR


ENDFOR ; END READING CERES FILES LOOP******************************************************
;

print,max(Fmap),MIN(Fmap)
;-------------------------------------------------------------
; Compute zonal means and area weighting
;-------------------------------------------------------------
;ZZ    = global_sum(Fmap)
Zonal = zonal_sum(Fmap);,mask=Pacific_mask);,/invert)

co    = cos(!pi*(findgen(180)-90)/180.)
Zonal_area = fltarr(180,Ndays_total )

For i=0,179 do Zonal_area(i,*) = Zonal(i,*)*co(i)
Fsum(*)    =    total(Zonal_area(30:105,*),1)/85.

FOR i=0,179 do Fsum_z(i,*) =    Zonal(i,*)

NFIG =0

;*************************************************************
; LOOP OVER NEGATIVE AND POSITIVE DECLINATION CASES
;*************************************************************
FOR jk = 0,1 DO BEGIN

  kp  = jk

  J0  =  julday(3,1,2000)
       
  Nfd = 5
Nsign = jk ;0 negativ declination.  1 positiv declination
;-----------------------------------------------------------
; Define Forbush Decrease dates
;-----------------------------------------------------------     
  JJd = [julday(10,31,2003),$;7 negativ dec -15
         julday(1, 19,2005),$;10             -20
         julday(11,10,2004),$;9             -18
         julday(11,25,2001),$;3             -21
         julday(12,15,2006)];14             -23


       
   IF Nsign EQ 1 THEN BEGIN      
      Nfd   = 6
;   positiv declination

       JJd = [$  ;positiv dec (all 6)
         julday(7, 16,2000),$;1  21
         julday(7, 17,2005),$;12 21
         julday(7, 27,2004),$;8  18
         julday(5, 31,2003),$;6  21
         julday(5, 15,2005),$;11 18
         julday(6, 25,2015) $;16 23. 
          ];12             
   ENDIF


 ; Initialize composite arrays
Fs   = 0    
Sum   = fltarr(nday)
Sum_z = fltarr(180,nday)
Sum_lon = fltarr(360,nday)
Sum_zs = fltarr(180,nday)
; Neutron monitor data
Sum_ou    = fltarr(nday)
Del_gcr = fltarr(Nfd)

;GLOBAL CERES DATA
print,JJd(0)-J0-15L-mdel,JJd(0)-J0+20L+DEL-mdel

;-----------------------------------------------------------
; Build composite centered on FD events
;-----------------------------------------------------------
    FOR I=0,Nfd-1 DO BEGIN
      FOR J=0,179 DO BEGIN                                                                        
        FOR k=0,359 DO SUM_map(k,j,*) = SUM_map(k,j,*) + Fmap(k,j,JJd(I+Fs)-J0-15L-mdel:JJd(I+Fs)-J0+20L+DEL-mdel)     
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

;-----------------------------------------------------------
; Remove linear trend at each grid point
;-----------------------------------------------------------
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
;***********************************************************
; Plot maps and zonal means
;***********************************************************

; Layout (2x3 figure)
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


; Average before/after FD
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

  
  Pacific_mask(*,0)=0. ;First row set to zero as it should be 
  A = WHERE(Pacific_mask EQ 0)
  Foo_map(A)     = -999
  Z_map_all(*,jm)= zonal_sum0(foo_map,360,180)
  A = WHERE(Z_map_all(*,jm) LT -900)
  Z_map_all(A,jm) = 0.

;colar tables
mkparula,0,255
    
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

 ;mkviridis
 mkparula,0,255




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
     
     XYOUTS,-85,7,Fig_num(kp) ,Color=cgcolor('black')
   
   ;READ FILES WITH DATA FROM FU-LIOU SIMULATION OF MONTH WITH NEGATIVE AND POSITIVE DECLINATION 
     close,1
     openr,1,DIR+'cloudresponce_neg_dec_5e2.dat' ; strength as solar cycle:5
     IF Nsign EQ 1 THEN BEGIN
      close,1
      openr,1,DIR+'cloudresponce_pos_dec_3e2.dat' ; 60% of solar cycle:3
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

mkparula,0,255

 fac = 0.06
 cgCOLORBAR, RANGE=[-8,8],POSITION=[0.3, 0.00+fac, 0.8,0.025+fac],title=bartit,$
  Tlocation='bottom',Color='charcoal',tickinterval=2
 
 device,/close
;OPEN AND DISPLAY THE FIGURE ON SCREEN
  SPAWN, 'open /Users/hesv/Documents/Work/CERES/GRLSvensmark/FIG2_GRL.ps'
 ; SPAWN, 'open /Users/hesv/Documents/Work/CERES/CERESmapFIG_GRL.pdf'
  STOP

 
end

FUNCTION zonal_sum0,x,nlon,nlat
  ;---------------------------------------------------------
  ; Computes zonal mean of 3D field X(lon,lat,time)
  ;
  ; Optional:
  ;  MASK   : 1/0 mask
  ;  INVERT : use complementary mask
  ;---------------------------------------------------------
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
  ;---------------------------------------------------------
  ; Computes zonal mean of a 2D field (lon,lat)
  ; Missing values assumed < -100
  ;---------------------------------------------------------
;PROGRAM TO MAKE ZONAL SUM ON A SPHERE
; X  : is the field that is summed X(lon,lat,time)
;mask: is a field used to map out certain regions 
;it has the values 1 or zero
;invert : Keyword makes it possible to take the complementary 
;area of mask
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

