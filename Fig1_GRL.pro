pro  Fig1_GRL
  ;**************************************
  ; Program: Fig1_GRL
  ; Purpose: Extract monthly radiative forcing files 
  ; and makes a figure + and data files with the forcing 
  ; for each month.
  ; for the Fu-Liou model using MODIS data
  ;**************************************

  ; Set plot output to PostScript format
  SET_PLOT,'ps'
  ; Output directory for plots
  dir = "Dir of the plot output" ;'/Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/'
  ; Input data directory
  datdir ="Dir of the data files produced by Fu-Liou model"; datdir+'/'
  
  ; Initialize device settings for PostScript output
  DEVICE,BITS=8,/PORTRAIT,file=dir+'MODIS_Parameters_FuLiou_MONTHS.ps'
  DEVICE,xsize=20,ysize=20
  DEVICE,xoffset=1,yoffset=2
  DEVICE,SET_FONT='Times-Roman'
  DEVICE,/Color
  !p.font=0
  !P.charsize=0.75
  !P.multi=[0,0,1]
  tek_color
  
  ; Kill existing Skim processes (PDF viewer) to avoid conflicts
  SPAWN, 'killall Skim'
 
   ; Initialize variables for reading and processing data 
  na = 0
  nb = 0
  degGCR = 0.
  swup   = 0.
  olr    = 0.
  a      = 0.
  b      = 0.
  c      = 0.
  d      = 0.
  e      = 0.
  f      = 0.
  g      = 0.
  h      = 0.
  foo    =''

; Define months and grid coordinates
  month = ['!7jan','!7feb','!7mar','!7apr','!7may','!7jun','!7jul','!7aug','!7sep','!7oct','!7nov','!7dec']
  lat = findgen(180)-89.5
  lon = findgen(360)- 179.5
  pos = cgLayout([3,4], OXMargin=[5, 5], OYMargin=[12, 12], XGap=3, YGap=3)
  CLOSE,/all
  ; Generate a Pacific-centered land-sea mask (0=sea, 1=land)
  Pacific_mask =land_sea_mask(360,180,/pacific) ;sea mask for centerd at the Pacific as the MODIS data are
 
  ; Input files from Fu-Liou model        
  FILENAME = [datdir+'/PT_MODIS_RESULT1e2.dat',$ ;str =5 NET =  1.55
              datdir+'/PT_MODIS_RESULT2e2.dat',$
              datdir+'/PT_MODIS_RESULT3e2.dat',$
              datdir+'/PT_MODIS_RESULT4e2.dat',$
              datdir+'/PT_MODIS_RESULT5e2.dat',$
              datdir+'/PT_MODIS_RESULT6e2.dat',$
              datdir+'/PT_MODIS_RESULT7e2.dat',$
              datdir+'/PT_MODIS_RESULT8e2.dat',$
              datdir+'/PT_MODIS_RESULT9e2.dat',$
              datdir+'/PT_MODIS_RESULT10e2.dat',$
              datdir+'/PT_MODIS_RESULT11e2.dat',$
              datdir+'/PT_MODIS_RESULT12e2.dat' ]     


  ;****************************************************************
  ; LOOP OVER MONTHS
  ;****************************************************************
   Nmonth    = 12
   Nstrength = 5 ;0-50. 5 correspond to 11% +-3% and is the FD reference.
  
   Mean_m        = fltarr(12)
   Month_forcing = fltarr(12,360,180)
   zonal         = FLTARR(180,12)
   mu_mean= FLTARR(12)
   
  FOR I=0,Nmonth-1 DO BEGIN
  
  ; Initialize data arrays for this month  
  DATA   = FLTARR(360,180,50,2,12)-999
  mu     = FLTARR(360,180) -999
  Smap   = FLTARR(360,180)-999
  Smap2  = FLTARR(360,180)-999
  SmapN  = FLTARR(360,180)-999
  Smap2N = FLTARR(360,180)-999 
  mask   = INTARR(360,180)
  
  
    
  CLOSE,/ALL  
  OPENR,1,FILENAME(I) ; Open input file for this month
  PRINT,FILENAME(I)
  Nline = FILE_LINES(FILENAME(I)) ;number of lines in file

    FOR J=0,Nline-1 DO BEGIN
    ; Read data line by line
     READF,1, na ,nb, delGCR,swup,olr,a,b,c,d,e,f,g,h,FORMAT='((2(I5,2x),11(F12.6,2x)))'
      Ngcr = FIX((delGCR+0.01)*5) ; Index for GCR strength
      DATA(na,nb,Ngcr,0,0) = swup
      DATA(na,nb,Ngcr,1,0) = olr
      mu(na,nb)            = a
      mask(na,nb)          = 1  
    ENDFOR
  
  ; Compute global mean of mu (solar zenit angle)
  print,       global_sum_mu(mu)
  mu_mean(I) = global_sum_mu(mu)

  ;Scaling of the radiative forcing so S^* <mu> = S_0/4
  fac = 0.41 ;MODIS mu .    Theoretical mu  fac = 0.41
  DATA       = DATA/(4.*mu_mean(I))   ;0.36 ;factor see above line mentioning scaling
  ; Extract reference and strength-dependent SW and LW maps
  Smap(*,*)   = DATA(*,*,        0,0,0) ;REference SW
  Smap2(*,*)  = DATA(*,*,Nstrength,0,0) ;SW at Nstrength GCR
  SmapN(*,*)  = DATA(*,*,        0,1,0) ;Reference LW
  Smap2N(*,*) = DATA(*,*,Nstrength,1,0) ;LW at Nstrength GCR
 ; Mask out land/invalid points
  AZ = WHERE(MASK EQ 0)
  Smap(Az)   = -999
  Smap2(Az)  = -999
  SmapN(Az)  = -999
  Smap2N(Az) = -999
  Forcing    = fltarr(12)

  Smap0  = fltarr(360,180)  
  Smap   = SHIFT(Smap,180,0)
  Smap2  = SHIFT(Smap2,180,0)
  SmapN  = SHIFT(SmapN,180,0)
  Smap2N = SHIFT(Smap2N,180,0)
  
  ;make the forcing***********************
  Smap0   = (Smap-Smap2) + (SmapN-Smap2N)
  
  print,'max SMAP ',i , max(SMAP0,Max_Subscript)
  print,max_subscript,Smap(Max_Subscript),Smap2(Max_Subscript),SmapN(Max_Subscript),Smap2N(Max_Subscript)

  
  
;  PRINT,WHERE(Smap gt 40.)

  foo_map = Smap0
  Am = WHERE(Pacific_mask EQ 0)
  foo_map(Am)     = !VALUES.F_NAN
  

  ;print,MEAN(Smap(A))
  print,'Global Average NET =  ',global_sum(foo_map,360,180)
  mean_M(I) =   global_sum(foo_map,360,180)
  print,'Global Average LW  =  ',global_sum(SmapN-Smap2N,360,180)  
  
  foo_map = Smap0
  Am = WHERE(Pacific_mask EQ 0)
  foo_map(Am)     = !VALUES.F_NAN
  
   ;FIND Outlayers
;   print,STDDEV(foo_map,/nan)
   Bm = where(foo_map gt 10)  
   foo_map(Bm) = !VALUES.F_NAN
;   print,foo_map(Bm)


  Month_forcing(I,*,*) = foo_map(*,*) 
  zonal(*,I) = zonal_sum(foo_map,360,180)
 
 ;plots forcing of each of the 12 month  
  ;loadcv, 72, /noqual,/reverse ;;loads the Matplotlib 'option B' colorbar 
  mkparula,0,255
 ytit=['!7lat','','','!7lat','','','!7lat','','','!7lat','','']
 xtit=['','','','','','','','','','!7lon','','!7lon']
  
  lat = findgen(180)-90 +0.5
  lon = findgen(360) +0.5
 ; loadcv, 72, /noqual,/reverse ;;loads the Matplotlib 'option B' colorbar
  cgplot,lon,lat,/nodata,yrange=[-90,90],xrange=[0,360],xstyle=1,ystyle=1,xtitle=xtit(I),ytitle=ytit(I)$
    ,ticklen=-0.02,yticks=6,xticks=6,NoErase= 1,Position=pos(*,I)

  cgMap_Set,0,180, /Cylindrical, /NoBorder, Limit=[-90, 0, 90, 360], Position=pos(*,I),noerase=1;
  cgimage,bytscl(Smap0,min=0,max=5),lon,lat,/overplot,noerase=1;/Save
 ; cgcontour,SHIFT(SMOOTH(Smap0,[10,10],/EDGE_WRAP),180,0),lon,lat, LEVELS=([0.5,1,2,4]),/overplot, /OnImage,C_Color=cgcolor('white'),charsize=0.7
  cgMap_Continents, Color='dark gray', /Fill
  ;cgMap_Continents, Color='black'
  ; cgMap_Grid,/Box_Axes;,BTHICK=0.0 
  xyouts,60,52,month(I),CHARSIZE=0.8,col=cgcolor('white')
  xyouts,32,32,'NET = '+STRING(global_sum(foo_map,360,180),FORMAT='(F3.1)')+'W/m!u2!n',CHARSIZE=0.6,col=cgcolor('white')
  
  
 
  CLOSE,/all
  
 ENDFOR
 fac =0.07
  cgcolorbar,POSITION=[0.4,0.02+fac,0.6,0.04+fac],title='!7Net absorbed radiation [Wm!u-2!n]',RANGE=[0,5],charsize=1.2 
 print,'Global Ocean Annual Average NET =  ',mean(mean_M) 

  PRINT,'mu_mean = ',MEAN(mu_mean)  
  CLOSE,/all
  DEVICE,/CLOSE
  PDFFIG ,dir, 'MODIS_Parameters_FuLiou_MONTHS.ps'
  
  ;Write a file with the zonal forcing. 
  close,1
  openw,1,dir+'zonal.dat'
  PRINTF,1,zonal
  close,1
 ;Write a file with the forcing shown in the plot. 
  close,1
  openw,1,dir+'Month_forcing4theory.dat'
  PRINTF,1,Month_forcing
  close,1
  
  SPAWN, 'open /Users/hesv/Documents/Work/Ed4_LaRC_FuLiou/idl/MODIS_Parameters_FuLiou_MONTHS.ps'

  STOP

END
FUNCTION global_sum,x,nlon,nlat
  ;print,'***GLOBAL_SUM********************************************
  ;EX = 0 latband around equtor
  ;EX = 1 Complementary (extratropics)
  z = 0.0d
  nom =0.0D0
  lat = (findgen(180)-89.5)/180.*!pi
  FOR i=0,nlon-1 DO BEGIN
    FOR j=0,nlat-1 DO BEGIN  ;Area correction
      IF FINITE(X(I,J)) THEN BEGIN
       z= z+x(i,J)*cos(lat(j))
       nom= nom +  cos(lat(j))
      ENDIF
    ENDFOR
  ENDFOR
  z = Z/nom
  RETURN,Z
  END
  FUNCTION zonal_sum,x,nlon,nlat
    ;=======================================================
    ; FUNCTION: global_sum
    ; Computes area-weighted global average of a 2D array
    ;=======================================================

  ;print,'***GLOBAL_SUM********************************************
  z = fltarr(nlat)

    FOR j=0,nlat-1 DO BEGIN  ;Area correction
      foo = x(*,j)
      A = WHERE(FINITE(foo),count)
       IF COUNT GT 0 THEN Z(J) = MEAN(foo(A))

    ENDFOR
 
  RETURN,Z
END
FUNCTION global_sum_mu,x
  ;=======================================================
  ; FUNCTION: global_sum_mu
  ; Computes area-weighted global mean for positive mu values
  ;=======================================================
  ;print,'***GLOBAL_SUM********************************************
  nlon = 360
  nlat = 180
  z = 0.0d
  nom =0.0D0
  lat = (findgen(180)-89.5)/180.*!pi
  FOR i=0,nlon-1 DO BEGIN
    FOR j=0,nlat-1 DO BEGIN  ;Area correction
      IF X(I,J) GT 0 THEN BEGIN
        z= z+x(i,J)*cos(lat(j))
        nom= nom +  cos(lat(j))
      ENDIF
    ENDFOR
  ENDFOR
  z = Z/nom
  RETURN, z
  END