pro FIG3_GRL
  ;=================================================
  ; FIG3_GRL: Plot Forbush Decrease (FD) Strength vs CERES SW TOA Radiation
  ;=================================================

  ;-----------------------------------------------
  ; SET UP PLOTTING DEVICE
  ;-----------------------------------------------
  SET_PLOT,'PS'
  DEVICE,bits=8,file='/Users/hesv/Documents/Work/CERES/GRLSvensmark/FIG3_GRL.ps'
  DEVICE,xsize=20,ysize=16
  DEVICE,xoffset=0.,yoffset=6
  DEVICE,/Color
  ;PLOTS NEW LINEAR FIGURE OF FD-STRENGTH AND CERES SW TOA RADIATION (2025)
  ; Close any previous PS viewer (Skim)
  SPAWN, 'killall Skim'
  ;-----------------------------------------------
  ; DEFINE DIRECTORIES
  ;-----------------------------------------------  
  dir = '/Users/hesv/Documents/Work/CERES/GRLSvensmark/'
  datdir = '/Users/hesv/Documents/Work/CERES/CERES2026/'
  
  
;***********************************************
;Read CERES DATA FILE NAMES
  CD,DATDIR
  name =  FILE_SEARCH('*.nc') 
;*********************************************** 
;-----------------------------------------------
; SET PLOTTING SYMBOL PARAMETERS
;-----------------------------------------------
  !p.charsize  = 1.
  !p.charthick = 1.0
  !x.thick     = 6.0
  !y.thick     = 6.0
  !p.thick     = 6.0
  !p.font      = 6
  !p.multi=[0,0,1]
  ;Make plot symbol (circle)
  amp = 1.5
  t = findgen(21)/20.*360.
  X =  amp*sin(t/180.*!pi)
  Y =  amp*cos(t/180.*!pi)
  USERSYM, X, Y ,/fill
  
  ;-----------------------------------------------
  ; SET PANEL LAYOUT
  ;-----------------------------------------------
  pos = cgLayout([3,1], OXMargin=[10, 10], OYMargin=[10, 3], XGap=4, YGap=3)
  DELx = 0.2
  del = pos(0,2)-pos(2,1)
  dx  = (pos(2,2)-pos(0,1)-delx -del)/2.
  pos(2,0)  = pos(2,0)+DELx                    ;X1 slut 1 fig
  pos(0,1)  = pos(0,1)+DELx                    ;X2 begynd 2 figur

  pos(2,1)  = pos(0,1)+dx       ;x3 slut 2 figur
  pos(0,2)  = pos(0,1)+dx+del           ;x4 begynd 3 figur
  ;x5 slut 3 figur

  print,pos(2,1)-pos(0,1),pos(2,2)-pos(0,2)

  ;-----------------------------------------------
  ; PARAMETERS FOR FORBUSH DECREASE (FD)
  ;-----------------------------------------------
  nst = 15 ;Day of minimum in FD
  ned = 25 ;End of interval where looking for a maximum
 
  ; Pre-allocate array for CERES SW TOA data 
  Fsum = fltarr(916+916+916+916+916+916+916+916+916+675+334)
  nt =[0,916,916,916,916,916,916,916,916,916,675,334]
  
  ntsum=FIX(total(nt,/cumulative))
  FOR JK=0,0 DO BEGIN

;READ CERES DATA FILES

    FOR K=0,10 DO BEGIN ;number of *.nc files
      ;-----------------------------------------------
      ; READ CERES NETCDF FILES
      ;-----------------------------------------------      
      ncid = NCDF_OPEN(datdir+name(k))

      filename=datdir+name(k)
      nstru = NCDF_INQUIRE(ncid) ;getting information of the number of variables
      print, nstru

      dims = lonarr(nstru.ndims)
      names = strarr(nstru.ndims)
      for i=0L,nstru.ndims-1 do begin  ;names of the dimensions
        NCDF_DIMINQ,ncid,i,dname,dsize
        names[i] = dname
        dims[i] = dsize
        print,'Dimension '+strtrim(i,2)+': ',dname,dsize
      endfor
      for i=0L,nstru.nvars-1 do begin  ;getting names of variables
        vardesc = NCDF_VARINQ(ncid, i)
        print,vardesc
      endfor



     varid3 = 'toa_sw_all_daily'

      NCDF_VARGET, ncid,  varid3, var_name3
    
      ;Close the netcdf file
      NCDF_CLOSE, ncid
      help,var_name3

      ;-----------------------------------------------
      ; APPLY CLOUD AND OCEAN MASKS
      ;-----------------------------------------------
      ;Which mask to use on the CERES map? Here using Liquid cloud fraction and ocean.
        Imask=0 
        nz =0
        FnameCloud_mask=['Liquid_cloud_fraction','High_cloud_fraction','Ice_cloud_fraction','global']
        
        noland = ['','_noland'] ;nz =   0 land included, 1 land excluded
        Fname = STRCOMPRESS(dir+'MASK_'+FnameCloud_mask(Imask)+noland(nz)+'.dat')
 
        cloud_mask = read_binary(Fname,data_type=2, data_dims=[360,180])
   
        Mask_sea = land_sea_mask(360,180,/pacific)
        cloud_mask=mask_sea*cloud_mask
        A = WHERE(cloud_mask eq 0)
        cloud_mask(A) = -999

        ZZ = -global_sum(var_name3,mask=cloud_mask)*0.42/0.71 ;0.42 fraction of oceans with liquid clouds (MODIS)
                                                              ;0.71 fraction of Earth with ocean

        Fsum(Ntsum(k):Ntsum(k)+Nt(k+1)-1) =  ZZ(*)

    ENDFOR ;READING CERES FILES
    
    ;-----------------------------------------------
    ; DEFINE FORBUSH DECREASE DATES (JJd)
    ;-----------------------------------------------
    ; Forbush dates
    J0  =  julday(3,1,2000)
    JJd = [julday(10,31,2003),$;1  16.2
           julday(1, 19,2005),$;2  12.3
           julday(9, 13,2005),$;3   8.8 
           julday(7, 16,2000),$;4   9.5
           julday(4, 12,2001),$;5   8.9
           julday(11,10,2004),$;6   8.0
           julday(9, 26,2001),$;7   4.9
           julday(7, 17,2005),$;8   6.3
           julday(7, 27,2004),$;9   6.8
           julday(5, 31,2003),$;10  6.1
           julday(11,25,2001),$;11  6.3
           julday(5, 15,2005),$;12  5.7      
           julday(8, 28,2001),$;13  4.4                    
           julday(12,15,2006),$;14  7.1 
           julday(3,  9,2012),$;15  9.7
           julday(6, 25,2015)] ;16  6.7

 
    Nfd = N_elements(JJd)     

    Sum = fltarr(36)
    

    ;************************************************************
    ;WRITE THE LONG TIMESERIES DATA TO FILE
    Jdays  =  J0 + findgen(N_elements(Fsum))
    caldat,Jdays,D,M,Y
    close,1
    openw,1,dir+STRcompress('long_'+varid3+'_'+FnameCloud_mask(Imask)+'HSV.dat')
    FOR i=0,N_elements(Fsum)-1 DO printf,1,Y(i),M(I),D(I), Fsum(I)
    CLOSE,1
    ;*************************************************************
    ;
    ;-----------------------------------------------
    ; ALLOCATE ARRAYS FOR FD ANALYSIS
    ;-----------------------------------------------
    Fdm = fltarr(Nfd)
    FdmV = fltarr(Nfd)
    
    ;-----------------------------------------------
    ; CALCULATE FD MAX RELATIVE TO REFERENCE
    ;-----------------------------------------------
    FOR I=0,Nfd-1 DO BEGIN
      Fdmax = FSUM(JJd(I)-J0-15L:JJd(I)-J0+20L)
      R = linfit(findgen(36)-15,Fdmax)
      Fdmax   = Fdmax -POLY(findgen(36)-15,R)
      Fdmax   = SMOOTH(Fdmax,5,/edge_truncate)
      FDm(I)  = MAX((FDmax(nst:ned))) -  MEAN((FDmax(5:15))) ;10 day finding max SW - 10 day reference period 
      IF JK EQ 1 THEN FDm(I)  = MIN((FDmax(nst:ned)))
      FDmV(I) = stddev(FDmax(0:10))
    ENDFOR
    ;-----------------------------------------------
    ; PLOT CERES FD STRENGTH VS SOLAR CYCLE
    ;-----------------------------------------------
    st = Total(FdmV)
    Fd_error = replicate(1.0,Nfd)

     NFDsub = findgen(Nfd)
    
    yrng = [[-.5,4.4],[-6,0]]
    

    ytit= ['!9D!7SW!ITOA!n [W/m!u2!n]','!9D!7SW!ITOA!n [W/m!u2!n]','!9D!7LW!ITOA!n [W/m!u2!n]']
cgplot,FD_strength(NFDsub),Fdm,psym=8,xrange=[30,130],xstyle=1,xtitle='Strength relative to solar cycle [%]',POSITION=pos(*,0)$
      ,ytitle=ytit(jk),xthick=5,NoErase=jk NE 0 ,yrange=yrng(*,jk),yTICKLEN=-.02,xTICKLEN=-.01,/nodata,yminor=5,xminor=2


    del =0.01
    POLYFILL,[30+del*10,130-del*10,130-del*10,30+del*10],[-0.5+del,-0.5+del,0,0] ,col=cgcolor('light gray') 
    POLYFILL,[100-5,100+5,100+5,100-5],[-0.5+del,-0.5+del,4.4-del,4.4-del] ,col=cgcolor('light gray') 
     xyouts,35,4.2,'Liquid Clouds'
 
    amp = 1.5
    t = findgen(21)/20.*360.
    X =  amp*sin(t/180.*!pi)
    Y =  amp*cos(t/180.*!pi)
    USERSYM, X, Y

    amp = 1.5
    t = findgen(21)/20.*360.
    X =  amp*sin(t/180.*!pi)
    Y =  amp*cos(t/180.*!pi)
    USERSYM, X, Y,/fill

    print,FD_strength(Nfdsub)
    
    AI =REVERSE(SORT(FD_strength(Nfdsub)))
    
    PRINT,FD_strength(AI)
    print,AI
    XX = FD_strength(AI)
    
    FOR I=0,Nfd-1 DO BEGIN
      caldat,JJd(AI(I)),D,M,Y
      print,I+1,Y,M,D,xx(i)
    ENDFOR


    oplot,[10,150],[-2.07,-2.07],li=3

    Foo1=0.0
    Foo2=0.0
    R = linfit(FD_strength(AI),Fdm(AI),sigma=Rs,CHISQR=foo2)
    foo1 = r(1) ;slope
    foo2 = r(0) ;Intersept
    Q1  = fltarr(3)
    Q2  = fltarr(3)
    Q1(0)=foo2 ;intersept
    Q2(0)=foo1 ;Slope
    oplot,findgen(200)+1,POLY((findgen(200)+1),R)


    oplot,findgen(200)+1, LIN_CONFIDENCE((FD_strength(AI)),Fdm(AI),(findgen(200)+1))+POLY((findgen(200)+1),R),li=1
    oplot,findgen(200)+1,-LIN_CONFIDENCE((FD_strength(AI)),Fdm(AI),(findgen(200)+1))+POLY((findgen(200)+1),R),li=1

    oplot,findgen(200)+1, LIN_CONFIDENCE((FD_strength(AI)),Fdm(AI),(findgen(200)+1),/meanp)+POLY((findgen(200)+1),R),li=2
    oplot,findgen(200)+1,-LIN_CONFIDENCE((FD_strength(AI)),Fdm(AI),(findgen(200)+1),/meanp)+POLY((findgen(200)+1),R),li=2
;  PLOT POINTS WITH NUMBERS
    oplot,FD_strength(AI),Fdm(AI),psym=8,col=cgcolor('blue')
    xyouts,FD_strength(AI)-1.2,Fdm(AI)-0.03,STRING(findgen(Nfd)+1,FORMAT='(I2)'),charsize=0.6,col=cgcolor('white')
  ENDFOR
  

  ;PLOTTING EMPIRICAL RESULTS FROM OCEAN HEAT CONTENT, CERES FD. 
   ;Fu_liouMODIS =[0.0,0.39,0.79,1.2,1.6,2.0]
   Fu_liouMODIS    =[0.30,0.61,0.92,1.22,1.55]  ;NET   W/m2 ;From Fu-Liou simulations using MODIS DATA
   Fu_liouMODIS_SW =[0.32,0.66,0.99,1.33,1.67]  ;NET   W/m2
   Fu_liouMODIS_STR=[20,40,60,80,100]           ;STRENGTH OF GCR reltive to solar cycle in %
   CLOUD11_TSI     =[1.3]                       ;W/m2
   Ocean_F         =[1.05,1.68,1.15,1.33]       ;Heatflux into the ocean over a solar cycle SHAVIV 2008
   delOcean        =[0.25,0.6,0.35,0.34]        ;Variance of ocean forcing

    fac = 1.

   OPLOT,Fu_liouMODIS_STR,Fu_liouMODIS_SW,psym=2,symsize=1.2,col=cgcolor('cyan')  
   OPLOT,Fu_liouMODIS_STR,Fu_liouMODIS,psym=2,symsize=1.2,col=cgcolor('red')  
   ;OPLOT,[100],CLOUD11_TSI,psym=7,symsize=1,col=cgcolor('cyan')
   cgLegend, Location=[40, 4.0], colors=[cgcolor('blue'),cgcolor('red'),cgcolor('cyan')],PSYMS=[8,2,2], $
     Titles=['CERES SW of 16 Forbush decreases','Fu-Luio (MODIS) NET','Fu-Luio (MODIS) SW'], $
     Length=0.0,charsize=1.,/data,VSpace=2.,symsize=1.2
;
;

;
;  ;Ocean heat content  
   OPLOT,replicate(100.05,1),[Ocean_F(0)]*fac,psym=4, col=cgcolor('orange'),symsize=1
   ERRplot,replicate(100.05,1),[Ocean_F(0)-delOcean(0)]*fac,[Ocean_F(0)+delOcean(0)]*fac,col=cgcolor('orange')
   OPLOT,replicate(102.05,1),[Ocean_F(1)]*fac,psym=5, col=cgcolor('Navy'),symsize=1
   ERRplot,replicate(102.05,1),[Ocean_F(1)-delOcean(1)]*fac,[Ocean_F(1)+delOcean(1)]*fac,col=cgcolor('Navy')
   OPLOT,replicate(100.05-3,1),[Ocean_F(2)]*fac,psym=6, col=cgcolor('Olive'),symsize=1
   ERRplot,replicate(100.05-3,1),[Ocean_F(2)-delOcean(2)]*fac,[Ocean_F(2)+delOcean(2)]*fac,col=cgcolor('Olive')
   OPLOT,replicate(100.05+4,1),[Ocean_F(3)]*fac,psym=7, col=cgcolor('Brown'),symsize=1
   ERRplot,replicate(100.05+4,1),[Ocean_F(3)-delOcean(3)]*fac,[Ocean_F(3)+delOcean(3)]*fac,col=cgcolor('Brown')
  
   cgLegend, Location=[40, 2.8], colors=[cgcolor('orange'),cgcolor('Navy'),cgcolor('Olive'),cgcolor('Brown')],PSYMS=[4,5,6,7], $
     Titles=['OHC (sc)','SST (sc)','SLR (sc)','Satellite (sc)' ], $
     Length=0.0,charsize=0.8,/data,VSpace=1.5
  
  
   ;-----------------------------------------------
   ; MONTE CARLO DISTRIBUTION FUNCTIONS
   ;-----------------------------------------------
  ;;;;;;;;;;;;;;;DISTRIBUTION FUNCTIONS;;;;;;;;;;;;
  !p.charsize  = 1.
  Nsamples  = 10000
  X = fltarr(Nsamples)
  Y = fltarr(Nsamples)
  ;find start
  Month_k  = 3
  Day_k    = 1
  Year_k   = 2000
  ;find end
  Month_e  = 10
  Day_e    = 10
  Year_e   = 2015

  Var_10day = fltarr(5000L*13)


  ;START and END DATE
  J0    =  julday(Month_k, Day_k,Year_k)
  Jend  =  julday(Month_e, Day_e,Year_e)
  ;MC SAMPLING LINEAR REGRESSION
  FOR K=0,Nsamples-1 DO BEGIN
    seed=k

    JJd = J0+16 + FIX((Jend-J0-40)*randomu(seed,Nfd))

    FOR I=0,Nfd-1 DO BEGIN
      Fdmax = FSUM(JJd(I)-J0-15L:JJd(I)-J0+20L)
      R = linfit(findgen(36)-15,Fdmax)
      Fdmax = Fdmax -POLY(findgen(36)-15,R)
      Fdmax = SMOOTH(Fdmax,5,/edge_truncate)
      FDm(I)  = MAX((FDmax(nst:ned))) -  MEAN((FDmax(5:15)))
    ENDFOR

    foo=0
    R = linfit(FD_strength(NFDsub),Fdm)
   
    X(k)=r(1) ;Slope
    Y(k)=r(0) ;Intersept
  ENDFOR
  Q1(1)=mean(y)
  Q2(1)=mean(x)
  Q1(2)=stddev(y)
  Q2(2)=stddev(x)
  print,'intercept ',foo2,mean(y),stddev(y)
  print, 'slope    ' ,foo1,mean(x),stddev(x)
  foo1= (foo1-mean(x))/stddev(x) ;slope in units of var
  X = (x-mean(x))/stddev(x)
  foo2= (foo2-mean(y))/stddev(y) ;intersept in units of var
  y = (y-mean(y))/stddev(y)

  ;SET_PLOT,'PS'
  panel = ['','','']
  tek_color
  H = HISTOGRAM(X,Binsize=0.2,MIN=-3,MAX=3)/float(Nsamples)
  cgplot,H,(findgen(N_elements(H))/(N_elements(H)+1)-0.5)*5 ,position=pos(*,1),NoErase=jk,xtickinterval =0.05,yrange=[-6,6],ytickinterval=1$
    ,/nodata,ystyle=1,xcharsize=1.

  xyouts,0.01,foo1+.1,'Slope'
  xyouts,0.01,-4,'A!iS!n= '+strcompress(number_formatter(Q2(0),decimals=2),/remove_all)
  xyouts,0.01,-4.5,'<A>!iMC!n= '+strcompress(number_formatter(0.0,decimals=1),/remove_all)
  xyouts,0.01,-5,'!9s!7!IMC!n= '+strcompress(number_formatter(Q2(2),decimals=2),/remove_all)

  POLYFILL,H,(findgen(N_elements(H))/(N_elements(H)+1)-0.5)*5  ,col=cgcolor('blue')

  oplot,H,(findgen(N_elements(H))/(N_elements(H)+1)-0.5)*5,thick=3,col=cgcolor('blue')

  oplot,[0,nsamples],[foo1,foo1],li=2,col=cgcolor('blue')

  H = HISTOGRAM(y,Binsize=.2,MIN=-3,MAX=3)/float(Nsamples)
  cgplot,H,(findgen(N_elements(H))/(N_elements(H)+1)-0.5)*5 ,position=pos(*,2),yrange=[-6,6],NoErase=jk,xtickinterval =0.05,ytickinterval=1$
    ,/nodata,ystyle=9,xcharsize=1.,$
    title='',xrange=[0,0.1],xstyle=1,$
    xtitle=panel(2),thick=3,$
    ytitle='',xmargin=[16,20],xticks=3,yTICKLEN=-.05
  AXIS,Yaxis=1,ytitle='!7Signal/variance',yRANGE=[-6,6]
  POLYFILL,H,(findgen(N_elements(H))/(N_elements(H)+1)-0.5)*5 ,col=cgcolor('blue')
  tek_color
  oplot,H,(findgen(N_elements(H))/(N_elements(H)+1)-0.5)*5,thick=3,col=cgcolor('blue')
  xyouts,0.01,foo2+.1,'Intercept'
  xyouts,0.01,5,'B!iS!n='+strcompress(number_formatter(Q1(0),decimals=2),/remove_all)
  xyouts,0.01,4.5,'<B>!iMC!n= '+strcompress(number_formatter(Q1(1),decimals=1),/remove_all)
  xyouts,0.01,4,'!c!9s!7!IMC!n= '+strcompress(number_formatter(Q1(2),decimals=2),/remove_all)

  xyouts,-0.12,6.4,'MC Distribution Functions'
  ;
  oplot,[0,nsamples],[foo2,foo2],li=2,col=cgcolor('blue')


  print,mean(Var_10day),stddev(Var_10day)
  
  ;-----------------------------------------------
  ; CLOSE DEVICE AND OPEN OUTPUT PS FILE
  ;-----------------------------------------------
  device,/close
  
  SPAWN, 'open /Users/hesv/Documents/Work/CERES/GRLSvensmark/FIG3_GRL.ps'  
  stop
END
  
  
FUNCTION global_sum,x,MASK=mask,Invert=invert
  ;-----------------------------------------------
  ; FUNCTION: global_sum
  ; Computes global sum on a sphere with optional mask
  ;-----------------------------------------------
  ;PROGRAM TO MAKE Global SUM ON A SPHERE
  ; X  : is the field that is summed X(lon,lat,time)
  ;mask: is a field used to map out certain regions
  ;with the values 1 or zero
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
  phi = ((FINDGEN(nlat)+0.5D0)*180.D0/nlat -90.0D0)*!PI/180.D0     ;LAT ARRAY
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
  ;-----------------------------------------------
  ; FUNCTION: FD_strength
  ; Returns estimated FD strength relative to GCR solar cycle
  ;-----------------------------------------------
;Estimated strength of FD relative to GCR solar cycle variation  
  list2=[126.99,115.60,80.54,88.50,96.16,54.49,41.33,49.51,61.01,55.54,57.76,36.55,50.75,57.05,105.56,55.03]; 2025 0.5-31 km
  IF N_elements(nfd) eq 1 THEN list2 = list2(0:nfd-1)
  IF N_elements(nfd) gt 1 THEN list2 = list2(nfd)
  RETURN,list2
END
FUNCTION LIN_CONFIDENCE, xi,yi,xc,meanp=mea
  ;-----------------------------------------------
  ; FUNCTION: LIN_CONFIDENCE
  ; Computes linear regression confidence intervals
  ;-----------------------------------------------
  ;Find linear confidence of fitted line
  R = LINFIT(xi,yi)
  delY =yi-POLY(xi,R)
  sigma = sqrt(1./(N_elements(xi)-2.)*total(delY^2))
  xmean = mean(xi)
  varxi2 = total((xi-xmean)^2)
  SEyc  = sigma*SQRT(1.+1./N_elements(xi) + (xc-xmean)^2/varxi2)
  IF Keyword_Set(mea) THEN  SEyc  = sigma*SQRT(1./N_elements(xi) + (xc-xmean)^2/varxi2)
  return,SEyc
end