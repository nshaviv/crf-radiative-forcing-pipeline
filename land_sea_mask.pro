FUNCTION LAND_SEA_MASK,nlon,nlat,PACIFIC=pacific,Land=land
;Function that gives a land_sea_mask 1 for sea and 0 for land
;The center for the map can be changed to the pacific basin by the pacific key_word
;The land keyword makes the map euqnl to 
    mask = intarr(360,180)
     mask=read_binary('/Users/hesv/Documents/Work/CERES/seamask.bin', data_type=2, data_dims=[360,180])
 
    Zero_opt = keyword_set(pacific)
    tlon     = keyword_set(nlon)
    tlat     = keyword_set(nlat)
    tland    = keyword_set(land)

 if Zero_opt then BEGIN
    pacific_mask = intarr(360,180)  
    FOR I=0,359 DO BEGIN
      Pacific_mask(i,*) = mask((180+i) MOD 360,*)
    ENDFOR
    mask = pacific_mask 
  endif
if (Nlon Ne 360) OR (Nlat NE 180) then mask = rebin(mask,NLON,NLAT)
if tland then                          mask = abs(mask-1)

return,mask
END