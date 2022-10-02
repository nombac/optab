PRO n_species, fname
  
  COMPILE_OPT IDL2
  !EXCEPT=0
  !QUIET=1

  file = H5F_OPEN(fname)
  n_layer = H5D_READ(h5D_OPEN(file, 'n_layer'))
  n_species= H5D_READ(h5D_OPEN(file, 'n_species'))
  ndens = H5D_READ(h5D_OPEN(file, 'ndens'))
  temp = H5D_READ(h5D_OPEN(file, 'temp'))
  pres = H5D_READ(h5D_OPEN(file, 'pres'))
  rho = H5D_READ(h5D_OPEN(file, 'rho'))
  species = H5D_READ(h5D_OPEN(file, 'species'))
  H5F_CLOSE, file
  
  fname0 =fname.Replace('.h5','')
  SPAWN, 'mkdir '+fname0

  ;; pres = pres / 1d6
  ;; ytitle = 'P [bar]'         & yrange = [1d-13,1d4 ] & ylog = 1
  ;; temp = 5040d0 / temp 
  ;; ytitle = '\theta [K^{-1}]' & yrange = [0d0  ,55d0] & ylog = 0
  xtitle = 'T [K]'      & xrange = [1d2 , 1d8] & xlog = 1
  ytitle = 'P [Ba]'     & yrange = [1d-1,1d10] & ylog = 1
  

;  title = 'Equation of state'
;  vrange = [-12,0]
  v = ALOG10(rho)

;  title = 'H1+1'
  FOR i = 0, n_species[0]-1 DO BEGIN
     title = species[i]
     IF(STRLEN(title) EQ 0) THEN CONTINUE
     FOR j = 0, n_layer[0]-1 DO BEGIN
        IF(ndens[i,j] * 0d0 NE 0d0) THEN BEGIN
           PRINT, ndens[i,j], species[i], j
           STOP
        ENDIF
     ENDFOR
     IF(MAX(ndens[i,*]) EQ 0d0) THEN CONTINUE
     v[*] = ALOG10(ndens[i,*])
     vrange = [0,25]

     set_plot,'ps'
     ps_name = fname0+'/'+title+'.eps'
     !p.font=0
     DEVICE, FILE=ps_name, /color, /encapsulated, /SCHOOLBOOK
     
     PRINT, ps_name
     LOADCT, 33
     !p.charsize=1.4
  
     CGPLOT, [0], [0], /NODATA, XR=xrange, YR=yrange, XST=1, YST=1, XTIT=TEXTOIDL(xtitle), YTIT=TEXTOIDL(ytitle), YLOG=ylog, XLOG=xlog, TITLE=TEXTOIDL(title) ;, XTICKF='exponent', YTICKF='exponent'
     FOR j = 0, n_layer[0]-1 DO BEGIN
        CGOPLOT, [temp[j]], [pres[j]], PSYM=15, SYMSIZE=3, COL=BYTSCL(FLOAT(v[j]), MIN=vrange[0], MAX=vrange[1])
     ENDFOR
     CGCOLORBAR, RANGE=vrange, TITLE=TEXTOIDL('log \rho [g cm^{-3}]'), POS=[0.77,0.20,0.80,0.80], CHARSIZE=1.5, /VERTICAL, /RIGHT

     CGTEXT, 0.8, 0.95, fname0, /NORM, SIZE=0.75, COL='DARK GRAY', FONT=1, TT_FONT='COURIER'

     DEVICE,/CLOSE
     SPAWN, 'convert "'+ps_name+'" -background white -mosaic +matte "'+ps_name+'.png"'
  ENDFOR

  ;; epstopdf = 'n'
  ;; READ, PROMPT='Do you have "epstopdf" (y/N)? ', epstopdf
  ;; IF(epstopdf EQ 'y') THEN BEGIN
  ;;    SPAWN, 'epstopdf '+ps_name+'.eps '+ps_name+'.pdf'
  ;; ENDIF
  ;; convert = 'n'
  ;; READ, PROMPT='Do you have "convert" (y/N)? ', convert
  ;; IF(convert EQ 'y') THEN BEGIN
  ;;    SPAWN, 'convert '+ps_name+'.eps -background white -mosaic +matte '+ps_name+'.png'
  ;; ENDIF

END
