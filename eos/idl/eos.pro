PRO eos, fname
  
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
  
  layer = LINDGEN(n_layer)

  ps_name =fname.Replace('.h5','')

  set_plot,'ps'
  !p.font=0
  DEVICE, FILE=ps_name+'.eps', /color, /encapsulated, /SCHOOLBOOK
  LOADCT, 10
  !p.charsize=1.4
  
  pres = pres / 1d6
  xtitle = 'P [bar]'         & xrange = [1d-13,1d10] & xlog = 1
  temp = 5040d0 / temp 
  ytitle = '\theta [K^{-1}]' & yrange = [0d0  ,55d0] & ylog = 0

  title = 'Equation of state'
  vrange = [-20,0]
  v = ALOG10(rho)
  
  CGPLOT, [0], [0], /NODATA, XR=xrange, YR=yrange, XST=1, YST=1, XTIT=TEXTOIDL(xtitle), YTIT=TEXTOIDL(ytitle), YLOG=ylog, XLOG=xlog, TITLE=TEXTOIDL(title)
  FOR j = 0, n_layer[0]-1 DO BEGIN
     CGOPLOT, [pres[j]], [temp[j]], PSYM=15, SYMSIZE=3.5, COL=BYTSCL(FLOAT(v[j]), MIN=vrange[0], MAX=vrange[1])
  ENDFOR
  CGCOLORBAR, RANGE=vrange, TITLE=TEXTOIDL('log \rho [g cm^{-3}]'), POS=[0.77,0.20,0.80,0.80], CHARSIZE=1.5, /VERTICAL, /RIGHT

  DEVICE,/CLOSE

  epstopdf = 'n'
  READ, PROMPT='Do you have "epstopdf" (y/N)? ', epstopdf
  IF(epstopdf EQ 'y') THEN BEGIN
     SPAWN, 'epstopdf '+ps_name+'.eps '+ps_name+'.pdf'
  ENDIF
  convert = 'n'
  READ, PROMPT='Do you have "convert" (y/N)? ', convert
  IF(convert EQ 'y') THEN BEGIN
     SPAWN, 'convert '+ps_name+'.eps -background white -mosaic +matte '+ps_name+'.png'
  ENDIF

END
