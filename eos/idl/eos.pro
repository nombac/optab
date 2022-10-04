PRO eos, fname=fname, arg=arg, syms=syms
  
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

  set_plot,'ps'
  !p.font=0
  DEVICE, FILE=fname0+'/eos.eps', /color, /encapsulated, /SCHOOLBOOK
  LOADCT, 33
  !p.charsize=1.4
  
  xtitle = 'T [K]'      & xrange = [1d2 ,1d8] & xlog = 1
  ytitle = 'P [Ba]'     & yrange = [1d-5,1d10] & ylog = 1

  IF(arg EQ 'rho') THEN BEGIN
     title = 'mass density'
     cb_title = 'log \rho [g cm^{-3}]'
     v = ALOG10(rho)
     vrange = [-16,-2]
  ENDIF ELSE IF(arg EQ 'mmw') THEN BEGIN
     title = 'mean molecular weight'
     cb_title = ''
     v = rho[*] / ndens[0,*] / !amu
     vrange = [0,3]
  ENDIF

  CGPLOT, [0], [0], /NODATA, XR=xrange, YR=yrange, XST=1, YST=1, XTIT=TEXTOIDL(xtitle), YTIT=TEXTOIDL(ytitle), YLOG=ylog, XLOG=xlog, TITLE=TEXTOIDL(title);, XTICKF='exponent', YTICKF='exponent'
  FOR j = 0, n_layer[0]-1 DO BEGIN
     IF(v[j]*0d0 NE 0d0) THEN BEGIN
        color = 'black'
     ENDIF ELSE BEGIN
        color = BYTSCL(FLOAT(v[j]), MIN=vrange[0], MAX=vrange[1])
     ENDELSE
     CGOPLOT, [temp[j]], [pres[j]], PSYM=15, SYMSIZE=syms, COL=color
  ENDFOR
  CGCOLORBAR, RANGE=vrange, TITLE=TEXTOIDL(cb_title), POS=[0.77,0.20,0.80,0.80], CHARSIZE=1.5, /VERTICAL, /RIGHT, COL='GRAY'

  CGTEXT, 0.8, 0.95, fname0, /NORM, SIZE=0.75, COL='DARK GRAY', FONT=1, TT_FONT='COURIER'

  DEVICE,/CLOSE

END
