PRO mono, dir=dir, layer=layer
  
  COMPILE_OPT IDL2
  !EXCEPT=0
  !QUIET=1

  SPAWN, 'ls -l ' + dir + '/output/mono_?????.h5 | wc -l', nmax
  nmax = LONG(nmax[0])
  IF(layer GT nmax-1) THEN BEGIN
     PRINT, 'ERROR: layer ',layer,' > layer max',nmax-1
     STOP
  ENDIF
  
  PRINT, FORMAT='(a,I5,a,$)', 'reading ', nmax, ' layers ... '
  n = 0
  READ_MONO_HDF5, dir, n+1, tmp, tmp2, grd, abs, sca, cnt, np
  IF(layer NE 0) THEN BEGIN
     prog = 0
     FOR n = 1, nmax-1 DO BEGIN
        READ_MONO_HDF5, dir, n+1, tmp, tmp2, grd, abs, sca, cnt, np
        IF((n - 1) / (nmax/10) EQ prog) THEN BEGIN
           PRINT, prog, FORMAT='(I01,$)'
           prog += 1
        ENDIF
        IF(n EQ layer) THEN BREAK
     ENDFOR
  ENDIF
  pre = np[0] * !k_bol * tmp

  grd0 = DINDGEN(101, START=1d0, INCREMENT=0.1d0)
  grd0 = 10d0^grd0
  planck = grd0^3 / (EXP(!c2 * grd0 / tmp[0]) - 1d0)
  
  ps_name = dir+'/output/mono_'+STRTRIM(STRING(layer,FORMAT='(I05)'),2)

  set_plot,'ps'
  !p.font=0
  DEVICE, FILE=ps_name+'.eps', /color, /encapsulated, /SCHOOLBOOK
  LOADCT, 10
  !P.charsize = 1.4
  
  xtitle = '\nu [cm^{-1}]' & xrange = [1d1,1d9] & xlog = 1
  ytitle = '\alpha [cm^{-1}]'  & yrange = [1d-25,1d10] & ylog = 1
  title = 'monochromatic opacity'

  planck /= MAX(planck)
  
  CGPLOT, [0], [0], XR=xrange, YR=yrange, XS=1, YS=1, XLOG=xlog, YLOG=ylog, XTIT=TEXTOIDL(xtitle), YTIT=TEXTOIDL(ytitle), TIT=TEXTOIDL(title)
  CGOPLOT, grd, abs, COL='RED'
  CGOPLOT, grd, cnt-sca, COL='GRAY'
  CGOPLOT, grd, sca, COL='BLUE'
  CGOPLOT, grd0, planck, LI=1
  CGTEXT, 0.25, 0.85, 'abs (cnt. + line)', COL='RED', /NORM, CHARSIZE=1
  CGTEXT, 0.25, 0.80, 'abs (cnt.)       ', COL='GRAY', /NORM, CHARSIZE=1
  CGTEXT, 0.25, 0.75, 'sca              ', COL='BLUE', /NORM, CHARSIZE=1
  ;; CGTEXT, 0.70, 0.80, TEXTOIDL('\theta = '+STRTRIM(STRING(5040d0/tmp,FORMAT='(F6.1)'),2)+' [K^{-1}]'), /NORM, CHARSIZE=1
  CGTEXT, 0.70, 0.80, TEXTOIDL('logT = '+STRTRIM(STRING(ALOG10(tmp),FORMAT='(F6.1)'),2)+' [K]'), /NORM, CHARSIZE=1
  ;; CGTEXT, 0.70, 0.85, TEXTOIDL('P = 10^{'+STRTRIM(STRING(ALOG10(pre/1d6), FORMAT='(F6.1)'),2)+'} [bar]'), /NORM, CHARSIZE=1
  CGTEXT, 0.70, 0.85, TEXTOIDL('logP = '+STRTRIM(STRING(ALOG10(pre), FORMAT='(F6.1)'),2)+' [Ba]'), /NORM, CHARSIZE=1

  DEVICE,/CLOSE
  PRINT, ''
  PRINT, 'output eps = ', ps_name+'.eps'
  
  SPAWN, 'epstopdf '+ps_name+'.eps '+ps_name+'.pdf'
  SPAWN, 'convert '+ps_name+'.eps -background white -mosaic +matte '+ps_name+'.png'

END
