PRO opac, dir=dir, mean=mean, syms=syms
  
  COMPILE_OPT IDL2
  !EXCEPT=0
  !QUIET=1
  DEFSYSV, '!K_BOL',    1.3806503d-16

  SPAWN, 'ls -l ' + dir + '/output/mono_?????.h5 | wc -l', nmax
  nmax = LONG(nmax[0])

  PRINT, FORMAT='(a,I5,a,$)', 'reading ', nmax, ' layers ... '
  n = 0
  READ_MONO_HDF5, dir, n+1, tmp, tmp2, grd, abs, sca, cnt, np, pla, pla2, ros, rho, plac, plal, plac2, plal2
  tmp_tot = DBLARR(nmax)
  rho_tot = DBLARR(nmax)
  pre_tot = DBLARR(nmax)
  ros_tot = DBLARR(nmax)
  pla_tot = DBLARR(nmax)
  pla2_tot = DBLARR(nmax)
  tmp_tot[n] = tmp
  rho_tot[n] = rho
  pre_tot[n] = np[0] * !k_bol * tmp
  ros_tot[n] = ros
  ;; pla_tot[n] = pla
  ;; pla2_tot[n] = pla2
  pla_tot[n] = plac + plal
  pla2_tot[n] = plac2 + plal2
  prog = 0
  FOR n = 1, nmax-1 DO BEGIN
     READ_MONO_HDF5, dir, n+1, tmp, tmp2, grd, abs, sca, cnt, np, pla, pla2, ros, rho, plac, plal, plac2, plal2
     IF((n - 1) / (nmax/10) EQ prog) THEN BEGIN
        PRINT, prog, FORMAT='(I01,$)'
        prog += 1
     ENDIF
     tmp_tot[n] = tmp
     rho_tot[n] = rho
     pre_tot[n] = np[0] * !k_bol * tmp
     ros_tot[n] = ros
     ;; pla_tot[n] = pla
     ;; pla2_tot[n] = pla2
     pla_tot[n] = plac + plal
     pla2_tot[n] = plac2 + plal2
  ENDFOR
  
  ros_tot[*] /= rho_tot[*]
  pla_tot[*] /= rho_tot[*]
  pla2_tot[*] /= rho_tot[*]
  
  ros_tot = ALOG10(ros_tot)
  pla_tot = ALOG10(pla_tot)
  pla2_tot = ALOG10(pla2_tot)

  FOR n = 0, nmax-1 DO BEGIN
     print, n, pla_tot[n], tmp_tot[n], rho_tot[n]
  ENDFOR
  vmax =  7d0
  vmin = -6d0
  ;; vmax =  4d0
  ;; vmin = -2d0
  
  color_ros = BYTE(((((ros_tot) > vmin) < vmax) - vmin) / (vmax - vmin) * 255)
  color_pla = BYTE(((((pla_tot) > vmin) < vmax) - vmin) / (vmax - vmin) * 255)
  color_pla2= BYTE(((((pla2_tot)> vmin) < vmax) - vmin) / (vmax - vmin) * 255)

  ps_name = dir+'/output/'+mean
  file = ps_name+'.eps'
  set_plot,'ps'
  !p.font=0
  DEVICE, FILE=file, /color, /encapsulated, /SCHOOLBOOK
;  LOADCT, 10
  JHCOLORS
  !P.CHARSIZE = 1.4

  IF(0) THEN BEGIN
     xtitle = 'P [bar]' & xrange = [1d-13,1d10] & xlog = 1
     ytitle = '\theta [K^{-1}]'           & yrange = [0d0  ,55d0] & ylog = 0
     xarr = pre_tot / 1d6
     yarr = 5040d0 / tmp_tot
  ENDIF

  ;; xtitle = 'T [K]'      & xrange = [1d2 ,1d8] & xlog = 1
  ;; ytitle = 'P [Ba]'     & yrange = [1d-1,1d10] & ylog = 1
  ;; yarr = pre_tot
  ytitle = '\rho [gcm^{-3}]' & yrange = [1d-21,1d0] & ylog = 1
  yarr = rho_tot
  xtitle = 'T [K]'  & xrange = [1d1,1d10] & xlog = 1
  xarr = tmp_tot




  CASE mean OF
     'ross': BEGIN
        title = 'Rosseland-mean opacity'
        btitle = 'log \kappa_{Ross} [cm^2/g]'
        col = color_ros
     END
     'pla': BEGIN
        title = 'Planck-mean opacity'
        btitle = 'log \kappa_{Pla} [cm^2/g]'
        col = color_pla
     END
     'pla2': BEGIN
        title = 'Planck-mean opacity at T_{rad}='+STRTRIM(STRING(FIX(tmp2),FORMAT='(I5)'),2)+'K'
        btitle = 'log \kappa_{Pla} [cm^2/g]'
        col = color_pla2
     END
     ELSE: BEGIN
        PRINT, 'ERROR'
        STOP
     END
  ENDCASE

  CGPLOT, [0], [0], XR=xrange, YR=yrange, XS=1, YS=1, XLOG=xlog, YLOG=ylog, XTIT=TEXTOIDL(xtitle), YTIT=TEXTOIDL(ytitle), TIT=TEXTOIDL(title)
  FOR n = 0, nmax-1 DO BEGIN
     CGOPLOT, [xarr[n]], [yarr[n]], PSYM=15, SYMSIZE=syms, COL=col[n]
  ENDFOR
  CGCOLORBAR, RANGE=[vmin, vmax], TITLE=TEXTOIDL(btitle), POS=[0.77,0.20,0.80,0.80], CHARSIZE=1.5, /VERTICAL, /RIGHT

  DEVICE,/CLOSE
  PRINT, ''
  PRINT, 'output eps = ', file
  
  SPAWN, 'epstopdf '+ps_name+'.eps '+ps_name+'.pdf'
  SPAWN, 'convert '+ps_name+'.eps -background white -mosaic +matte '+ps_name+'.png'

END
