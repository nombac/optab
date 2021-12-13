PRO read_mono_hdf5, dir, n, temp, temp2, grd, abs, sca, cnt, np, pla, pla2, ros, rho, plac, plal, plac2, plal2

  file = dir+'/output/mono_'+STRTRIM(STRING(n, FORMAT='(I05)'),2)+'.h5'

  CATCH, error
  stat = 0
  IF(error NE 0) THEN BEGIN
     stat = 1
     RETURN 
     CATCH, /cancel
  END
  
  f = H5F_OPEN(file)
  d = h5D_OPEN(f, 'temp')
  temp = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'temp2')
  temp2 = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'rho')
  rho = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'grd')
  grd = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'abs')
  abs = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'sca')
  sca = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'cnt')
  cnt = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'nden')
  np = H5D_read(d)
  H5D_CLOSE, d
  ngas = np[0]
  d = h5D_OPEN(f, 'pla')
  pla = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'pla2')
  pla2 = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'ros')
  ros = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'plac')
  plac = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'plal')
  plal = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'pla2')
  pla2 = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'plac2')
  plac2 = H5D_read(d)
  H5D_CLOSE, d
  d = h5D_OPEN(f, 'plal2')
  plal2 = H5D_read(d)
  H5D_CLOSE, d
  
  H5F_CLOSE, f

END
