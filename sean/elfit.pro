pro elfit

  file = 'snr_grid_npca=2,nfwhm=150,drsub=2,meansubmeanadd_collapsed.fitssnrmap.fits'

  img = readfits(file)
  xs = intarr((size(img))[1:2])
  ys = xs

  for i=0, (size(img))[1]-1 do begin ; assumes images are square
     xs[i,*] = i
     ys[*,i] = i
  endfor

  img[finite(img, /nan)] = 0

  p = mpfitellipse(xs, ys, /tilt, weights=img)

  tvscl, img, 0, 0, /nan

  phi = dindgen(101)*2D*!dpi/100
  ; New parameter P[4] gives tilt of ellipse w.r.t. coordinate axes
  ; We must rotate a standard ellipse to this new orientation
  xm = p[2] + p[0]*cos(phi)*cos(p[4]) + p[1]*sin(phi)*sin(p[4])
  ym = p[3] - p[0]*cos(phi)*sin(p[4]) + p[1]*sin(phi)*cos(p[4])
  oplot, xm, ym, color='ff'xl 

  stop

end

