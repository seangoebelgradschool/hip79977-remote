pro calc_disk_sb,image_in,theta,dy=dy,r_o=r_o,r_f=r_f,dr=dr,fwhm=fwhm_o,extcube=extcube,bothsides=bothsides,snrlim=snrlim,otherside=otherside,offset=offset,pick=pick,siglosscor=siglosscor,selsiglossfile=selsiglossfile,$
setspineheight=setspineheight,spineheightfile=spineheightfile,$
debug=debug

;***to do: force a crash if the user tries to calc sb on an uncalibrated cube (i.e. raw counts)

;Calculates the Disk Position Surface Brightness Profile Over a Range of Stellocentric Distances
; - For now, does this for three bins: J, H, and K
;- Method: rotates the disk to disk PA, finds the spine at each r separation (or user can override this and set a fixed value)
;          does a moving-box median to the disk in J, H, and K (box size = 1 FWHM in J, H, and K)
;          calculates the "signal to noise" of this median-filtered disk in J, H, and K
;          reports the disk intensity along the spine (intensity vs. separation on east and west sides both if /bothsides switch thrown)
;          divides by pixel scale to get flux density/square arc-sec --> Surface Brightness
;          reports uncertainty in SB using the "SNR" and the above SB map
;          transforms to units of magnitudes/square arc-sec
;          plots results, outputs results to .txt files
;
;****requires***
;- a data cube (rotated north-up)
;- the PA of the disk
;
;***can accept keywords to calculate ...
; - dy: the vertical extent of the disk over which to search for the spine (default = 61 pixels)
;- dr: the step size in r for spine calculation (default = 1)
;- fwhm: the FWHM of the image in pixels (default = 3)
;- snrlim: the limit in SNR to consider pixels for the spine fitting (default = 1)
;- r_o and r_f: the inner and out radii for spine fitting (default = 20 and 55)

;**other keywords
; - /extcube: switch to throw if the data cube has a primary and secondary extension, a la GPI and CHARIS
; - /bothsides: to calculate this for both sides of the image (e.g. no misalignment)

setcolors,/system_variables
if ~keyword_set(r_o) then begin
r_o=10
r_f=43
endif

if ~keyword_set(theta) then theta=112.4

if keyword_set(otherside) then theta+=180

; does simple analysis for edge-on disks
;...computes PA of disk, overall PA, FWHM profile and SB

if ~keyword_set(dr) then dr =3 
if ~keyword_set(fwhm) then fwhm_o=3.

;*******Throughput Correction*****
;If you have a signal loss data cube (ONLY WORKS FOR CHARIS RIGHT NOW!!!!!!!!!) ...
; then input it.   You can just name it or you can select it via gui


if keyword_set(siglosscor) then begin

if keyword_set(selsiglossfile) then begin
;siglosscube_in=dialog_pickfile(Title="Pick the Cube Mapping Disk Attenuation")
print, "sean overrode this pickfile command."
siglosscube_in = '/home/sgoebel/thayne/pipeline/biggerr/reduc/proc/model_atten/egrater_0.600000,1.00000,4,-3.50000,112.400,0,0,112.400,84.6000,70,2,0_attenfactcube.fits'
;***
if ~keyword_set(extcube) then begin
siglosscube=readfits(siglosscube_in,h1l)
endif else begin
siglosscube=readfits(siglosscube_in,ext=1,h1l)
endelse
;***

endif else begin

;***
if ~keyword_set(extcube) then begin
siglosscube=readfits(siglosscorr,h1l)
endif else begin
siglosscube=readfits(siglosscorr,h1l,ext=1)
endelse
;***

endelse

endif



;*******Input Data Cube*****
;Input
if ~keyword_set(pick) then begin

if ~keyword_set(extcube) then begin
image=readfits(image_in,h1)
endif else begin
image=readfits(image_in,ext=1,h1)
endelse

;image=filter_image(image,fwhm=fwhm_o)
;snrmap=readfits('snrmap.fits')
endif else begin
;image_in=dialog_pickfile(Title="Select Input Cube")
print, "sean overrode this pickfile command."
image_in = '/home/sgoebel/thayne/pipeline/biggerr/reduc/proc/biggerr,2,6,1,5,70_cal.fits'
if ~keyword_set(extcube) then begin
image=readfits(image_in,h1)
endif else begin
image=readfits(image_in,ext=1,h1)
endelse

endelse

s=size(image)

;if the program thinks you have a data cube, then assume CHARIS data.   If an image, assume HiCIAO data.
if s[0] eq 3 then begin
 dimx=s[1] 
 dimy=s[2] 
 dimcube=s[3]
get_charis_wvlh,h1,wvlh
lambda=1d-3*wvlh
Dtel=7.9d0
pixscale=0.0164
fwhm=1.*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale
fwhm_j=median(fwhm[0:4],/even)
fwhm_h=median(fwhm[7:13],/even)
fwhm_k=median(fwhm[15:20],/even)
endif else begin
lambda=1.65
dimx=s[1]
dimy=s[2]
pixscale=0.0083
fwhm=0.206*lambda/7.9/0.0083
endelse

xc=dimx/2 & yc=dimy/2

;vertical extent of the search strip
if ~keyword_set(dy) then dy=21.

if ~keyword_set(snrlim) then snrlim=1

;x and y indices

if ~keyword_set(bothsides) then begin
nr=(r_f-r_o)/dr+1
rrange=findgen(nr)*dr+r_o
endif else begin
nr=(r_f-r_o)/dr
rrange=findgen(nr)*dr+r_o
rrange2=-1*rrange
rrange=[rrange,rrange2]
nr=2*nr
endelse
print,rrange,nr
;stop

if ~keyword_set(offset) then offset=0
height=findgen(dy)-dy/2.+offset

if s[0] eq 3 then begin
image_rot=rotat_cube(image,-90+theta)
image_rot_col=median(image_rot,dimension=3,/even)
image_rot_j=median(image_rot[*,*,0:4],dimension=3,/even)
image_rot_h=median(image_rot[*,*,7:13],dimension=3,/even)
image_rot_k=median(image_rot[*,*,15:20],dimension=3,/even)

if keyword_set(siglosscor) then begin
siglosscube=rotat_cube(siglosscube,-90+theta)
sigloss_j=median(siglosscube[*,*,0:4],dimension=3,/even)
sigloss_h=median(siglosscube[*,*,7:13],dimension=3,/even)
sigloss_k=median(siglosscube[*,*,15:20],dimension=3,/even)
endif



;endif


writefits,'image_rot.fits',image_rot_col
writefits,'image_rot_j.fits',image_rot_j
writefits,'image_rot_h.fits',image_rot_h
;stop

sb_j=fltarr(nr)
sb_h=fltarr(nr)
sb_k=fltarr(nr)

esb_j=fltarr(nr)
esb_h=fltarr(nr)
esb_k=fltarr(nr)

;moving-box median filter
 filt_j=filter_image(image_rot_j,median=fwhm_j)
 filt_h=filter_image(image_rot_h,median=fwhm_h)
 filt_k=filter_image(image_rot_k,median=fwhm_k)
filtcol=filter_image(image_rot_col,median=median(lambda*pixscale))

if keyword_set(siglosscor) then begin
filt_siglossj=filter_image(sigloss_j,median=fwhm_j)
filt_siglossh=filter_image(sigloss_h,median=fwhm_h)
filt_siglossk=filter_image(sigloss_k,median=fwhm_k)
endif

;now do SB "snr map"
 snratio_exten,filtcol,fwhm=median(lambda*pixscale),/zero,/finite,snrmap=snrmap_col,noisemap=noisemap_col
 snratio_exten,filt_j,fwhm=fwhm_j,/zero,/finite,snrmap=snrmap_j,noisemap=noisemap_j
 snratio_exten,filt_h,fwhm=fwhm_h,/zero,/finite,snrmap=snrmap_h,noisemap=noisemap_h
 snratio_exten,filt_k,fwhm=fwhm_k,/zero,/finite,snrmap=snrmap_k,noisemap=noisemap_k
 writefits,'snrmap_j.fits',snrmap_j
 writefits,'snrmap_h.fits',snrmap_h
 writefits,'snrmap_k.fits',snrmap_k
 writefits,'snrmap_col.fits',snrmap_col

;if keyword_set(siglosscor) then begin
;writefits,'siglossj.fits',filt_siglossj
;writefits,'siglossj.fits',sigloss_j
;stop
;endif
 
 
;for i=0L,nr-1 do print, rrange[i],i
 for i=0L,nr-1 do begin
; image_range_col=image_rot_col[xc+rrange[i]-dr:xc+rrange[i]+dr,yc+height]
image_range_col=image_rot_col[xc+rrange[i],yc+height]
;image_snrmap_col=snrmap_col[xc+rrange[i]-dr:xc+rrange[i]+dr,yc+height]
image_snrmap_col=snrmap_col[xc+rrange[i],yc+height]
 image_range_j=image_rot_j[xc+rrange[i],yc+height]
 image_range_h=image_rot_h[xc+rrange[i],yc+height]
 image_range_k=image_rot_k[xc+rrange[i],yc+height]

;help,image_range_col
;good=where(finite(image_range_col) eq 1,ngood)
;if ngood eq 0 then continue
;image_range_col=image_range_col[good]
;height=height[good]
plot,height,image_range_col,psym=4
;stop
;g=mpfitpeak(height,image_range_col,apar,/lorentzian,perror=aparerror,weights=image_snrmap_col,nterms=4)
;g=gaussfit(height, image_range_col, apar, nterms=5)
;print,apar
;stop
;spine_height=apar[1]
if keyword_set(setspineheight) then begin
	spine_height=setspineheight
endif else if keyword_set(spineheightfile) then begin
	readcol, 'spinefit.txt', myx, myy
	result = poly_fit(myx, myy, 4)
	spine_height = result[4]*(xc+rrange[i])^4 + $
			result[3]*(xc+rrange[i])^3 + $ 
			result[2]*(xc+rrange[i])^2 + $ 
			result[1]*(xc+rrange[i]) + $ 
			result[0] - yc
	;spine_height = myy[where(myx eq xc+rrange[i])] - yc
endif
;spine_height=-3
;print,spine_height
;stop
;plot,height,image_range_col,psym=4
;oplot,height,g,linestyle=1
vline, spine_height
;****note: better method would be to robustly calculate this: OK for simple disk geometry, better method needed for more complex geometries.
;spine_height=-2 ;for HIP 79977
;spine_height=7 ;for HD 15115
 ;trace_pos=image_rot[xc+rrange[i],yc+round(spine_height)]

 print,'spine height is',spine_height,rrange[i],xc+rrange[i],yc+round(spine_height)
;stop
;help,filt_j
 sb_j[i]=filt_j[xc+rrange[i],yc+round(spine_height)]/pixscale^2.
 sb_h[i]=filt_h[xc+rrange[i],yc+round(spine_height)]/pixscale^2.
 sb_k[i]=filt_k[xc+rrange[i],yc+round(spine_height)]/pixscale^2.
 esb_j[i]=noisemap_j[xc+rrange[i],yc+round(spine_height)]/pixscale^2.
 esb_h[i]=noisemap_h[xc+rrange[i],yc+round(spine_height)]/pixscale^2.
 esb_k[i]=noisemap_k[xc+rrange[i],yc+round(spine_height)]/pixscale^2.


fract_errorj=esb_j[i]/sb_j[i]
fract_errorh=esb_h[i]/sb_h[i]
fract_errork=esb_k[i]/sb_k[i]

if keyword_set(siglosscor) then begin
attenj=filt_siglossj[xc+rrange[i],yc+round(spine_height)]
sb_j[i]/=attenj
esb_j[i]=sb_j[i]*fract_errorj
attenh=filt_siglossh[xc+rrange[i],yc+round(spine_height)]
sb_h[i]/=attenj
esb_h[i]=sb_h[i]*fract_errorh
attenk=filt_siglossk[xc+rrange[i],yc+round(spine_height)]
sb_k[i]/=attenk
esb_k[i]=sb_k[i]*fract_errork



endif



 writefits,'filt_j.fits',filt_j
 writefits,'filt_h.fits',filt_h
 if keyword_set(debug) then wait,1
 endfor

writecol,'sb_j.txt',rrange*pixscale,-2.5*alog10(sb_j*1d-3/1560.),1.0857*esb_j/sb_j
writecol,'sb_h.txt',rrange*pixscale,-2.5*alog10(sb_h*1d-3/1024.),1.0857*esb_h/sb_h
writecol,'sb_k.txt',rrange*pixscale,-2.5*alog10(sb_k*1d-3/666.7),1.0857*esb_k/sb_k

sbj=-2.5*alog10(sb_j*1d-3/1560.)
sbh=-2.5*alog10(sb_h*1d-3/1024.)
sbk=-2.5*alog10(sb_k*1d-3/666.7)
esbj=1.0857*esb_j/sb_j
esbh=1.0857*esb_h/sb_h
esbk=1.0857*esb_k/sb_k

rarc=rrange*pixscale
west=where(rarc gt 0,complement=east)
rarc[east]*=-1.

set_plot,'ps'
setcolors,/system_variables
device,filename='sb.eps',/encapsulated,/color
plot,rarc,sbj,/nodata,xrange=[0.2,1.0],yrange=[15,11],xthick=5,ythick=5,$
;plot,rarc,sbj,/nodata,xrange=[-0.8,0.8],yrange=[14,11],xthick=5,ythick=5,$
xtitle='R (arc-sec)',ytitle=textoidl('mag arcsec^{-2}'),charsize=1.25,charthick=3
;oploterror,rarc,sbj,esbj,color=!blue,psym=4
;oploterror,rarc,sbh,esbh,color=!green,psym=4
;oploterror,rarc,sbk,esbk,color=!red,psym=4
oploterror,rarc[west],sbj[west],esbj[west],color=!blue,psym=-4
oploterror,rarc[east],sbj[east],esbj[west],color=!blue,psym=-4,linestyle=1
oploterror,rarc[west],sbh[west],esbh[west],color=!green,psym=-4
oploterror,rarc[east],sbh[east],esbh[east],color=!green,psym=-4,linestyle=1
oploterror,rarc[west],sbk[west],esbk[west],color=!red,psym=-4
oploterror,rarc[east],sbk[east],esbk[east],color=!red,psym=-4,linestyle=1
device,/close

set_plot,'ps'
setcolors,/system_variables
device,filename='sbj.eps',/encapsulated,/color
plot,rarc,sbj,/nodata,xrange=[0.2,1.0],yrange=[16,11],xthick=5,ythick=5,$
;plot,rarc,sbj,/nodata,xrange=[-0.8,0.8],yrange=[14,11],xthick=5,ythick=5,$
xtitle='R (arc-sec)',ytitle=textoidl('mag/arcsec^{2}'),charsize=1.25,charthick=3
;oploterror,rarc,sbj,esbj,color=!blue,psym=4
;oploterror,rarc,sbh,esbh,color=!green,psym=4
;oploterror,rarc,sbk,esbk,color=!red,psym=4
oploterror,rarc[west],sbj[west],esbj[west],color=!blue,psym=-4
oploterror,rarc[east],sbj[east],esbj[west],color=!blue,psym=-4,linestyle=1
device,/close

set_plot,'ps'
setcolors,/system_variables
device,filename='sbh.eps',/encapsulated,/color
plot,rarc,sbh,/nodata,xrange=[0.2,1.0],yrange=[16,11],xthick=5,ythick=5,$
;plot,rarc,sbj,/nodata,xrange=[-0.8,0.8],yrange=[14,11],xthick=5,ythick=5,$
xtitle='R (arc-sec)',ytitle=textoidl('mag/arcsec^{2}'),charsize=1.25,charthick=3
;oploterror,rarc,sbj,esbj,color=!blue,psym=4
;oploterror,rarc,sbh,esbh,color=!green,psym=4
;oploterror,rarc,sbk,esbk,color=!red,psym=4
oploterror,rarc[west],sbh[west],esbj[west],color=!blue,psym=-4
oploterror,rarc[east],sbh[east],esbj[west],color=!blue,psym=-4,linestyle=1
device,/close

set_plot,'ps'
setcolors,/system_variables
device,filename='sbk.eps',/encapsulated,/color
plot,rarc,sbh,/nodata,xrange=[0.2,1.0],yrange=[16,11],xthick=5,ythick=5,$
;plot,rarc,sbj,/nodata,xrange=[-0.8,0.8],yrange=[14,11],xthick=5,ythick=5,$
xtitle='R (arc-sec)',ytitle=textoidl('mag/arcsec^{2}'),charsize=1.25,charthick=3
;oploterror,rarc,sbj,esbj,color=!blue,psym=4
;oploterror,rarc,sbh,esbh,color=!green,psym=4
;oploterror,rarc,sbk,esbk,color=!red,psym=4
oploterror,rarc[west],sbk[west],esbj[west],color=!blue,psym=-4
oploterror,rarc[east],sbk[east],esbj[west],color=!blue,psym=-4,linestyle=1
device,/close

device,filename='colorsjh.eps',/encapsulated,/color
plot,rarc,sbj-sbh,/nodata,xrange=[0.2,1.0],yrange=[-2,2],xthick=5,ythick=5,$
xtitle='R (arc-sec)',ytitle=textoidl('J-H'),charsize=1.25,charthick=3
oploterror,rarc[west],sbj[west]-sbh[west],sqrt(esbj[west]^2.+esbh[west]^2.),color=!blue,psym=4
oploterror,rarc[east],sbj[east]-sbh[east],sqrt(esbj[east]^2.+esbh[east]^2.),color=!green,psym=2
;oploterror,rarc,sbh-sbk,sqrt(esbh^2.+esbk^2.),color=!green,psym=4
device,/close

device,filename='colorshk.eps',/encapsulated,/color
plot,rarc,sbh-sbk,/nodata,xrange=[0.2,1.0],yrange=[-2,2],xthick=5,ythick=5,$
xtitle='R (arc-sec)',ytitle=textoidl('H-K'),charsize=1.25,charthick=3
oploterror,rarc[west],sbh[west]-sbk[west],sqrt(esbk[west]^2.+esbh[west]^2.),color=!blue,psym=4
oploterror,rarc[east],sbh[east]-sbk[east],sqrt(esbk[east]^2.+esbh[east]^2.),color=!green,psym=2
;oploterror,rarc,sbh-sbk,sqrt(esbh^2.+esbk^2.),color=!green,psym=4
device,/close

device,filename='colorsjk.eps',/encapsulated,/color
plot,rarc,sbj-sbk,/nodata,xrange=[0.2,1.0],yrange=[-2,2],xthick=5,ythick=5,$
xtitle='R (arc-sec)',ytitle=textoidl('J-K'),charsize=1.25,charthick=3
oploterror,rarc[west],sbj[west]-sbk[west],sqrt(esbk[west]^2.+esbj[west]^2.),color=!blue,psym=4
oploterror,rarc[east],sbj[east]-sbk[east],sqrt(esbk[east]^2.+esbj[east]^2.),color=!green,psym=2
;oploterror,rarc,sbh-sbk,sqrt(esbh^2.+esbk^2.),color=!green,psym=4
device,/close

endif else begin

sb_h=fltarr(nr)
esb_h=fltarr(nr)
image_rot=rotat(image,-90+theta)
filt_h=filter_image(image_rot,median=fwhm)
snratio_exten,filt_h,fwhm=fwhm_h,/zero,/finite,snrmap=snrmap_h,noisemap=noisemap_h,rmax=300
writefits,'image_rot.fits',image_rot
writefits,'filt_h.fits',filt_h
 writefits,'snrmap_h.fits',snrmap_h
;stop
for i=0L,nr-1 do begin

image_range_h=image_rot[xc+rrange[i],yc+height]
g=mpfitpeak(height,image_range_h,apar,/gaussian,perror=aparerror,weights=image_range_h^2.,nterms=5)
spine_height=apar[1]
spine_height=-3
plot,height,image_range_h,psym=4
oplot,height,g,linestyle=0
;wait,1
sb_h[i]=filt_h[xc+rrange[i],yc+round(spine_height)]/pixscale^2.
esb_h[i]=noisemap_h[xc+rrange[i],yc+round(spine_height)]/pixscale^2.

endfor

writecol,'sb_h.txt',rrange*pixscale,-2.5*alog10(sb_h*1d-3/1024.),1.0857*esb_h/sb_h
sbh=-2.5*alog10(sb_h*1d-3/1024.)+13
esbh=1.0857*esb_h/sb_h
rarc=rrange*pixscale
west=where(rarc gt 0,complement=east)
rarc[east]*=-1.
set_plot,'ps'
setcolors,/system_variables
device,filename='sbh.eps',/encapsulated,/color
;plot,rarc,sbh,/nodata,xrange=[0.3,1.8],yrange=[18,11],xthick=5,ythick=5,$
plot,rarc,sbh,/nodata,xrange=[0.3,1.],yrange=[16,11],xthick=5,ythick=5,$
;plot,rarc,sbj,/nodata,xrange=[-0.8,0.8],yrange=[14,11],xthick=5,ythick=5,$
xtitle='R (arc-sec)',ytitle=textoidl('mag/arcsec^{2}'),charsize=1.25,charthick=3
oploterror,rarc[west],sbh[west],esbh[west],color=!green,psym=-4
oploterror,rarc[east],sbh[east],esbh[east],color=!green,psym=-4,linestyle=1
device,/close

endelse



end
