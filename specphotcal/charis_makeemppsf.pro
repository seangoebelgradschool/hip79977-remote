pro charis_makeemppsf,pfname,pick=pick,filt=filt,fsize=fsize,prad=prad,psfsize=psfsize,nonorm=nonorm,nozero=nozero,$
outname=outname

;make an empirical PSF from a set of registered images
;02-05-2018 - now puts empirical PSF in its own directory for later use.
;10-11-2017 - added 'nonorm' switch to not normalize the PSF.  Useful for a quick and dirty flux calibration estimate from an already-processed data set.
;added 'zero' switch to ensure that the background is zero'd

if ~keyword_set(outname) then outname='psfcube.fits'

datadir='./reduc/'
datadir1=datadir+'reg/'
;create directory for PSF model if it is not there ...
reducdir='./psfmodel/'
file_mkdir,reducdir

if ~keyword_set(fsize) then fsize=21
dist_circle,subarray_dist,fsize

if ~keyword_set(pick) then begin

prefname='n'
suffname='reg'
param,'fnum_sat',flist,/get,pfname=pfname
filenum=nbrlist(flist)
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)

test=readfits(datadir1+files[0],/exten,h1)
h0=headfits(datadir1+files[0],ext=0)

endif else begin

files=dialog_pickfile(Title="Pick Files",/multiple)
test=readfits(files[0],/exten,h1)
h0=headfits(files[0],ext=0)

endelse

;Now Get the Wavelength Vector
get_charis_wvlh,h0,wvlhs
lambda=wvlhs*1d-3

Dtel = 7.9d0 ;effective Subaru pupil diameter
fwhm=0.206*lambda/Dtel
fwhm/=0.0164  ;charis pixel scale

nwvlh=n_elements(lambda)

nfiles=n_elements(files)

dim=sxpar(h1,'naxis1')
xc=dim/2 & yc=dim/2

;imfourcube=fltarr(dim,dim,nwvlh,nfiles)

dpsf=21
if keyword_set(psfsize) then dpsf=psfsize

nspot=4

psf_spot=fltarr(dpsf,dpsf,nspot)
psf_imagelam=fltarr(dpsf,dpsf,nwvlh,nfiles)
psf_image=fltarr(dpsf,dpsf,nfiles)


for il=0L,nwvlh-1 do begin

print,'analyzing wavelength slice ',il+1

for nf=0L,nfiles-1 do begin
if ~keyword_set(pick) then begin
h0=headfits(datadir1+files[nf],/silent,ext=0)
im=(readfits(datadir1+files[nf],h1,/exten,/silent))[*,*,il]
endif else begin
h0=headfits(files[nf],/silent,ext=0)
im=(readfits(files[nf],h1,/exten,/silent))[*,*,il]
endelse

profrad_tc,im,p2d=pr
im-=pr

if keyword_set(filt) then begin
im-=filter_image(im,median=fsize)
endif

;sat spot
;psf_spot=fltarr(dim,dim,4)
for spot =0,3 do begin
 hdrval=double(strsplit(sxpar(h1,'SATS'+strtrim(il,2)+'_'+strtrim(spot,2)),' ',/extract))
 shifts=hdrval-floor(hdrval)

;print,'nfile and nlambda',nf,il
;print,'stuff1'
; help,im,dpsf
;print,'stuff2'
;print,floor([hdrval[0],hdrval[1]])

;break out of computation if NaN value for fits header
if (hdrval[0] gt 1000. or hdrval[1] gt 1000. $
or finite(hdrval[0]) eq 0 or finite(hdrval[1]) eq 0) then begin
psf_spot[*,*,spot]=!values.f_nan
goto,breakoutspot
endif

 s=subarr(im,dpsf,floor([hdrval[0],hdrval[1]]))

 psf_spot[*,*,spot]=shift_sub(s,-1*shifts[0],-1*shifts[1])
 psf_spot[*,*,spot]-=median(psf_spot[*,*,spot])

 if ~keyword_set(nozero) then begin
 dum=psf_spot[*,*,spot]
 bckgd=median(dum[where(subarray_dist gt 2*fwhm[il])],/even)
 dum-=bckgd
 psf_spot[*,*,spot]=dum
 endif

if ~keyword_set(nonorm) then begin
psf_spot[*,*,spot]/=total(psf_spot[*,*,spot],/nan)
endif
breakoutspot:
endfor

;stuff=median(psf_spot,dimension=3,/even)
psf_imagelam[*,*,il,nf]=median(psf_spot,dimension=3,/even)
;if keyword_set(zero) then begin
;psf_imagelam[*,*,il,nf]=stuff-median(stuff,/even)
;endif



endfor
endfor

if nfiles gt 1 then begin
psf_image=median(psf_imagelam,/even,dimension=4)
endif else begin
for il=0L,nwvlh-1 do begin
psf_imagelam[*,*,il]=fixpix_rs(psf_imagelam[*,*,il],iter=3)
endfor
psf_image=psf_imagelam
endelse

if ~keyword_set(nonorm) then begin
for il=0L,nwvlh-1 do psf_image[*,*,il]/=total(psf_image[*,*,il],/nan)
;for il=0L,nwvlh-1 do psf_image[*,*,il]-=median(psf_image[*,*,il],/even)
endif

if ~keyword_set(nonorm) then begin
writefits,reducdir+'psfcube_med'+'.fits',0
writefits,reducdir+'psfcube_med'+'.fits',psf_image,/append
endif else begin
writefits,reducdir+'psfcube_med_nonorm'+'.fits',0
writefits,reducdir+'psfcube_med_nonorm'+'.fits',psf_image,/append
endelse

;writefits,'psfimage_tot'+'.fits',psf_image
;writefits,'psfimage_'+strtrim(nf,2)+'.fits',psf_image[*,*,nf]


cube=median(psf_image,dimension=3,/even)
if ~keyword_set(nonorm) then begin
writefits,reducdir+'psf_medcol.fits',cube
endif else begin
writefits,reducdir+'psf_medcolnonorm.fits',cube
endelse


;cube=median(imfourcube,dimension=4,/even)
;writefits,'cube.fits',cube


end
