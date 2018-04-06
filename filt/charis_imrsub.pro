pro charis_imrsub,pfname,$
mask=mask,rmax=rmax,psfsub=psfsub,prad=prad,pang=pang,sdi=sdi,$
medbox=medbox,prefname=prefname,suffname=suffname

;Spatially filters a CHARIS data cube.   Limited functionality right now
;*** 5/25/2017 - quick and dirty fix to adapt to CHARIS data


;A. Spatial Filtering
;1. radial profile subtraction
;2. unsharp masking

;B. PSF subtraction
;1. If "psfsub" switch is thrown, does median combination of images --> reference PSF and subtracts



;box size for median filtering
;if ~keyword_set(medbox) then medbox=11

; setupdir,reducdir=reducdir
reducdir='./reduc/'

;data directory
datadir=reducdir+'reg/'

;determine reduction subdirectory
subdir='rsub/'
reducdir+=subdir

;create list of filenames
;param,'obsdate',date,/get,pfname=pfname & date=strtrim(date,2)
param,'fnum_sat',flist,/get,pfname=pfname
param,'fwhm',fwhm,/get,pfname=pfname


;**** Prefix names for your files (Make sure to change these with different data!!!)**
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'

filenum=nbrlist(flist)

filesout=filelist(filenum,prefix=prefname,suffix='rsub')
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)

;Reads in parallactic angle
readcol,'reduc.log',filenum,allxc,allyc,allrsat,allha,allpa


paranglef=fltarr(nfiles)

param,'spang*',spang,/get,pfname=pfname
param,'cenmask',cenmask,/get,pfname=pfname
cenmasks=strsplit(cenmask,',',/extract)

;defines width of spider mask ... in case you want to use it
param,'spmask',spmask,/get,pfname=pfname

;get dim from first image header
htest=headfits(datadir+files[0])
fimagetest=readfits(datadir+files[0],htest1,ext=1)

;assume square array
dim=(size(fimagetest,/dim))[0]
sz=size(fimagetest,/dim)
d0=sz[0]
if ~keyword_set(xc0) then xc0=dim/2
if ~keyword_set(yc0) then yc0=dim/2

;Maximum radial distance for radial profile subtraction
if ~keyword_set(rmax) then rmax=dim/2


;mask the sat spots?
if keyword_set(mask) then begin
spi=mkspider(d0,spmask,spang,/justind)
endif

;Loop to subtract

imcubeout=fltarr(sz[0],sz[1],sz[2])
imgiantcube=fltarr(sz[0],sz[1],sz[2],nfiles)

for i=0L,nfiles - 1 do begin
print,'reading in file number ',i+1
imcube=readfits(datadir+files[i],h1cube,ext=1)
h0cube=headfits(datadir+files[i])

paranglef[i]=allpa[i]
;loop on cube element

;*************FILTERING*****
;You have the choice of doing one or more of the following
;1. subtract reference PSF (not implemented)
;2. subtract 2D radial profile
;3. Fourier filter
;4. unsharp masking


for ii=0L,sz[2]-1 do begin

 imslice=imcube[*,*,ii]
 ;***2D Radial profile subtraction

if keyword_set(prad) then begin
profrad_tc,imslice,1,allrsat[i],rmax,p2d=pr,p1d=p1d,rayon=r
;profrad_tc,imslice,1,allrsat[i],rmax,p2d=pr,p1d=p1d,rayon=r,/med4q
imslice-=pr
endif


if keyword_set(pang) then begin
profang_tc,imslice,1,allrsat[i],rmax,p2d=pr,rayon=r
;writefits,'pr.fits',pr
;writefits,'imslice.fits',imslice
imslice-=pr
;writefits,'imslice2.fits',imslice
;stop
endif

;***filtering
if keyword_set(lowpass) then begin
;im=filter_image(im,median=0.5*fwhm)
;fwhm=2*fwhm
imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax]=filter_image(imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax],median=0.5*fwhm)
endif


;****''unsharp masking''
if keyword_set(medbox) then begin

;im=im-filter_image(im,median=medbox)
imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax]=imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax] $
- filter_image(imslice[d0/2-rmax:d0/2+rmax,d0/2-rmax:d0/2+rmax],median=medbox)

endif

;**** Fourier Filtering
;not well understood right now - 1/18/2016
;if keyword_set(fft) then begin
;imslice=fft_filt(imslice,/highpass,boxsize=fft)
;endif

if keyword_set(mask) then imslice[spi]=!values.f_nan
imcubeout[*,*,ii]=imslice
endfor

writefits,reducdir+filesout[i],0,h0cube
writefits,reducdir+filesout[i],imcubeout,h1cube,/append
imgiantcube[*,*,*,i]=imcubeout
endfor

;Now, for combining and derotating
imgiantcuberot=imgiantcube
imgiantavgcube=imgiantcube

if keyword_set(psfsub) then begin
imcuberef=median(imgiantavgcube,dimension=4,/even)
endif
;writefits,'imcuberef.fits',imcuberef
;stop

;SDI
if keyword_set(sdi) then begin 
imgiantsdicubeavg=imgiantcube
imgiantsdicuberot=imgiantcube
endif

for i=0L,nfiles-1 do begin
imt=readfits(datadir+files[i],ext=1,h1cubeold)
h0cubeold=headfits(datadir+files[i])

print,'Subtracting File Number ',long(i)
;derotating

pardiff=paranglef[0]-paranglef[i]
print,pardiff
imuse =imgiantcuberot[*,*,*,i]
if keyword_set(psfsub) then imuse-=imcuberef
;writefits,'imuse.fits',imuse
;stop

;SDI
if keyword_set(sdi) then begin
imsdi=imuse
nslices=(size(imsdi,/dim))[2]
charis_alignspeckle,imsdi,h0cubeold,h1cubeold,refslice=nslices-1,/nolocs
imsdiavg=median(imsdi,dimension=3,/even)
;writefits,'alignspeckle.fits',imsdi
;writefits,'sdiavg.fits',imsdiavg

dist_circle,roi,sz[1]
good=where(roi le 70)
get_charis_wvlh,h0cubeold,wvlhinput
wvlhratio=(wvlhinput/wvlhinput[11])^(-2.)
for ij=0L,nslices-1 do begin 
;medslice=median((imsdi[*,*,ij])[good],/even)
;medavg=median(imsdiavg[good],/even)
weightpsf=total(imsdi[*,*,ij]*imsdiavg,/nan)/total(imsdiavg*imsdiavg,/nan)
imsdi[*,*,ij]-=weightpsf*imsdiavg
;imsdi[*,*,ij]-=imsdiavg*wvlhratio[ij]
;imsdi[*,*,ij]-=imsdiavg*medslice/medavg
;imsdi[*,*,ij]-=imsdiavg*(median((imsdi[*,*,ij])[good],/even))/(median(imsdiavg[good],/even))
;imsdi[*,*,ij]-=imsdiavg*((median((imsdi[*,*,ij])[good]),/even))/(median(imsdiavg[good],/even))
endfor
;writefits,'sdisub.fits',imsdi
charis_alignspeckle,imsdi,h0cubeold,h1cubeold,refslice=nslices-1,/reverse,/nolocs
imuse=imsdi
;writefits,'ralignsdi.fits',imuse
;stop
endif


;imgiantavgcube[*,*,*,i]=imuse
imusef=rotat_cube(imuse,1*pardiff,missing=!values.f_nan)
;imusef=imuse
imgiantcuberot[*,*,*,i]=imusef
;writefits,'gah.fits',imusef
;stop
endfor

imcubefinal=median(imgiantcuberot,dimension=4,/even)
writefits,'imcubefinal_derot.fits',0,htest
writefits,'imcubefinal_derot.fits',imcubefinal,htest1,/append
imfinal=median(imcubefinal,dimension=3,/even)
writefits,'imcubefinal_derot_collapsed.fits',0,htest
writefits,'imcubefinal_derot_collapsed.fits',imfinal,htest1,/append

imcubeavg=median(imgiantavgcube,dimension=4,/even)
writefits,'imcubefinal_avg.fits',0,htest
writefits,'imcubefinal_avg.fits',imcubeavg,htest1,/append
imfinalavg=median(imcubeavg,dimension=3,/even)
writefits,'imcubefinal_avg_collapsed.fits',0,htest
writefits,'imcubefinal_avg_collapsed.fits',imfinalavg,htest1,/append


end



endfor
