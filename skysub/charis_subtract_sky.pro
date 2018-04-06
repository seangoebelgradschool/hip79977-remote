pro charis_subtract_sky,pfname,medbox=medbox,prefname=prefname,suffname=suffname

;****Performs sky subtraction on CHARIS data cubes.  Useful for removing detector noise (and sky noise at K band)
;08/16/17 - First version of sky subtraction routine.  Asks user to identify sky frames, subtracts them.  Puts in suffix of 'skysub'.

;Procedure:
;*1. Asks user to identify sky frames, creates a master sky cube from median at each pixel
;*2. Loops the data cubes, performs sky subtraction
;*3. (Optional) Subtracts off a moving-box median filter from residuals

; setupdir,reducdir=reducdir
reducdir='./reduc/'

;data directory
datadir=reducdir+'expand/'

;determine reduction subdirectory
subdir='expand/'
reducdir+=subdir

;create list of filenames
param,'fnum_sat',flist,/get,pfname=pfname
param,'fwhm',fwhm,/get,pfname=pfname

;**** Prefix names for your files (Make sure to change these with different data!!!)**
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='e'

filenum=nbrlist(flist)
filesout=filelist(filenum,prefix=prefname,suffix=suffname+'_skysub')
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)


;;get dim from first image header
htest=headfits(datadir+files[0])
fimagetest=readfits(datadir+files[0],htest1,ext=1)

;assume square array
dim=(size(fimagetest,/dim))[0]
sz=size(fimagetest,/dim)
d0=sz[0]

;****Select sky frames
skyframes=dialog_pickfile(Title="Select Sky Frames",/multiple_files)
nsky=n_elements(skyframes)
skycube=fltarr(sz[0],sz[1],sz[2],nsky)
skytime=fltarr(nsky)

;Loop to define the sky cube
for i=0L,nsky-1 do begin
a=readfits(skyframes[i],ext=1,h1sky)
h0sky=headfits(skyframes[i])
exp1time=sxpar(h0sky,'exp1time')
coadds=sxpar(h0sky,'coadds')
skytime[i]=exp1time*coadds
skycube[*,*,*,i]=a
endfor
;in case you have variable integration times for your sky frames
for i=0L,nsky-1 do skycube[*,*,*,i]*=(max(skytime)/skytime[i])

;now collapse to get a master sky
if n_elements(skycube gt 1) then begin
master_skycube=median(skycube,dimension=4,/even)
endif else begin
master_skycube=skycube
endelse

;****Loop to perform sky subtraction

for i=0L,nfiles-1 do begin

imcube=readfits(datadir+files[i],h1cube,ext=1)
h0cube=headfits(datadir+files[i])
exp1time=sxpar(h0cube,'exp1time')
coadds=sxpar(h0cube,'coadds')
exptime_sci=exp1time*coadds


;subtract sky.   For now do a simple subtraction since the PSF halo contaminates an estimate of the sky background
imcube-=master_skycube*exptime_sci/max(skytime)

writefits,reducdir+filesout[i],0,h0cube
writefits,reducdir+filesout[i],float(imcube),h1cube,/append
endfor

end

