pro charis_imregister2,pfname,prefname=prefname,suffname=suffname,$
nosats=nosats,simple=simple,smallsteps=smallsteps,bandaid=bandaid,$
refslice=refslice,$
check=check,$
;satfit=satfit,$
;fluxnorm=fluxnorm,$
;fratio=fratio,$
;rcore=rcore,$
;rhalo=rhalo,$
;fluxout=fluxout,$
rewrite=rewrite,verbose=verbose

;08/14/2017 - Simpler program, renamed to charis_imregister2.pro
;07/10/2017 - Renamed to charis_imregister.pro
;07/10/2017 - Now at least for the sat-spot = true case, saves the newsatspot positions in fits header.
;04/13/2017 - Added switch to do x-corr fitting of the halo instead of satellite spots
;04/08/2017 - Redone for CHARIS to find the satellite positions
;03/09/2014 - Quick hack-version for GPI image registration across cube
; components.  

;1. Uses an initial guess for the header positions from the reference slice that is hardwired for now.
;2. Computes centroid position for this reference slice.  
;3. Goes on to compute satellite spot positions at other wavelengths.
;4. shift channels based on median offset.

;puts the star at the dimx/2 dimy/2 position. 

; setupdir,reducdir=reducdir
reducdir='./reduc/'

;data directory
datadir=reducdir+'expand/'
;datadir+=subdir

;determine reduction subdirectory
subdir='reg/'
reducdir+=subdir

;define a temporary directory
tmpdir=reducdir+'tmp/'
file_mkdir,tmpdir


;create list of filenames
param,'obsdate',date,/get,pfname=pfname & date=strtrim(date,2)
param,'fnum_sat',flist,/get,pfname=pfname
param,'rsat',rsat,/get,pfname=pfname
param,'raper',r_aper,/get,pfname=pfname


;**** Prefix names for your files (Make sure to change these with different data!!!)**
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='e'
if ~keyword_set(rhalo) then rhalo=100
if ~keyword_set(rcore) then rcore=r_aper

filenum=nbrlist(flist)
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
filesout=filelist(filenum,prefix=prefname,suffix='reg')
;*************

paranglef=fltarr(nfiles)

if keyword_set(fluxnorm) then begin
fluxval=fltarr(nfiles)
endif

param,'spang*',spang,/get,pfname=pfname
param,'cenmask',cenmask,/get,pfname=pfname
param,'spmask',spmask,/get,pfname=pfname
param,'RA',ra2000,/get,pfname=pfname
param,'fwhm',fwhmi,/get,pfname=pfname
ra2000=(24/360.)*ra2000
cenmask=strsplit(cenmask,',',/extract)
cenmasks=cenmask

;get dim from first image header
h=headfits(datadir+files[0],/exten)
h0=headfits(datadir+files[0])
print,files[0]
dim=sxpar(h,'naxis1')
if ~keyword_set(xc0) then xc0=dim/2
if ~keyword_set(yc0) then yc0=dim/2

;**********Auto-Corr Registration of first image if no sat spots********
;**if no sat spots, do x-corr registration of first image
if keyword_set(nosats) then begin 

fimage=median(readfits(datadir+files[0],/exten),dimension=3,/even)
cenmasks=cenmask
cenmasks[2]=cenmask[2]+0

; if you know you are already close to the centroid position, this switch prevents 
; - the x-corr search from moving too far away from initial guess.


if keyword_set(smallsteps) then begin

if ~keyword_set(subreg) then begin
;--x-corr of image with 180 deg rotation of itself.  assumes speckle symmetry to zeroth order
fidcenter=getcenter(fimage,xc0,yc0,spang,cenmasks,/smallsteps)
endif else begin
;-- self-subtracts halo from rotated halo.  Works well with ACS, not so much with others.
fidcenter=getcentersub(fimage,xc0,yc0,spang,cenmasks,/smallsteps)
endelse

endif else begin

if ~keyword_set(subreg) then begin
fidcenter=getcenter(fimage,xc0,yc0,spang,cenmasks)
endif else begin
fidcenter=getcentersub(fimage,xc0,yc0,spang,cenmasks)
endelse

endelse

xc1=fidcenter[0] & yc1=fidcenter[1]

;In case you need to put a bandaid on the registration, use with caution!!!
if keyword_set(bandaid) then begin
xc1=xc1+bandaid[0]
yc1=yc1+bandaid[1]
endif else begin
if keyword_set(bandaid) then bandaid[*,*]=0
endelse

print,'centroid',xc1,yc1
s=size(fimage) & d0=s[1]

;define initial delta(x,y) offsets from image slice dimension center
deltax=d0/2-xc1
deltay=d0/2-yc1
deltaxc=deltax
deltayc=deltay

print,fimage[round(xc1),round(yc1)]

if ~keyword_set(noshift) then begin
fimage=shift_sub(fimage,1.*deltax,1.*deltay)
endif else begin

if keyword_set(bandaid) then begin
;if you just want to manually set the centroid position from the starting position.
print,'shifting by',bandaid[*]
fimage=shift_sub(fimage,bandaid[0],bandaid[1])
deltaxc=bandaid[0]
deltayc=bandaid[1]
endif

endelse

fimage_ref=fimage
print,'delta x y',deltax, deltay,d0/2.

;***mask the spider***
spi=mkspider(d0,10,spang,/justind)
;spiref=spi

;*** Define region of interest to determine image registration***
;This is flexible to include all of the image, a circular region of interest,
;spidermask
spi=mkspider(d0,cenmasks[2],spang,/justind)

;aind=get_cind2(d0,d0,cenmasks[0]+cenmasks[2],cenmasks[0])
aind=get_aind(d0,d0,cenmasks[0],cenmasks[1]-cenmasks[0])

;Intersection of Spider and radial extent = region used for Image Registration
rois=intersect(spi,aind)

endif
;*******************************

;***Make file saving the centroid positions (just for accounting purposes), the inner saturation radius (don't worry),
;the hour angle (need it for LOCI/ADI), and the Parallactic Angle (*Definitely* need for LOCI/ADI)

;Writing out the reduc.log file header and first entry
openw,1,'reduc.log'
printf,1,'File_no,','XC,','YC,','Rsat','HA','Par.Angle'

;to save output flux withi some aperture
if keyword_set(fluxout) then begin
openw,2,'fluxout.dat'
endif

;to save core to halo flux ratio
if keyword_set(fratio) then begin
openw,22,'fluxratio.dat'
endif

;take initial frame and grab the filter name
filtname=strtrim(sxpar(h0,'CAL_BAND'))
if  filtname eq 'H' then begin

sat_xguess=[64,112,137,89]
sat_yguess=[114,138,89,64]

endif else begin
;if not H then default to low res spots for now
sat_xguess=[70,108,128,90]
sat_yguess=[110,129,90,71]
endelse

;You probably want to set for array element 0 for now.
if ~keyword_set(refslice) then refslice=0

for i=0L,nfiles-1 do begin

print,'Registering File Number ',i+1,' ',files[i]
;read in file
a=readfits(datadir+files[i],h1,/exten,/silent)
h0=headfits(datadir+files[i],ext=0,/silent)
a_in=a
sz=size(a)


parangle0=float(sxpar(h0,'PA'))
ha0=float(sxpar(h0,'HA'))

if ~keyword_set(nosats) then begin
get_charis_wvlh,h0,wvlhs
scl=wvlhs[refslice]/wvlhs
;do an initial centroid measurement from channel 1
a_firstchan=a[*,*,0]
a_firstchan-=filter_image(a_firstchan,median=15)
 gcntrd,a_firstchan,sat_xguess[0],sat_yguess[0],satcx0,satcy0,3
 gcntrd,a_firstchan,sat_xguess[1],sat_yguess[1],satcx1,satcy1,3
 gcntrd,a_firstchan,sat_xguess[2],sat_yguess[2],satcx2,satcy2,3
 gcntrd,a_firstchan,sat_xguess[3],sat_yguess[3],satcx3,satcy3,3

 cens0=[[satcx0,satcy0],[satcx1,satcy1],[satcx2,satcy2],[satcx3,satcy3]]
; cens0x=median([satcx0,satcx1,satcx2,satcx3],/even)
; cens0y=median([satcy0,satcy1,satcy2,satcy3],/even)
; cens0=[cens0x,cens0y]
 
 cens=dblarr(2,4,sz[3])
 cens[*,*,0]=cens0
 c0=(total(cens0,2)/4) # (fltarr(4)+1.)
 cens0p=cv_coord(from_rect=cens0-c0,/to_polar)
 
;sat_centroids=fltarr(4,2,n_elements(wvlh))
for ii=0L,sz[3]-1 do begin
 cens[*,*,ii]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[ii]],/to_rect)+c0

 a_test=a[*,*,ii]
medbox=15
 a_test-=filter_image(a_test,median=medbox)  ; this is roughly 5 lambda/D = safe
 a_test=smooth(a_test,3)
;do gaussian centroid estimate for satellite spot positions.  
 gcntrd,a_test,cens[0,0,ii],cens[1,0,ii],satcx0,satcy0,3
 gcntrd,a_test,cens[0,1,ii],cens[1,1,ii],satcx1,satcy1,3
 gcntrd,a_test,cens[0,2,ii],cens[1,2,ii],satcx2,satcy2,3
 gcntrd,a_test,cens[0,3,ii],cens[1,3,ii],satcx3,satcy3,3
writefits,'test.fits',a_test

 cens[0,0,ii]=satcx0
 cens[1,0,ii]=satcy0
 cens[0,1,ii]=satcx1
 cens[1,1,ii]=satcy1
 cens[0,2,ii]=satcx2
 cens[1,2,ii]=satcy2
 cens[0,3,ii]=satcx3
 cens[1,3,ii]=satcy3
endfor

PSFcens=fltarr(2,sz[3])
deltax=fltarr(sz[3])
deltay=fltarr(sz[3])

;initial stuff for satellite spot locations
for s=0L,sz[3]-1 do begin
; for j=0L,3 do begin
;  sxaddpar,h1,'SATS'+strtrim(s,2)+'_'+strtrim(j,2),$
;   string(strtrim(cens[*,j,s],2),format='(F7.3," ",F7.3)'),$
;   'Location of sat. spot '+strtrim(j,2)+' of slice '+strtrim(s,2)
; endfor

;set centroid position for each slice to be (median(x_positions),median(y_positions))
PSFcens[*,s]=[(median(cens[0,*,s],/even)),(median(cens[1,*,s],/even))]
;sxaddpar,h1,'PSFC_'+strtrim(s,2),$
; string(strtrim(PSFcens[*,s],2),format='(F7.3," ",F7.3)'),$
;  'PSF Center of slice '+strtrim(s,2)
deltax[s]=xc0-PSFcens[0,s]
deltay[s]=yc0-PSFcens[1,s]
endfor

slices=findgen(sz[3])-refslice
fitx=poly_fit(slices[2:sz[3]-1-2],PSFcens[0,2:sz[3]-1-2]-PSFcens[0,refslice],2)
fity=poly_fit(slices[2:sz[3]-1-2],PSFcens[1,2:sz[3]-1-2]-PSFcens[1,refslice],2)

refoffsetx=PSFcens[0,refslice]-xc0
refoffsety=PSFcens[1,refslice]-yc0

;for s=0L,sz[3]-1 do begin
;cubic=-0.5 interpolation to shift centroid to d0/2, d0/2

deltax=fitx[0]+fitx[1]*(slices-refslice)+fitx[2]*(slices-refslice)^2.+refoffsetx
deltay=fity[0]+fity[1]*(slices-refslice)+fity[2]*(slices-refslice)^2.+refoffsety

for s=0L,sz[3]-1 do begin
a_in[*,*,s]=shift_sub(a_in[*,*,s],-1*deltax[s],-1*deltay[s])
;***new stuff
cens[0,*,s]-=deltax[s]
cens[1,*,s]-=deltay[s]
PSFcens[*,s]=[(median(cens[0,*,s],/even)),(median(cens[1,*,s],/even))]
if keyword_set(verbose) then print,'PSFCENS is ',i,PSFcens[*,s]
;stop
 for j=0L,3 do begin
  sxaddpar,h1,'SATS'+strtrim(s,2)+'_'+strtrim(j,2),$
   string(strtrim(cens[*,j,s],2),format='(F7.3," ",F7.3)'),$
    'Location of sat. spot '+strtrim(j,2)+' of slice '+strtrim(s,2)
 endfor

sxaddpar,h1,'PSFC_'+strtrim(s,2),$
 string(strtrim(PSFcens[*,s],2),format='(F7.3," ",F7.3)'),$
  'PSF Center of slice '+strtrim(s,2)
endfor


;endfor

psfcenx=median(PSFcens[0,*],/even)
psfceny=median(PSFcens[1,*],/even)
sxaddpar,h1,"PSFCENTX",psfcenx,'Mean PSF center X' 
sxaddpar,h1,"PSFCENTY",psfceny,'Mean PSF center Y' 
if keyword_set(verbose) then print,'Median PSF Center is ',psfcenx,psfceny
a=a_in

endif else begin

;redefine the image to be a dummy variable

;if no satellite spots, collapse the cube along wavelength dimension to better illuminate halo.
;- really should do this on speckle-aligned images.  But hard to do if you don't know the centroid in the first place.
;--  Do this without alignment for now.

for ii=0L,sz[3]-1 do begin
im_in=a[*,*,ii]
if ~keyword_set(simple) then begin
offset=getrelshift1(fimage_ref,im_in,rois,0,0)
endif else begin
offset=getrelshift2(fimage_ref,im_in,rois,0,0)
endelse

deltax=offset[0] & deltay=offset[1]

;print,'xc yc', xc0-offset[0],yc0-offset[1]
;for ii=0L,sz[3]-1 do begin
 a_in[*,*,ii]=shift_sub(a_in[*,*,ii],deltax,deltay)
endfor

psfcenx=xc0-deltax
psfceny=yc0-deltay
a=a_in

endelse


writefits,reducdir+filesout[i],0,h0
writefits,reducdir+filesout[i],a,h1,/append
printf,1,filenum[i],psfcenx,psfceny,rsat,ha0,parangle0

endfor
;sxaddpar,h1,'SATS'+strtrim(long(ii),2)+'_'+strtrim(j,2),strtrim(satpos[0,j,ii]+deltax,2)+' '+strtrim(satpos[1,j,ii]+deltay,2)
;print,'SATS'+strtrim(long(ii),2)+'_'+strtrim(j,2),satpos[0,j,ii],satpos[0,j,ii]+deltax,satpos[1,j,ii],satpos[1,j,ii]+deltay
;print,deltax,deltay

close,1

end
