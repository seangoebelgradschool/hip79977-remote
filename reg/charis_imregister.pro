pro charis_imregister,pfname,prefname=prefname,suffname=suffname,$
simple=simple,smallsteps=smallsteps,bandaid=bandaid,$
guessoffsets=guessoffsets,$
refslice=refslice,$
nosats=nosats,$
xcorr=xcorr,$
check=check,$
subreg=subreg,$
smask=smask,$
;satfit=satfit,$
;fluxnorm=fluxnorm,$
;fratio=fratio,$
;rcore=rcore,$
;rhalo=rhalo,$
;fluxout=fluxout,$
rewrite=rewrite,verbose=verbose

;09/27/2017 - Switch to estimate spot positions in registered images where spot positions cannot be calculated
;09/07/2017 - Switch to escape the spot placements to -9999 in case of failure.  Allows to you to see whether some frames are just off
;	or whether there was a centroid adjustment later.
;08/13/2017 - Option to find satellite spot positions from x-corr of region of interests instead of gaussian centroiding.
;07/10/2017 - Renamed to charis_imregister.pro
;07/10/2017 - Now at least for the sat-spot = true case, saves the newsatspot positions in fits header.
;04/13/2017 - Added switch to do x-corr fitting of the halo instead of satellite spots
;04/08/2017 - Redone for CHARIS to find the satellite positions
;03/09/2014 - Quick hack-version for GPI image registration across cube
; components.  

;*****Nominal Program (no switches)*****
;1. Uses an initial guess for the header positions from the reference slice that is hardwired for now.
;2. Computes centroid position for this reference slice.  
;3. Goes on to compute satellite spot positions at other wavelengths.
;4. Computes the centroid from the spot positions in each slice
;5. Fits a 3rd-order polynominal to the centroid vs. slice
;6. Shifts image slice to put the centroid at dimx/2 dimy/2

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

;***source mask file****
if keyword_set(smask) then begin
readcol,'smask.dat',xsource,ysource
endif

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
;filtnamelist=['
case filtname of 
'H': begin
       sat_xguess=[64,112,137,89]
       sat_yguess=[114,138,89,64]
     end
'K': begin
  
       sat_xguess=[50,116,150,84]
       sat_yguess=[117,150,83,49]
     end
 else: begin

;endif else begin
;if not highres then default to low res spots for now
         sat_xguess=[70,108,128,90]
         sat_yguess=[110,129,90,71]
       end
endcase

if keyword_set(guessoffsets) then begin
sat_xguess+=guessoffsets[0]
sat_yguess+=guessoffsets[1]
endif
;endelse

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

;center the pa value on zero
if parangle0 gt 180 then parangle0-=360

;*****Image Registration Loop.
;**- If you have satellite spots then you have a choice of getting the centroid from gaussian-fitting or from joint cross-correlation

if ~keyword_set(nosats) then begin
get_charis_wvlh,h0,wvlhs
scl=wvlhs[refslice]/wvlhs


;****do an initial centroid measurement from channel 1
;***assumes that the reference channel sat spots are visible.

a_firstchan=a[*,*,0]
a_firstchan-=filter_image(a_firstchan,median=15)
 gcntrd,a_firstchan,sat_xguess[0],sat_yguess[0],satcx0,satcy0,3
 gcntrd,a_firstchan,sat_xguess[1],sat_yguess[1],satcx1,satcy1,3
 gcntrd,a_firstchan,sat_xguess[2],sat_yguess[2],satcx2,satcy2,3
 gcntrd,a_firstchan,sat_xguess[3],sat_yguess[3],satcx3,satcy3,3

 
if (satcx0 lt 0 or satcy0 lt 0) then begin 
gcntrd,a_firstchan,sat_xguess[0],sat_yguess[0],satcx0,satcy0,3,/keepcenter
endif
if (satcx1 lt 0 or satcy1 lt 0) then begin 
print,'v1',satcx1,satcy1
gcntrd,a_firstchan,sat_xguess[1],sat_yguess[1],satcx1,satcy1,3,/keepcenter
print,'v2',satcx1,satcy1
endif
if (satcx2 eq -1 or satcy2 eq -1) then begin 
gcntrd,a_firstchan,sat_xguess[2],sat_yguess[2],satcx2,satcy2,3,/keepcenter
endif
if (satcx3 eq -1 or satcy3 eq -1) then begin 
gcntrd,a_firstchan,sat_xguess[3],sat_yguess[3],satcx3,satcy3,3,/keepcenter
endif

 cens0=[[satcx0,satcy0],[satcx1,satcy1],[satcx2,satcy2],[satcx3,satcy3]]
 ;cens0x=median([satcx0,satcx1,satcx2,satcx3],/even)
 ;cens0y=median([satcy0,satcy1,satcy2,satcy3],/even)
 ;cens0f=[cens0x,cens0y]
;print,cens0f
 
 cens=dblarr(2,4,sz[3])
 cens[*,*,0]=cens0
 c0=(total(cens0,2)/4) # (fltarr(4)+1.)
 cens0p=cv_coord(from_rect=cens0-c0,/to_polar)
 ;print,c0
 ;stop
initial_centroid=median(cens0,dimension=2,/even)

;sat_centroids=fltarr(4,2,n_elements(wvlh))
dist_circle,g,[sz[1],sz[2]]

PSFcens=fltarr(2,sz[3])
deltax=fltarr(sz[3])
deltay=fltarr(sz[3])

for ii=0L,sz[3]-1 do begin
 cens[*,*,ii]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[ii]],/to_rect)+c0

 a_test=a[*,*,ii]
medbox=15
 a_test-=filter_image(a_test,median=medbox)  ; this is roughly 5 lambda/D = safe
 a_test=smooth(a_test,3)

;****if you have a very bright pt source that is skewing the sat spot determination, apply an NAN mask here.
if keyword_set(smask) then begin
for sm=0L,n_elements(x)-1 do begin
gcntrd,a_test,xsource[sm],ysource[sm],xsourceout,ysourceout
dist_circle,sloc,[sz[1],sz[2]],xsourceout,ysourceout
maskme=where(sloc le 2)
a_test[maskme]=!values.f_nan
endfor
endif

;***
;****1. Cross-Correlation Fitting to find the star center. 2. Gaussian Centroiding.
;do gaussian centroid estimate for satellite spot positions.  

if keyword_set(xcorr) then begin
;***Cross-Correlation Fitting


;calculating satellite spot locations
;for s=0L,sz[3]-1 do begin

;*Define region of interest
;sat 1
dist_circle,sat1,[sz[1],sz[2]],cens[0,0,ii],cens[1,0,ii]

;sat 2
dist_circle,sat2,[sz[1],sz[2]],cens[0,1,ii],cens[1,1,ii]

;sat 3
dist_circle,sat3,[sz[1],sz[2]],cens[0,2,ii],cens[1,2,ii]

;sat 4
dist_circle,sat4,[sz[1],sz[2]],cens[0,3,ii],cens[1,3,ii]

search_radius=2.
roi=where((sat1 le search_radius) or (sat2 le search_radius) or (sat3 le search_radius) or (sat4 le search_radius),nroi)
;roi=where(sat1 le search_radius)

;gah=fltarr(sz[0],sz[1])
;gah[roi]=1.
;writefits,'roi.fits',sat1
;stop
;print,PSFcens[*,ii],xc0,yc0
;initial_centroid=[sz[1]/2,sz[2]/2]
print,initial_centroid
;stop

if ~keyword_set(smallsteps) then begin
fidcenter=charis_getcenter(a_test,initial_centroid[0],initial_centroid[1],roi)

endif else begin
fidcenter=charis_getcenter(a_test,initial_centroid[0],initial_centroid[1],roi,/smallsteps)
endelse
;print,fidcenter[0],fidcenter[1]
;stop

PSFcens[0,ii]=fidcenter[0]
PSFcens[1,ii]=fidcenter[1]

deltax[ii]=xc0-PSFcens[0,ii]
deltay[ii]=yc0-PSFcens[1,ii]
;endfor
;stop
print,'Center from X-Corr Fitting is ',PSFcens[*,ii],' at slice',ii
endif else begin

;*Gaussian Centroiding

 gcntrd,a_test,cens[0,0,ii],cens[1,0,ii],satcx0,satcy0,3
 gcntrd,a_test,cens[0,1,ii],cens[1,1,ii],satcx1,satcy1,3
 gcntrd,a_test,cens[0,2,ii],cens[1,2,ii],satcx2,satcy2,3
 gcntrd,a_test,cens[0,3,ii],cens[1,3,ii],satcx3,satcy3,3
;writefits,'test.fits',a_test

if (satcx0 lt 0 or satcy0 lt 0) then begin 
gcntrd,a_test,cens[0,0,ii],cens[1,0,ii],satcx0,satcy0,3,/keepcenter
endif
if (satcx1 lt 0 or satcy1 lt 0) then begin 
;print,'v1',satcx1,satcy1
cntrd,a_test,cens[0,1,ii],cens[1,1,ii],satcx1,satcy1,3,/keepcenter
print,'v2',satcx1,satcy1
;stop
endif
if (satcx2 lt 0 or satcy2 lt 0) then begin 
gcntrd,a_test,cens[0,2,ii],cens[1,2,ii],satcx2,satcy2,3,/keepcenter
endif
if (satcx3 lt 0 or satcy3 lt 0) then begin 
gcntrd,a_test,cens[0,3,ii],cens[1,3,ii],satcx3,satcy3,3,/keepcenter
endif


 cens[0,0,ii]=satcx0
 cens[1,0,ii]=satcy0
 cens[0,1,ii]=satcx1
 cens[1,1,ii]=satcy1
 cens[0,2,ii]=satcx2
 cens[1,2,ii]=satcy2
 cens[0,3,ii]=satcx3
 cens[1,3,ii]=satcy3

;kill switch
if (satcx0 lt 0 or satcy0 lt 0 or satcx1 lt 0 or satcy1 lt 0 or satcx2 lt 0 or satcy2 lt 0 or satcx3 lt 0 or satcy3 lt 0) then begin
cens[*,*,ii]=!values.f_nan
endif
; print,satcx0,satcy0,satcx1,satcy1,satcx2,satcy2,satcx3,satcy3
;endfor

;PSFcens=fltarr(2,sz[3])
;deltax=fltarr(sz[3])
;deltay=fltarr(sz[3])

;initial stuff for satellite spot locations
;for s=0L,sz[3]-1 do begin
; for j=0L,3 do begin
;  sxaddpar,h1,'SATS'+strtrim(s,2)+'_'+strtrim(j,2),$
;   string(strtrim(cens[*,j,s],2),format='(F7.3," ",F7.3)'),$
;   'Location of sat. spot '+strtrim(j,2)+' of slice '+strtrim(s,2)
; endfor

;set centroid position for each slice to be (median(x_positions),median(y_positions))
PSFcens[*,ii]=[(median(cens[0,*,ii],/even)),(median(cens[1,*,ii],/even))]
;sxaddpar,h1,'PSFC_'+strtrim(s,2),$
; string(strtrim(PSFcens[*,s],2),format='(F7.3," ",F7.3)'),$
;  'PSF Center of slice '+strtrim(s,2)
deltax[ii]=xc0-PSFcens[0,ii]
deltay[ii]=yc0-PSFcens[1,ii]
;endfor

print,'Center from Gaussian Fitting is ',PSFcens[*,ii],' at slice',ii

endelse

endfor

;****Polynomial fit to individual slice centroid measurements****
slices=findgen(sz[3])-refslice
;good=where(finite(PSFcens[0,*]) ne 0 and finite(PSFcens[1,*]) ne 0,ngood)
good=where(finite(PSFcens[0,2:sz[3]-1-2]) ne 0 and finite(PSFcens[1,2:sz[3]-1-2]) ne 0,ngood)
print,'good is', ngood

;fitx=poly_fit(slices[2:sz[3]-1-4],PSFcens[0,2:sz[3]-1-4]-PSFcens[0,refslice],2)
;fity=poly_fit(slices[2:sz[3]-1-4],PSFcens[1,2:sz[3]-1-4]-PSFcens[1,refslice],2)
;fitx=robust_poly_fit(slices[2:sz[3]-1-2],PSFcens[0,2:sz[3]-1-2]-PSFcens[0,refslice],1,numit=15)
;fity=robust_poly_fit(slices[2:sz[3]-1-2],PSFcens[1,2:sz[3]-1-2]-PSFcens[1,refslice],1,numit=15)

;fitx=robust_poly_fit(slices[2:sz[3]-1-2],PSFcens[0,2:sz[3]-1-2]-PSFcens[0,refslice],2,numit=15)
;fity=robust_poly_fit(slices[2:sz[3]-1-2],PSFcens[1,2:sz[3]-1-2]-PSFcens[1,refslice],2,numit=15)

;****for now put in break to exit in case the program cannot find any centroids
if ngood gt 0 then begin
fitx=robust_poly_fit((slices[2:sz[3]-1-2])[good],(PSFcens[0,2:sz[3]-1-2])[good]-PSFcens[0,refslice],2,numit=15)
fity=robust_poly_fit((slices[2:sz[3]-1-2])[good],(PSFcens[1,2:sz[3]-1-2])[good]-PSFcens[1,refslice],2,numit=15)

;fitx=poly_fit(slices[2:sz[3]-1-2],PSFcens[0,2:sz[3]-1-2]-PSFcens[0,refslice],2)
;fity=poly_fit(slices[2:sz[3]-1-2],PSFcens[1,2:sz[3]-1-2]-PSFcens[1,refslice],2)

refoffsetx=PSFcens[0,refslice]-xc0
refoffsety=PSFcens[1,refslice]-yc0

;for s=0L,sz[3]-1 do begin
;cubic=-0.5 interpolation to shift centroid to d0/2, d0/2

deltax=fitx[0]+fitx[1]*(slices-refslice)+fitx[2]*(slices-refslice)^2.+refoffsetx
deltay=fity[0]+fity[1]*(slices-refslice)+fity[2]*(slices-refslice)^2.+refoffsety

print,refoffsetx,xc0,PSFcens[0,refslice]
;stop
;plot,slices,PSFcens[0,*],yrange=[min(PSFcens[0,*])-2.5,max(PSFcens[0,*]+2.5)],ystyle=1
plot,slices,PSFcens[0,*],yrange=[sz[1]/2-7,sz[1]/2+7],ystyle=1
oplot,slices,deltax+xc0,linestyle=1

endif else begin
deltax[*]=-9999
deltay[*]=-9999
refoffsetx=-9999
refoffsety=-9999
endelse


;stop

;if you cannot find the sat spot position in some slice then use first channel sat position and extrapolate
; cens[*,*,ii]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[ii]],/to_rect)+c0
cenxref=cens[0,*,0]-deltax[0]
cenyref=cens[1,*,0]-deltay[0]

for s=0L,sz[3]-1 do begin
a_in[*,*,s]=shift_sub(a_in[*,*,s],-1*deltax[s],-1*deltay[s])
;***new stuff
cens[0,*,s]-=deltax[s]
cens[1,*,s]-=deltay[s]
PSFcens[*,s]=[(median(cens[0,*,s],/even)),(median(cens[1,*,s],/even))]



;switch if sat spot position is NaN ...
centotal=total(cens[*,*,s])
;cenref=([(median(cens[0,*,0],/even)),(median(cens[1,*,0],/even))]) # (fltarr(4)+1.)
; c0=(total(cens0,2)/4) # (fltarr(4)+1.)
;print,cenref
cenref=(total(cens[*,*,0],2)/4) # (fltarr(4)+1.)

if finite(centotal) eq 0 then begin
cens0=cens[*,*,0]
;print,centotal
;print,cenref
;print,'cenx',cenref[0,*]
;print,'ceny',cenref[1,*]
 cens0p=cv_coord(from_rect=cens0-cenref,/to_polar)
 ;cens[*,*,s]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[s]],/to_rect)+cenref
 cens[*,*,s]=cv_coord(from_polar=[cens0p[0,*],cens0p[1,*]/scl[s]],/to_rect)+cenref
;print,cens[*,*,s]
;print,c0
;stop
PSFcens[*,s]=[(median(cens[0,*,s],/even)),(median(cens[1,*,s],/even))]
endif

 for j=0L,3 do begin

  sxaddpar,h1,'SATS'+strtrim(s,2)+'_'+strtrim(j,2),$
   string(strtrim(cens[*,j,s],2),format='(F7.3," ",F7.3)'),$
    'Location of sat. spot '+strtrim(j,2)+' of slice '+strtrim(s,2)
 endfor

sxaddpar,h1,'PSFC_'+strtrim(s,2),$
 string(strtrim(PSFcens[*,s],2),format='(F7.3," ",F7.3)'),$
  'PSF Center of slice '+strtrim(s,2)

if keyword_set(verbose) then print,'PSFCENS is ',i,PSFcens[*,s],' shifts are ',deltax[s],deltay[s]

endfor



psfcenx=median(PSFcens[0,*],/even)
psfceny=median(PSFcens[1,*],/even)
sxaddpar,h1,"PSFCENTX",psfcenx,'Mean PSF center X' 
sxaddpar,h1,"PSFCENTY",psfceny,'Mean PSF center Y' 
if keyword_set(verbose) then print,'Median PSF Center is ',psfcenx,psfceny
a=a_in

endif else begin

PSFcens=fltarr(2,sz[3])
deltax=fltarr(sz[3])
deltay=fltarr(sz[3])

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

;deltax=offset[0] & deltay=offset[1]

deltax[ii]=offset[0]
deltay[ii]=offset[1]
PSFcens[0,ii]=xc0-deltax[ii]
PSFcens[1,ii]=yc0-deltay[ii]
;print,'xc yc', xc0-offset[0],yc0-offset[1]
;for ii=0L,sz[3]-1 do begin
; a_in[*,*,ii]=shift_sub(a_in[*,*,ii],deltax,deltay)
endfor

slices=findgen(sz[3]-refslice)
fitx=robust_poly_fit(slices[2:sz[3]-1-2],PSFcens[0,2:sz[3]-1-2]-PSFcens[0,refslice],2)
fity=robust_poly_fit(slices[2:sz[3]-1-2],PSFcens[1,2:sz[3]-1-2]-PSFcens[1,refslice],2)

refoffsetx=PSFcens[0,refslice]-xc0
refoffsety=PSFcens[1,refslice]-yc0

deltax=fitx[0]+fitx[1]*(slices-refslice)+fitx[2]*(slices-refslice)^2.+refoffsetx
deltay=fity[0]+fity[1]*(slices-refslice)+fity[2]*(slices-refslice)^2.+refoffsety
;print,deltax,deltay

;cens=dblarr(2,4,sz[3])

for ii=0L,sz[3]-1 do begin
a_in[*,*,ii]=shift_sub(a_in[*,*,ii],-1*deltax[ii],-1*deltay[ii])
;***new stuff
;cens[0,*,ii]-=deltax[ii]
;cens[1,*,ii]-=deltay[ii]
PSFcens[0,ii]=xc0-deltax[ii]
PSFcens[1,ii]=yc0-deltay[ii]

;PSFcens[*,ii]=[(median(cens[0,*,ii],/even)),(median(cens[1,*,ii],/even))]
if keyword_set(verbose) then print,'PSFCENS is ',i,PSFcens[0,ii],PSFcens[1,ii],' shifts are ',deltax[ii],deltay[ii]

sxaddpar,h1,'PSFC_'+strtrim(ii,2),$
 string(strtrim(PSFcens[*,ii],2),format='(F7.3," ",F7.3)'),$
  'PSF Center of slice '+strtrim(ii,2)
endfor 

print,PSFcens[0,*]
;stop
plot,slices,PSFcens[0,*],yrange=[sz[1]/2-5,sz[1]/2+5],ystyle=1
oplot,slices,deltax+xc0,linestyle=1

psfcenx=median(PSFcens[0,*],/even)
psfceny=median(PSFcens[1,*],/even)
sxaddpar,h1,"PSFCENTX",psfcenx,'Mean PSF center X'
sxaddpar,h1,"PSFCENTY",psfceny,'Mean PSF center Y'
if keyword_set(verbose) then print,'Median PSF Center is ',psfcenx,psfceny
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
