pro charis_specphot_cal,pfname,pick=pick,r_ap=r_ap,datacube=datacube,calcube=calcube,modamp=modamp,$
ap_factor=ap_factor,$
;attenfac=attenfac,ap_factor=ap_factor,$
filter=fname,starname=starname,mag=mag,$
prefname=prefname,suffname=suffname,test=test,$
verbose=verbose,$
fluxunits=fluxunits,outfilename=outfilename,help=help,guide=guide

;Spectro-Photometric Calibration program (work in progress)
; - uses images with satellite spots, info file with spectral type and brightness of star in H band
; - feeds this information into 'calculation' subroutine to get a flux-calibrated spectrum of your star
; - computes satellite spot attenuation based on the amplitude of modulation and the MJD date
; - flux calibrates data cubes based on model spectrum and on sat spots
;-...'calculation.pro subroutine done -7/12/2017
;-assumes same exposure time right now -7/12/2017
;
;*to implement ...
; switches on 'attenfac', 'starname','filter','mag', etc ... to ...
; override the choice of star name, filter, brightness

if (N_PARAMS() eq 0 or keyword_set(help) or keyword_set(guide)) then begin

Print,'Calibrate Photometric Flux subroutine'
Print,'This program does flux calibration for a CHARIS datacube or sequence of data cubes'
Print,' using data of a known star taken with a similar setup'
Print,'  '
Print,'****Calling Sequence ****'
Print,'  '
Print,'charis_specphot_cal,pfname,pick=pick,r_ap=r_ap,datacube=datacube,calcube=calcube,modamp=modamp,attenfac=attenfac,ap_factor=ap_factor'
Print,'filter=fname,starname=starname,mag=mag,prefname=prefname,suffname=suffname,test=test,verbose=verbose,fluxunits=fluxunits,outfilename=outfilename,help=help'
Print,' '
Print,'Parameters ...'
Print,'pfname - input info file for applying correction to pre-determined sequence of registered images'
Print,'pick - keyword to instead manually select images to be calibrated'
Print,'r_ap - aperture radius for calibration'
Print,'datacube - Data to Calibrate'
Print,'calcube - Datacube used for Flux Calibration'
Print,'pref/suffname - Manually override the prefix and suffix of file names used (useful for the pfname switch but applied to images other than "reg" images)'
Print,'test - just do a test, do not override file name, stop after one'
print,'funit - What Flux Density units? 0=mJy,1=Jy,2=W/m^2/um,3=ergs/s/cm^2/A,4=ergs/s/cm^2/Hz,5=contrast'
print,''
print,'Limitations ...'
print,'Need to implement override switches on attenfac, starname,filtername, star mag ...'
goto,skiptoend
endif

;applies a photometric calibration either to an entire sequence of files or to one or more identified image

;default is to apply to sequence of registered but not PSF subtracted images
;allow for override by choosing 'pick'


;***Preliminary***
;define data directory
reducdir='./reduc/'

;data directory
datadir=reducdir+'reg/'
subdir='reg/'
reducdir+=subdir

;information about the target properties.
param,'HMAG',hmag,/get,pfname=pfname
param,'EHMAG',ehmag,/get,pfname=pfname
param,'SPTYPE',spt,/get,pfname=pfname

;identify the aperture radius if you haven't already
;correction*: just use the FWHM of the image
;if ~keyword_set(r_ap) then param,'raper',r_ap,/get,pfname=pfname

;if you want to keep the aperture equal to the image slice FWHM then do nothing.  If you want to change the scaling set the ap_factor keyword
if ~keyword_set(ap_factor) then ap_factor = 1



;*****Identify calibration datacube.
;
;This must be a datacube with visible satellite spots.
;You can opt to choose the same science datacubes or some subset of them as well.

if ~keyword_set(calcube) then $
calcube=dialog_pickfile(Title="Select Your Datacube Used for Flux Calibration")
;***for now just allow for one cube. change this later!!!

if ~keyword_set(modsize) then begin
;if look for the ASTROGRID modulation amplitude in the fits header
test=readfits(calcube,h1test,ext=1)
h0test=headfits(calcube)

modsize1=sxpar(h0test,'X_GRDAMP',count=modcount1)
modsize2=sxpar(h1test,'X_GRDAMP',count=modcount2)

if (modcount1 eq 0 and modcount2 eq 0) then begin
read,"Manually Enter the ASTROGRID modulation amplitidue in nanometers ",modsize
endif else begin
if modcount1 eq 0 then  modsize=modsize2*1d3
if modcount2 eq 0 then  modsize=modsize1*1d3
endelse

endif 

;***star to satellite contrast ratio*********
;*Assumptions: the wavelength range for the cal file is the same as the science data you are wanting to calibrate.
;;             the calibration hasn't changed.  If it has you will have to enter the right value manually.
;;This has been determined for the 1550nm spot using the latest values for modulation. 
;;;the precision is good to ~2% over a factor of 2 mod size (~factor of 4 brightness)
;;;;this strictly speaking has only been determined for the low-res mode.  Need to test high res.

;First, get the CHARIS wavelengths
get_charis_wvlh,h0test,wvlh_test

;Now get the time of observation in modified Julien Days
mjd_in=sxpar(h0test,'MJD',count=mjdcount)
if mjdcount eq 0 then begin  ;if for some reason MJD is not printed, then enter this manually
Read,'Are these data taken after August 28 2017 (if yes, type 0)?',mjdanswer
if mjdanswer eq 0 then mjd_in=58000.0
endif

;Now from the modulation size, get the CHARIS attenuation factor/channel
attenfac=get_charis_sat_contrast(modsize,wvlh_test*1d-3,mjd_in)




;*****Identify which datacubes you want to Flux-Calibrate

;*If you want to manually set it then do this

;if you are calling the to-be-flux-cal'd file from command line or calling this program from the 1d_spectrum program ...
if keyword_set(datacube) then begin
files_full_path=datacube
nfiles=n_elements(files_full_path)
files=files_full_path
filesbase=files_full_path

;pick out the directory under which your files reside.  WARNING: assumes that all files are in same subdirectory!!!!!
your_path=strpos(files_full_path[0],'/',/reverse_search)+1

;and put the new files exactly where they were before; when you open the file later, look in the same place.
reducdir=strmid((files_full_path[0]),0,your_path)
datadir=reducdir


;otherwise, select the datacubes from the list in the .info file or select them via GUI
;if ~keyword_set(pick) then begin
;if ~keyword_set(suffname) then suffname='reg'

;;***trimming off path
;for i=0L,nfiles-1 do begin
;z=strsplit(files_full_path[i],'/',/extract)
;files[i]=z[n_elements(z)-1]
;zbase=strsplit(files[i],'.fits',/extract)
;filesbase[i]=zbase[0]  
;endfor

;***trimming off path
for i=0L,nfiles-1 do begin
z=strsplit(files_full_path[i],'/',/extract)
files[i]=z[n_elements(z)-1]
suffposition=strpos(files[i],'.fits',/reverse_search)
filesbase[i]=strmid(files[i],0,suffposition)
endfor


suffname2='_cal'
goto,breakout
endif

if ~keyword_set(pick) then begin
;assume we are working in the 'reg' subdirectory
param,'fnum_sat',flist,/get,pfname=pfname
filenum=nbrlist(flist)

if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then suffname='reg'

;input files
suffname2=suffname+'_cal'
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)


;assume a simple overwrite
filesout=filelist(filenum,nfiles,prefix=prefname,suffix=suffname2)

goto,breakout
endif else begin

;******Hack-tastic inelegant solution until I figure out a cleaner way of parsing the IDL path.
files_full_path=dialog_pickfile(Title="Select Datacubes to Flux Calibrate",/multiple)
nfiles=n_elements(files_full_path)
files=files_full_path
filesbase=files_full_path

;pick out the directory under which your files reside.  WARNING: assumes that all files are in same subdirectory!!!!!
your_path=strpos(files_full_path[0],'/',/reverse_search)+1

;and put the new files exactly where they were before; when you open the file later, look in the same place.
reducdir=strmid((files_full_path[0]),0,your_path)
datadir=reducdir
;print,'REDUCDIR ',reducdir

;if ~keyword_set(suffname) then suffname='reg'

;***trimming off path
for i=0L,nfiles-1 do begin
z=strsplit(files_full_path[i],'/',/extract)
files[i]=z[n_elements(z)-1]
suffposition=strpos(files[i],'.fits',/reverse_search)
filesbase[i]=strmid(files[i],0,suffposition)
endfor



suffname2='_cal'

goto,breakout

;assume a simple overwrite
;filesout=filelist(filenum,nfiles,prefix=prefname,suffix=suffname2)


endelse

breakout:

;end
;*****Now perform photometry on the calimage, record its integration time
;***open the file, pull the fits header from the file, measure the integrated signal, measure the noise, do snrcalc/wavelength slice.
ncal=n_elements(calcube)
;test calimage
imcaltest=readfits(calcube[0],hcaltest,/exten)
dimcal=(size(imcaltest,/dim))

;not doing the 'goodcode hex2bin stuff with GPI. assume all sat spots are good

satflux=fltarr(ncal,dimcal[2])
esatflux=fltarr(ncal,dimcal[2])

for ical=0L,ncal-1 do begin
h0cal=headfits(calcube[ical])
imcal=readfits(calcube[ical],h1cal,ext=1)

texp=sxpar(h0cal,'exp1time')
ncoadd=sxpar(h0cal,'coadds')
tint_cal=texp*float(ncoadd)

;***Get the CHARIS wavelengths***
get_charis_wvlh,h0cal,wvlh_charis
;fwhm=0.206*wvlh_charis*1d-3/7.9d0
Dtel=7.9d0
pixscale=0.0164
fwhm=1.0*(1.d-6*(wvlh_charis*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale
if keyword_set(verbose) then print,'fwhm is',fwhm

;Now nominally assume that the aperture radius is 1/2 the FWHM
aperrad=double(fwhm*0.5)
if keyword_set(verbose) then print,aperrad[0],fwhm[0]
;top

;sat=sxpar(h1,'SATMASK')
;goodcode=hex2bin(sat,(size(image,/dim))[2])
;good=long(where(goodcode eq 1))
cens=fltarr(2,4,dimcal[2])

;loop on wavelength
 for s=0L,dimcal[2]-1 do begin
;for s=0L,n_elements(good)-1 do begin

  for j=0,3 do begin  ;assume 4 sat spots
   tmp=fltarr(2)+!values.f_nan
    sf=string(s)
    tmp= sxpar(h1cal,'SATS'+strtrim(long(sf),2)+'_'+strtrim(j,2))  
;**total, ugly hack-tastic workaround since I can't seem to get IDL to treat the string properly otherwise...
    tmpz=strsplit(tmp,' ',/extract)
    tmp=[tmpz[0],tmpz[1]]
    cens[*,j,s]=tmp
  endfor
 endfor

;**now you have the positions of the satellite spots.  time to extract the photometry.
;assume that you don't have any errors triggered like can happen with the GPI pipeline

; extract the flux of the satellite spots
  sat1flux = fltarr(n_elements(cens[0,0,*]))    ;;top left
  sat2flux = fltarr(n_elements(cens[0,0,*]))    ;;bottom left
  sat3flux = fltarr(n_elements(cens[0,0,*]))    ;;top right
  sat4flux = fltarr(n_elements(cens[0,0,*]))    ;;bottom right
  sat1noise = fltarr(n_elements(cens[0,0,*]))    ;;top left
  sat2noise = fltarr(n_elements(cens[0,0,*]))    ;;bottom left
  sat3noise = fltarr(n_elements(cens[0,0,*]))    ;;top right
  sat4noise = fltarr(n_elements(cens[0,0,*]))    ;;bottom right
  ;mean_sat_flux = fltarr(n_elements(cens[0,0,*]))
  ;stdev_sat_flux=fltarr(n_elements(cens[0,0,*]))

for il=0L,dimcal[2]-1 do begin
;do a radial profile subtraction of the image.  Should get close to zero'ing the background.
;assume for now that the satellite spot signal does not affect the median profile calculation.
imslice=imcal[*,*,il]
profrad_tc,imslice,1,1,dimcal[0]/2,p2d=pr
imslice-=pr
;writefits,'cal.fits',imslice
;stop

aperradf=aperrad[il]*ap_factor
if keyword_set(verbose) then print,aperradf,fwhm[il],fwhm[il]/2.,ap_factor,aperrad[il]*ap_factor,aperrad[il]

;aperrad=0.5*fwhm[il]*ap_factor
;aperture photometry at each slice assuming a zero background
phpadu=1.0
sat_skyrad=[2,6]*aperradf

if keyword_set(verbose) then  print,cens[0,0,il],cens[1,0,il],imslice[cens[0,0,il],cens[1,0,il]]
  ;stop
aper, imslice,cens[0,0,il],cens[1,0,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat1flux[il]=flux
 ;print,'flux and sky1 is',flux,sky
aper, imslice,cens[0,1,il],cens[1,1,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat2flux[il]=flux
 ;print,'flux and sky2 is',flux,sky
aper, imslice,cens[0,2,il],cens[1,2,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat3flux[il]=flux
 ;print,'flux and sky3 is',flux,sky
aper, imslice,cens[0,3,il],cens[1,3,il],flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
     sat4flux[il]=flux
 ;print,'flux and sky4 is',flux,sky
;stop
       imslicemask=imslice
;SNR calculation per slice to get the background rms within an aperture, perform this on images with the sat spots masked out.
       dist_circle,mask1,dimcal[0],[cens[0,0,il],cens[1,0,il]]
       dist_circle,mask2,dimcal[0],[cens[0,1,il],cens[1,1,il]]
       dist_circle,mask3,dimcal[0],[cens[0,2,il],cens[1,2,il]]
       dist_circle,mask4,dimcal[0],[cens[0,3,il],cens[1,3,il]]
       masksat=where(mask1 le aperradf or mask2 le aperradf or mask3 le aperradf or mask4 le aperradf,nmasksat)
       imslicemask[masksat]=!values.f_nan
       ;help,imslicemask
       snratio_sub,imslicemask,fwhm=aperradf*2,/finite,noisemap=noisemap_slice,/silent

       sat1noise[il]=noisemap_slice(cens[0,0,il],cens[1,0,il])
       sat2noise[il]=noisemap_slice(cens[0,1,il],cens[1,2,il])
       sat3noise[il]=noisemap_slice(cens[0,2,il],cens[1,2,il])
       sat4noise[il]=noisemap_slice(cens[0,3,il],cens[1,3,il])


  ;      stddev_sat_flux=  fltarr(n_elements(cens[0,0,*]))

       stdev_sat_in=0.25*sqrt(sat1noise[il]^2.+sat2noise[il]^2+sat3noise[il]^2.+sat4noise[il]^2.)
       stdev_sat_sys=stdev([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]])
       ;stdev_sat_flux[ical,il]=sqrt(stdev_sat_in^2.+stdev_sat_sys^2.)
       ;mean_sat_flux[ical,il]=mean([sat1flux[il]],[sat2flux[il]],[sat3flux[il]],[sat4flux[il]])

       esatflux[ical,il]=sqrt(stdev_sat_in^2.+stdev_sat_sys^2.)
       ;esatflux[ical,il]=sqrt(stdev_sat_in^2.+stdev_sat_sys^2.)

       satflux[ical,il]=mean([sat1flux[il],sat2flux[il],sat3flux[il],sat4flux[il]])
       ;print,mean_sat_flux[ical,il],stdev_sat_flux[ical,il]
       if keyword_set(verbose) then print,satflux[ical,il],esatflux[ical,il],satflux[ical,il]/esatflux[ical,il]
endfor
endfor

if ncal gt 1 then begin

satfluxf=reform(mean(satflux,dimension=1))
esatfluxf=reform(mean(esatflux,dimension=1))
endif else begin
satfluxf=reform(satflux)
esatfluxf=reform(esatflux)
endelse

snratio_sats=satfluxf/esatfluxf
if keyword_set(verbose) then print,snratio_sats

if keyword_set(verbose) then begin
window,1
;ploterror,wvlh_charis*1.d-3,satflux,esatflux,xrange=[1,3]
ploterror,wvlh_charis*1.d-3,satflux*(wvlh_charis/wvlh_charis[0])^(-1*2.),esatflux,xrange=[1,3]
oplot,wvlh_charis*1.d-3,satflux[0]*(wvlh_charis/wvlh_charis[0])^(-1*4.),linestyle=1
endif

      


funitlist=[0,1,2,3,4,5]   ;which corresponds to ...
funitname=['mJy','Jy','W/m^2/um','ergs/s/cm^2/A','ergs/s/cm^2/Hz','contrast']

if ~keyword_set(fluxunits) then begin
fluxunits_cube= funitlist[0]
funitname_output=funitname[0]
endif
if keyword_set(fluxunits) then begin
fluxunits_cube=funitlist[fluxunits]
funitname_output=funitname[fluxunits]
endif

;***Flux-Calibrated Model Spectrum
if fluxunits_cube le 4 then begin
 ;Physical Units
;speccal=charis_photometric_calibration_calculation(h0cal,h1cal,wvlh_charis*1d-3,spectype=spt,star_mag=hmag,units=fluxunits_cube,/diagnostic)
speccal=charis_photometric_calibration_calculation(h0cal,h1cal,wvlh_charis*1d-3,spectype=spt,star_mag=hmag,units=fluxunits_cube)
endif else begin
 ;Contrast Units
speccal=findgen(n_elements(wvlh_charis))*0+1
endelse

;print,satflux
;plot,wvlh_charis*1d-3,satfluxf
;stop
;*****Flux-Calibration in Physical Units
;;;conversion factor/channel: multiply the cube slices by this
conv_fact=(1.0/satfluxf)*speccal*attenfac
if keyword_set(verbose) then print,'satfluxf ', satfluxf
print,'speccal',speccal
print,'conv',conv_fact


;;;flux-cal-ing the data cubes

for icube=0L,nfiles-1 do begin

 a=readfits(datadir+files[icube],ext=1,h1)
 h0=headfits(datadir+files[icube])
 print,'the file name is',files[icube]
 print,'the data directory is ',datadir
 for ichannel=0L,n_elements(wvlh_charis)-1 do begin
 a[*,*,ichannel]*=conv_fact[ichannel]
 ;print,'hi2!',median(a[*,*,ichannel])
 ;sxaddpar,h1,'SatSpot SNR',snratio_sats[ichannel]
 endfor

 ;funitname_output=strtrim(funitname_output)

sxaddpar,h1,'FLUXUNIT',strtrim(funitname_output,1),"Flux Density Units in Cube"

;aperture
for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'R_AP'+strtrim(ichannel,2),aperrad[ichannel]*ap_factor,"Aperture Radius (Pixels)"

;sky annulus
;for ichannel=0L,n_elements(wvlh_charis)-1 $
; do sxaddpar,h1,'R_AP'+strtrim(ichannel,2),aperrad[ichannel]*ap_factor,"Aperture Radius (Pixels)"
sxaddpar,h1,'Sky_In',strtrim(strc(2)),"Inner Radius for Sky Annulus (in Units of Aperture Radius)"
sxaddpar,h1,'Sky_Out',strtrim(strc(4)),"Outer Radius for Sky Annulus (in Units of Aperture Radius)"

;scaling
for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'FSCALE'+strc(ichannel),conv_fact[ichannel],"scale to convert counts to "+strc(funitname_output)

;error in spectrophotometric calibration/slice
for ichannel=0L,n_elements(wvlh_charis)-1 $
 do sxaddpar,h1,'CERR'+strc(ichannel),1/snratio_sats[ichannel],"Cal Fractional error for slice "+strc(ichannel)

;temporary hack of file path 
 if (keyword_set(pick) or keyword_set(datacube)) then begin
 print,'file! ',reducdir+filesbase[icube]+suffname2+'.fits'
 writefits,reducdir+filesbase[icube]+suffname2+'.fits',0,h0
 writefits,reducdir+filesbase[icube]+suffname2+'.fits',a,h1,/append
 
 if keyword_set(datacube) then outfilename=reducdir+filesbase[icube]+suffname2+'.fits'
 
 endif else begin
 writefits,reducdir+filesout[icube],0,h0
 writefits,reducdir+filesout[icube],a,h1,/append
 endelse

endfor
skiptoend:
close,/all
 
end 
