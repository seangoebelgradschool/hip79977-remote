pro charis_calc_final_contrast,pfname,pick=pick,r_ap=r_ap,datacube=datacube,calcube=calcube,modsize=modsize,nosub=nosub,$
attenfac=attenfac,ap_factor=ap_factor,$
prefname=prefname,suffname=suffname,outputfile=outputfile

;*****TO FINISH!!!
;****WARNING!!!!!!!!!!!!
;This version of the code is procedurally hard-wired!!!!!!
;  Assumes that spectrophotometric calibration has not been performed.  Need to re-write to be flexible and ensure someone cannot break code.
;v3.0 - Oct 11 2017 - copied from charis_calc_raw_contrast and charis_specphot_cal, rewritten to focus specifically on PSF subtracted images.
;v2.0 - Oct 4 2017 - Allows switch to compute the "contrast" as 5*halo brightness normalized by the satellite spots.  Sends output to diff directory
;v1.0 - Oct 1 2017 - reports contrast as rms of the noise (probably a bit dangerous since the residuals in a non-PSF subtracted, speckle filled image are non-gaussian)

file_mkdir,'final_contrasts'
outputdir='final_contrasts'


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


;if you want to keep the aperture equal to the image slice FWHM then do nothing.  If you want to change the scaling set the ap_factor keyword
if ~keyword_set(ap_factor) then ap_factor = 1

;*****Identify calibration datacube.
;
;This must be a datacube with visible satellite spots. Best to use sequence-combined, spatially-filtered cube 

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


;*****Identify files to estimate contrasts****
files_full_path=dialog_pickfile(Title="Select Datacubes For Contrast Computation",/multiple)
nfiles=n_elements(files_full_path)
files=files_full_path
filesbase=files_full_path

;pick out the directory under which your files reside.  WARNING: assumes that all files are in same subdirectory!!!!!
your_path=strpos(files_full_path[0],'/',/reverse_search)+1

;and put the new files exactly where they were before; when you open the file later, look in the same place.
reducdir=strmid((files_full_path[0]),0,your_path)
;`datadir=reducdir
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

;goto,breakout



;*********Satellite Flux Measurements********
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
sat_skyrad=[2,4]*aperradf

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
       masksat=where(mask1 le aperrad or mask2 le aperrad or mask3 le aperrad or mask4 le aperrad,nmasksat)
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

;*****Flux-Calibration in units of contrast
satfluxf=satfluxf/attenfac

;test calimage, getting dimensions.  Assumes that the dimensions of the cal cube and the PSF-sub cube are the same.
imcaltest=test
dimcal=(size(imcaltest,/dim))
profrad_tc,imcaltest[*,*,0],1,1,dimcal[0]/2,p1d=p1test,rayon=rayontest


;output contrasts
median_contrast=fltarr(nfiles,n_elements(p1test))
contrast_per_channel=fltarr(n_elements(p1test),dimcal[2],nfiles)
giantcube=fltarr(dimcal[0],dimcal[1],dimcal[2],nfiles)

noisemap_slicef=fltarr(n_elements(p1test),dimcal[2],nfiles)  ;hardwired for now, assume square arrays

for isci=0L,nfiles-1 do begin
print,'Calculating Final Contrast for File ',files[isci]
h0sci=headfits(reducdir+files[isci])
imsci=readfits(reducdir+files[isci],h1sci,ext=1)

;***Get the CHARIS wavelengths***
get_charis_wvlh,h0sci,wvlh_charis
Dtel=7.9d0
pixscale=0.0164
fwhm=1.0*(1.d-6*(wvlh_charis*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale
if keyword_set(verbose) then print,'fwhm is',fwhm

;Now nominally assume that the aperture radius is 1/2 the FWHM
aperrad=double(fwhm*0.5)

if ~keyword_set(ap_factor) then ap_factor = 1.0



for il=0L,dimcal[2]-1 do begin
;do a radial profile subtraction of the image.  Should get close to zero'ing the background.
;assume for now that the satellite spot signal does not affect the median profile calculation.
imslice=imsci[*,*,il]
profrad_tc,imslice,1,1,dimcal[0]/2,p2d=pr
imslice-=pr

aperradf=aperrad[il]*ap_factor


imslicemask=imslice

snratio_sub,imslicemask,fwhm=aperradf*2,/finite,noisemap=noisemap_slice,/silent,/zero,/filt

profrad_tc,noisemap_slice,1,1,100,rayon=rayon,p1d=dum

noisemap_slicef[*,il,isci]=dum

imsci[*,*,il]=imslicemask
giantcube[*,*,il,isci]=imslicemask

contrast_per_channel[*,il,isci]=5*noisemap_slicef[*,il,isci]/satfluxf[il]
print,'contrast and satflux are',contrast_per_channel[*,il,isci],satfluxf[il]
endfor

;calculate the contrast in the median_collapsed cube
imcol=median(imsci,dimension=3,/even)
medsatflux=median(satfluxf,/even)
snratio_sub,imcol,fwhm=median(aperrad*ap_factor,/even)*2.,/finite,noisemap=noisemap_col,/silent,/filt,/zero

;noise map of wavelength-collapsed image
profrad_tc,noisemap_col,1,1,100,p1d=noisemap_colf,rayon=rayon


median_contrast[isci,*]=5*noisemap_colf/medsatflux


writecol,'./'+outputdir+'/'+files[isci]+'_finalcontrast.txt',rayon*pixscale,5*noisemap_colf/medsatflux,$
5*noisemap_slicef[*,4,isci]/satfluxf[4],5*noisemap_slicef[*,11,isci]/satfluxf[11],5*noisemap_slicef[*,18,isci]/satfluxf[18],fmt='(f,f,f,f,f)'


endfor


end

