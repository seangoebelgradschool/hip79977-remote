pro charis_extract_1d_spectrum,pfname,datacube=datacube,coords=coords,fitcollapse=fitcollapse,fitslice=fitslice,planefit=planefit,zerosky=zerosky,noprad=noprad,filt_slice=filt_slice,$
pick=pick,r_ap=r_ap,calcube=calcube,$
attenfac=attenfac,ap_factor=ap_factor,$
filter=fname,starname=starname,mag=mag,$
ndfilt=ndfilt,$
prefname=prefname,suffname=suffname,test=test,$
fluxunits=fluxunits,help=help
;pro charis_extract_1d_spectrum,pfname,datacube=datacube,pick=pick,r_ap=r_ap,calcube=calcube,$
;attenfac=attenfac,ap_factor=ap_factor,$
;filter=fname,starname=starname,mag=mag,$
;ndfilt=ndfilt,$
;prefname=prefname,suffname=suffname,test=test,$
;fluxunits=fluxunits,help=help

;Built off of GPI's spectral extraction program but probably more robust.

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

if keyword_set(datacube) then begin
extractcube=datacube
endif else begin
extractcube=dialog_pickfile(Title="Select the Data Cube From Which You Want to Extract a Spectrum")
endelse

;read in data cube, if the cube has been flux-calibrated then proceed, if not then do flux calibration

test=readfits(extractcube,h1,ext=1)

caltest=sxpar(h1,'FLUXUNIT',count=ct)

;If 'FLUXUNIT' has been found, then your data cube is flux-calibrated. proceed with extracting a spectrum.  If it has not been found, then do flux calibration first.


if ct eq 0 then begin    
;********porting charis_specphot_cal.pro here*****
;*************************************************

charis_specphot_cal,pfname,pick=pick,r_ap=r_ap,datacube=extractcube,calcube=calcube,$
attenfac=attenfac,ap_factor=ap_factor,$
filter=fname,starname=starname,mag=mag,$
ndfilt=ndfilt,$
prefname=prefname,suffname=suffname,test=test,$
fluxunits=fluxunits,help=help,outfilename=outfilename
print,'outfilename is',outfilename

extractube=outfilename
print,'extractcube is ..',extractcube
;stop
;the outfile is now the file from which you will extract a cube

endif 


;****Extract a Spectrum*********
;*******************************

;1. Define position
if ~keyword_set(coords) then begin
read,"Input X Coordinate ",xpos
read,"Input Y Coordinate ",ypos
endif else begin
xpos=coords[0]  ;x position
ypos=coords[1]  ;y position
endelse
xposf=xpos
yposf=ypos

;2. Read in the datacube
data=readfits(extractcube,h1,ext=1)
h0=headfits(extractcube)
get_charis_wvlh,h0,wvlh_charis ;pull wavelengths
Dtel=7.9d0
pixscale=0.0164
fwhm=1.0*(1.d-6*(wvlh_charis*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale

nlambda=(size(data,/dim))[2]
aperrad=fltarr(nlambda)

;3. Define extraction radius
 for ir=0L,nlambda-1 do begin
  aperrad[ir]=sxpar(h1,'R_AP'+strtrim(string(ir),2))
 endfor
  ;print,aperrad
  ;stop

;4.If you want to fit the position of the companion in the collapsed cube, then throw the fitcollapse switch
if keyword_set(fitcollapse) then begin
data_col=median(data,dimension=3,/even)
med_wvlh=median(wvlh_charis,/even)
fwhm_col=1.0*(1.d-6*(med_wvlh*1d-3)/Dtel)*(180.*3600./!dpi)/pixscale
gcntrd,data_col,xpos,ypos,xposf,yposf,fwhm_col

endif


fluxdens_spectrum=fltarr(nlambda)
efluxdens_spectrum=fltarr(nlambda)

for il=0L,nlambda-1 do begin
fwhm_slice=fwhm[il]
;aperradf=0.5*fwhm_slice
aperradf=aperrad[il]
sat_skyrad=aperradf*[2,4]

if keyword_set(filt_slice) then begin
noisyslice=data[*,*,il]
noisyslice-=filter_image(noisyslice,median=10*fwhm_slice)
data[*,*,il]=noisyslice
endif

if ~keyword_set(noprad) then begin
profrad_tc,data[*,*,il],p2d=pr
data[*,*,il]-=pr
endif


;if you want to fit the companion position per wavelength slice then throw the fitslice switch
if keyword_set(fitslice) then begin
gcntrd,data[*,*,il],xpos,ypos,xposf,yposf,fwhm_slice
endif

;Aperture Photometry
phpadu=1.0 ;dummy number
if ~keyword_set(zerosky) then begin
aper,data[*,*,il],xposf,yposf,flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent
endif else begin
aper,data[*,*,il],xposf,yposf,flux,eflux,sky,skyerr,phpadu,aperradf,sat_skyrad,[0,0],/flux,/exact,/nan,/silent,setskyval=0.
endelse

;now compute the SNR properly
snratio_sub,data[*,*,il],fwhm=2*aperradf,coord=[xposf,yposf],snrval=snrval,/finite,/zero,/silent

fluxdens_spectrum[il]=flux
efluxdens_spectrum[il]=fluxdens_spectrum[il]/snrval

endfor

writecol,'spectrum.dat',wvlh_charis*1d-3,fluxdens_spectrum,efluxdens_spectrum

set_plot,'ps'
device,filename='spectrum.eps',bits=8,/encapsulated
plot,wvlh_charis*1d-3,fluxdens_spectrum,xrange=[1,2.5],/nodata,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)'
oploterror,wvlh_charis*1d-3,fluxdens_spectrum,efluxdens_spectrum,errstyle=1,linestyle=0,thick=4
device,/close

set_plot,'x'
plot,wvlh_charis*1d-3,fluxdens_spectrum,xrange=[1,2.5],/nodata,xthick=5,ythick=5,xtitle='Wavelength (Microns)',ytitle='Flux Density (mJy)'
oploterror,wvlh_charis*1d-3,fluxdens_spectrum,efluxdens_spectrum,errstyle=1,linestyle=0,thick=4


end
