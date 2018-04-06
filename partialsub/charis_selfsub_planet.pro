pro charis_selfsub_planet,pfname,reducname=reducname,pixmask=pixmask,corr=corr,prefname=prefname,psf=psf,empflux=empflux,$
magfile=magfile

;**09/04/2017**
;****INCOMPLETE CODE****
;Version 1.0 - for CHARIS data, check line 473 about HA and PARANG.  
;**08/24/2017**
;Version 0.1 - for CHARIS data, simplified

;**07/25/2015**
;Version 0.1 - using empty data cube to calibrate throughput in GPI data pipeline.  
;Right now, it is simple.
;Assume ...
;1. single assumption for planet PSF (e.g. no protoplanets)
;2. one planet only.

;**07/22/2015**
;Version 1.0 - Throughput loss/calibration
;This program will supplement the code that calculates throughput loss corrections with fake point sources.
;It uses the engine from selfsub_locidisk, where we save the LOCI coefficients in each region and apply them
; to an empty data cube to figure out throughput loss. Compare to real data to figure out throughput loss/photometry.
;it will work provided that you perform local masking on the subtraction zone, otherwise negative dark speckles get 
;put on planet to remove it: end up underestimating planet brightness by some amount.
;removes ...
;'the big overhead stuff' for IRCS.  assume you'll never need this again

;**05/4/2014**
;Version 0.1.1 - Code to mimick Marois' SOSIE PSF phot/astro calibration.  E.g. remove the planets from the reference PSF.
;For now, do just for a single planet... -- TC
;Assumes all files are of the same dimension with a common center and are square arrays
;update: good through line 377

;*****************************

;setupdir,reducdir=reducdir

;determine reduction subdirectory
;pos1=strpos(pfname,'/')+1
;pos2=strpos(pfname,'.',/reverse_search)

reducdir='./reduc/'
;subdir='procsub/'
subdir='proc/'
reducdir1=reducdir+'reg/'
datadir=reducdir1
reducdir+=subdir

reducdirorg=reducdir
;*****

if ~keyword_set(prefname) then begin
;***Prefixes and such***
prefname='n'
endif

suffname='reg'
suffname2='rsub'

if ~keyword_set(reducname) then begin

;file name, with full path
reducname_full_path=dialog_pickfile(Title="Select Your Processed Image")

your_path=strpos(reducname_full_path,'/',/reverse_search)+1

;determining the long-form of the reduction directory
reducdirorg=strmid((reducname_full_path),0,your_path)

;now determining the base name of your reduced file
z=strsplit(reducname_full_path,'/',/extract)
reducname=z[n_elements(z)-1]
endif

;***stuff for CHARIS***
;Telescope Diameter for Subaru
Dtel=7.9d0  ;effective aperture of Subaru
;pixel scale 16.4 mas/pixel
pixscale=0.0164

;****the data cubes
    hrefcube=headfits(reducdirorg+reducname,ext=0)
    refcube=readfits(reducdirorg+reducname,ext=1,h1)
    reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
    hrefcol=headfits(reducdirorg+reducname_col,ext=0)
    refcol=readfits(reducdirorg+reducname_col,ext=1,h1col)

dimx=sxpar(h1,'naxis1')
dimy=sxpar(h1,'naxis2')
xc=dimx/2 & yc=dimy/2

;Now Get the Wavelength Vector
get_charis_wvlh,hrefcube,lambda
lambda*=1d-3
nlambda=n_elements(lambda)
;determine array of FWHM
fwhm=1.0*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale

;;***Now get basic image properties

header=headfits(reducdirorg+reducname,ext=1)
header0=headfits(reducdirorg+reducname)
dimx=sxpar(header,'naxis1')
dimy=sxpar(header,'naxis2')
xc=dimx/2 & yc=dimy/2

;ALOCI parameters from fits header in output file
znfwhm=sxpar(header,'loci_nfw')
zna=sxpar(header,'loci_na')
zgeom=sxpar(header,'loci_geo')
zdrsub=sxpar(header,'loci_drs')
zrmin=sxpar(header,'loci_rmi')
zrmax=sxpar(header,'loci_rma')

zcorrlim=sxpar(header,'corr_lim',count=countcor)
zcutoff=sxpar(header,'svd',count=countcutoff)
znref=sxpar(header,'nref',count=countnref)

;how much of the image do you want to reduce?
;**for now, don't set these. 
;if ~keyword_set(redmax) then redmax=zrmax
;if ~keyword_set(redmin) then redmin=zrmin

;define the temporary directory
tmpdir=reducdir+'tmp/'

;parameters of the sequence
param,'obsdate',date,/get,pfname=pfname & date=strtrim(date,2)
param,'fnum_sat',flist,/get,pfname=pfname

filenum=nbrlist(flist)
files=filelistf(filenum,nfiles,prefix=prefname,suffix=suffname)
filesfc=filelistf(filenum,nfiles,prefix=prefname,suffix=suffname+'_fc')


lat=1.*sxpar(header0,'lat')
lng=1.*sxpar(header0,'lng')

;reading hour angle and parallactique angle at the beginning of exposures

readcol,'reduc.log',filenum,allxc,allyc,allrsat,allha,allpa



;******Define the Position of the Planet****
;******

if ~keyword_set(magfile) then begin
file_mag=dialog_pickfile(Title="Select Position and Photometry File")
readcol,file_mag,xp,yp,mag_p
endif else begin
readcol,mag,xp,yp,mag_p
endelse
xp_o=xp
yp_o=yp

        ;Define number of planets you have to simulate
        nplanet=n_elements(xp)

	;Brightness of the planet: an array of length n_lambda
        flux_planet=fltarr(nlambda,nplanet)

;;****** Brightness of the Planet ****

        ;loop on planets
        for iplanet=0L,nplanet-1 do begin
	;Option 1 - Just read whatever value you have from the file (do nothing), assume constant brightness across filters, assume pseudo-magnitude units
        if ~keyword_set(empflux) then begin
          flux_planet[*,iplanet]=10^(0.4*(25-mag_p[iplanet]))

	;Option 2 - Look in the reduced image, measure the brightness of the planet, assume some throughput loss, get predicted planet brightness
	endif else begin
         ;loop on wavelength to measure brightness via aperture photometry
            phpadu=1.0
          for ilin=0L,nlambda-1 do begin
            imslice=refcube[*,*,ilin]

             ;for now, size the aperture to be equal to the point source FWHM
            aper,imslice,xp[iplanet],yp[iplanet],flux,eflux,sky,skyerr,phpadu,0.5*fwhm[ilin],[2,4]*fwhm[ilin],setskyval=0,/flux,/exact,/nan,/silent
            ;aper,imslice,xp,yp,flux,eflux,sky,skyerr,phpadu,0.5*fwhm[ilin],[2,4]*fwhm[ilin],setskyval=0,/flux,/exact,/nan,/silent
            flux_planet[ilin,iplanet]=flux
            print,'magplanet!',lambda[ilin],flux
          endfor ;ilin

           ;for now, assume 75% throughput as a fiducial value
           attenfac=0.75 
           flux_planet/=attenfac
         
        endelse  
             plot,lambda,flux_planet[*,iplanet],xtitle='Wavelength (Microns)',ytitle='Flux Density'
        endfor  ;iplanet

;opening to use
param,'raper',raper,/get,pfname=pfname


;loop over input magnitudes/input guessed magnitudes
header=headfits(reducdirorg+reducname,ext=1)
testfile=readfits(reducdirorg+reducname,ext=1)
image_dim=(size(testfile))[1]

dx =fltarr(n_elements(xp))
dy=fltarr(n_elements(yp))
rc=fltarr(n_elements(xp))
for i=0L,n_elements(xp)-1 do begin
;get radius position of planet
dx[i]=(xp[i]-image_dim/2)
dy[i]=(yp[i]-image_dim/2)

print,'coords!',[dx,dy]+image_dim/2
rc[i]=sqrt(dx[i]^2.+dy[i]^2.)
endfor


;***now define PSF size 
;***PSF
dpsf=21
intr_psf=fltarr(dpsf,dpsf,nlambda)
;***
print,'loop size is ',nplanet

;****Define Your Prior: The Intrinsic PSF****
;****

;1.
;use an empirical PSF

if keyword_set(psf) then begin

;-automatically
fname=findfile('psfcube_med.fits',count=c)

if c eq 0 then begin
;-manually
fname=dialog_pickfile(Title="Pick Reference PSF Files",/MULTIPLE_FILES)
c=n_elements(fname)
endif

;loop in wavelength
for il=0,nlambda-1 do begin
;Since you have coronographic data and/or saturated data

; best to just use a model PSF, e.g. c=0
;if c gt 0 then begin
    fpeak=0.
    for n=0,c-1 do begin
        im=readfits(fname[n],ext=1)
        im=im[*,*,il]
        im=subarr(im,dpsf)
;make the PSF a unit PSF: signal of 1 within a FWHM-sized aperture
        im/=myaper_tc(im,dpsf/2,dpsf/2,0.5*fwhm[il])
        if n eq 0 then psf=im/c else psf+=im/c
    endfor
   intr_psf[*,*,il]=psf
endfor

;2.
;use a simple gaussian
endif else begin

for il=0L,nlambda -1 do begin
;print,il,il,lambda[il]


;approximate the PSF by a gaussian; should be okay for high Strehl and simple PSF structure (e.g. Subaru or VLT, not Keck)
a=psf_gaussian(npixel=dpsf,fwhm=fwhm[il],/double,/normalize)
a/=myaper(a,dpsf/2,dpsf/2,0.5*fwhm[il])

intr_psf[*,*,il]=a
endfor
endelse
;stop

writefits,'psf.fits',intr_psf



;Now you have flux defined, PSF shape defined.

;x,y coordinates of pixels of the PSF centered on pixel 0.0

xpsf=lindgen(dpsf)#replicate(1l,dpsf)-dpsf/2
ypsf=replicate(1l,dpsf)#lindgen(dpsf)-dpsf/2

;indices of pixels of the psf centered on pixel 0.0 in the big picture
ipsf=xpsf+ypsf*dimx

;Now estimate your "real" input PSF scaled to match the flux of your object

;***diagnostic
input_planets=fltarr(dimx,dimy,nlambda,nplanet)

;***for now, this solves a stupid numerical problem.  So the 'input flux' is not perfectly passed on.  With this fix some input flux is to get attenuation right.

for irc=0L,nplanet -1 do begin


for il=0L,nlambda-1 do begin
imgpsf=fltarr(dimx,dimy)
imgpsf[*]=0

;in_mag=mag_p[irc]
;in_flux=10^(0.4*(25.-in_mag))
in_flux=flux_planet[il,irc]
;print,in_flux
xe=floor(xp_o[irc]) & ye=floor(yp_o[irc])
dxf=xp_o[irc]-xe & dyf=yp_o[irc]-ye

psfs=shift_sub(intr_psf[*,*,il],dxf,dyf)
;psfs=translate(intr_psf[*,*,il],dxf,dyf)

imgpsf[ipsf+xe+ye*dimx]=psfs*in_flux
input_planets[*,*,il,irc]+=imgpsf
aper,input_planets[*,*,il],xp[irc],yp[irc],fluxo,efluxo,sky,skyerr,1,0.5*fwhm[il],[2,6]*fwhm[il],setskyval=0,/flux,/exact,/nan,/silent
;gcntrd,input_planets[*,*,il],xp[irc],yp[irc],xout,yout,fwhm[il]
;aper,input_planets[*,*,il],xout,yout,fluxo,efluxo,sky,skyerr,1,0.5*fwhm[il],[2,4]*fwhm[il],setskyval=0,/flux,/exact,/nan,/silent
;print,'input and output fluxes are ',flux_planet[il],fluxo,fluxo/flux_planet[il],xp[irc],yp[irc],xout,yout
;print,'input and output positions are ',xp[irc],yp[irc],xout,yout
;print,fluxo/flux_planet[il,irc],fluxo,flux_planet[il,irc]
flux_planet[il,irc]=fluxo
endfor
endfor
;stop
writefits,'inputplanet.fits',input_planets
;stop

;+**fabrique la psf etiree a rayon=rc
;stretched to the psf radius rc=

imt=fltarr(dpsf,dpsf,nplanet)
nsc=21 ;nombre de sous-compagnon

;for n=0,nfiles-1 do begin

;*****Keywords for RA,DEC, individual exposure time (exp1time), and coadds (so exp1time*ncoadds = exptime)*****
    ;h=headfits(datadir+files[n],ext=1)
    ;radeg=sxpar(h,'ra')
    ;decdeg=sxpar(h,'dec')
    param,'RA',radeg,/get,pfname=pfname
    param,'DEC',decdeg,/get,pfname=pfname

    ;exptime=float(sxpar(h,'exp1time'))
    ;exp1time=exptime
    ;ncoadds=float(sxpar(h,'coadds'))


;lines 376 to 
;for now, just keep these as is.

res_flux=fltarr(nlambda,nplanet)  ;residual flux (1-attenuation)
diff_r=fltarr(nlambda,nplanet)  ;astrometric bias in angular separation
diff_ang=fltarr(nlambda,nplanet) ;astrometric bias in position angle

;loop calculates the subtraction residuals for fake planets of a given brightness

;define the brightness steps

flux=fltarr(nplanet)
mag_flux=fltarr(nplanet)

;Calculate the nominal companion flux within aperture radius,raper

    print,' Adding planet PSF to empty cubes...'
;  "adding companions.."

for n=0,nfiles-1 do begin

        h0=headfits(datadir+files[n],ext=0,/silent)
        im=readfits(datadir+files[n],h,/silent,ext=1)

;number of wavelength channels
        nwvlh=(size(im,/dim))[2]

;set to empty data cube
        imff=im
        im[*]=1.d-15
        ;parallactic angle of sub companions

        exptime=float(sxpar(h0,'exp1time'))
        ncoadds=float(sxpar(h0,'coadds'))
        ha=allha[n]

;***check these things...9/5/2017-TC
        ;x=ha+(findgen(nsc)/(nsc-1.))*(exptime*ncoadds/3600.) ;SI HA DEBUT POSE
        dha=(findgen(nsc)/(nsc-1.))*(exptime*ncoadds/3600.) ;SI HA DEBUT POSE
        dha-=max(dha)/2   ;center the array of hour angles upon the mean hour angle
        x=ha+dha
        pa0=parangle(x,decdeg,lat)
        ;dpa=pa0-pa0[0]
         dpa=pa0-median(pa0,/even)

;*** charis records the mean parallactic angle.  so dpa is a range centered on allpa[n]
         ;print,'dpa',dpa
         ;print,'decdeg',decdeg
         pa=allpa[n]+dpa
         ;print,'pa is',pa
         ;stop

;loop over wavelength
 for il=0L,nwvlh-1 do begin

  for irc=0,nplanet-1 do begin
      ;use the planet flux input computed well above

            ;determine angles des sous-compagnons

            plan_coords=cv_coord(/double,from_rect=[dx[irc],dy[irc]],/to_polar,/degrees)
  
;for normal data 
            ;asc=(northpa+(plan_coords[0]+dpa-max(dpa)/2)+(pa-allpa[0]))*!dtor  

            ;asc=(plan_coords[0]+dpa-max(dpa)/2)*!dtor
            ;asc=(plan_coords[0]+max(dpa)/2-dpa-(pa-allpa[0]))*!dtor
            asc=(plan_coords[0]-(pa-allpa[0]))*!dtor


            ;determine x,y positions of sub-companions

            xsc=rc[irc]*cos(asc)+xc
            ysc=rc[irc]*sin(asc)+yc

                img=fltarr(dimx,dimy)
                img[*]=0.
                ;writefits,'img.fits',img
              ;stop
            for isc=0,nsc-1 do begin
                xe=floor(xsc[isc]) & ye=floor(ysc[isc])
                dxf=xsc[isc]-xe    & dyf=ysc[isc]-ye

                ;psfs=translate(intr_psf[*,*,il],dxf,dyf)

                psfs=shift_sub(intr_psf[*,*,il],dxf,dyf)
                ;psfs-=min(psfs)
                 
; for some reason, you have to do this in two steps with IFS data
                ;img[ipsf+xe+ye*dimx]=psfs*flux_planet[il,irc]/smear_fac[irc]/nsc
                img[ipsf+xe+ye*dimx]=psfs*flux_planet[il,irc]/nsc
                im[*,*,il]+=img
                
            endfor   ;isc
            ;writefits,'imm.fits',im[*,*,il]
      endfor ;irc
 endfor  ;il
;register image
        writefits,'imm.fits',im
   ;     stop
        writefits,datadir+filesfc[n],0,h0
        writefits,datadir+filesfc[n],im,h,/append
        im[*]=0
 endfor
;stop
    ;*** ADI treatment


;****ADI****
goto,skipstuf
;LOCI keywords**
header=headfits(reducdir+reducname)
znfwhm=sxpar(header,'loci_nfw')
zna=sxpar(header,'loci_na')
zgeom=sxpar(header,'loci_geo')
zdrsub=sxpar(header,'loci_drs')
zrmin=sxpar(header,'loci_rmi')
zrmax=sxpar(header,'loci_rma')
zcorrlim=sxpar(header,'corrlim',count=countcor)
;endif
skipstuf:

;good to here, 8/25/2017, 1828pm
if countcutoff eq 0 then zcutoff=1.d-99
print,'LOCI parameters are ',znfwhm,zdrsub,zna,zgeom,zrmin,zrmax
print,'countcor',countcor

itc=0
    charis_sublocirx,pfname,prefname=prefname,$
;nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,/svd,nref=znref,cutoff=zcutoff,outfile='res'+string(itc,format='(i1)')+'.fits',/usecoeff,/fc
nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,/svd,nref=znref,cutoff=zcutoff,outfile='res'+string(itc,format='(i1)')+'.fits',/usecoeff,/fc
;nfwhm=znfwhm,drsub=zdrsub,na=zna,geom=zgeom,rmin=zrmin,rmax=zrmax,/svd,cutoff=zcutoff,outfile='res'+string(itc,format='(i1)')+'.fits',/usecoeff,/fc
      

;Now read back in the files, rescale them to minimize chi-squared, add fits header information, and save them with unique input file name        
   
;your original file 
 ;  remember ...
    ;hrefcube=headfits(reducdir+reducname,ext=0)
    ;refcube=readfits(reducdir+reducname,ext=1,h1)
    ;reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
    ;hrefcol=headfits(reducdir+reducname_col,ext=0)
    ;refcol=readfits(reducdir+reducname_col,ext=1,h1col)


;the empty cube, after processing
    h0cube=headfits(reducdir+'res'+string(itc,format='(i1)')+'.fits')
    imcube=readfits(reducdir+'res'+string(itc,format='(i1)')+'.fits',h1cube,ext=1)

    modelname='res'+string(itc,format='(i1)')+'.fits'
    modelbasename=(strsplit(modelname,'.',/extract,count=modcount))[0]

;wavelength-collapsed version.
    h0col=headfits(reducdir+'res'+string(itc,format='(i1)')+'_collapsed.fits')
    imcol=readfits(reducdir+'res'+string(itc,format='(i1)')+'_collapsed.fits',h1col,ext=1)

;***now Loop on Wavelength and then Planet to get Attenuation vs. Channel 
   for ilambda=0L,nlambda-1 do begin 
     imslice=imcube[*,*,ilambda]
    for iplanet=0L,nplanet-1 do begin 
    
     ;compute the flux at the original position 
     aper,imslice,xp[iplanet],yp[iplanet],flux,eflux,sky,skyerr,1,0.5*fwhm[ilambda],[2,6]*fwhm[ilambda],setskyval=0,/flux,/exact,/nan,/silent

     ;now compare to the original, input flux to get the attenuation.
     res_flux[ilambda,iplanet]=flux/flux_planet[ilambda,iplanet]
     print,'Throughput for Planet ',string(iplanet+1),' at Wavelength ',string(lambda[ilambda]),' is ...',res_flux[ilambda,iplanet]
    endfor
   endfor

;***Throughput for wavelength collapsed image
    for iplanet=0L,nplanet-1 do begin
     aper,imcol,xp[iplanet],yp[iplanet],flux,eflux,sky,skyerr,1,0.5*median(fwhm,/even),[2,4]*median(fwhm,/even),setskyval=0,/flux,/exact,/nan,/silent 
    endfor
 

;this is the raw difference.
    writefits,'diff.fits',refcol-imcol
    writefits,'new.fits',imcol
    ;stop
    help,im2,im

end
