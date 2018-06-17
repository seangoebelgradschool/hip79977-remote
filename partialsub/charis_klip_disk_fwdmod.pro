pro charis_klip_disk_fwdmod,pfname,reducname=reducname,prefname=prefname,$
lrmin=lrmin,lrmax=lrmax,snrlim=snrlim,roirange=roirange,roimanual=roimanual,$
psf=psf,gausspsf=gausspsf,manscale=manscale,calcattenfact=calcattenfact,$
xext_fwhm=xext_fwhm,yext_fwhm=yext_fwhm,$
ffac=ffac,ntc=ntc,test=test,$
n_threads=n_threads, codeno=codeno, hardcodeparams=hardcodeparams, $
noreprocess=noreprocess ;sean's additions

;****2/5/2018**
;Version 2.0 - CHARIS/KLIP Forward-Modeling for Disks (complete).   
;		code is cleaned up, new structure using PSF models (rewrite of charis_makeemppsf)
;		new output of attenuated model
;		new switch (calcattenfact) to only output the model
;cube and attenuation factor if you throw it
; n_threads and codeno should be used if this code is being run
;multiple times simultaneously. They decide how the disk models are
;divided. If this is being run twice, then for the first instance,
;n_threads=2, codeno=1. For the second instance, n_threads=2, codeno=2.
;if the hardcodeparams flag isn't set, then the code prompts
;the user to select the model PSF and model disk.

;****1/29/2018**
;Version 1.1 - CHARIS/KLIP Forward-Modeling for Disks

;***10/13/2015
;Version 1.1
;include switch to use empirical PSF
;**

;*** 7/11/2015
;Version 1.0.1
;cleaned up code a bit.

;*** 3/14/2015
;Version 0.4
;edit roi to focus on the regions with high SNR detections from the disk.
;
;*** 3/11/2015
;Version 0.3
;forward-model grater image, compare with real image, compute chisq/min(chisq)scaling and best fit/flux-calibrated input image
;now edit to just do grater models!!!!!!!!
;
;*** 2/19/2015
;Version 0.2
;now modified to compare with a GRaTeR synthetic image.
;right now just does one single image

;**** 2/17/2015
;Version 0.1
;modified to do forward modeling of disks with GPI data
;where 'disk' equals just some gaussian you put in.
;hardwire this for GPI for right now
;
;**** 3/20/2012
;modified to model the self subtraction of a disk.
;*****PARAMETERS HARDWIRED for now

;if ~keyword_set(snrlim) then snrlim = 2.

if ~keyword_set(prefname) then begin
;***Prefixes and such***
prefname='n'
endif

suffname='reg'
suffname2='regcal'
;suffname2='rsub'

 
;*****************************

; calculates the attenuation from KLIP
; ntc = number of theta companion, default = 15
; / test to stop at the end of the loop to see
;And the residual flux for testing ntc

;setupdir,reducdir=reducdir

reducdir='./reduc/'
subdir='proc/'
reducdir1=reducdir+'reg/'
datadir=reducdir1
reducdir+=subdir

reducdirorg=reducdir

reducdir_model1=reducdir+'/model_in/'
reducdir_model2=reducdir+'/model_psfsub/'
reducdir_modelstat=reducdir+'/model_stat/'
reducdir_modelatten=reducdir+'/model_atten/'
file_mkdir,reducdir_model1
file_mkdir,reducdir_model2
file_mkdir,reducdir_modelstat
if keyword_set(calcattenfact) then file_mkdir,reducdir_modelatten
stat_filename = reducdir_modelstat+'fit_outcomes'+timestamp()+'.dat'
openw,unit,stat_filename, /get_lun ;Sean's note: the old version would fail to save the file if
; the code was stopped/crashed. This appends data, saving each time.
printf,unit,'MODEL_NAME','G','KSI0','ALP_I','ALP_O','BETA','XDO','YDO','E','R0','THETA0','CHISQ/dof',$
   format='(a91,a5,1x,9(a6,1x),a9)'
free_lun, unit

;******


;***stuff for CHARIS***
;Telescope Diameter for Subaru
Dtel=7.9d0  ;effective aperture of Subaru
;pixel scale 16.4 mas/pixel
pixscale=0.0164


if ~keyword_set(reducname) then begin
reducname='final.fits'
endif

;****the data cubes
    hrefcube=headfits(reducdirorg+reducname,ext=0)
    refcube=readfits(reducdirorg+reducname,ext=1,h1)
    reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
    hrefcol=headfits(reducdirorg+reducname_col,ext=0)
    refcol=readfits(reducdirorg+reducname_col,ext=1,h1col)


imx=sxpar(h1,'naxis1')
dimx=sxpar(h1,'naxis2')
dimy=dimx   ;assume square arrays
xc=dimx/2 & yc=dimy/2

;Now Get the Wavelength Vector
get_charis_wvlh,hrefcube,lambda
lambda*=1d-3
nlambda=n_elements(lambda)
;determine array of FWHM
fwhm=1.0*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale
medfwhm=median(fwhm,/even)

;nominally, we use medfwhm for the aperture in snrmap calculations and evaluating chi-sq.   
;option to change medfwhm to 2*aperture radius if need be but will require additional code changes where applicable.
;param,'raper',raper,/get,pfname=pfname


;since you're just a disk model at a single location.
nrc=1

;FOR DISK keep centered on central pixel!
rc=0.

;number of angles to position


header=headfits(reducdir+reducname,ext=1)

;***********

znfwhm=sxpar(header,'klip_nfw')
zdrsub=sxpar(header,'klip_drs')
zrmin=sxpar(header,'klip_rmi')
zrmax=sxpar(header,'klip_rma')
znpca=sxpar(header,'klip_pca')


;define the temporary directory
tmpdir=reducdir+'tmp/'


;parameters of the sequence
param,'obsdate',date,/get,pfname=pfname & date=strtrim(date,2)
param,'fnum_sat',flist,/get,pfname=pfname

filenum=nbrlist(flist)
files=filelistf(filenum,nfiles,prefix=prefname,suffix=suffname)
filesfc=filelistf(filenum,nfiles,prefix=prefname,suffix=suffname+'_fc')
;filesfc2=filelistf(filenum,nfiles,prefix=prefname,suffix=suffname2+'_fc')


;header of the first image to determine dimension

h=headfits(datadir+files[0],ext=1)
dimx=sxpar(h,'naxis1')
dimy=sxpar(h,'naxis2')
xc=dimx/2 & yc=dimy/2
;print,'xc yc',xc,yc

intr_psf=fltarr(dimx,dimy,nlambda)
;if ~keyword_set(psf) then begin
if keyword_set(gausspsf) then begin
for il=0L,nlambda -1 do begin
;print,il,il,lambda[il]

;determine FWHM
fwhm_channel=(1.d-6*lambda[il]/Dtel)*(180.*3600./!dpi)/pixscale

;approximate the PSF by a gaussian; should be okay for high Strehl
a=psf_gaussian(npixel=201,fwhm=fwhm_channel,centroid=[xc,yc],/normalize,/double)
intr_psf[*,*,il]=a
endfor
goto,breakoutpsf
endif 
if keyword_set(psf) then begin

if ~keyword_set(hardcodeparams) then begin
   psffile=dialog_pickfile(Title="Choose Empirical PSF")
endif else begin
   psffile='/home/sgoebel/thayne/pipeline/charisred/psfcube_med.fits'
endelse

intr_psf=readfits(psffile,ext=1)
goto,breakoutpsf
endif

psfdir='./psfmodel/'
print,'Using model PSF in "./psfmodel/" directory: psfcube_med.fits'
intr_psf=readfits(psfdir+'psfcube_med.fits',ext=1)

breakoutpsf:

lat=1.*sxpar(h,'lat')
lng=1.*sxpar(h,'lng')

;reading hour angle and parallactic angle at the beginning of exposures

readcol,'reduc.log',filenum,allxc,allyc,allrsat,allha,allpa

;radii for which to calculate the self subtraction

if ~keyword_set(lrmax) then lrmax=dimx/2
if ~keyword_set(lrmin) then lrmin=5

;***SNR Map***
;-self-consistently calculate the SNR map of the collapsed image
;- save SNR map, the aperture-summed image, and the uncertainty map
;snrmap=readfits('snrmap.fits')
;sigma for the image
;sigdata=abs(refcol)/snrmap
snratio_sub,refcol,fwhm=medfwhm,rmax=lrmax,/zero,/finite,snrmap=snrmap,noisemap=sigdata,imcol=refcol_con
;refcol_con=sumaper_im(refcol,medfwhm*0.5,lrmax,/nan)
;sigdata=abs(refcol_con)/snrmap

;writefits,'refcol.fits',refcol_con
;writefits,'sigdata.fits',sigdata
;writefits,'snr.fits',snrmap
;;stop


;*****defining region of interest for evaluating the disk model...
;-0. manually (based on geometry)
;-1. pre-defined region of interest
;-2. snrlimit
;-3. combination of 1 and 2

;array of radial distances
dist_circle,rarray,[dimx,dimy],[xc,yc]

if keyword_set(roimanual) or (~keyword_set(snrlim) and ~keyword_set(roirange)) then begin

read,'Set the Width of the Evaluation Region (set to 201 if no restriction)',gah1
read,'Set the Height of the Evaluation Region (set to 201 if no restriction)',gah2
read,'Set the Position Angle of the Evaluation Region (set to 0 if no restriction)',gah3
roirange=[gah1,gah2,gah3]
charis_generate_roi,output=roi_geom,roidim=[dimx,dimy],roirange=roirange
roitotal=roi_geom
goto,breakoutroi
endif

if keyword_set(snrlim) and ~keyword_set(roirange) then begin
;roi_snrmap=fltarr(dimx,dimy)
;roi_snrmap[*]=1
charis_generate_roi,output=roi_snrmap,roidim=[dimx,dimy],snrmap=snrmap
roitotal=roi_snrmap
goto,breakoutroi
endif 

if keyword_set(roirange) and ~keyword_set(snrlim) then begin
charis_generate_roi,output=roi_geom,roidim=[dimx,dimy],roirange=roirange
roitotal=roi_geom
goto,breakoutroi
endif

if keyword_set(roirange) and keyword_set(snrlim) then begin
charis_generate_roi,output=roi_snrmap,roidim=[dimx,dimy],snrmap=snrmap
charis_generate_roi,output=roi_geom,roidim=[dimx,dimy],roirange=roirange
;roitotal=intersect(roi_snrmap,roi_geom)
roitotal=long(roi_snrmap*roi_geom)
goto,breakoutroi
endif

breakoutroi:
roi_rad=where(rarray ge lrmin and rarray le lrmax)
;roi_snrmap=where(snrmap ge snrlim)
;roitotal=intersect(roitotal,roi_rad)
;roitotal_ind=(roitotal)
roitotal=where(roitotal eq 1)
roitotal=intersect(roi_rad,roitotal)
roiregion=fltarr(dimx,dimy)
roiregion[roitotal]=1
writefits,'roifit.fits',roiregion
;stop

;roitotal[roitotal_ind]=1
;writefits,'roitotal.fits',roitotal
;stop

;***PSF
;****modify for disk!
dpsf=dimx
;h=headfits(datadir+files[0])

if ~keyword_set(xext_fwhm) then xext_fwhm=120.
if ~keyword_set(yext_fwhm) then yext_fwhm=10.

;Here is where you can put in a disk model from GRaTeR ...
;find all model files
if keyword_set(hardcodeparams) then begin
   result=file_search('/home/sgoebel/thayne/pipeline/charisred/data/disk_models/*fits')
endif else begin
   result=dialog_pickfile(Title="Select Your Model Disks",/multiple)
endelse

if keyword_set(noreprocess) then begin
   ;Move previously-processed files to main directory

   spawn, 'mv /home/sgoebel/thayne/pipeline/copy2/reduc/proc/model_in/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_in/'
   spawn, 'mv /home/sgoebel/thayne/pipeline/copy2/reduc/proc/model_psfsub/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_psfsub/'

   spawn, 'mv /home/sgoebel/thayne/pipeline/copy3/reduc/proc/model_in/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_in/'
   spawn, 'mv /home/sgoebel/thayne/pipeline/copy3/reduc/proc/model_psfsub/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_psfsub/'

   spawn, 'mv /home/sgoebel/thayne/pipeline/copy4/reduc/proc/model_in/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_in/'
   spawn, 'mv /home/sgoebel/thayne/pipeline/copy4/reduc/proc/model_psfsub/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_psfsub/'

   spawn, 'mv /home/sgoebel/thayne/pipeline/copy5/reduc/proc/model_in/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_in/'
   spawn, 'mv /home/sgoebel/thayne/pipeline/copy5/reduc/proc/model_psfsub/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_psfsub/'

   spawn, 'mv /home/sgoebel/thayne/pipeline/copy6/reduc/proc/model_in/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_in/'
   spawn, 'mv /home/sgoebel/thayne/pipeline/copy6/reduc/proc/model_psfsub/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_psfsub/'

spawn, 'mv /home/sgoebel/thayne/pipeline/copy7/reduc/proc/model_in/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_in/'
   spawn, 'mv /home/sgoebel/thayne/pipeline/copy7/reduc/proc/model_psfsub/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_psfsub/'

spawn, 'mv /home/sgoebel/thayne/pipeline/copy8/reduc/proc/model_in/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_in/'
   spawn, 'mv /home/sgoebel/thayne/pipeline/copy8/reduc/proc/model_psfsub/*fits /home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_psfsub/'



   ;remove previously-processed files from list
   ;result=file_search('reduc/proc/model_in/*fits')
   done=file_search('/home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_in/*fits')
   mncopy = result
   for i=0, n_elements(mncopy)-1 do begin ;trim off directory info
      mncopy[i] = strmid(mncopy[i], max(strsplit(mncopy[i], '/')), strlen(mncopy[i]))
   endfor
   for i=0, n_elements(done)-1 do begin ;trim off directory and suffix info
      done[i] = repstr(done[i], '_inputcolsc', '')
      done[i] = strmid(done[i], max(strsplit(done[i], '/')), strlen(done[i]))
      
      loc = where(mncopy eq done[i]) ; has the file already been analyzed
      if loc ne -1 then begin
         print, "skipping", mncopy[loc]
         ;print
         remove, loc, result, mncopy
      endif
   endfor
endif

if ~keyword_set(n_threads) or ~keyword_set(n_threads)  then begin
   n_threads=1
   codeno=1
endif 
count=n_elements(result)
print, count, "new models found total."
modelname=result[count*(codeno-1)/n_threads : count*codeno/n_threads-1]

ntc=n_elements(modelname)
;print,modelname

;Now, loop on the model grid!

for itc=0,ntc-1 do begin ;loop sur les angles
   print,"Working on model", strcompress(itc+1), " of", strcompress(ntc)
;the input disk model
psf=readfits(modelname[itc],h_psf)

;now read in the important parameters
model_g=sxpar(h_psf,'G')
model_ksi=sxpar(h_psf,'KSI0')
model_alphain=sxpar(h_psf,'ALPHAIN')
model_alphaout=sxpar(h_psf,'ALPHAOUT')
model_beta=sxpar(h_psf,'BETA')
model_xdo=sxpar(h_psf,'XDO')
model_ydo=sxpar(h_psf,'YDO')
model_ecc=sxpar(h_psf,'E')
model_theta=sxpar(h_psf,'THETA0')
model_itilt=sxpar(h_psf,'ITILT')
model_r0=sxpar(h_psf,'r0')

;scaled to match image (a guess)
synthcube_scale_guess=1.*total(psf[roitotal]*refcol[roitotal],/nan)/total(psf[roitotal]*psf[roitotal],/nan)
print,'scale is ',synthcube_scale_guess
psf*=synthcube_scale_guess
psf+=1.d-10
psf_input=psf


;treat the PSF as a cube.
psf_chan=fltarr(dimx,dimy,nlambda)
psf_chanconv=fltarr(dimx,dimy,nlambda)
for il=0L,nlambda-1 do begin

;convolve the input disk model with the PSF
psf_chan[*,*,il]=convolve(psf,intr_psf[*,*,il])
psf_chanconv[*,*,il]=psf_chan[*,*,il]

endfor

writefits,'psf_chan.fits',psf_chanconv
writefits,'psf.fits',median(psf_chanconv,dimension=3,/even)

;x,y coordinates of pixels of the PSF centered on pixel 0.0

xpsf=lindgen(dpsf)#replicate(1l,dpsf)-dpsf/2
ypsf=replicate(1l,dpsf)#lindgen(dpsf)-dpsf/2

;indices of pixels of the psf centered on pixel 0.0 in the big picture
ipsf=xpsf+ypsf*dimx

;stretched to the psf radius rc=

imt=fltarr(dpsf,dpsf,nrc)
nsc=21 ;number of sub-images for each image slice

for n=0,nfiles-1 do begin

;Comment: here are useful all the images to seq
;working out the smearing even if one removes for
;the median(badim). Minor effect, negligible...


;*****Keywords for RA,DEC, individual exposure time (exp1time), and coadds (so exp1time*ncoadds = exptime)*****
    h=headfits(datadir+files[n],ext=0)
    radeg=sxpar(h,'ra')
    decdeg=sxpar(h,'dec')

    exptime=float(sxpar(h,'exp1time'))
    ncoadds=float(sxpar(h,'coadds'))
endfor

;*************************************************************
;***Change back to original version, which was the right way to calculate this!
;***smearing effect is negligible.  eliminated for a cleaner code.

smear_fac=fltarr(nrc)
;for irc=0,nrc-1 do smear_fac[irc]=myaper_tc(imt[*,*,irc],dpsf/2,dpsf/2,raper)
smear_fac[*]=1.
imt=0
;-**


;loop calculates the angles and the self-subtraction
print,'ntc is ',ntc

print,'NEXT!',itc
;add in fake disks here
    ;constructed images false companions
    print,'itc='+strtrim(itc+1,2)+'/'+strtrim(ntc,2)

    print,' adding model disks ...'
;  "adding companions.."
    for n=0,nfiles-1 do begin

        ;im=readfits(reducdir+files[n],h,/silent)
        h0=headfits(datadir+files[n],ext=0,/silent)
        im=readfits(datadir+files[n],h,/silent,ext=1)

;now get the north PA offset
        getrot,h,northpa,cdelt
        print,northpa,allpa[n]

        nwvlh=(size(im,/dim))[2]
;empty cube
        im[*]=0
        ;parallactic angle of sub companions

        exptime=float(sxpar(h0,'exp1time'))
        ncoadds=float(sxpar(h0,'coadds'))
        ha=allha[n]
        x=ha+(findgen(nsc)/(nsc-1.))*(exptime*ncoadds/3600.) ;SI HA DEBUT POSE
        pa0=parangle(x,decdeg,lat)
        dpa=pa0-pa0[0]
        ;dpa[*]=0
        pa=allpa[n]+dpa
        ;print,pa
        ;stop

;loop over wavelength channels
  for il=0L,nwvlh-1 do begin

        for irc=0,nrc-1 do begin
            ;determine angles des sous-compagnons

            ;asc=(replicate(tcomp[itc],nsc)-(pa-allpa[0]))*!dtor
            ;asc=(replicate(tcomp[itc]+(irc mod 10)*36.+(irc mod 2)*(180.-36.),nsc)+(pa-allpa[0]))*!dtor

            ;asc=(replicate(tcomp[itc]+(irc mod 10)*36.+(irc mod 2)*(180.-36.),nsc)-(pa-allpa[0]))*!dtor

            ;determine x,y positions of sub-companions

;            xsc=rc[irc]*cos(asc)+xc
;            ysc=rc[irc]*sin(asc)+yc

            ;xsc=rc[irc]*cos(-1*asc)+xc
            ;ysc=rc[irc]*sin(-1*asc)+yc
           
            for isc=0,nsc-1 do begin
                ;xe=floor(xsc[isc]) & ye=floor(ysc[isc])
                ;dxf=xsc[isc]-xe    & dyf=ysc[isc]-ye

; so basically don't shift the model disk at all.  very small effect.
                dxf=0 & dyf=0
                psfs=shift_sub(psf_chan[*,*,il],dxf,dyf)
                

;If you have an IFS data format written by normal human beings then the first line is right.
                im[*,*,il]+=(rotat(psfs,(-1*northpa))/smear_fac[irc]/nsc)
                ;im[*,*,il]+=(rotat(psfs,(pa[isc]-allpa[0]+northpa))/smear_fac[irc]/nsc)
;If instead you have GPI/NICI data this is flipped
                ;im[*,*,il]+=(rotat(psfs,-1*(pa[isc]-allpa[n]+northpa))/smear_fac[irc]/nsc)
                ;endif

            endfor
        endfor
  endfor
        ;register image
        writefits,datadir+filesfc[n],0,h0
        writefits,datadir+filesfc[n],im,h,/append
        im[*]=0
    endfor

    ;*** ADI treatment

;    imrsub,pfname,prefname=prefname,/nomask,/skiprefme,suffname=suffname,rmax=lrmax,/prad,/fc
    ;imrsub,pfname,fft=21,prefname=prefname,/nomask,/skiprefme,suffname=suffname,rmax=lrmax,/prad,/fc

;****ADI****
znfwhm=sxpar(header,'klip_nfw')
zdrsub=sxpar(header,'klip_drs')
zrmin=sxpar(header,'klip_rmi')
zrmax=sxpar(header,'klip_rma')
znpca=sxpar(header,'klip_pca')

print,'settings',znfwhm,zdrsub,zrmin,zrmax,znpca

;PCA/KLIP with Forward-Modeling

charis_subsklip,pfname,prefname=prefname,$
/meansub,$
nfwhm=znfwhm,drsub=zdrsub,npca=znpca,rmin=zrmin,rmax=zrmax,/usecoeff,/fc,/fwdmod

;Now read back in the files, rescale them to minimize chi-squared, add fits header information, and save them with unique input file name
    h0cube=headfits(reducdir+'final_fc.fits') 
    imcube=readfits(reducdir+'final_fc.fits',ext=1,h1cube)
    h0col=headfits(reducdir+'final_fc_collapsed.fits')
    imcol=readfits(reducdir+'final_fc_collapsed.fits',ext=1,h1col)

;model base name, extract string here and in line below to determine how to write out scaled input model and scaled PSF sub model
    modelbasename=(strsplit(modelname[itc],'/',/extract,count=modcount))[modcount-1]
   ;modelbasename0=(strsplit(modelbasename,'.',/extract,count=modcount))[0]
    modelbasename0=strmid(modelbasename, 0, max(strsplit(modelbasename, '.'))-1)
;and read back in the original collapsed PSF and PSF/vs. channel
    psf_input=readfits('psf.fits')
    psf_input_spec=readfits('psf_chan.fits')

   ;the collapsed cube 
;****
    ;define earlier, up top ;refcol_con=sumaper_im(refcol,raper,lrmax,/nan)
    imcol_con=sumaper_im(imcol,medfwhm*0.5,lrmax,/nan)
    synthcol_scale=total((refcol_con*imcol_con/sigdata^2)[roitotal],/nan)/total((imcol_con*imcol_con/sigdata^2)[roitotal],/nan)

;if you want to skip scaling in case it messes up
    if keyword_set(manscale) then synthcol_scale=1
    imcol_con*=synthcol_scale
    imcol*=synthcol_scale
    psf_input*=synthcol_scale


    diffcube=refcol_con-imcol_con
    diffcubeu=refcol-imcol

;now define the convolved real data and model flux for comparison
    ;chisq_col=total(abs(diffcube[roitotal]/sigdata[roitotal])^2.,/nan)/(n_elements(roitotal)-1)
    ;chisqmap=total(abs(diffcube[roitotal]/sigdata[roitotal])^2.,/nan)
    chisqmap=fltarr(dimx,dimy)
    chisqmap[roitotal]=abs(diffcube[roitotal]/sigdata[roitotal])^2.

    writefits,'imcol_con.fits',imcol_con
    writefits,'refcol_con.fits',refcol_con
    writefits,'diffcube.fits',diffcube
    writefits,'diffcubeu.fits',diffcubeu
    writefits,'sigdata.fits',sigdata
    writefits,'chisqmap.fits',chisqmap

;Now, bin the image to resolution, take the total chisquared value of the binned image, find number of bins, and 
; then compute the reduced chi-squared in the binned image.
     
    charis_calc_chisqbin,chisqmap,medfwhm,chisq_bin_out,nbins,chisqmap_bin=chisqmap_bin
    chisq_col=chisq_bin_out
    nres=nbins
    writefits,'chisqmapbin.fits',chisqmap_bin

    sxaddpar,h1col,'sfac',synthcol_scale
    sxaddpar,h1col,'csq_col (binned image)',chisq_col
    sxaddpar,h1col,'npts (binned image)',nres
    sxaddpar,h1col,'G',model_g
    sxaddpar,h1col,'ksi0',model_ksi
    sxaddpar,h1col,'alphain',model_alphain
    sxaddpar,h1col,'alphaout',model_alphaout
    sxaddpar,h1col,'beta',model_beta
    sxaddpar,h1col,'r0',model_r0
    sxaddpar,h1col,'xdo',model_xdo
    sxaddpar,h1col,'ydo',model_ydo
    sxaddpar,h1col,'e',model_ecc
    sxaddpar,h1col,'theta0',model_theta
    sxaddpar,h1col,'itilt',model_itilt

;PSF-subbed model
    writefits,reducdir_model2+modelbasename0+'_psfsubcol.fits',0,h0col
    writefits,reducdir_model2+modelbasename0+'_psfsubcol.fits',imcol,h1col,/append
;input model (with scaling applied)
    writefits,reducdir_model1+modelbasename0+'_inputcolsc.fits',0,h0col
    writefits,reducdir_model1+modelbasename0+'_inputcolsc.fits',psf_input,h1col,/append
;***  

if keyword_set(calcattenfact) then begin 
  ;this, it turns out, is very easy assuming that the forward-model accurately predicts signal loss ...
    ;the input data cube

    ;The file is outputed to the /model_psfsub directory
    
    input_datacube=psf_chanconv
    annealed_datacube=imcube
    
    
    attenfact_cube=annealed_datacube/input_datacube
    badattenfact=where(abs(attenfact_cube) gt 1.01 or finite(attenfact_cube) eq 0,nbadattenfact) ;this will always be outside the disk trace.
    if nbadattenfact ge 1 then attenfact_cube[badattenfact]=1
    writefits,reducdir_model1+modelbasename0+'_attenfactcube.fits',0,h0col
    writefits,reducdir_model1+modelbasename0+'_attenfactcube.fits',attenfact_cube,h1col,/append

endif
    

;    print,modelbasename,' ','r-chisq is ',total(abs(diffcube[roi_rad])^2.,/nan)/(n_elements(roi_rad)-1)
    print,modelbasename,' ','r-chisq is ',chisq_col
;total(abs(diffcube[roi_rad])^2.,/nan)/(n_elements(roi_rad)-1)

;just a diagnostic
    writefits,'diff.fits',0,hrefcol
    writefits,'diff.fits',diffcubeu,h1col,/append

;now that you have everything done, save the model parameters in a file with their chisq values
    openu,unit,stat_filename, /get_lun
    printf,unit,modelbasename,model_g,model_ksi,model_alphain,model_alphaout,model_beta,$
              model_xdo,model_ydo,model_ecc,model_r0,model_theta,chisq_col,format='(a91,9(f6.3,1x),f8.0,1x,f6.3)'
    free_lun, unit

    ;update text file with status
    spawn, 'pwd', dir
    dirsub = strmid(dir, strlen(dir)-1, strlen(dir))
    if dirsub eq 'd' then dirsub='1' ;if in charisred then print 1 to file
    if ~file_test('/home/sgoebel/thayne/pipeline/charisred/status.txt') then begin
       spawn, "echo >/home/sgoebel/thayne/pipeline/charisred/status.txt"
    endif
    spawn, "echo "+systime() + ' ' + dirsub  + ' '+ modelbasename+'>>/home/sgoebel/thayne/pipeline/charisred/status.txt'

endfor

free_lun, 31

end 
