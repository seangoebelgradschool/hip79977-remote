pro recalc_chisq;, 
  ;stuff that would normally be passed in
  ;reducname=reducname, lrmin=lrmin, lrmax=lrmax, roirange=roirange
reducname='settled2_npca=2,drsub=6,nfwhm=1,rmin=5,rmax=50,meansub,meanadd.fits'
lrmin=10
lrmax=50
roirange=[100,25,-22]

;read in actual data reduction
reducdir='./reduc/'
subdir='proc/'
reducdir1=reducdir+'reg/'
datadir=reducdir1
reducdir+=subdir

reducdirorg=reducdir

psf_input=readfits('psf.fits')
psf_input_spec=readfits('psf_chan.fits')
refcube=readfits(reducdirorg+reducname,ext=1,h1)
reducname_col=(strsplit(reducname,'.',/extract))[0]+'_collapsed.fits'
;imcol=readfits(reducdir+'final_fc_collapsed.fits',ext=1,h1col)
refcol=readfits(reducdirorg+reducname_col,ext=1,h1col)

;***stuff for CHARIS***
;Telescope Diameter for Subaru
Dtel=7.9d0  ;effective aperture of Subaru
;pixel scale 16.4 mas/pixel
pixscale=0.0164

;Now Get the Wavelength Vector
get_charis_wvlh,hrefcube,lambda
lambda*=1d-3
nlambda=n_elements(lambda)
;determine array of FWHM
fwhm=1.0*(1.d-6*lambda/Dtel)*(180.*3600./!dpi)/pixscale
medfwhm=median(fwhm,/even)

imx=sxpar(h1,'naxis1')
dimx=sxpar(h1,'naxis2')
dimy=dimx   ;assume square arrays
xc=dimx/2 & yc=dimy/2

charis_generate_roi,output=roi_geom,roidim=[dimx,dimy],roirange=roirange
roitotal=roi_geom
roi_rad=where(rarray ge lrmin and rarray le lrmax)
roitotal=where(roitotal eq 1)
roitotal=intersect(roi_rad,roitotal)
roiregion=fltarr(dimx,dimy)
roiregion[roitotal]=1

;READ IN MODEL DISK INFORMATION
;result=file_search('/home/sgoebel/thayne/pipeline/charisred/data/disk_models/*fits')
result=file_search('/home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_psfsub/', count=count)
modelname=result
ntc = count

itc=0 ; FOR LOOP HERE
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


snratio_sub,refcol,fwhm=medfwhm,rmax=lrmax,/zero,/finite,snrmap=snrmap,noisemap=sigdata,imcol=refcol_con

imcol_con=sumaper_im(imcol,medfwhm*0.5,lrmax,/nan)
synthcol_scale=total((refcol_con*imcol_con/sigdata^2)[roitotal],/nan)/total((imcol_con*imcol_con/sigdata^2)[roitotal],/nan)

imcol_con*=synthcol_scale
;imcol*=synthcol_scale
;psf_input*=synthcol_scale

diffcube=refcol_con-imcol_con
;diffcubeu=refcol-imcol

chisqmap=fltarr(dimx,dimy)
chisqmap[roitotal]=abs(diffcube[roitotal]/sigdata[roitotal])^2.

charis_calc_chisqbin,chisqmap,medfwhm,chisq_bin_out,nbins,chisqmap_bin=chisqmap_bin
chisq_col = chisq_bin_out
print,modelbasename,' ','r-chisq is ',chisq_col


end
