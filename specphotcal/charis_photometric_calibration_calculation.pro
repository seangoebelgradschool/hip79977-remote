function charis_photometric_calibration_calculation,pri_header,ext_header,lambda,spectype=spectype,star_mag=star_mag,units=units,$
diagnostic=diagnostic

;***Spectrophotometric calibration subroutine, following GPI DRP v1.4 methods
;-reads in header, wavelength array, spectral type, and star brightness (in magnitudes)
;-searches the Pickles library for a match to the spectral type string
;-takes Pickles spectrum at resolution of CHARIS, applies color correction to get to zero mag in H band 
;-returns a calibrated Pickles spectrum at CHARIS' resolution and the star's brightness in integrated flux

;**Limitations**
;- currently has trouble if you cannot find a match in the Pickles library.  **to fix!!!**

;****Filter Responses for color corrections
;standard filter response functions subdirectory
;we are going to use the MKO H band filter for color corrections.
filtresponsedir='~/Research/Planets/DI/filter_response/'
charis_std_filter='Hband'
readcol,filtresponsedir+charis_std_filter+'.dat',charis_std_filter_wvlh,charis_std_filter_response

good=where(charis_std_filter_response gt 10)

charis_std_filter_wvlh=charis_std_filter_wvlh[good]
charis_std_filter_response=charis_std_filter_response[good]

charis_std_filter_response=charis_std_filter_response(sort(charis_std_filter_wvlh))
charis_std_filter_wvlh=charis_std_filter_wvlh(sort(charis_std_filter_wvlh))


;****change the model directory to your path!!!****
modelspecdir='~/idl_tools/gpi_pipeline/pipeline/config/pickles/'
readcol,modelspecdir+'AA_README',pickles_fnames,pickles_sptypes,pickles_temps,skipline=113,numline=79,format=('A,A,F'),/silent

 for i=0L,n_elements(pickles_fnames)-1 do begin
  result=strmatch(spectype,pickles_sptypes[i],/FOLD_CASE)
   if result ne 0 then begin
    ref_model_spectrum=modelspecdir+pickles_fnames[i]+'.fits'
    print,'refmodelname!',ref_model_spectrum
    ;stop
    if file_test(ref_model_spectrum) eq 0 then return,error ('CANNOT FIND SPECTRUM FILE in Pickles Library')

        ;load the spectrum
        print,'refmodelspec',ref_model_spectrum
        pickles=mrdfits(ref_model_spectrum,1)
        model_wavelengths=pickles.wavelength
        model_flux=pickles.flux     ;units are in ergs/s/cm^2/A with V=0

        ;the first match you find, stay with it and break out of this loop
        break
    endif
 endfor

;Bin the spectrum from full resolution to matching CHARIS's resolution
 filter=sxpar(pri_header,'FILTNAME')
 dloglam=sxpar(pri_header,'DLOGLAM')
 specresolution=1./dloglam
  ;case filter of 
  ; 'broadband': specresolution=30.
  ; 'J': specresolution=100.
  ; 'H': specresolution=100.
  ; 'K': specresolution=100.
  ;endcase

  fwhmloc=VALUE_LOCATE(model_wavelengths/1e4,[(lambda[0]),(lambda[0]+dloglam)])
  fwhm=float(fwhmloc[1]-fwhmloc[0])
  gaus=PSF_GAUSSIAN(Npixel=3.*fwhm, FWHM=fwhm, NDIMEN =1, /NORMAL )
  charis_model_flux0=CONVOL( reform(model_flux), gaus , /EDGE_TRUNCATE )

  ;interpolation

  charis_model_flux=interpol(charis_model_flux0,model_wavelengths/1e4,lambda)   
  charis_model_flux_colorcorr=interpol(charis_model_flux0,model_wavelengths/1e4,charis_std_filter_wvlh*1d-3)   

;goto,skipdiag 
 if keyword_set(diagnostic) then begin
 setcolors,/system_variables,/silent
  ;plot,model_wavelengths/1e4,charis_model_flux0,color=!green,xrange=[1,2.5]
  ;oplot,model_wavelengths/1e4,reform(model_flux),thick=10
  plot,model_wavelengths/1e4,charis_model_flux0*model_wavelengths^2.,color=!green,xrange=[1,2.5]
  oplot,model_wavelengths/1e4,reform(model_flux)*model_wavelengths^2.,thick=10
;  oplot,model_wavelengths/1e4,charis_model_flux0,color=!green,thick=5
  oplot,lambda,charis_model_flux,color=!blue,thick=2.5 
 endif
;stop
skipdiag:

;for now, assume that no filter correction.  It's slight if anything and broadband mode covers JHK simultaneously.

;****Derive Color Correction for Model Spectrum since each model is defined as V=0.  
;;;;;use vega spectrum for this

 vega=mrdfits(modelspecdir+'alpha_lyr_stis_005.fits',1)
 charis_vega_flux0=convol(reform(vega.flux),gaus,/edge_truncate)


 ;interpolate model to filter wavelengths
 model_vega_flux=interpol(charis_vega_flux0,vega.wavelength/1e4,charis_std_filter_wvlh*1d-3)

if keyword_set(diagnostic) then begin
setcolors,/system_variables,/silent
plot,vega.wavelength/1e4,reform(vega.flux),psym=4,xrange=[1.0,2.5],color=!blue,thick=5
oplot,model_wavelengths/1e4,reform(model_flux),linestyle=1,color=!green,thick=5
;stop
  
setcolors,/system_variables,/silent
;plot,charis_std_filter_wvlh*1d-3,model_vega_flux,color=!blue,thick=5
;oplot,charis_std_filter_wvlh*1d-3,charis_model_flux_colorcorr,color=!green,thick=3 
;print,min(charis_model_flux_colorcorr),mean(charis_model_flux_colorcorr),max(charis_model_flux_colorcorr)
;print,max(charis_std_filter_wvlh*1d-3)
;stop
;help,charis_std_filter_wvlh,charis_model_flux_colorcorr,model_vega_flux 

window,1
 ;plot,charis_std_filter_wvlh*1d-3, model_vega_flux,xrange=[0,3],color=!blue,thick=5,yrange=[min(charis_model_flux_colorcorr),max(charis_model_flux_colorcorr)]
 ;oplot,charis_std_filter_wvlh*1d-3,charis_model_flux_colorcorr,color=!green,thick=3
;print,int_tabulated(charis_std_filter_wvlh*1d-3,charis_model_flux_colorcorr,/double,/sort)
;print,int_tabulated(charis_std_filter_wvlh*1d-3,model_vega_flux,/double,/sort)
;stop
endif

;***Deriving the color correction.
 ;star_color_correction=-2.5*alog10(int_tabulated(charis_std_filter_wvlh*1d-3,charis_model_flux_colorcorr,/double,/sort)/int_tabulated(charis_std_filter_wvlh*1d-3,model_vega_flux,/double,/sort))
 star_color_correction=-2.5*alog10(int_tabulated(charis_std_filter_wvlh*1d-3,charis_std_filter_response*charis_model_flux_colorcorr,/double,/sort)/int_tabulated(charis_std_filter_wvlh*1d-3,model_vega_flux*charis_std_filter_response,/double,/sort))

if keyword_set(diagnostic) then begin
 plot,charis_std_filter_wvlh*1d-3,charis_model_flux_colorcorr*charis_std_filter_wvlh^2.
 oplot,charis_std_filter_wvlh*1d-3,model_vega_flux*charis_std_filter_wvlh^2.,linestyle=1
;stop
endif

 print,'color correction factor is',star_color_correction 

;***Absolute flux calibration.  Going from V=0, H=-colorcorr, to V=whatever, H=the H band magnitude of the star
 charis_model_flux*=10.0^(-(star_mag-star_color_correction)/2.5)   ;units still in ergs/s/cm^2/A
 model_flux_colorcorr=model_flux*10.0^(-(star_mag-star_color_correction)/2.5)

if keyword_set(diagnostic) then begin
setcolors,/system_variables,/silent
 ;print,(charis_std_filter_wvlh*1d-3)[0]
 plot,charis_std_filter_wvlh*1d-3,charis_model_flux_colorcorr,color=!green
 oplot,lambda,charis_model_flux,color=!blue
 ;charis_model_flux_colorcorr=interpol(charis_model_flux0,model_wavelengths/1e4,charis_std_filter_wvlh*1d-3)   
 bah=-2.5*alog10(int_tabulated(charis_std_filter_wvlh*1d-3,charis_model_flux_colorcorr,/double,/sort)/int_tabulated(charis_std_filter_wvlh*1d-3,model_vega_flux,/double,/sort))
 ;print,'BAHHH',bah,star_mag
endif


;******Now convert units
   unitslist=['mJy','Jy','W/m^2/um','ergs/s/cm^2/A','ergs/s/cm^2/Hz']
   case units of 
    0: begin ;'mJy'
          unitname=unitslist[0]
           ; 1Jy = 10^-23 erg/s/Hz/cm2
                                ;c/l=f
        ;c/l^2 dl= df  -- so dl/df= l^2 / c
                          c=2.99792458d14                     ; um / s
                                conv_fact=((lambda[*]^2)/c) ; this is in um*s
                                conv_fact*=1e4 ; now in A*s
                                conv_fact*=(1.0/10.0^(-23.0))
                                conv_fact*=1d3  ;from Jy to mJy
                                 end
    1: begin ;'Jy'
          unitname=unitslist[0]
           ; 1Jy = 10^-23 erg/s/Hz/cm2
                                ;c/l=f
        ;c/l^2 dl= df  -- so dl/df= l^2 / c
                          c=2.99792458d14                     ; um / s
                                conv_fact=((lambda[*]^2)/c) ; this is in um*s
                                conv_fact*=1e4 ; now in A*s
                                conv_fact*=(1.0/10.0^(-23.0))
                                 end
    2: begin ;'W/m^2/um'
          unitname=unitslist[1]
                                conv_fact=fltarr(N_ELEMENTS(lambda))+1
                                conv_fact*=10.0^(4) ; from  ergs/s/cm^2/A to ergs/s/cm^2/um
                                conv_fact*=(100.0^2.0) ; from ergs/s/cm^2/um to ergs/s/m^2/um
                                conv_fact*=1e-7 ; from ergs/s/m^2/um to J/s/m2/um = W/m2/um
     
       end

    3: begin ;'ergs/s/cm^2/A'
          unitname=unitslist[2]
                                conv_fact=1.0
       end
    4: begin ;'ergs/s/cm^2/Hz'
          unitname=unitslist[3]
                                ;c/l=f
        ;c/l^2 dl= df  -- so dl/df= l^2 / c
                          c=2.99792458d14                     ; um / s
                                conv_fact=((lambda^2)/c) ; this is in um*s
                                conv_fact*=1e4 ; now in A*s

       end
   endcase
  print,'Units are in ',unitname

if keyword_set(diagnostic) then begin
setcolors,/system_variables,/silent
;window,2 
; plot,model_wavelengths/1e4,model_flux_colorcorr,xthick=3,ythick=3,/nodata,xtitle='wavelength',ytitle='Flux Density
  oplot,model_wavelengths/1e4,model_flux*10.0^(-(star_mag-0*star_color_correction)/2.5),color=!green
  oplot,model_wavelengths/1e4,model_flux_colorcorr,color=!blue
;stop
window,3
; plot,lambda,conv_fact*charis_model_flux,xthick=3,ythick=3,/nodata,xtitle='wavelength',ytitle='Flux Density
  oplot,lambda,conv_fact*charis_model_flux,color=!green
  ;oplot,lambda,model_flux_colorcorr,color=!blue

 ;stop
endif

charis_model_flux*=conv_fact
if keyword_set(diagnostic) then $
 plot,lambda,charis_model_flux,color=!blue

;stop
return, charis_model_flux
  
;  stop
stop
end
