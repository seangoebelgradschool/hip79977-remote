function get_charis_sat_contrast,modulation_size,wavelengths,mjd,errors=errors

;Given the modulation amplitude, returns the attenuation of the satellite spots/wavelength.  Formally only directly tested for broadband data but probably fine for highres.
;****Based on internal source testing 07/20/2017 and earlier in summer 2017

;Decision tree.

;--If the modulation size is 50nm, then you know it precisely
;--If the modulation size is 25nm, then you know it precisely
;--If the modulation size is something else, then extrapolate, knowing that the precision is good to ~a few percent or so.

;set to 'long' value.
 modulation_size=long(round(modulation_size)) 


atten_fac=fltarr(n_elements(wavelengths))
central_wvlh=1.550

modname=-99
if modulation_size eq 50 then modname=0
if modulation_size eq 25 then modname=1

if mjd lt 57995.0 then begin
case modulation_size of
    50: begin
      ;use measured attenuation of 5.991 +/-0.060 mags.
       atten_fac=10^(-0.4*5.991)
         conv_fac=(central_wvlh/wavelengths)^2.
       atten_fac*=conv_fac
       errors=0.06/1.0857*atten_fac
       end
    25: begin
      ;use measured attenuation of 7.504 +/- 0.060 mags
       atten_fac=10^(-0.4*7.5037)
         conv_fac=(central_wvlh/wavelengths)^2.
       atten_fac*=conv_fac
       errors=0.06/1.0857*atten_fac
       end
   else: begin
      ;rough use scaling.
      atten_fac=(modulation_size/25.)^2.*10^(-0.4*7.5037)
         conv_fac=(central_wvlh/wavelengths)^2.
      atten_fac*=conv_fac
      errors=0.06/1.0857*atten_fac
      end
endcase

endif else begin


case modulation_size of
    50: begin
      ;use measured attenuation of 4.92 +/-0.050 mags.
       atten_fac=10^(-0.4*4.92)
         conv_fac=(central_wvlh/wavelengths)^2.
       atten_fac*=conv_fac
       errors=0.05/1.0857*atten_fac
       end
    25: begin
      ;use measured attenuation of 6.412 +/- 0.050 mags
       atten_fac=10^(-0.4*6.412)
         conv_fac=(central_wvlh/wavelengths)^2.
       atten_fac*=conv_fac
       errors=0.05/1.0857*atten_fac
       end
   else: begin
      ;rough use scaling.
      atten_fac=(modulation_size/25.)^2.*10^(-0.4*6.412)
         conv_fac=(central_wvlh/wavelengths)^2.
      atten_fac*=conv_fac
      errors=0.05/1.0857*atten_fac
      end
endcase


endelse

return,atten_fac
end
