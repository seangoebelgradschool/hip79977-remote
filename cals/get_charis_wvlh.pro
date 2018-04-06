pro get_charis_wvlh,header,wavelengths,wavecalfile=wavecalfile
;v2.

;Takes the fits header from Charis and returns the wavelengths
;Note: must change wavecal path in line 9

filtname=strtrim(sxpar(header,'CAL_BAND',count=filtcount))  ;for right now, using the calibration band to define the mode since FILTNAME gets unpopulated.

;SEAN'S EDIT: PATH NEEDS TO BE UPDATED
wavecalpath='/home/sgoebel/thayne/pipeline/charisred/cals/'
;~/idl_tools/ADI_dl/charisred/cals/'

if filtcount eq 0 then read,'Select the filter to use (lowres,J,H,K)',filtname

if ~keyword_set(wavecalfile) then begin

;broadband
if filtname eq 'lowres' then begin
wavecalfile='lowres_wvlh.txt'

endif

if filtname eq 'J' then begin
;jband
wavecal='highres_jwvlh.txt'
endif

if filtname eq 'H' then begin
;h band
wavecalfile='highres_hwvlh.txt'
endif

if filtname eq 'K' then begin
;kband
wavecalfile='highres_kwvlh.txt'
endif


;wavecalfile='polychromekeyR30.fits'
;you can override the wavecalfile
endif

;override
;wavelengths=readfits(wavecalpath+wavecalfile,h1)
readcol,wavecalpath+wavecalfile,wavelengths,/silent



end
