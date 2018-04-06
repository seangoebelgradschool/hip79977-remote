pro klipgrid

; calls charis_subsklip with a variety of parameters.

; Why is idl such a flaming pile of shit?
CATCH, Error_status 
IF Error_status NE 0 THEN BEGIN ;In case the code tries to crash for no reason in particular
   PRINT, 'Error index: ', Error_status
   PRINT, 'Error message: ', !ERROR_STATE.MSG

   CATCH, /CANCEL ;restart code instead of exiting
ENDIF

count = 0
for npca_val = 2,2 do begin ; 3 is default
   for nfwhm_val = 1.0, 1. do begin ; 1 is defult
      for drsub_val = 5, 5 do begin ; 5 is default
         count += 1
         print, "Working on file", strcompress(count)

         filename = 'grid_npca='+ strcompress(npca_val, /remove)+ $
                    ',nfwhm='+strmid(repstr(strcompress(nfwhm_val),'.',''),1,3)+$
                    ',drsub='+ strcompress(drsub_val, /remove)+ $
                    ',zeromeansub.fits'
         
         filename='deleteme2.fits'

         if file_test('reduc/proc/'+filename) then begin
            print, "Skipping "+filename
            print
            continue ;don't process again
         endif

         drsub_value = drsub_val ; IDL PIECE OF CRAP BUG
         npca_value = npca_val ; IDL PIECE OF CRAP BUG

         spawn, '/bin/rm reduc/proc/*adisub*fits'

         charis_subsklip, 'HIP79977h.info', npca=npca_value, $
                          nfwhm=nfwhm_val, rmin=5, rmax=50, $
                          drsub=drsub_value, outfile=filename, $
                          /meansub, /meanadd

         ;print, npca_val
         ;print, nfwhm_val
         ;print, drsub_valu
         print, "Finished ", filename, "."
         print
      endfor
   endfor
endfor

end
