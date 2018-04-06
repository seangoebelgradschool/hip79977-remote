pro batchsnr

;searches for klip-generated files and then calculates the SNR for them.

spawn, 'ls /home/group18/sgoebel/thayne/pipeline/charisred/reduc/proc/grid*fits', mylist
;remove 'collapsed' files from list
mylist = mylist[where(strmatch(mylist, '*collapsed*', /fold_case) eq 0)]

for i=0, (size(mylist))[1]-1 do begin
   filename = repstr(mylist[i], 'reduc/proc/', 'reduc/proc/snr_')

   print
   if file_test(filename+'snrmap.fits') then begin
      print, "Skipping "+filename
      print
      continue                  ;don't process again
   endif

   print, 'Working on file', strcompress(i+1), ' of', $
          strcompress((size(mylist))[1])
   print

   snratio, mylist[i], outname=filename
endfor

end
