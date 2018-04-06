pro scrapefitsinfo
; Scrapes disk parameter information from fits headers and prints it to a txt document.

  dirs = ['/home/sgoebel/thayne/pipeline/charisred/reduc/proc/model_psfsub/',$
          '/home/sgoebel/thayne/pipeline/copy2/reduc/proc/model_psfsub/',$
          '/home/sgoebel/thayne/pipeline/copy3/reduc/proc/model_psfsub/']

  openw,unit, '/home/sgoebel/thayne/pipeline/charisred/grid_search_stats.txt', /get_lun

  printf,unit,'MODEL_NAME','G','KSI0','ALP_I','ALP_O','BETA','XDO','YDO','E','R0','THETA0','CHISQ/dof', format='(a91,a5,1x,9(a6,1x),a9)'

  for i=0, n_elements(dirs)-1 do begin
     print, "Working on directory ", strcompress(i+1, /remove), " of ", strcompress(n_elements(dirs),/remove), "."
     
     result=file_search(dirs[i]+'*fits', count=count)
     ;spawn, 'ls -t '+dirs[i]+'*.fits', result
     ;result = result[0:3]
     ;count = 4

     for j=0, count-1 do begin
        hdr = headfits(result[j], ext=1)

        filename=strmid(result[j], max(strsplit(result[j], '/')), strlen(result[j]))
        filename = repstr(filename, '_psfsubcol', '')
        
        chi_sq = sxpar(hdr, 'CSQ_COL')
        g = sxpar(hdr, 'G')
        ksi_0 = sxpar(hdr, 'KSI0')
        alpha_in = sxpar(hdr, 'ALPHAIN')
        alpha_out = sxpar(hdr, 'ALPHAOUT')
        beta = sxpar(hdr, 'BETA')
        r0 = sxpar(hdr, 'R0')
        xdo = sxpar(hdr, 'XDO')
        ydo = sxpar(hdr, 'YDO')
        e = sxpar(hdr, 'E')
        theta_0 = sxpar(hdr, 'THETA0')
        itilt = sxpar(hdr, 'ITILT')

        printf,unit,filename,g,ksi_0,alpha_in,alpha_out,beta,$
              xdo,ydo,e,r0,theta_0,chi_sq,format='(a91,9(f6.3,1x),f8.0,1x,f6.3)'
     endfor
  endfor

  free_lun, unit

end
