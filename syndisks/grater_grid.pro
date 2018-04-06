pro grater_grid
;generates a grid of grater files
;Sean changed the filenames. Also e wasn't being passed in.

;set these
itilt_in=[84.5]
beta_in=2.2

r0_range= [73, 78, 83, 88]
gscat_range=[0.5, 0.55, 0.6, 0.65, 0.7]
ksi_range=[0.5, 1, 1.5]
alphain_range= [5, 7.5]
alphaout_range=[-2, -2.5, -3]
pa_range=[112.5]
ydo_range=[0]
xdo_range=[0]
e_range=[0]
theta_range=[112.5]
openw,18,'modellist'
openw,19,'model_params'
printf,18,'MODEL_LIST'
printf,19,'MODEL_NAME','G','KSI0','ALP_I','ALP_O','BETA','XDO','YDO','INC','THETA0',$
   format='(a25,a6,1x,8(a7,1x))'

print,n_elements(gscat_range)*n_elements(ksi_range)*n_elements(alphain_range)*n_elements(alphaout_range)*n_elements(pa_range)*n_elements(r0_range)*n_elements(xdo_range)*n_elements(ydo_range)*n_elements(itilt_in)
for q=0L,n_elements(r0_range) - 1 do begin
for i=0L,n_elements(gscat_range)-1 do begin
 for ii=0L,n_elements(ksi_range) -1 do begin
  for iii=0L,n_elements(alphain_range)-1 do begin
   for iv=0L,n_elements(alphaout_range)-1 do begin
    for j=0L,n_elements(pa_range)-1 do begin
     for ji=0L,n_elements(xdo_range)-1 do begin  
      for jii=0L,n_elements(ydo_range)-1 do begin
       for jiii=0L,n_elements(e_range)-1 do begin
        for jiv=0L,n_elements(itilt_in)-1 do begin
        ;for jiv=0L,n_elements(theta_range)-1 do begin
;outfilename='data/disk_models/egrater'+strtrim(i,1)+strtrim(ii,1)+strtrim(iii,1)+strtrim(iv,1)+strtrim(j,1)+strtrim(ji,1)+strtrim(jii,1)+strtrim(jiii,1)+strtrim(jiv,1)+strtrim(q,1)+'_model.fits'

outfilename='data/disk_models/egrater_'+$
            strtrim(gscat_range[i],1)+','+$
            strtrim(ksi_range[ii],1)+','+$
            strtrim(alphain_range[iii],1)+','+$
            strtrim(alphaout_range[iv],1)+','+$
            strtrim(pa_range[j],1)+','+$
            strtrim(xdo_range[ji],1)+','+$
            strtrim(ydo_range[jii],1)+','+$
            strtrim(theta_range[0],1)+','+$
            strtrim(itilt_in[jiv],1)+','+$
            strtrim(r0_range[q],1)+','+$
            strtrim(e_range[jiii],1)+'.fits'

print,outfilename

call_grater,g=gscat_range[i],ksi0=ksi_range[ii],alphain=alphain_range[iii],$
  alphaout=alphaout_range[iv],pa=pa_range[j],xdo=xdo_range[ji],ydo=ydo_range[jii],$
   e=e_range[jiii],theta0=theta_range[0],outfile=outfilename,$
   itilt=itilt_in[jiv],r0=r0_range[q],$
dstar=123.,pfov=0.0164,nx=201,ny=201
;and now the hardwired stuff
  
skipoverme:
printf,18,outfilename,format='(a20)'
printf,19,outfilename,gscat_range[i],ksi_range[ii],alphain_range[iii],alphaout_range[iv],pa_range[j],$
              xdo_range[ji],ydo_range[jii],itilt_in[jiv],theta_range[0],format='(a25,9(f7.3,1x))'
         endfor
        endfor
       endfor
      endfor
     endfor
    endfor
   endfor
  endfor
 endfor
 endfor

close, 18, 19

end
