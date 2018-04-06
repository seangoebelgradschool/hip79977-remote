pro charis_calc_analysis,compare=compare,avg=avg,outfile=outfile,datatype=datatype,haloprof=haloprof

;simple utility to evaluate contrasts and compare contrasts
;Use wisely and read/understand the source code.  do NOT assume that it is plotting what you imagine it to plot
;Again, Use at your own risk!!!!

;select CHARIS output text files with broadband contrast, J/H/K contrast 
; compute an averaged contrast 
; compare contrasts

;updates
;10/12/2017 - added switches 'datatype' that can be used to write 'lyot' or 'vortex' or 'spc' on the plots to identify the type of data plotted.
;10/12/2017 - added switch 'haloprof' to tell code that you are plotting the normalized halo profile, not really the 5-sigma contrast curve.  
;             so have to divide by 5 since input code (...raw_contrast) thinks you are wanting the halo intensity profile

;writes results to working directory

outputcontrasts=dialog_pickfile(Title="Select Output Contrast Curves to Analyze",/multiple_files)
noutputs=n_elements(outputcontrasts)

;test file
readcol,outputcontrasts[0],rtest,btest,jtest,htest,ktest,format='(f,f,f,f,f)'
nr=n_elements(rtest)
print,nr

;loop the output files,compute a median contrast curve in broadband, J, H, and K
bcontrast=fltarr(nr,noutputs)
jcontrast=fltarr(nr,noutputs)
hcontrast=fltarr(nr,noutputs)
kcontrast=fltarr(nr,noutputs)


for i=0L,noutputs-1 do begin

readcol,outputcontrasts[i],rarcsec,bcon,jcon,hcon,kcon,format='(f,f,f,f,f)'
if ~keyword_set(haloprof) then begin
bcontrast[*,i]=bcon
jcontrast[*,i]=jcon
hcontrast[*,i]=hcon
kcontrast[*,i]=kcon
endif else begin

bcontrast[*,i]=bcon/5.
jcontrast[*,i]=jcon/5.
hcontrast[*,i]=hcon/5.
kcontrast[*,i]=kcon/5.

endelse

endfor

if keyword_set(avg) then begin
;now compute the averaged contrast
bcontrastf=median(bcontrast,dimension=2,/even)
jcontrastf=median(jcontrast,dimension=2,/even)
hcontrastf=median(hcontrast,dimension=2,/even)
kcontrastf=median(kcontrast,dimension=2,/even)

set_plot,'x'
setcolors,/system_variables,/silent
plot,rtest,bcontrastf,/nodata,xrange=[0,1],yrange=[5d-6,5d-3],/ylog,$ 
xtitle='Separation (arc-sec)',ytitle=textoidl('Avg. 5-\sigma Contrast ( \Delta F)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
loadct,13,/silent
oplot,rtest,jcontrastf,color=(4/21.)*256,thick=3
oplot,rtest,hcontrastf,color=(11/21.)*256,thick=3
oplot,rtest,kcontrastf,color=(18/21.)*256,thick=3

set_plot,'ps'
if keyword_set(outfile) then begin
output=outfile
endif else begin
output='avgcontrast'
endelse
device,filename=output+'.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent

plot,rtest,bcontrastf,/nodata,xrange=[0,1],yrange=[5d-6,5d-3],/ylog,$
xtitle='Separation (arc-sec)',ytitle=textoidl('Avg. 5-\sigma Contrast ( \Delta F)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
loadct,13,/silent
oplot,rtest,jcontrastf,color=(4/21.)*256,thick=3
oplot,rtest,hcontrastf,color=(11/21.)*256,thick=3
oplot,rtest,kcontrastf,color=(18/21.)*256,thick=3

setcolors,/system_variables,/silent
oplot,rtest,bcontrastf,color=!cyan,thick=10

setcolors,/system_variables,/silent
al_legend,color=[(5./22.)*256,(12/22.)*256,(19./22.)*256,!cyan],['J','H','Ks','Broadband'],/bottom,$
box=0,linestyle=0,thick=[3,3,3,10],charsize=1.25,charthick=3
device,/close

writecol,outfile+'.txt',rtest,bcontrastf,jcontrastf,hcontrastf,kcontrastf,fmt='(f,f,f,f,f)'
endif

nsmooth=3
if keyword_set(compare) then begin

if noutputs gt 3 or noutputs lt 3 then begin
print,'can only compare three contrast curves right now'
goto,breakmeout
endif 

set_plot,'x'
setcolors,/system_variables,/silent
plot,rtest,bcontrast[*,0],/nodata,xrange=[0,1],yrange=[5d-6,5d-3],/ylog,$ 
xtitle='Separation (arc-sec)',ytitle=textoidl('Avg. 5-\sigma Contrast ( \Delta F)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
loadct,13,/silent
oplot,rtest,jcontrast[*,0],color=(4/21.)*256,thick=3
oplot,rtest,hcontrast[*,0],color=(11/21.)*256,thick=3
oplot,rtest,kcontrast[*,0],color=(18/21.)*256,thick=3

setcolors,/system_variables,/silent
oplot,rtest,bcontrast[*,0],color=!cyan,thick=3
oplot,rtest,bcontrast[*,1],color=!cyan,thick=3,linestyle=1
oplot,rtest,bcontrast[*,2],color=!cyan,thick=3,linestyle=2

set_plot,'ps'
if keyword_set(outfile) then begin
output=outfile
endif else begin
output='comp'
endelse
device,filename=output+'bcontrast'+'.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent

plot,rtest,bcontrast[*,0],/nodata,xrange=[0,1],yrange=[1d-5,5d-3],/ylog,$
xtitle='Separation (arc-sec)',ytitle=textoidl('5-\sigma Contrast ( \Delta F)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
loadct,13,/silent
setcolors,/system_variables,/silent
oplot,rtest,bcontrast[*,0],color=!black,thick=3,linestyle=2
oplot,rtest,bcontrast[*,1],color=!black,thick=3,linestyle=1
oplot,rtest,bcontrast[*,2],color=!black,thick=3,linestyle=0
device,/close


device,filename=output+'jcontrast'+'.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent

plot,rtest,bcontrast[*,0],/nodata,xrange=[0,1],yrange=[1d-5,1d-3],/ylog,$
xtitle='Separation (arc-sec)',ytitle=textoidl('5-\sigma Contrast ( \Delta F) (J)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
loadct,13,/silent
setcolors,/system_variables,/silent
oplot,rtest,jcontrast[*,0],color=!blue,thick=3,linestyle=2
oplot,rtest,jcontrast[*,1],color=!blue,thick=3,linestyle=1
oplot,rtest,jcontrast[*,2],color=!blue,thick=3,linestyle=0
device,/close

device,filename=output+'hcontrast'+'.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent

plot,rtest,bcontrast[*,0],/nodata,xrange=[0,1],yrange=[1d-5,1d-3],/ylog,$
xtitle='Separation (arc-sec)',ytitle=textoidl('5-\sigma Contrast ( \Delta F) (H)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
loadct,13,/silent
setcolors,/system_variables,/silent
oplot,rtest,hcontrast[*,0],color=!green,thick=3,linestyle=2
oplot,rtest,hcontrast[*,1],color=!green,thick=3,linestyle=1
oplot,rtest,hcontrast[*,2],color=!green,thick=3,linestyle=0
device,/close


device,filename=output+'kcontrast'+'.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent

plot,rtest,bcontrast[*,0],/nodata,xrange=[0,1],yrange=[1d-5,1d-3],/ylog,$
xtitle='Separation (arc-sec)',ytitle=textoidl('5-\sigma Contrast ( \Delta F) (K)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
loadct,13,/silent
setcolors,/system_variables,/silent
oplot,rtest,kcontrast[*,0],color=!red,thick=3,linestyle=2
oplot,rtest,kcontrast[*,1],color=!red,thick=3,linestyle=1
oplot,rtest,kcontrast[*,2],color=!red,thick=3,linestyle=0
device,/close


device,filename=output+'allcontrast'+'.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent
if ~keyword_set(haloprof) then begin
plot,rtest,bcontrast[*,0],/nodata,xrange=[0.1,1],yrange=[5d-6,1.5d-3],/ylog,$
xtitle='Separation (arc-sec)',ytitle=textoidl('5-\sigma Contrast'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
endif else begin
plot,rtest,bcontrast[*,0],/nodata,xrange=[0.1,1],yrange=[1d-5,5d-3],/ylog,$
xtitle='Separation (arc-sec)',ytitle='Halo Intensity Profile',xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
endelse
loadct,13,/silent
setcolors,/system_variables,/silent
oplot,rtest,smooth(bcontrast[*,2],nsmooth),color=!black,thick=5,linestyle=0
oplot,rtest,smooth(bcontrast[*,0],nsmooth),color=!black,thick=5,linestyle=2
oplot,rtest,smooth(bcontrast[*,1],nsmooth),color=!black,thick=5,linestyle=1
;goto,skipover
;oplot,rtest,smooth(jcontrast[*,0],nsmooth),color=!blue,thick=3,linestyle=2
;oplot,rtest,smooth(jcontrast[*,1],nsmooth),color=!blue,thick=3,linestyle=1
;oplot,rtest,smooth(jcontrast[*,2],nsmooth),color=!blue,thick=3,linestyle=0
oplot,rtest,smooth(hcontrast[*,0],nsmooth),color=!green,thick=3,linestyle=2
oplot,rtest,smooth(hcontrast[*,1],nsmooth),color=!green,thick=3,linestyle=1
oplot,rtest,smooth(hcontrast[*,2],nsmooth),color=!green,thick=3,linestyle=0
;oplot,rtest,smooth(kcontrast[*,0],nsmooth),color=!red,thick=3,linestyle=2
;oplot,rtest,smooth(kcontrast[*,1],nsmooth),color=!red,thick=3,linestyle=1
;oplot,rtest,smooth(kcontrast[*,2],nsmooth),color=!red,thick=3,linestyle=0
skipover:
setcolors,/system_variables
al_legend,['SR = 0.80','SR = 0.75','SR = 0.70'],linestyle=[0,1,2],thick=5,box=0,charsize=1.2,charthick=3,/right
al_legend,color=[!black,!green],['Broadband','H band Slice'],box=0,/bottom,linestyle=0,thick=3,$
charsize=1.2,charthick=3
if keyword_set(datatype) then begin
al_legend,datatype,charthick=3,charsize=1.2,box=0,/top
endif
device,/close

;ratios of contrast


device,filename='ratio'+'allcontrast'+'.eps',/encapsulated,/color,bits=8
setcolors,/system_variables,/silent

plot,rtest,bcontrast[*,0],/nodata,xrange=[0,1],yrange=[0,8],$
xtitle='Separation (arc-sec)',ytitle=textoidl('5-\sigma Contrast ( \Delta F) (K)'),xthick=5,ythick=5,xstyle=1,ystyle=1,$
charsize=1.5,charthick=3
loadct,13,/silent
setcolors,/system_variables,/silent
oplot,rtest,bcontrast[*,1]/bcontrast[*,2],color=!black,thick=3
oplot,rtest,bcontrast[*,0]/bcontrast[*,2],color=!black,thick=3,linestyle=2
goto,skipover2
oplot,rtest,jcontrast[*,1]/jcontrast[*,2],color=!blue,thick=3
oplot,rtest,jcontrast[*,0]/jcontrast[*,2],color=!blue,thick=3,linestyle=2
oplot,rtest,hcontrast[*,1]/hcontrast[*,2],color=!green,thick=3
oplot,rtest,hcontrast[*,0]/hcontrast[*,2],color=!green,thick=3,linestyle=2
oplot,rtest,kcontrast[*,1]/kcontrast[*,2],color=!red,thick=3
oplot,rtest,kcontrast[*,0]/kcontrast[*,2],color=!red,thick=3,linestyle=2
skipover2:
oplot,rtest,bcontrast[*,1]*0+1,linestyle=2,color=!gray
device,/close




writecol,output+'compare'+'.txt',rtest,bcontrast[*,1]/bcontrast[*,2],bcontrast[*,0]/bcontrast[*,2],$
jcontrast[*,1]/jcontrast[*,2],jcontrast[*,0]/jcontrast[*,2],$
hcontrast[*,1]/hcontrast[*,2],hcontrast[*,0]/hcontrast[*,2],$
kcontrast[*,1]/kcontrast[*,2],kcontrast[*,0]/kcontrast[*,2],$
fmt='(9(f8.3,1x))'
;,f,f,f,f,f,f)'



endif


breakmeout:

end
