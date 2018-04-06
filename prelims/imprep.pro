pro imprep,zip=zip,gz=gz,filelist=filelist,reducdir=reducdir,sigclip=sigclip,boxsize=boxsize,$
robust=robust,$
pfname=pfname,$
setsky=setsky,$
sort=sort,$
rename=rename,$
preproc=preproc,$
prefname=prefname,$
noextniri=noextniri,$
channel=channel,$
fixheader=fixheader,$
cal=cal,$
obslog=obslog,$
linearize=linearize,$
charismanual=charismanual,$
noexten=noexten,$
collapse=collapse,dither=dither,instrument=instrument
;Establishes uniform fits header information for all science frames you want to reduce
;individual integration time, coadds (so int.*coadds = total int), altitude, azimuth, LST, PA, beam (if dither), latitude, longitude, HA, WCS

;07/24/2017 -- Modified HiCIAO section to offer a robust calculation of parallactic angle if keyword is missing
;01/15/2016 -- cleaned up code.  
;11/26/2015 -- modified HiCIAO section to redefine the HA to be the average HA (thus the PA is the avg PA of the sequence)
;03/?/2014 -- modified to include GPI
;01/02/2013 -- modified to include Subaru/HiCIAO
;12/14/2012 -- modified to include Gemini/NICI data (desktop version modified for HiCIAO (get that too!)
;05/03/2012 -- 'setsky' switch to turn all frames into sky frames.  ***WORKS ONLY FOR KECK!***
;7/21/2011 -- Only does beam for VLT/NaCo and MMT/Clio data right now.  Need to change this

;7/27/2011 -- Right now, allows dither switch for ...
;nothing

if ~keyword_set(channel) then channel= 1
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(boxsize) then boxsize=150
datadir='./data/raw/'

reducdir1='./data/'
if ~keyword_set(reducdir) then begin
subdir='prep/'
reducdir=reducdir1+subdir
endif else begin
reducdir=reducdir1
endelse

;If the data are already reduced (e.g. GPI pipeline products), then just put in reg file
;for now, this only works with GPI and only works with data you rename!!!!
if keyword_set(preproc) then begin
;reducdir='./reduc/expand/'
reducdir='./reduc/reg/'
endif

if keyword_set(cal) then begin
reducdir='./data/cal/'
endif

file_mkdir,reducdir

print,'reducdir is',reducdir

if keyword_set(filelist) then begin
;List of Science Frames from text file
readcol,filelist,files,format='a'
nfiles=n_elements(files)

endif else begin

;assuming all files you are interested in, have a .fit or .fits string somewhere in the name!
if ~keyword_set(cal) then begin
files=file_search(datadir+'*fit*')
endif else begin
files=file_search(datadir+'/../cal/flat/'+'*fit*') 
files2= file_search(datadir+'/../cal/dark/'+'*fit*')
files=[files,files2]
print,files
endelse

if keyword_set(sort) then begin
times=dblarr(n_elements(files))
for i=0L,n_elements(files)-1 do begin
h1=headfits(files[i])

;****THIS IS A HACK RIGHT NOW!!!!!

;times[i]=double(sxpar(h1,'mjd'))
times[i]=ten(sxpar(h1,'lst'))
;times[i]=ten(sxpar(h1,'p_hststr'))

if times[i] lt 20 then times[i]+=24
print,'times',times[i],' ',files[i]
endfor

files=files(sort(times))
times=times(sort(times))

endif



nfiles=n_elements(files)

endelse
print,'nfiles is',nfiles


;**Keck/NIRC2**
;*DONE*
;**except that beam not set**

if instrument eq 'nirc2' then begin

 for i=0L,nfiles-1 do begin
 
  bb=readfits(files[i],h1,/silent)
  
  exp1time=sxpar(h1,'itime')
  sxaddpar,h1,'exp1time',exp1time,after='itime',format="f9.5"
  sxdelpar,h1,'itime'

  altitude=sxpar(h1,'el')
  sxaddpar,h1,'altitude',altitude,after='el',format="f9.5"
  sxdelpar,h1,'el'

  azimuth=sxpar(h1,'az') 
  sxaddpar,h1,'azimuth',azimuth,after='az',format="f9.5"
  sxdelpar,h1,'az'

  pa=sxpar(h1,'parang')
  sxaddpar,h1,'PA',pa,after='parang',format="f12.6"
  sxdelpar,h1,'parang'

  ;latitude and longitude
   
  lat=19.825d0
  lng=-155.4802d0

;do these steps because your code otherwise has serious trouble with 
; producing images that display properly in DS9 or fitsview
;***
  bb=double(bb)
  sxdelpar,h1,'bscale'
  sxdelpar,h1,'bzero'
  sxdelpar,h1,'bitpix'
;***

  sxaddpar,h1,'lat',lat,after='telescop',format="f10.6"
  sxaddpar,h1,'lng',lng,after='lat',format="f11.6"

;Are these predetermined to be sky frames?  If so, let the fits header know!
if keyword_set(setsky) then begin
sxdelpar,h1,'object'
sxaddpar,h1,'object','sky',after='observer',format="a3"
endif

; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]
print,fileoutname
a=bb

if keyword_set(linearize) then begin
linearize_nirc2,a,h1
endif

if keyword_set(gz) then begin 
;If you're grabbing data from the archive, it's already compressed and has a '.gz' suffix. 
;In such a case, do not add '.gz' to the filename

is_archive=sxpar(h1,'KOAID',count=arch_count)

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

if (arch_count) eq 0 then begin
fileoutname=fileoutname+'.gz'
endif

   
  writefits,reducdir+fileoutname,a,h1,/compress
endif else begin

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

  writefits,reducdir+fileoutname,a,h1
endelse

  ;writefits,reducdir+files[i],a,h1
 endfor

goto,breakout
endif 


;**Gemini/NIRI**

;*need to do PA*
;**beam not set**

if instrument eq 'niri' then begin

 for i=0L,nfiles-1 do begin

;first, deal with NIRI's multifits extension nonsense
  if ~keyword_set(noextniri) then begin
  bb=readfits(files[i],h2,ext=1)
  h1=headfits(files[i],ext=0)
  endif else begin
  bb=readfits(files[i],h2)
  h1=headfits(files[i])
  endelse
  naxis1=sxpar(h2,'naxis1')
  naxis2=sxpar(h2,'naxis2')
  naxis=sxpar(h2,'naxis')
  extend=sxpar(h1,'extend')
  
  sxdelpar,h1,'extend'
  sxaddpar,h1,'extend','F'

  sxaddpar,h1,'naxis',naxis
  sxaddpar,h1,'naxis1',naxis1
  sxaddpar,h1,'naxis2',naxis2
  
  exp1time=sxpar(h1,'exptime')
  sxaddpar,h1,'exp1time',exp1time,after='exptime',format="f9.5"
  
  altitude=sxpar(h1,'elevatio')
  sxaddpar,h1,'altitude',altitude,after='elevatio',format="f9.5"
  sxdelpar,h1,'elevatio'

  lst=sxpar(h1,'st')
  sxaddpar,h1,'lst',lst,after='st',format="f9.5"
  sxdelpar,h1,'st'

  ;latitude and longitude
   
  lat=19.825d0
  lng=-155.4802d0

  sxaddpar,h1,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h1,'lng',lng,after='lat',format="f12.7"

  if keyword_set(cal) then goto, skippa
  sxdelpar,h1,'PA'
  pa=getposang(h1)

  ;pa=sxpar(h1,'parang')
  sxaddpar,h1,'PA',pa,after='crpa',format="f12.6"
  ;sxdelpar,h1,'parang'
  skippa:
 
; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]
a=bb

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

if keyword_set(gz) then begin 
;fileoutname=fileoutname+'.gz'
writefits,reducdir+fileoutname,a,h1,/compress
endif else begin
writefits,reducdir+fileoutname,a,h1
endelse

  ;writefits,reducdir+files[i],a,h1
 endfor

goto,breakout
endif 

;**Gemini/NICI**

;*need to do PA*
;**beam not set**

if instrument eq 'nici' then begin

 for i=0L,nfiles-1 do begin

;first, deal with NICI's multifits extension nonsense
;*decide if you want to do the red channel (0) or blue channel (1)
  ;if ~keyword_set(noextniri) then begin
  if channel eq 1 then begin
  if ~keyword_set(noextniri) then begin
  bb=readfits(files[i],h2,ext=1)
  endif else begin
  bb=readfits(files[i],h1)
  endelse
  h1=headfits(files[i],ext=0)
  
  exp1time=sxpar(h1,'itime_r')
  coadds=sxpar(h1,'ncoadd_r')
  endif
  if channel ge 2 then begin
  
  bb=readfits(files[i],h2,ext=2)
  h1=headfits(files[i],ext=0)
  
  exp1time=sxpar(h1,'itime_b')
  coadds=sxpar(h1,'ncoadd_b')
  endif

  sxaddpar,h1,'exp1time',exp1time,after='exptime',format="f9.5"
  sxaddpar,h1,'coadds',coadds,format="f9.5"

  if ~keyword_set(noextniri) then begin
  naxis1=sxpar(h2,'naxis1')
  naxis2=sxpar(h2,'naxis2')
  naxis=sxpar(h2,'naxis')
  endif else begin
  naxis1=sxpar(h1,'naxis1')
  naxis2=sxpar(h1,'naxis2')
  naxis=sxpar(h1,'naxis')
  endelse
  extend=sxpar(h1,'extend')
  
  sxdelpar,h1,'extend'
  sxaddpar,h1,'extend','F'

  sxaddpar,h1,'naxis',naxis
  sxaddpar,h1,'naxis1',naxis1
  sxaddpar,h1,'naxis2',naxis2
  
  
  altitude=sxpar(h1,'el')
  sxaddpar,h1,'altitude',altitude,after='el',format="f9.5"
  sxdelpar,h1,'elevatio'

;***WCS information
   ;goto,skipwcs
   if ~keyword_set(noextniri) then begin
   radecsys=sxpar(h2,'RADECSYS')
   sxaddpar,h1,'RADECSYS',radecsys
   ctype1=sxpar(h2,'CTYPE1')
   sxaddpar,h1,'CTYPE1',ctype1

   cd11=sxpar(h2,'CD1_1')
   sxaddpar,h1,'CD1_1',cd11


   cd22=sxpar(h2,'CD2_2')
   sxaddpar,h1,'CD2_2',cd22

   crval2=sxpar(h2,'CRVAL2')
   sxaddpar,h1,'CRVAL2',crval2

   crpix1=sxpar(h2,'CRPIX1')
   sxaddpar,h1,'CRPIX1',crpix1

   ctype2=sxpar(h2,'CTYPE2')
   sxaddpar,h1,'CTYPE2',ctype2

   cd12=sxpar(h2,'CD1_2')
   sxaddpar,h1,'CD1_2',cd12

   cd21=sxpar(h2,'CD2_1')
   sxaddpar,h1,'CD2_1',cd21

   crval1=sxpar(h2,'CRVAL1')
   sxaddpar,h1,'CRVAL1',crval1

   endif else begin

   radecsys=sxpar(h1,'RADECSYS')
   sxaddpar,h1,'RADECSYS',radecsys
   ctype1=sxpar(h1,'CTYPE1')
   sxaddpar,h1,'CTYPE1',ctype1

   cd11=sxpar(h1,'CD1_1')
   sxaddpar,h1,'CD1_1',cd11


   cd22=sxpar(h1,'CD2_2')
   sxaddpar,h1,'CD2_2',cd22

   crval2=sxpar(h1,'CRVAL2')
   sxaddpar,h1,'CRVAL2',crval2

   crpix1=sxpar(h1,'CRPIX1')
   sxaddpar,h1,'CRPIX1',crpix1

   ctype2=sxpar(h1,'CTYPE2')
   sxaddpar,h1,'CTYPE2',ctype2

   cd12=sxpar(h1,'CD1_2')
   sxaddpar,h1,'CD1_2',cd12

   cd21=sxpar(h1,'CD2_1')
   sxaddpar,h1,'CD2_1',cd21

   crval1=sxpar(h1,'CRVAL1')
   sxaddpar,h1,'CRVAL1',crval1


   endelse
   skipwcs:

  ;latitude and longitude
   
  lat=-30.24075d0
  lng=-70.73669d0

  sxaddpar,h1,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h1,'lng',lng,after='lat',format="f12.7"

  if keyword_set(cal) then goto, skippa2
  sxdelpar,h1,'PA'
  print,'files',files[i]

; do you need to fix the header?  This is what happens when others partially 
; process the data.  If so, fix the entries for Hour Angle and Local Sidereal Time

  if keyword_set(fixheader) then begin
  fixheader,h1,fitskey='HA '
  fixheader,h1,fitskey='ST '
  
  endif

  lst=sxpar(h1,'st')
  sxaddpar,h1,'lst',lst,after='st',format="f9.5"
  sxdelpar,h1,'st'

  pa=getposang(h1)

  sxaddpar,h1,'PA',pa,after='crpa',format="f12.6"
  skippa2:
 
; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]
a=bb

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

if keyword_set(gz) then begin 
;fileoutname=fileoutname+'.gz'
writefits,reducdir+fileoutname,a,h1,/compress
endif else begin
writefits,reducdir+fileoutname,a,h1
endelse

  ;writefits,reducdir+files[i],a,h1
 endfor

goto,breakout
endif 


;**Subaru/IRCS**
;*do PA, beam, hour angle*

;***For new version, for old version use IRCSo
if instrument eq 'ircs' then begin

 for i=0L,nfiles-1 do begin
 
  bb=readfits(files[i],h1,/silent)
  

  ;Are these predetermined to be sky frames?  If so, let the fits header know!
if keyword_set(setsky) then begin
sxdelpar,h1,'object'
sxaddpar,h1,'object','sky',after='observer',format="a3"
endif
   
  ;latitude and longitude
  lat=19.825d0
  lng=-155.4802d0

  sxaddpar,h1,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h1,'lng',lng,after='lat',format="f12.7"

  lst=sxpar(h1,'lst')
  lst=ten(lst) 
 
  ra=sxpar(h1,'ra')
  ra=ten(ra)
  
  ha=double(lst)-double(ra)
  sxaddpar,h1,'HA',ha,after='altitude',format="f9.5"
  
  sx=sxpar(h1,'D_IMRPAD',count=ctpa)

  if ctpa gt 0 then begin
   pa=sx
  endif else begin
  pa=getposang(h1,/split)
  endelse
  sxaddpar,h1,'PA',pa,after='parang',format="f12.6"

;Dithering?
  if keyword_set(dither) then begin
  beam=sxpar(h1,'I_DTHPOS')
   sxaddpar,h1,'beam',beam,after='PA',format='i5'
  endif

;Was this taken in data-cube mode?
;if so, you need to reset the number of coadds
sz=size(bb)
if sz[0] eq 3 then begin
sxdelpar,h1,'coadd'
sxdelpar,h1,'coadds'

coadds=long(sz[3])
exp1time=sxpar(h1,'exp1time')
exptime=coadds*exp1time
print,'coadds',coadds
sxaddpar,h1,'exptime',exptime,after='exp1time',format="f12.6"
sxaddpar,h1,'coadds',coadds,after='exptime',format="f12.6"
endif

;do you want to collapse the cube for a quick reduction?
if sz[0] eq 3 and keyword_set(collapse) then begin

af=total(bb,3)

bb=af
endif
  
 
; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

a=bb

if keyword_set(gz) then begin 
fileoutname=fileoutname+'.gz'
writefits,reducdir+fileoutname,a,h1,/compress
endif else begin
writefits,reducdir+fileoutname,a,h1
endelse
;  writefits,reducdir+files[i],a,h1
 endfor

goto,breakout
endif 

;***For old version, for old version use IRCSo
if instrument eq 'ircso' then begin

 for i=0L,nfiles-1 do begin
 
  bb=readfits(files[i],h1,/silent)
  
  ;latitude and longitude
   
  lat=19.825d0
  lng=-155.4802d0

  sxaddpar,h1,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h1,'lng',lng,after='lat',format="f12.7"

  lst=sxpar(h1,'lst')
  lst=ten(lst) 
 
  ra=sxpar(h1,'ra')
  ra=ten(ra)
  
  ha=double(lst)-double(ra)
  sxaddpar,h1,'HA',ha,after='altitude',format="f9.5"
  
  sx=sxpar(h1,'D_IMRPAD',count=ctpa)
  if ctpa gt 0 then begin
   pa=sx
  endif else begin
  pa=getposang(h1,/split)
  endelse
  sxaddpar,h1,'PA',pa,after='parang',format="f12.6"

;Dithering?
  if keyword_set(dither) then begin
  beam=sxpar(h1,'I_DTHPOS')
   sxaddpar,h1,'beam',beam,after='PA',format='i5'
  endif
  
 
; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

a=bb
if keyword_set(gz) then begin 
fileoutname=fileoutname+'.gz'
writefits,reducdir+fileoutname,a,h1,/compress
endif else begin
writefits,reducdir+fileoutname,a,h1,/compress
endelse
;  writefits,reducdir+files[i],a,h1
 endfor

goto,breakout
endif 

;**Subaru/HiCIAO**

if instrument eq 'hiciao' then begin

if keyword_set(robust) then begin
if  ~keyword_set(pfname) then read,"Select Info File for Precise PA Calculation ",pfname

param,'RA',ra2000,/get,pfname=pfname
param,'DEC',dec2000,/get,pfname=pfname

endif

if keyword_set(obslog) then begin
 openw,1,'obslog.txt'
endif

 for i=0L,nfiles-1 do begin
 
  bb=readfits(files[i],h1,/silent)
  
  ;latitude and longitude
  ;Are these predetermined to be sky frames?  If so, let the fits header know!

if keyword_set(setsky) then begin
sxdelpar,h1,'object'
sxaddpar,h1,'object','sky',after='observer',format="a3"
endif
   
  lat=19.825504d0
  lat_obs=lat
  lng=-155.47602d0

  ;set ra and dec at beginning of loop
  if ~keyword_set(robust) then begin
  ra0=ra2000
  dec0=dec2000
  endif else begin
  ra0=sxpar(h1,'crval1')
  dec0=sxpar(h1,'crval2')
  endelse

  sxaddpar,h1,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h1,'lng',lng,after='lat',format="f12.7"

;****Exposure Time and Clock Time during Exposure
;if the time interval of the exposure is given, then use that instead to figure out avg. HA
; and avg PA

exp1time=sxpar(h1,'exp1time')
coadds=sxpar(h1,'coadd')

hastart=sxpar(h1,'p_hststr',count=startcount)
haend=sxpar(h1,'p_hstend',count=endcount)

if (startcount gt 0 and endcount gt 0) then begin
;redo estimate of HA
hastart=ten(hastart)
haend=ten(haend)
if haend lt hastart then haend+=24
clocktime=3600*(haend-hastart)
endif else begin
clocktime=exp1time*coadd
endelse

;get the UT start
utstart=ten(sxpar(h1,'UT'))
utend=utstart+clocktime
hms0=sixty(utstart)
hms1=sixty(utend)

 ;;get the date
dateobs=sxpar(h1,'DATE')
     ymd0 = double(strsplit(dateobs,'-',/extract))
jd0=julday(ymd0[1], ymd0[2], ymd0[0],hms0[0],hms0[1],hms0[2])

epoch0 = (jd0 - 2451545d0)/365.25d0 + 2000d0


;***Determine Hour Angle***
  lst=sxpar(h1,'lst')

  lst=(strsplit(lst,':',/extract))
  lst0=ten(lst[0],lst[1],lst[2])
  lst1=lst0+clocktime/3600.


  ;precess coordinates to current epoch
precess,ra0,dec0,2000d0,epoch0
;print,lst0,ra0/15.,lst0-ra0/15.

ha0=lst0-ra0/15.
ha1=lst1-ra0/15.

if ha0 lt -12 then ha0+=24
if ha1 lt -12 then ha1+=24

; a hack in case of weirdness with the LST fits header
  if lst0 gt 24 then lst0-=24
  if lst1 gt 24 then lst1-=24
;  print,'ha ',files[i],'',lst

  ;ra=sxpar(h1,'crval1')
  ;dec=sxpar(h1,'crval2')
 
;  ha=lst-ra/15. 

sxaddpar,h1,'HA',0.5*(ha0+ha1),after='altitude',format="f9.5"
;sxaddpar,h1,'HA',ha,after='altitude',format="f9.5"
  
;  sxaddpar,h1,'PA',pa,after='parang',format="f12.6"

;Dithering?
  if keyword_set(dither) then begin
  beam=sxpar(h1,'I_DTHPOS')
   sxaddpar,h1,'beam',beam,after='PA',format='i5'
  endif



;Was this taken in data-cube mode?
;if so, you need to reset the number of coadds
sz=size(bb)
if sz[0] eq 3 then begin
sxdelpar,h1,'coadd'
sxdelpar,h1,'coadds'

coadds=long(sz[3])
exp1time=sxpar(h1,'exp1time')
exptime=coadds*exp1time



sxaddpar,h1,'exptime',exptime,after='exp1time',format="f12.6"
sxaddpar,h1,'coadds',coadds,after='exptime',format="f12.6"
endif

;***Determine Parallactic Angle***

x=ha0+(findgen(1000)/999.)*(clocktime/3600.) ;if is the HA at the beginning of the exposure
;print,avg(x),'avg',0.5*(ha0+ha1),'beginning'

panew=calc_avparang(ha0,ha1,dec0,lat=lat_obs)
if panew lt 0 then panew+=360

pa=int_tabulated(x,parangle(x,dec2000,lat),/double)/(clocktime/3600.)
;print,'stuff',ha0,pa,panew,exp1time,coadds,i
print,'parang_output',long(i)+1,panew,pa,0.5*(ha0+ha1),ha0,ha1

 sxaddpar,h1,'coadds',coadds,after='exp1time',format="i5"
 sxaddpar,h1,'PA',panew,after='parang',format="f12.6"
 ;sxaddpar,h1,'PA',pa,after='parang',format="f12.6"

;print,'PA is',PA
;print,'coadds',coadds

;do you want to collapse the cube for a quick reduction?
if sz[0] eq 3 and keyword_set(collapse) then begin

af=total(bb,3)

bb=af
endif
  
 
; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
;print,'new file name is ',fileoutname
endif

a=bb

if keyword_set(gz) then begin 
fileoutname=fileoutname+'.gz'
writefits,reducdir+fileoutname,a,h1,/compress
endif else begin
writefits,reducdir+fileoutname,a,h1
endelse
;  writefits,reducdir+files[i],a,h1

if keyword_set(obslog) then begin
filt=sxpar(h1,'FILTER01')
print,'FILT!',filt
;printf,1,files[i],fileoutname,exp1time*coadds,' ',filt,0.5*(ha0+ha1),pa,panew
printf,1,files[i],fileoutname,exp1time*coadds,filt,0.5*(ha0+ha1),pa,panew,format='(a,1x,a,1x,f6.2,1x,a3,1x,f9.5,1x,f9.5,1x,f9.5)'
endif
 endfor

goto,breakout
endif 


;**VLT/NaCo**
;*
if instrument eq 'naco' then begin

 for i=0L,nfiles-1 do begin
 
  bb=readfits(files[i],h1,/silent)
  print,'reading in ...',files[i]
  sz=size(bb)
 
;what kind of image is it? 
  imtype=sxpar(h1,'IMAGETYP')
   
  hgrep_extract,h1,'HIERARCH ESO DET DIT ',exp1time,/linenum
  hd=strsplit(h1[exp1time],'=',/extract)
  hd2=strsplit(hd[1],' ',/extract)
  sxaddpar,h1,'exp1time',hd2[0],after='IMAGETYP',format="f9.5"
  
  hgrep_extract,h1,'HIERARCH ESO DET NDIT ',coadds,/linenum
  hd=strsplit(h1[coadds],'=',/extract)
  hd2=strsplit(hd[1],' ',/extract)
  sxaddpar,h1,'coadds',hd2[0],after='exp1time',format="f9.5"

  hgrep_extract,h1,'HIERARCH ESO TEL ALT',altitude,/linenum
  hd=strsplit(h1[altitude],'=',/extract)
  hd2=strsplit(hd[1],' ',/extract)
  sxaddpar,h1,'altitude',hd2[0],after='coadds',format="f9.5"


  hgrep_extract,h1,'HIERARCH ESO TEL AZ',azimuth,/linenum
  hd=strsplit(h1[azimuth],'=',/extract)
  hd2=strsplit(hd[1],' ',/extract)
  sxaddpar,h1,'azimuth',hd2[0],after='altitude',format="f9.5"
 
  hgrep_extract,h1,'ESO ADA POSANG',posang_pos,/linenum
  if n_elements(posang_pos) gt 1 then posang_pos=posang_pos[0]
  hd=strsplit(h1[posang_pos],'=',/extract)
  hd2=strsplit(hd[1],' ',/extract)
  sxaddpar,h1,'PA',hd2[0],after='azimuth',format="f12.6"

;if image_type = object then you can dither and stuff.  
;if image_type /= object then no dithering and the following doesn't apply 

imtype=strsplit(imtype,' ',/extract)
imtype=imtype[0]
print,'imtype ',imtype

;if the frame is labeled as a sky frame then no dither
if imtype eq 'SKY' then goto,skipdither

;if the frame has no label but you know it's an acquisition exposure,
;then no dither
orr=sxpar(h1,'ORIGFILE')
ora=strsplit(orr,'_',/extract)
prefix=string(ora[1],format='(a3)')

if prefix eq 'ACQ' then goto,skipdither
;dither?
;assume a four-point dither pattern for NaCo, that's what they always do
;change if necessary

  beam=0
  if keyword_set(dither) then begin
 
   hgrep_extract,h1,'HIERARCH ESO SEQ CUMOFFSETX',xoffset,/linenum
    hd=strsplit(h1[xoffset],'=',/extract)
    hd2=strsplit(hd[1],' ',/extract)
    xoff=double(hd2[0])
   
   hgrep_extract,h1,'HIERARCH ESO SEQ CUMOFFSETY',yoffset,/linenum
    hd=strsplit(h1[yoffset],'=',/extract)
    hd2=strsplit(hd[1],' ',/extract)
    yoff=double(hd2[0])

;where are you in the dither pattern?
   if xoff gt 0 and yoff gt 0 then beam = 0
   if xoff lt 0 and yoff gt 0 then beam = 1
   if xoff lt 0 and yoff lt 0 then beam = 2
   if xoff gt 0 and yoff lt 0 then beam = 3

   sxaddpar,h1,'beam',beam,after='PA',format='i5'
  endif

;do you want to collapse the cube???
  if keyword_set(collapse) then begin


;do you want to sigma-clip or is this a quick and dirty reduction?
   if keyword_set(sigclip) then begin
a=fltarr(sz[1],sz[2])
;approximate centroid positions from fits header
    xcen=sz[1]/2.+xoff
    ycen=sz[2]/2.+yoff
    ;a=total(bb,3)/double(sz[3])
     a=bb[*,*,sz[3]-1]

    print,'xycen',xcen,ycen
;boxsize for sigma-clipping.  Change if necessary
    if (xcen-boxsize/2 lt 0 or ycen-boxsize/2 lt 0) then goto,skipclip
for j=xcen-boxsize/2,xcen+boxsize/2 do begin
 for jj=ycen-boxsize/2,ycen+boxsize/2 do begin
   ;medianclip,bb[j,jj,0:sz[3]-2],pixval,5,maxiter=3
   meanclip,bb[j,jj,0:sz[3]-2],pixval,5,maxiter=3
   a[j,jj]=pixval
 endfor
endfor
skipclip:
   endif else begin

  ;a=total(bb,3)/double(sz[3])
   a=bb[*,*,sz[3]-1]

   endelse 

  endif else begin

a=bb
  endelse

skipdither:
 
  ;latitude and longitude
   
  lat=-24.62706d0
  lng=-70.404194d0

  sxaddpar,h1,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h1,'lng',lng,after='lat',format="f9.5"


;use the file name embedded in the fits header
fileoutname=sxpar(h1,'origfile') 
if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif
if keyword_set(gz) then begin 
fileoutname=fileoutname+'.gz'
writefits,reducdir+fileoutname,a,h1,/compress
endif else begin
;if keyword_set(zip) then fileoutname=fileoutname+'.Z'
writefits,reducdir+fileoutname,a,h1
endelse
;  writefits,reducdir+files[i],a,h1
 endfor

goto,breakout

endif 

;**MMT/Clio**
;* get hour angle

if instrument eq 'clio2' then begin

 for i=0L,nfiles-1 do begin
 
  a=readfits(files[i],h1,/silent)
 
  sz=size(a)
  d0=sz[0]

  exp1time=sxpar(h1,'int')
  exp1time=exp1time/1000.
  sxaddpar,h1,'exp1time',exp1time,after='int',format="f9.5"
  sxdelpar,h1,'int'
  
  altitude=sxpar(h1,'alt')
  sxaddpar,h1,'altitude',altitude,after='alt',format="f9.5"
  sxdelpar,h1,'alt'
  
  azimuth=sxpar(h1,'az') 
  sxaddpar,h1,'azimuth',after='az',format="f9.5"
  sxdelpar,h1,'az'
 
  ;latitude and longitude

  lat=ten(31,41,19.6) ;latitude of MMT
lng=-1*ten(110,53,4.4) ;longitude of MMT
lng=1.d0*lng & lat=1.d0*lat



  sxaddpar,h1,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h1,'lng',lng,after='lat',format="f11.6"

altaz2hadec,double(altitude),double(azimuth),lat,ha,dumbdec
sxaddpar,h1,'HA',ha,after='azimuth',format="f9.5"

radeg=sxpar(h1,'cra')
decdeg=sxpar(h1,'cdec')
pixscale_deg=0.029917*3600.

;Astrometry header information
sxaddpar,h,'crpix1',d0/2.+1
sxaddpar,h,'crpix2',d0/2.+1
sxaddpar,h,'crval1',radeg
sxaddpar,h,'crval2',decdeg
sxaddpar,h,'cdelt1',pixscale_deg
sxaddpar,h,'cdelt2',pixscale_deg
sxaddpar,h,'ctype1','RA---TAN'
sxaddpar,h,'ctype2','DEC--TAN'
sxaddpar,h,'cunit1','degree'
sxaddpar,h,'cunit2','degree'
sxaddpar,h,'PC001001',0.00000
sxaddpar,h,'PC001002',-1.0000
sxaddpar,h,'PC002001',-1.0000
sxaddpar,h,'PC002002',0.0000


; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

if keyword_set(gz) then begin 
fileoutname=fileoutname+'.gz'
writefits,reducdir+fileoutname,a,h1,/compress
endif else begin
writefits,reducdir+fileoutname,a,h1
endelse

;  writefits,reducdir+files[i],a,h1
endfor
endif


if instrument eq 'clio1' then begin

 for i=0L,nfiles-1 do begin
 
  a=readfits(files[i],h1,/silent)

  sz=size(a)
  d0=sz[1] & d1=sz[2]
  dz=sz[0]
  
  if dz eq 2 then coadds = 1.
  if dz gt 2 then coadds = sz[3]
  exp1time=sxpar(h1,'int')
  exp1time=exp1time/1000.
  sxaddpar,h1,'exp1time',exp1time,after='int',format="f9.5"
  sxdelpar,h1,'int'
  sxaddpar,h1,'coadds',coadds,after='exp1time',format="f9.5"
  
  altitude=sxpar(h1,'alt')
  sxaddpar,h1,'altitude',altitude,after='alt',format="f9.5"
  sxdelpar,h1,'alt'

  azimuth=sxpar(h1,'az') 
  sxaddpar,h1,'azimuth',azimuth,after='az',format="f9.5"
  sxdelpar,h1,'az'
 
  ;pa=sxpar(h1,'parang')
  ;sxaddpar,h1,'PA',pa,after='parang',format="f9.5"
  ;sxdelpar,h1,'parang'

  ;latitude and longitude

  lat=ten(31,41,19.6) ;latitude of MMT
lng=-1*ten(110,53,4.4) ;longitude of MMT
lng=1.d0*lng & lat=1.d0*lat

  sxaddpar,h1,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h1,'lng',lng,after='lat',format="f11.6"

altaz2hadec,double(altitude),double(azimuth),lat,ha,dumbdec
sxaddpar,h1,'HA',ha,after='azimuth',format="f9.5"

;get local sidereal time from hour angle and right ascension

radeg=sxpar(h1,'cra')
decdeg=sxpar(h1,'cdec')

lst=radeg*24./360.+ha

sxaddpar,h1,'LST',lst,after='azimuth',format="f9.5"
pixscale_deg=0.04857/3600.

;Astrometry header information
sxaddpar,h1,'crpix1',d0/2.+1
sxaddpar,h1,'crpix2',d1/2.+1
sxaddpar,h1,'crval1',(360/24.)*radeg
sxaddpar,h1,'crval2',decdeg
sxaddpar,h1,'cdelt1',pixscale_deg
sxaddpar,h1,'cdelt2',pixscale_deg
sxaddpar,h1,'ctype1','RA---TAN'
sxaddpar,h1,'ctype2','DEC--TAN'
sxaddpar,h1,'cunit1','degree'
sxaddpar,h1,'cunit2','degree'
sxaddpar,h1,'PC001001',0.00000
sxaddpar,h1,'PC001002',-1.0000
sxaddpar,h1,'PC002001',-1.0000
sxaddpar,h1,'PC002002',0.0000

sxdelpar,h1,'bzero'
sxdelpar,h1,'bscale'
sxdelpar,h1,'o_bscale'
 
; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

;workaround to get rid of 'bzero'
;writefits,'bah.fits',a,h1
;a-=median(a,/even)
a=float(a)
if keyword_set(gz) then begin 
fileoutname=fileoutname+'s'+'.gz'
writefits,reducdir+fileoutname,a,h1,/compress
endif else begin
writefits,reducdir+fileoutname+'s',a,h1
endelse
;  writefits,reducdir+files[i],a,h1
 endfor

goto,breakout
endif 




;**HST/ACS Data**

;*need to do PA, set HA = 0, change exptimes*
;*for now, assume that the data are dark subtracted and flatfielded already --> set exptime = 0!!!!
;**beam not set**

if instrument eq 'acs' then begin

 for i=0L,nfiles-1 do begin

;first, deal with ACS's multifits extension nonsense

  bb=readfits(files[i],h2,ext=1)
  h1=headfits(files[i],ext=0)

;remove extname
  sxdelpar,h2,'extname'  
  sxaddpar,h2,'orig_name',files[i],format="a14"
  exp1time=sxpar(h1,'exptime')
  sxaddpar,h2,'exptime',exp1time,after='bunit',format="f20.3"
  sxaddpar,h2,'exp1time',exp1time,after='exptime',format="f20.3"
  sxaddpar,h2,'coadds',1,after='exp1time',format="f9.5"
  ;sxdelpar,h1,'exptime'
  
  ;altitude=sxpar(h1,'elevatio')
  ;since obviously there is no 'altitude' in space
  sxaddpar,h2,'altitude',0,after='coadds',format="f9.5"

  sxaddpar,h2,'lst',0,after='altitude',format="f9.5"

; date of obs.  just for bookkeeping
 date=sxpar(h1,'date-obs')

  ;latitude and longitude
   
  lat=0.d0
  lng=0.d0
  sxaddpar,h2,'telescop','HST',after='altitude',format="a3"
  sxaddpar,h2,'instrument','ACS',after='telescop',format="a10"
  sxaddpar,h2,'lat',lat,after='lst',format="f9.5"
  sxaddpar,h2,'lng',lng,after='lat',format="f12.7"

  pa=sxpar(h2,'orientat') 
  ;pa=sxpar(h2,'pa_aper')

  if pa eq 0 then begin
   pa=sxpar(h2,'orientat')
  endif
  sxaddpar,h2,'PA',pa,after='lat',format="f20.8"
  sxaddpar,h2,'Date-Obs',date,after='PA',format="a10"

; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]
a=bb

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

if keyword_set(gz) then begin 
;fileoutname=fileoutname+'.gz'
writefits,reducdir+fileoutname,a,h2,/compress
endif else begin
writefits,reducdir+fileoutname,a,h2
endelse

  ;writefits,reducdir+files[i],a,h2
 endfor

goto,breakout
endif 


;**HST/WFC3 Data**

;*need to do PA, set HA = 0, change exptimes*
;*for now, assume that the data are dark subtracted and flatfielded already --> set exptime = 0!!!!
;**beam not set**

if instrument eq 'wfc3' then begin

 for i=0L,nfiles-1 do begin

;first, deal with WFC3's multifits extension nonsense

  bb=readfits(files[i],h2,ext=1)
  h1=headfits(files[i],ext=0)


;remove extname
  sxdelpar,h2,'extname'  
  sxaddpar,h2,'orig_name',files[i],format="a14"
  exp1time=sxpar(h1,'exptime')
  sxaddpar,h2,'exptime',exp1time,after='bunit',format="f20.3"
  sxaddpar,h2,'exp1time',exp1time,after='exptime',format="f20.3"
  sxaddpar,h2,'coadds',1,after='exp1time',format="f9.5"
  ;sxdelpar,h1,'exptime'
  
  ;since obviously there is no 'altitude' in space
  sxaddpar,h2,'altitude',0,after='coadds',format="f9.5"

  ;lst=sxpar(h1,'st')
  sxaddpar,h2,'lst',0,after='altitude',format="f9.5"

; date of obs.  just for bookkeeping
 date=sxpar(h1,'date-obs')

  ;latitude and longitude
   
  lat=0.d0
  lng=0.d0
  sxaddpar,h2,'telescop','HST',after='altitude',format="a3"
  sxaddpar,h2,'instrument','WFC3',after='telescop',format="a10"
  sxaddpar,h2,'lat',lat,after='lst',format="f9.5"
  sxaddpar,h2,'lng',lng,after='lat',format="f12.7"

  pa=sxpar(h2,'oorienta') 
  ;pa=sxpar(h2,'pa_aper')
  print,'pa is',pa
  sxaddpar,h2,'PA',pa,after='lat',format="f20.8"
  sxaddpar,h2,'Date-Obs',date,after='PA',format="a10"

; just save the last bit
fname=strsplit(files[i],'/',/extract)
nfname=n_elements(fname)
fileoutname=fname[nfname-1]
a=bb

if keyword_set(rename) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif

if keyword_set(gz) then begin 
;fileoutname=fileoutname+'.gz'
writefits,reducdir+fileoutname,a,h2,/compress
endif else begin
writefits,reducdir+fileoutname,a,h2
endelse

  ;writefits,reducdir+files[i],a,h2
 endfor

goto,breakout
endif 

if instrument eq 'gpi' then begin
openw,1,'reduc.log'
printf,1,'File_no,','XC,','YC,','Rsat','HA','Par.Angle'

for i=0L,nfiles-1 do begin
bb=readfits(files[i],h1,/exten)
h0=headfits(files[i],ext=0)
sz=size(bb)

  exp1time=sxpar(h1,'itime')
  coadds=sxpar(h1,'coadds')

  lst=sxpar(h0,'st')
  sxaddpar,h0,'lst',lst,after='st',format="f9.5"
  ;sxdelpar,h1,'st'

 ;Hour Angle
  ha=ten(sxpar(h0,'HA'))
  ;ha=ten(sxpar(h0,'HA'))*360/24.
  sxdelpar,h0,'HA'
  ;ha=hms2deg(ha)
  sxaddpar,h0,'HA',ha,after='lst',format="f9.5"


;for now, we're also at Cerro Pachon
  lat=-30.24075d0
  lng=-70.73669d0

  altitude=sxpar(h0,'elevatio')
  sxaddpar,h0,'altitude',altitude,after='elevatio',format="f9.5"

;now put in integration time
  sxaddpar,h0,'exp1time',exp1time,after='altitude',format="f9.5"
  sxaddpar,h0,'coadds',coadds,after='exp1time',format="f9.5"

  sxaddpar,h0,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h0,'lng',lng,after='lat',format="f12.7"

;And, parallactic angle
  pa0=sxpar(h0,'PA')
  sxaddpar,h0,'PAoff',pa0,after='PA',format="f12.6"
  sxdelpar,h0,'PA'
  pa=sxpar(h0,'PAR_ANG')
  sxaddpar,h0,'PA',pa,after='PAoff',format="f12.6"
  

;For now, just assume GPI will preprocess stuff
if keyword_set(rename) then begin

if ~keyword_set(preproc) then begin
fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
endif else begin
xcen=sxpar(h1,'PSFCENTX')
ycen=sxpar(h1,'PSFCENTY')
printf,1,float(i)+1,xcen,ycen,0,ha,pa
fileoutname=prefname+nbr2txt(i+1,4)+'reg.fits'
;fileoutname=prefname+nbr2txt(i+1,4)+'e.fits'
endelse

endif

if keyword_set(gz) then begin
writefits,reducdir+fileoutname+'.gz',0,h0,/compress
writefits,reducdir+fileoutname+'.gz',bb,h1,/compress,/append
endif else begin
writefits,reducdir+fileoutname,0,h0
writefits,reducdir+fileoutname,bb,h1,/append
endelse
endfor



endif

if instrument eq 'charis' then begin
;for now, we're at Maunakea
  lat=19.825504d0
  lng=-155.47602d0

;test the first image, see if it has the fits header info you need.

;if ~keyword_set(noexten) then begin
bbtest=readfits(files[0],h1test,ext=1)
htest=headfits(files[0])
;endif else begin
;bbtest=readfits(files[0],h1test)
;h0test=h1test
;endelse
;writefits,'ha.fits',bbtest
;stop

;if you want to manually set the exposure time and coadds then input them here...
if keyword_set(charismanual) then begin

read,'**Enter exposure time (seconds) ',exp1time_manual

read,'**Enter number of coadds ',coadds_manual
endif

;openw,1,'reduc.log'
;printf,1,'File_no,','XC,','YC,','Rsat','HA','Par.Angle'

for i=0L,nfiles-1 do begin
;if ~keyword_set(noexten) then begin
bb=readfits(files[i],h1,/exten)
h0=headfits(files[i],ext=0)
;endif else begin
;bb=readfits(files[i],h1)
;h0=h1
;endelse

sz=size(bb)

;if you don't manually set the exposure time and coadds, then use the number of reads from the file to determine exptime
  if keyword_set(charismanual) then begin
  exp1time=exp1time_manual
  endif else begin
  exp1time=1.475
  endelse
   
  if keyword_set(charismanual) then begin
  coadds=coadds_manual
  endif else begin
  firstread=sxpar(h1,'firstrd')
  lastread=sxpar(h1,'lastrd')
  coadds=long((lastread-firstread)+1)
  ;coadds=sxpar(h1,'coadds')
  endelse


  lst1=sxpar(h0,'st',count=lst1count)
  lst2=sxpar(h0,'lst',count=lst2count)
  if lst1count gt 0 then lst=lst1
  if lst2count gt 0 then lst=lst2

 ;Hour Angle, probably is not entered in here.
  ha=sxpar(h0,'HA',count=hacount)

  if hacount gt 0 then begin
  ha=ten(ha)
  sxdelpar,h0,'HA'
  ;ha=hms2deg(ha)
  sxaddpar,h0,'HA',ha,after='lst',format="f9.5"
  endif

 ;;get the date, then calculate the epoch 
  dateobs=sxpar(h0,'UTC-DATE')
  utcdate=ten(sxpar(h0,'UTC-TIME'))
  hms0=sixty(utcdate)


  ymd0 = double(strsplit(dateobs,'-',/extract))
  jd0=double(julday(ymd0[1], ymd0[2], ymd0[0],hms0[0],hms0[1],hms0[2]))

  epoch0 = (jd0 - 2451545d0)/365.25d0 + 2000d0
   ;print,jd0,epoch0

  ;***if the LST is not entered, then you need to compute it based on the MJD and longitude of Maunakea
  if lst1count eq 0 and lst2count eq 0 then begin
   jd=2400000.5+double(sxpar(h0,'mjd'))
   ct2lst,lst,lng,blah,jd

  ;now pull the right ascension and compute the difference between LST and RA
   ;put ra and dec in degrees for now.

  ;ra_string=strsplit(sxpar(h0,'ra'),':',/extract)
  ;ra=ten(ra_string[0],ra_string[1],ra_string[2])

  ra=15*ten(sxpar_charis(h0,'ra',count=racount,/justfirst))
  dec=ten(sxpar_charis(h0,'dec',count=racount,/justfirst))

   ;ha=lst-ra/15.
   ;  print,'ra is',ra,files[i]
  endif

  ;precess coordinates to current epoch
  rap=ra
  decp=dec
  precess,rap,decp,2000d0,epoch0

  ;compute the hour angle for current epoch, ra2000, and lst
  ha=lst-rap/15.


  if hacount eq 0 then begin
  sxaddpar,h0,'HA',ha,after='SHUTTER',format="f9.5"
  endif
  
  sxaddpar,h0,'lst',lst,after='HA',format="f9.5"
  ;sxdelpar,h1,'st'

  ;altitude=sxpar(h0,'elevatio')
  ;sxaddpar,h0,'altitude',altitude,after='elevatio',format="f9.5"

;now put in integration time
  sxaddpar,h0,'exp1time',exp1time,after='altitude',format="f9.5"
  sxaddpar,h0,'coadds',coadds,after='exp1time',format="f9.5"
  
  sxaddpar,h0,'telescop','Subaru',after='maskivar',format="a10"
  sxaddpar,h0,'lat',lat,after='telescop',format="f9.5"
  sxaddpar,h0,'lng',lng,after='lat',format="f12.7"

;And, parallactic angle
  pa0=sxpar(h0,'PARANG')
  sxaddpar,h0,'PA',pa0,after='PARANG',format="f12.6"


;****wavelength information****
;...will add later
sxaddpar,h1,'CRPIX3',1,after='CRPIX2',format="i4"
sxaddpar,h1,'CUNIT3','microns',after='CRPIX3',format="a8"
;sxaddpar,h1,'CRVAL3',crvalue3
  

;For now, just assume CHARIS will create the data cubes
if keyword_set(rename) then begin

;if ~keyword_set(preproc) then begin
;fileoutname=prefname+nbr2txt(i+1,4)+'.fits'
;endif else begin
;xcen=sxpar(h1,'PSFCENTX')
;ycen=sxpar(h1,'PSFCENTY')
;printf,1,float(i)+1,xcen,ycen,0,ha,pa
;fileoutname=prefname+nbr2txt(i+1,4)+'reg.fits'

fileoutname=prefname+nbr2txt(i+1,4)+'e.fits'
;endelse

endif
reducdir='./reduc/expand/'
;if keyword_set(noexten) then sxaddpar,h1,'XTENSION','IMAGE'
if keyword_set(gz) then begin
writefits,reducdir+fileoutname+'.gz',0,h0,/compress
writefits,reducdir+fileoutname+'.gz',bb,h1,/compress,/append
endif else begin
writefits,reducdir+fileoutname,0,h0
writefits,reducdir+fileoutname,bb,h1,/append
endelse
endfor



endif




breakout:

end
