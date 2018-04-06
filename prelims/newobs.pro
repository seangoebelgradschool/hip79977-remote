; Procedure written by David Lafreniere
; modified by Markus Janson
; modified for simplicity and speed by Thayne Currie
;v 2.0 - redone 

pro newobs, date=date,object=object,filter=filter,suffix=suffix,skipmot=skipmot

;setupdir,datadir=datadir,reducdir=reducdir,logfile=logfile

if ~keyword_set(date) then begin
date=''
read,'Enter date: ',date
endif
if ~keyword_set(object) then begin
object=''
read,'Enter object name: ',object
endif
if ~keyword_set(filter) then begin
filter=''
read,'Enter filter name: ',filter
endif
if ~keyword_set(suffix) then begin
suffix=''
read,'Enter suffix for .info file: ',suffix
endif

obj=strcompress(object,/remove_all)

obj_main = obj
;print,'here2',obj_main
;file_mkdir,reducdir+obj+suffix
;file_mkdir,obj+suffix
;pfname=reducdir+obj+suffix+'.info'
pfname=obj+suffix+'.info'
param,'name',object,/set,pfname=pfname

;Right now, we don't use a target database even though SEEDS did.  This is legacy code. Kept for future re-use.
if keyword_set(targetdatabase) then begin
   ; read in the targetdatabase file
   readcol, targetdatabase, obj, ra, dec, radeg, decdeg, simbad, hip, spt, vmag, plx, plx_error, $
            prop_ra, prop_ra_err, prop_dec, prop_dec_err, twomass, hmag, hmag_err, $
            format='a,a,a,f,f,a,a,a,f,f,f,f,f,f,f,a,f,f', skipline=3, delimiter='|'

   ; look for the object in that file
   obj_ind=where(strtrim(obj,2) eq object)
   
   if (obj_ind ne -1) then begin
      print, 'Resolved data for ', strtrim(object,2), ' in ', strtrim(targetdatabase), '.'

      ; copy the information into the same variables as querysimbad and queryvizier
      ra=(radeg[obj_ind])[0]
      dec=(decdeg[obj_ind])[0]

      s={pmra:(prop_ra[obj_ind])[0],e_pmra:0.0,pmde:(prop_dec[obj_ind])[0],e_pmde:0.0,plx:(plx[obj_ind])[0],$
         e_plx:(plx_error[obj_ind])[0],hmag:(hmag[obj_ind])[0],e_hmag:(hmag_err[obj_ind])[0]}

   endif else begin
   print, 'Sorry, but we did not find ', strtrim(object,2), ' in ', strtrim(targetdatabase), '.'

   endelse

;***Now begin the loop
endif else begin
   ;retrieve info on target
   ;@#$% Removed /ep2000 keyword from vizierquery -- not recognized.
   ;@#$% Enforcing the SIMBAD query to make up for the broken Vizier query.
   ;s=queryvizier('I/239/hip_main',object,20./60.,/allcolumns)
   ;if size(s,/type) eq 3 then begin

   querysimbad,object,ra,dec,id

;2mass database for JHK
   flag_s2mass=0
   s2mass=queryvizier('2MASS-PSC',[ra,dec],20./60.,/allcolumns)

;HIP database for parallax, vmag, stype etc.
   print,ra,dec

   flag_ship =0
   ship=queryvizier3('I/239/hip_main',[ra,dec],20./60.,/allcolumns)
   ;ship=queryvizier2('I/239/hip_main',[ra,dec],20./60.,/allcolumns)
   if size(ship,/type) ne 8 then flag_ship=1
   ;print,ship.e_plx
endelse

;Put parameters into parameter file

;from querysimbad
param,'ra',ra,'RA, equinox J2000, epoch J2000',/set,pfname=pfname
param,'dec',dec,'DEC, equinox J2000, epoch J2000',/set,pfname=pfname
param,'ra_hms',deg2hms(ra),'RA, equinox J2000, epoch J2000',/set,pfname=pfname
param,'dec_dms',deg2dms(dec),'DEC, equinox J2000, epoch J2000',/set,pfname=pfname


if flag_ship eq 1 then begin
plx=0
eplx=0
pmra=0
epmra=0
pmde=0
epmde=0
sptype='-999'
distance=plx
endif else begin

;Parallax and error on parallax might not have been found
plxset = where(strmatch( tag_names(ship), 'plx', /fold_case) eq 1)
eplxset = where(strmatch( tag_names(ship), 'e_plx', /fold_case) eq 1)
pmraset = where(strmatch( tag_names(ship), 'pmra', /fold_case) eq 1)
epmraset = where(strmatch( tag_names(ship), 'e_pmra', /fold_case) eq 1)
pmdecset = where(strmatch( tag_names(ship), 'pmde', /fold_case) eq 1)
epmdecset = where(strmatch( tag_names(ship), 'e_pmde', /fold_case) eq 1)
sptypeset=where(strmatch( tag_names(ship),'sptype', /fold_case) eq 1)


if (plxset gt 0) then plx=ship.plx else begin
    plx=0
    print, 'Parallax not found. Setting parallax=0.'
endelse
if (eplxset gt 0) then eplx=ship.e_plx else begin 
    eplx=0
    print, 'Error on parallax not found. Setting error on parallax=0.'
endelse

if plx gt 0 then distance = 1.d3/plx else distance =-999

if (pmraset gt 0) then begin
pmra=ship.pmra
epmra=ship.e_pmra
endif else begin
;if (pmraset gt 0) then (pmra=ship.pmra) and (epmra=ship.epmra) else begin
    pmra=0
    epmra=0
    print, 'PM RA not found. Setting (error on) PMRA=0.'
endelse
if (pmdecset gt 0) then begin
pmde=ship.pmde
epmde=ship.e_pmde
endif else begin 
    pmde=0
    epmde=0
    print, 'PM DEC not found. Setting (error on) PMDEC=0.'
endelse

if (sptypeset gt 0) then sptype=ship.sptype else begin
sptype = '-999'
Print, 'Spectral type not found, setting to dummy variable.'
endelse

endelse

;Write parallax and proper motion + associated errors
param,'plx',plx,' Parallax, mas',/set,pfname=pfname
param,'distance',distance,' Distance (pc)',/set,pfname=pfname
param,'eplx',eplx,' Error on parallax, mas',/set,pfname=pfname
param,'pmra',pmra,' Proper motion RA [*cos(DEC)], mas/yr',/set,pfname=pfname
param,'epmra',epmra,' Error on PMRA',/set,pfname=pfname
param,'pmdec',pmde,' Proper motion DEC, mas/yr',/set,pfname=pfname
param,'epmdec',epmde,' Error on PMDEC',/set,pfname=pfname
param,'sptype',sptype,' Spectral Type',/set,pfname=pfname


;@#$% Replaced the following line to use SIMBAD coordinates instead of
;     the nonexistent Vizier ones
;s=queryvizier('II/246/out',[s._raj2000,s._dej2000],20./60.,/allcolumns)

;Presumably, your source is in the 2MASS point source catalog!
if ~keyword_set(targetdatabase) then begin
  print,'radec',ra,dec
  s=queryvizier3('II/246/out',[ra,dec],20./60.,/allcolumns)
  ;s=queryvizier3('II/246/out',[ra,dec],20./60.,/allcolumns,/canada)
  ind=sort(s._r)
  s=s[ind[0]]
endif

;Put in spectral type if you have it


param,'hmag',s.hmag,'H magnitude, 2MASS',/set,pfname=pfname
param,'ehmag',s.e_hmag,'Error on HMAG',/set,pfname=pfname


goto,skipover
;retrieve data for this target
;read log file
readcol,reducdir+logfile,$
  fname,obj,objtyp,objclass,f1,f2,f3,exptime,coadds,obsdate,hwpangle, $
  format='a,a,a,a,a,a,a,f,f,a,f', skipline=5 
fname=strtrim(fname,2)
obj=strtrim(obj,2)
objtyp=strtrim(objtyp,2)
objclass=strtrim(objclass,2)
f1=strtrim(f1,2)
f2=strtrim(f2,2)
f3=strtrim(f3,2)
obsdate=strtrim(obsdate,2)
filters=strarr(n_elements(f1))
for n=0,n_elements(f1)-1 do filters[n]=filtername2(f1[n],f2[n],f3[n])

good=where(obsdate eq date and obj eq object and (filters eq filter) and objtyp eq 'OBJECT' and ((objclass eq 'SCIENCE') or $
	(objclass eq 'ENGINEER') or $
	(objclass eq 'SUBROUTINE') or $
	(objclass eq 'TOOL') or $
	(objclass eq 'DI') or (objclass eq 'PDI') or $
	(objclass eq 'LAUNCHER')), $
c,complement=bad)

if c eq 0 then stop
if (bad[0] ne -1) then begin
    remove,bad,fname,obj,objtyp,objclass,exptime,coadds,obsdate
    remove,bad,f1,f2,f3,filters,hwpangle
endif

;Get rid of '_ds' tag on destriped filenames
;and remove duplicates from list
wds = strpos(fname, '_ds')
wnotds = where(wds EQ -1)
if (wnotds[0] ne -1) then begin
    wds[wnotds] = strlen(wds[wnotds])
endif
;strmid has an odd syntax for extracting strings from an array of
;strings. All this command does is chop off the _ds.
fname= strmid(fname, 0, transpose(wds))

;Figure out which files aren't duplicates
ord = uniq(fname, sort(fname))

;Keep only the unique files
fname = fname[ord]
obj = obj[ord]
objtyp = objtyp[ord]
objclass = objclass[ord]
f1 = f1[ord]
f2 = f2[ord]
f3 = f3[ord]
obsdate = obsdate[ord]
filters = filters[ord]
exptime = exptime[ord]
coadds = coadds[ord]
hwpangle = hwpangle[ord]

;extract all exptimes used
uexptime=exptime[uniq(exptime,sort(exptime))]

;determine the identity of saturated exposures
indsat=where(exptime eq max(uexptime),nsat)
writecol, reducdir+obj_main+suffix+'_sat.txt', fname(indsat)
;determine the identity of non-saturated exposures
indsat=where(exptime eq min(uexptime),nnsat)
writecol, reducdir+obj_main+suffix+'_nsat.txt', fname(indsat)

;determine the HWP state of the frames (for polarimetry)
if keyword_set(polarimetry) then begin
	indpQ = where(hwpangle eq 0, npQ)
	if (npQ gt 0) then output = fname(indpQ) else output = ''
	writecol, reducdir+obj_main+suffix+'_pQ.txt', output, indpQ
	indmQ = where(hwpangle eq 45, nmQ)
	if (nmQ gt 0) then output = fname(indmQ) else output = ''
	writecol, reducdir+obj_main+suffix+'_mQ.txt', output, indmQ
	indpU = where(hwpangle eq 22.5, npU)
	if (npU gt 0) then output = fname(indpU) else output = ''
	writecol, reducdir+obj_main+suffix+'_pU.txt', output, indpU
	indmU = where(hwpangle eq 67.5, nmU)
	if (nmU gt 0) then output = fname(indmU) else output = ''
	writecol, reducdir+obj_main+suffix+'_mU.txt', output, indmU
endif

skipover:

param,'obsdate',date,'Observation date',/set,pfname=pfname
param,'fnum_sat','1-100','File #s, sat, main seq',/set,pfname=pfname
param,'cenref',0,'Number of file to use as ref for centering',/set,pfname=pfname

;default, canned values for other parameters
;other algorithm parameters
param,'procdir','./proc/','Subdirectory to put processed files',/set,pfname=pfname
param,'SPANG1',120,'Angle of spider diffraction spike 1',/set,pfname=pfname
param,'SPANG2',240,'Angle of spider diffraction spike 2',/set,pfname=pfname
param,'SPANG3',300,'Angle of spider diffraction spike 3',/set,pfname=pfname
param,'SPANG4',60,'Angle of spider diffraction spike 4',/set,pfname=pfname
param,'CENMASK','0,5,100','R1, R2 and width of mask used for centering',/set,pfname=pfname
param,'SPMASK',10,'Width of spider mask',/set,pfname=pfname
param,'rsat',-1.,'PSF saturation radius',/set,pfname=pfname
param,'fwhm',4.,'PSF FWHM in pixels',/set,pfname=pfname
param,'ffwhm',4.,'PSF FWHM of fake point sources',/set,pfname=pfname
param,'nfwhm',0.75,'Min displacement for sub, in # of PSF FWHM',/set,pfname=pfname
param,'deltar','5.,50.,120.,1.','Width of subtraction annuli, LOCI',/set,pfname=pfname
param,'na',300.,'Area of optimization regions, LOCI',/set,pfname=pfname
param,'geom',1.,'Geometry of optimization regions, LOCI',/set,pfname=pfname
param,'raper',2,'Photometry aperture radius in pixels',/set,pfname=pfname
param,'fpeak',-1.,'Peak PSF flux in aperture of radius raper',/set,pfname=pfname

end
