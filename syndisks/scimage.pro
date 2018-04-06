function scimage, $
   g=g, $                       ; anisotropic scattering factor, -1<g<1 (0=isotropic)
   xdo=xdo,ydo=ydo, $           ; disk offset in the disk frame, in AU
   e=e,theta0=theta0, $         ; eccentricity and argument of pericenter
   itilt=itilt, $               ; disk inclination, in degrees (0=pole-on)
   PA = PA, $                   ; position angle of the disk major axis (East of North), in degrees
   nx=nx,ny=ny, $               ; numbers of pixels along the x- and y-axis
   pfov=pfov, $                 ; pixel field of view, in arcsec
   dens = dens, $               ; density structure
   dstar=dstar                  ; stellar distance, in parsec


;+
;
; WARNING: This version does not work for edge-on disks (itilt = +90deg modulo 180deg)
;
; Author: J.-C. Augereau, Grenoble.
;         Part of the GRaTer code for optically thin disks (Augereau et al. 1999)
;         + Simplified version, June 6, 2010, for E. Buenzli
;         + Case of an eccentric ring added, July 8, 2010
;-

;-- Default values:
if not keyword_set(g) then g = 0.d0
if not keyword_set(xdo) then xdo = 0.d0
if not keyword_set(ydo) then ydo = 0.d0
if not keyword_set(e) then e = 0.d0
if not keyword_set(theta0) then theta0= 0.d0
if not keyword_set(itilt) then itilt = 0.d0
if not keyword_set(PA) then PA = 0.d0
if not keyword_set(nx) then nx = 350
if not keyword_set(ny) then ny = 350
if not keyword_set(pfov) then pfov = 0.01325d0
if not keyword_set(dstar) then dstar = 35.d0
if not keyword_set(dens) then $
   dens = {type:'2powerlaw', density0:1.d0, r0:55.d0, alphain:10.d0, alphaout:-4.d0, $
           beta:1.d0, gamma:2.d0, ksi0:2.d0, rmin:0.d0, rmax:0.d0}


;-- Initiate the calculation:
pixAU = pfov*dstar                               ; 1 pixel in AU
itilt_rad = itilt *!DPi / 180.d0                 ; radians
PA_rad = (90.d0+PA) *!DPi / 180.d0               ; radians
e2 = e*e                                         ; squared eccentricity
twopi = 2.d0*!Dpi                                ; 2 times Pi
theta0_rad = theta0 * !DPi / 180.d0              ; radians
cs = cos(itilt_rad)                              ; cosinus of inclination angle
sn = sin(itilt_rad)                              ; sinus of inclination angle
xc = nx / 2                                      ; star center along the x-axis, in pixels
yc = ny / 2                                      ; star center along the y-axis, in pixels
n2D = float(nx) * float(ny)                      ; total number of pixels in an image
idxy = findgen(n2D)                              ; ids of all pixels
xtmp = (idxy mod nx) - xc                        ; x values for all pixels in the sky plane, in pixels
y = (floor(idxy / nx)) - yc                      ; y values for all pixels in the sky plane, in pixels
x = (cos(PA_rad)*xtmp + sin(PA_rad)*y)*pixAU     ; rotation to get the disk major axis properly oriented, x in AU
y = (-sin(PA_rad)*xtmp + cos(PA_rad)*y)*pixAU    ; same for y, in AU


;-- Compute rmin, rmax and zmax:
if dens.rmin eq 0 then dens.rmin = sqrt(xdo^2+ydo^2)+pixAU ; in AU
p = 0.005d0                     ; percentage
case dens.type of    
   '2powerlaw': begin
      if dens.rmax eq 0.d0 then dens.rmax = dens.r0 * p^(1.d0/dens.alphaout) ; maximum distance of integration, AU
      r = sqrt((x-xdo)^2+(y/cs-ydo)^2)                                       ; miplane distance to the disk center, AU
      zmax = dens.ksi0 *(r/dens.r0)^dens.beta * (-alog(p))^(1.d0/dens.gamma) ; max altitude, depends on pixel position, AU

      apeak = dens.r0 * (-(dens.alphain+dens.beta)/(dens.alphaout+dens.beta)) $
              ^(1.d0/(2.d0*(dens.alphain+dens.beta-2.d0*dens.alphaout+dens.beta)))
      print,'   Density peak semi-major axis = '+strtrim(string(apeak),2)+' AU. Eccentricity = '+strtrim(string(e),2)
      print,'   Density peak position between '+strtrim(string(apeak*(1.d0-e)),2)+' AU and '+ $
            strtrim(string(apeak*(1.d0+e)),2)+' AU'
      print,'   Disk vertical FWHM = '+strtrim(string(2.d0*dens.ksi0*alog(2.d0)^(1.d0/dens.gamma)),2)+ $
                                             ' AU at a semi-major axis of '+strtrim(string(dens.r0),2)+' AU'
   end
   'hd141569': begin
      r = sqrt((x-xdo)^2+(y/cs-ydo)^2)                                       ; miplane distance to the disk center, AU
      zmax = dens.ksi0 *(r/dens.r0)^dens.beta * (-alog(p))^(1.d0/dens.gamma) ; max altitude, depends on pixel position, AU
   end
endcase
print,'   Minimum distance of integration rmin = '+strtrim(string(dens.rmin),2)+' AU'
print,'   Maximum distance of integration rmax = '+strtrim(string(dens.rmax),2)+' AU'


;-- Line of sight distances:
lz0 =  y * tan(itilt_rad)                    ; distance to the line of sight to reach the disk midplane (z_D=0), AU
lzp = zmax/cs + lz0                          ; distance to reach +zmax, AU
lzm = -zmax/cs + lz0                         ; distance to reach -zmax, AU
dl = abs(lzp-lzm)                            ; l range, in AU
lmax2 = dens.rmax^2 - (x^2+y^2)              ; squared maximum l value to reach the outer disk radius, in AU^2
  ;idok = where(lmax2 ge 0.,complement=idex)    ; exclude regions where there is no disk to accelerate the computation
  ;dl(idok) = dl(idok) < (2.*sqrt(lmax2(idok))) ; l range cannot be larger than sqrt(lmax2)
idok = idxy

;-- Loop on the distances along the line of sight (l):
Nl = 25                                                        ; half number of l values
lwidth = 100.0d0                                               ; control the distribution of distances along l
ll = (exp(indgen(Nl)*alog(lwidth+1.d0)/(Nl-1.d0))-1.d0)/lwidth ; between 0 and 1
ll = [-reverse(ll(1:*)),ll]                                    ; between -1 and 1
Nl = 2*Nl-1                                                    ; total number of distances along the line of sight
ycs = cs*y(idok)                                               ; pre-calculated values, AU
zsn = -sn*y(idok)                                              ; pre-calculated values, AU
xd = x(idok)                                                   ; x_disk, in AU
limage = dblarr(Nl,nx,ny)+0.d0
image = dblarr(nx,ny)+0.d0

for il=0,Nl-1 do begin

   ;-- Distance along the line of sight:
   l = lz0(idok) + ll(il) * dl(idok) ; AU

   ;-- Rotation about x-axis => cartesian coordinates in the disk frame:
   yd = ycs + sn * l            ; y_Disk, in AU
   zd = zsn + cs * l            ; z_Disk, in AU

   ;-- Distances and polar angles in the frame centered on the star position:
   d2star = xd^2+yd^2+zd^2                                      ; squared distance to the star, in AU^2
   dstar = sqrt(d2star)                                         ; distance to the star, in AU
   rstar = sqrt(xd^2+yd^2)                                      ; midplane distance to the star (r coordinate), in AU
   thetastar = atan(yd,xd)                                      ; polar angle in radians
   n2pi = ceil(min(thetastar) / twopi)                          ; number of 2pi values
   thetastar = (thetastar + (n2pi+1)*twopi) mod twopi           ; to get the angle between 0 and 2Pi

   ;-- Phase angles:
   phi = acos((rstar*sn*sin(thetastar)+zd*cs)/dstar) ; in radians

   ;-- Polar coordinates in the disk frame, and semi-major axis:
   r = sqrt((xd-xdo)^2+(yd-ydo)^2)                      ; midplane distance to the disk center (r coordinate), in AU
   theta = atan(yd-xdo,xd-xdo)                          ; polar angle in radians
   n2pi = ceil(min(theta) / twopi)                      ; number of 2pi values
   theta = (theta + (n2pi+1)*twopi) mod twopi           ; to get the angle between 0 et 2Pi
   a = r*(1.d0+e*cos(theta-theta0_rad)) / (1.d0-e2)     ; semi-major axis, in AU

   ;-- Scattered light:
   rho = density(a,thetastar,zd,dens)                              ; volume density
   Phase = (1.d0-g*g)/(4.d0*!DPi*(1.d0+g*g-2.d0*g*cos(phi))^(1.5)) ; phase function
   image(idok) = rho * Phase / d2star                              ; image array used here as a tmp array
   limage(il,*,*) = image

endfor
id = where(finite(limage,/nan) eq 1)       ; remove NANs ...
if (size(id))(0) ne 0 then limage(id) = 0. ; ... if any

;-- Integration along the line of sight:
image = dblarr(nx,ny)+0.d0
for il=1,Nl-1 do image = image + (ll(il) - ll(il-1)) * (limage(il-1,*,*)+limage(il,*,*))
image = image * dl / 2.d0 * pfov^2
print,'>> Total flux : '+ strtrim(string(total(image,/nan)),2) + ' [arbitrary units]'

return,image
end
