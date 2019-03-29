
;+
; ap_inter_sol: interpolated solar wind hourly averages
; this auroral power version uses a weighted average over last 4 hours
; returns imf hourly averages
; data is in ascii files in omni2/omni2_yyyy.dat
; 
; :history:
; Patrick Newell, April 2006.  Updated to omni2 October 2009.
; 2012-02-13 R. Redmon: minor cleanup; no technical changes.
; 
; 2019-03 Christian Moestl IWF Graz - added input from PREDSTORM
; sample file
;/home/cmoestl/python/PREDSTORM/results/savefiles/predstorm_v1_realtime_stereo_a_save_2019-03-29-08_00.txt
;-
pro ap_inter_sol_predstorm, yearcall, day, sod, mag, bx, by, bz, ni, v, Ec
common ap_inter_sol,last_year,nrec,maga,bxa,bya,bza,nia,va


nh = 4       ;hours previous to integrate over
wh = 0.65     ;reduce weighting by factor of wh each hour back
ncoup = 33             ;number of coupling functions
;print,sod
;plot,[-1,1],[2,0]
hour = fix(sod/3600)

res = (float(sod) - 3600.*hour)/3600.
if( n_elements(last_year) eq 0 )then begin
  last_year = 0                        ;first call to routine
  nrec = 24*366
  maga = fltarr(nrec)
  bxa = fltarr(nrec)
  bya = fltarr(nrec)
  bza = fltarr(nrec)
  nia = fltarr(nrec)
  va = fltarr(nrec)
endif

; yearcall, day, hour are integer.  year = two digit version of yearcall
mag = 0.
bx = 0.
by = 0.
bz = 0.
ni = 0.
v = 0.

;integer*4 yr, dy, hr
r = 0
i = 0

year = yearcall

;if( year ne last_year )then begin
;close,23
;nrec = 24*365
;if ( 4*fix(year/4) eq year )then nrec  = 24*366                      ;leap year check
;sol_file = !ov_config.dir_omni2+'omni2_' + string(year, format='(i4)') + '.dat'

;*************** add in config file
sol_file='/home/cmoestl/python/predstorm/results/savefiles/predstorm_v1_realtime_stereo_a_save_2019-03-29-08_00.txt'



;############################
openr,23,sol_file, ERROR=err_open
if (err_open NE 0) then print_exit, 'ap_inter_sol_predstorm: error opening '+sol_file
close,23

;file format
;#         time      matplotlib_time B[nT] Bx   By     Bz   N[ccm-3] V[km/s] Dst[nT]   Kp   aurora [GW]
;2019  3 23  0  0  0 737141.000000   3.2  -0.3   2.4   2.1         5     303       0   0.7   0.0
;2019  3 23  1  0  0 737141.041667   3.2  -1.8   1.8   2.0         4     301     -15   0.7   0.0


data = READ_ASCII(sol_file,data_start=1)

;number of columns
col=16
datsize=n_elements(data.field01)
;number of rows
row=datsize/col
data1=fltarr(row, col)

for p=0, row-1 do begin 
    data1[p,0:15]=data.field01[0+p*col:15+p*col]
endfor

;write data into arrays

yearcall=data1[*,0]
month=data1[*,1]
dayinmonth=data1[*,2]
hour=data1[*,3]
minute=data1[*,4]
sec=data1[*,5]

mag = data1[*,7]
bx = data1[*,8]
by = data1[*,9]
bz = data1[*,10]
ni = data1[*,11]
v = data1[*,12]
w = fltarr(row)
Ec = fltarr(row)



mag_arr = fltarr(nh)
bx_arr = fltarr(nh)
by_arr = fltarr(nh)
bz_arr = fltarr(nh)
ni_arr = fltarr(nh)
v_arr = fltarr(nh)
Ec_arr = fltarr(ncoup,nh)
w_arr = fltarr(nh)



;yearcall, day, sod, mag, bx, by, bz, ni, v, Ec

stop;



;#################

;f101 = '(2i4,i3,i5,2i3,2i4,14f6.1,f9.0,f6.1,f6.0,2f6.1,f6.3,f6.2,'
;f101 = f101 + 'f9.0,f6.1,f6.0,2f6.1,f6.3,2f7.2,f6.1,i3,i4,i6,i5,f10.2,'
;f101 = f101 + '5f9.2,i3,i4,2f6.1,2i6)'
;yr = 0
;dy = 0
;hr = 0
;for r = 0, nrec-1 do begin
;readf,23,format=f101,yr,dy,hr,brn,sc,plid,nfine,npla,fmag,mag,dumlat,dumlong,$
;         bx,bygse,bzgse,by,bz,dumsig1,dumsig2,dumsig3,dumsig4,dumsig5,$
;         Tp,ni,v,vlong,vlat,nanp,pdum,sigt,sig_n,sigv,sigphi,sigtheta,$
;         sigratio,dumE,beta,machn,kp10,rspot,dstdum,aedum,prot_dum1,$
;         prot_dum2,prot_dum3,prot_dum4,prot_dum5,prot_dum6,dumflux,$
;         ap_dum,f107_dum,pcn_dum,al_dum,au_dum
;testr = hr + 24*(dy-1)
;if( r ne testr )then begin
;   print,'Discrepancy between expected time and solwind file time:'
;   print,'hr,dy,r,testr=',hr,dy,r,testr
;endif


;if( mag ge 999. )then begin
;   bx = 0.
;   by = 0.
;   bz = 0.
;   mag = 0.
;endif
;if( v ge 999. )then begin
;   v = 0.
;   ni = 0.
;endif

;maga[r] = mag
;bxa[r] = bx
;bya[r] = by
;bza[r] = bz
;nia[r] = ni
;va[r] = v
;endfor
;last_year = year
;endif





;r1 = hour + 24*(day-1)
;if( r1 ge 8784 ) then begin
;  r1 = 8783
;endif
;if( (4*(year/4)) ne year )then begin
;  if( r1 ge 8760 )then begin
;    r1 = 8759
;  endif
;endif





;################### interpolate last nh hours




w_arr = fltarr(nh)
mag_arr = fltarr(nh)
bx_arr = fltarr(nh)
by_arr = fltarr(nh)
bz_arr = fltarr(nh)
ni_arr = fltarr(nh)
v_arr = fltarr(nh)
Ec_arr = fltarr(ncoup,nh)
ng = 0


for i = 0, nh-1 do begin
  ri = r1 - i
  if( ri lt 0 )then begin
     ri = 0
  endif
  mag_arr[i] = maga[ri]
  bx_arr[i] = bxa[ri]
  by_arr[i] = bya[ri]
  bz_arr[i] = bza[ri]
  ni_arr[i] = nia[ri]
  v_arr[i] = va[ri]
  sol_coup,bx_arr[i],by_arr[i],bz_arr[i],v_arr[i],ni_arr[i],Ec
  Ec_arr[*,i] = Ec
  w_arr[i] = wh^i
  ng = ng + 1
  if( (mag_arr[i] eq 0.) or (v_arr[i] eq 0.) or (ni_arr[i] eq 0.) )then begin
    w_arr[i] = 0.
    ng = ng - 1
  endif
endfor



a = (float(sod) - 3600.*hour)/3600.
w_arr[0] = sqrt(a)*w_arr[0]

if( ng lt 3 )then begin
   mag = 0.
   bx = 0.
   by = 0.
   bz = 0.
   v = 0.
   ni = 0.
   return
endif

wt = total(w_arr)
w_arr = w_arr/wt
mag = total(w_arr*mag_arr)
bx = total(w_arr*bx_arr)
by = total(w_arr*by_arr)
bz = total(w_arr*bz_arr)
v = total(w_arr*v_arr)
ni = total(w_arr*ni_arr)
for i=0,ncoup-1 do begin
  Ec[i] = total(w_arr*Ec_arr[i,*])
endfor

return
end






