@utilities_season.pro

;+
; season_epoch.pro:  IDL program to make a precipitation map appropriate to a
; specified date and time -- with seasonal variations incorporated
; The focus is on the combined (two hemisphere) result, which then
; applies to the Northern Hemisphere.  The southern hemisphere can be
; taken to be 182 days out of phase with the called time
;
; :history:
; Patrick Newell, November 2009.
;     Writes a ps file (plot) and a text file giving flux in each bin
; Janet Machol - 12 Jan 2012 - replaced loops with array math in several places
; Janet Machol - 12 Jan 2012 - replaced hardwired dir names with ones from !ov_config
; Rob Redmon - 16 Feb 2012 - 
;     Added keyword for supplying own Solar Wind data.  
;           Implemented for real time though could be used for injecting synthetic solar wind data on the fly.
;     Added keyword GRID_VALUES to pass flux back to caller.  Note that this variable is locally called 'je'.
; Janet Machol - 3 March 2012 
;      Moved read array statements to fill_arrays.pro
;      Added electron auroral type as a first element to b1a, b1p, b2a, b2p and added '_all' to array names
;        e.g., b1p_all = fltarr(3,4,nmlt,nmlat) was b1p = fltarr(4,nmlt,nmlat)
;-
pro season_epoch,out_file,atype,jtype,yyyy,doy,sod, $
  b1a_all, b2a_all, b1p_all, b2p_all, jtype3_or_4, Prob_all, jtime, SW_DATA=SW_DATA, GRID_VALUES=je
common prob_array, Prob


; find dFdt (solar wind coupling function) via realtime (ACE) or retrospective (OMNI)
if ( n_elements( SW_DATA ) gt 0 ) then begin
    ap_inter_sol_realtime, SW_DATA, mag, bx, by, bz, ni, v, Ec    ; very similar to ap_inter_sol.  uses input Solar Wind 
                  ;instead of reading a file.  Also does not de-weight the current hour (should possibly reconsider this).
endif $
else begin
  ;  ap_inter_sol,yyyy,doy,sod,mag,bx,by,bz,ni,v,Ec
     ap_inter_sol_predstorm,yyyy,doy,sod,mag,bx,by,bz,ni,v,Ec

endelse

dFdt = Ec(19)

if( mag eq 0. )then begin
  print,'No IMF data available',yyyy,doy,sod
  return
endif


;call with atype=0 for diffuse,  atype=1=mono,  2=wave,  3=ions
;jtype=1 for electron energy flux; 2=ion energy flux; 3=e- number flux
;jtype=4=ion number flux;  5=e- average energy  6=ion average energy
;jtype=7 = correlation coefficient, 8=corr coeff ions
;afile = regression coefficients for auroral energy (or number) flux
;pfile = regression coefficients for probabilities
;yyyy = four digit year
;doy = day of year
;sod = second of day
;n_or_s=3 for combine north and south.  In effect this is the only
;option that works.  The result is appropriate to the northern
;hemisphere.  To get a result appropriate to the southern hemisphere,
;call with doy = 365 - actual doy
;out_file is the ps file output

;;;grid size of model (fixed)
atype_string = ['diff','mono','wave','ions']
nmlt = 96                           ;number of mag local times in arrays
nmlat = 160                         ;number of mag latitudes in arrays
ndF = 12
n_or_s = 3
je_all = fltarr(4,nmlt,nmlat)       ;unweighted fluxes

;fill probability arrays if first time in season_epoch or if jtype
;  has switched to or from number fluxes (jtype is 3 or 4) 
if (n_elements(Prob_all) EQ 0) then $    
  fill_arrays, b1a_all, b2a_all, b1p_all, b2p_all, $
    jtype, jtype3_or_4, ndF, nmlat, nmlt, Prob_all

if (((jtype3_or_4 EQ 0) AND ((jtype EQ 3) OR (jtype EQ 4))) $
  OR ((jtype3_or_4 EQ 1) AND (jtype NE 3) AND (jtype NE 4))) then $
  fill_arrays, b1a_all, b2a_all, b1p_all, b2p_all, $
    jtype, jtype3_or_4, ndF, nmlat, nmlt, Prob_all

sf0 = 0

for iseason=0,3 do begin
;Prob = Prob_all(iseason,0:2,0:nmlt-1,0:nmlat-1)
  Prob = fltarr(3,ndF,nmlt,nmlat)

  Prob=reform(Prob_all[iseason, *,*,*,*])  
  
  
  prob_curr = fltarr(nmlt,nmlat)
  if( atype le 2 )then begin
    for i=0,nmlt-1 do begin
      for j=0,nmlat-1 do begin
        b1t = b1p_all(atype, iseason,i,j)
        b2t = b2p_all(atype, iseason,i,j)
        prob_curr(i,j) = prob_estimate(b1t,b2t,dFdt,atype,i,j)
      endfor
    endfor
  endif
  if( atype ge 3 )then prob_curr[*,*] = 1.0            

  je = fltarr(nmlt,nmlat)
  for i=0,nmlt-1 do begin
    for j=0,nmlat-1 do begin
      if( prob_curr(i,j) gt 0. )then begin
        je(i,j)=(dFdt*b2a_all(atype, iseason,i,j)+ b1a_all(atype, iseason,i,j))*prob_curr(i,j)
        if( je(i,j) lt 0. )then je(i,j) = 0.
     
        if( (atype le 2) and (jtype eq 1) )then begin
          if( je(i,j) gt 10. )then je(i,j) = 0.5
          if( je(i,j) gt 5. )then je(i,j) = 5.
        endif
 
        if( (atype le 2) and (jtype eq 3) )then begin
          if( je(i,j) lt 0. )then je(i,j) = 0.
          if( je(i,j) gt 2.0e10 )then je(i,j) = 0.
          if( je(i,j) gt 2.0e9 )then je(i,j) = 1.0e9
        endif

        if( (atype eq 3) and (jtype eq 2) )then begin
          if( je(i,j) lt 0. )then je(i,j) = 0.
          if( je(i,j) gt 4. )then je(i,j) = 0.25
          if( je(i,j) gt 2. )then je(i,j) = 2.
        endif

        if( (atype eq 3) and (jtype eq 4) )then begin
          if( je(i,j) lt 0. )then je(i,j) = 0.
          if( je(i,j) gt 5.0e8 )then je(i,j) = 0.
          if( je(i,j) gt 1.0e8 )then je(i,j) = 1.0e8
        endif
      endif
    endfor
  endfor
  je_all(iseason,*,*) = je
endfor

season_weights,doy,w0,w1,w2,w3
je[*,*] =w0*je_all[0,*,*] + w1*je_all[1,*,*] + w2*je_all[2,*,*] + $  
         w3*je_all[3,*,*]

time_string = string(yyyy,format='(i4.4)')
time_string = time_string + ' ' + string(doy,format='(i3.3)')
time_string = time_string + ' ' + string(sod,format='(i5.5)')

jmin = j_plot_min(atype,jtype)
jmax = j_plot_max(atype,jtype)

filename=out_file
jtitle = atype_string(atype) + time_string

if( (jtype eq 3) or (jtype eq 4) )then begin
 jtitle = atype_string(atype) + ' num ' + time_string
endif


sf0=0                          ;set season for printing in draw_je
if (w1 GE 0.5) then sf0=1
if (w2 GE 0.5) then sf0=2
if (w3 GE 0.5) then sf0=3

draw_je,filename,je,n_or_s,jmin,jmax,jtitle,jtype,dS,sf0,jtime

close,/all                            ;prevents IDL from running out of LUs

return
end


;---------------------------------------------------------------------------------
pro fill_arrays, b1a_all, b2a_all, b1p_all, b2p_all, $
jtype, jtype3_or_4, ndF, nmlat, nmlt, Prob_all

common data_interval,y0,d0,yend,dend,files_done
 
jtype3_or_4=0                                     ;not number fluxes
if ((jtype EQ 3) OR (jtype EQ 4)) then jtype3_or_4=1  ;number fluxes

;Files giving the seasonal probability regression files have to exist
season = ['winter','spring','summer','fall']
atype_string = ['diff','mono','wave','ions']

;When probability fits fail, this array is used as a backup
Prob_all = fltarr(4,3,ndF,nmlt,nmlat)
b1p_all = fltarr(4, 4,nmlt,nmlat)
b2p_all = fltarr(4, 4,nmlt,nmlat)
b1a_all = fltarr(4, 4,nmlt,nmlat)
b2a_all = fltarr(4, 4,nmlt,nmlat)

for atype=0,3 do begin           ;loop on electron auroral types
  for iseason=0,3 do begin
    afile = !ov_config.dir_premodel + season(iseason) + '_' + atype_string(atype)
    if( (jtype eq 3) or (jtype eq 4) )then afile = afile + '_n'
    afile = afile + '.txt'

    pfile = !ov_config.dir_premodel + season(iseason) + '_' $
          + 'prob_b_' + atype_string(atype) + '.txt'

;need to calculate probability functions from linear regressions
;read those in
    if( atype le 2 )then begin       
      b1 = 0. 
      b2 = 0.
      yend = 1900
      dend = 1
      y0 = 1900
      d0 = 1
      files_done = 0
      sf0 = 0
      openr,20,pfile, ERROR=err_open
      if (err_open NE 0) then print_exit, 'season_epoch: error opening '+pfile
    
      readf,20,y0,d0,yend,dend,files_done,sf0
      for i =0, nmlt-1 do begin
        for j=0, nmlat-1 do begin
          readf,20,b1,b2
          b1p_all(atype, iseason,i,j) = b1
          b2p_all(atype, iseason,i,j) = b2
        endfor
      endfor
      for i=0,nmlt-1 do begin
        for j=0,nmlat-1 do begin
          readf,20, Prob_all(iseason,atype,0:ndF-1,i,j)
        endfor
      endfor
      close,20
    endif

;now read in regression coefficients for auroral flux
    openr,20,afile, ERROR=err_open 
    if (err_open NE 0) then print_exit, 'season_epoch: error opening '+afile
    readf,20,y0,d0,yend,dend,files_done,sf0
    i = 0
    j = 0
    b1 = 0.
    b2 = 0.
    while (not eof(20) )do begin
      readf,20,i,j,b1,b2,rF
      b1a_all(atype, iseason,i,j) = b1
      b2a_all(atype, iseason,i,j) = b2
    endwhile
    close,20
  endfor                          ;end loop on iseason
endfor                          ;end loop on auroral types

return
end











