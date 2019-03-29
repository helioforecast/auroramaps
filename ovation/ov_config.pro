;+
; Fill Ovation Prime Configuration System Structure
; 
; :author: Rob Redmon, NOAA/NGDC Jan-2012
;
; :history:
; 2012-01-09 Initial version.
; 2012-01-19 Additions by JM.
; 2012-05-24 Add use of a ov_config_xxx.txt files to hold config parameters. JM
; 
; :examples:
; ov_config
;-
pro ov_config
    compile_opt idl2 ; DEFINT32 strictarr
    close, /all
 
    ;---------------------------------------------------------------------------------------------
    ;                          USER ADJUSTABLE PARAMETERS
    ;---------------------------------------------------------------------------------------------
    ;Enter directory and name of configuration file here (examples below).
    ;config_path='C:\AIDLWorkspace\ovation_pr_sourceforge\trunk\config\'
    ;config_path='C:\ovation-prime\config\"
    ;config_file='ov_config.txt'
    
    config_path='config/'
    config_file='ov_config_1.txt'
    print, 'config_path:', config_path
    print, 'config_file:', config_file
    
    
    if ( n_elements( config_path ) eq 0 or n_elements( config_file ) eq 0 ) then begin
        junk = dialog_message( "Please: 1) Edit ov_config.pro and set the variables config_path and config_file.  2) Edit the config text file." )
        print, 'ov_config: Please: 1) Edit ov_config.pro and set the variables config_path and config_file.  2) Edit the config text file.'
        stop
    endif

    ;---------------------------------------------------------------------------------------------
    ;                          END OF USER ADJUSTABLE PARAMETERS
    ;---------------------------------------------------------------------------------------------

    ov_create_config_struct           ;create structure filled with dummy variables

    ;----------- fill config array from config file --------------------------------------------
    openr, lun0, config_path+config_file, /get_lun, ERROR=err_open  ;
    if (err_open NE 0) then print_exit, 'ov_config: error opening '+config_file, err_num=err_open

    str=' '           
    n_elements_ov=n_tags(!ov_config)       ;number of integer elements in config structure
    filled= intarr(n_elements_ov)          ;used to keep track if structure is filled

    for ii=0, n_elements_ov-1 do begin     ;read in numeric values
    
      repeat begin                         ;skip blank lines, comments and lines missing value
        if EOF(lun0) then print_exit, 'ov_config: '+config_file + ' missing elements' 
        readf, lun0, str
        result=strsplit(str, count=count, /extract)   
      endrep until NOT((count EQ 0) OR (strcmp(strmid(result[0],0,1), ';')EQ 1))
 
      if (count EQ 1) then print_exit, 'ov_config: error in '+config_file     ;exit if single element is not a comment 
      
      n_values=count-1                                    ;number of values before comments
      for jj=count-1, 1, -1 do begin                      ;exclude comments from count of values
        if (strcmp(strmid(result[jj],0,1), ';') EQ 1) $
           then n_values=jj-1  
      endfor
      
      case result[0] of

        'do_test:':            !ov_config.do_test=fix(result[1])    ;run tests (or not)
        'mode:':               !ov_config.mode=fix(result[1])       ; 0: historical, 1: real-time, 2: SWPC real-time                               

        'epoch_start:':        begin                                         ;read in as: year, month, day, hour, min
                                 if (n_values NE 5) then print_exit, $
                                   'ov_config: epoch_start has too few values'
                                 !ov_config.epoch_start=julday(result[2], $     ;julday(mon, day, yr, hour, min, sec)
                                   result[3], result[1], result[4], result[5], 0)
                               end

        'epoch_end:':          begin                                                ;read in as: year, month, day, hour, min
                                 if (n_values NE 5) then print_exit, $
                                   'ov_config: epoch_end has too few values'
                                 !ov_config.epoch_end=julday(result[2], $         ;julday(mon, day, yr, hour, min, sec)
                                    result[3], result[1], result[4], result[5], 0)
                               end

        'epoch_delta_minute:': !ov_config.epoch_delta_minute=fix(result[1])        ;time spacing between model runs
        'do_psfile:':          !ov_config.do_psfile=fix(result[1])                 ;(don't) create plot file(s)

        'do_txtfile:':         !ov_config.do_txtfile=fix(result[1])                ;create/don't createtext file(s)
        'do_interpolation:':   !ov_config.do_interpolation=fix(result[1])          ;(don't) interpolate near midnight
        'n_types:':            !ov_config.n_types=fix(result[1])                   ;number of types of output files

        'atype:':              !ov_config.atype=fix(result[1:n_values])            ;array of auroral type(s)        
        'jtype:':              !ov_config.jtype=fix(result[1:n_values])            ;array of output type(s)    
        'show_model_stats:':   !ov_config.show_model_stats=fix(result[1])          ;(don't) print model stats on plots 

        'do_database:':        !ov_config.do_database=fix(result[1])               ;(don't) update MySQL database 
        'proxy_host:':         if (n_values GT 0) then !ov_config.proxy_hostname=result[1] ;proxy host name if behind firewall                     
        'proxy_port:':         if (n_values GT 0) then !ov_config.proxy_port=result[1]     ;proxy name if behind firewall                          

        'dir_col:':       !ov_config.dir_col=result[1]            ;color file directory
        'dir_log:':       !ov_config.dir_log=result[1]            ;log file directory
        'dir_omni2:':     !ov_config.dir_omni2=result[1]          ;omni2 data directory
        'dir_output:':    !ov_config.dir_output=result[1]         ;output directory
        'dir_premodel:':  !ov_config.dir_premodel=result[1]       ;premodel directory
        'dir_testdata:':  !ov_config.dir_testdata=result[1]       ;test_data directory
        
         else: print_exit, 'ov_config error. file: '+ config_file + ', line: '+ str
      endcase

      struct_names=tag_names(!ov_config)
      for jj=0, n_elements_ov-1 do begin
        if (strcmp(result[0], struct_names[jj]+':', /FOLD_CASE) EQ 1) then filled[jj]=1
      endfor

    endfor
    close, lun0

    index=where(filled EQ 0, count)
    if (count GT 0) then print_exit, 'ov_config: config file missing '+struct_names[index[0]]
   
    ;---------------------------------------------------------------------------------------------
    if (!ov_config.do_test GT 0) then ov_config_for_test $         ;change config for test case
    else if (!ov_config.mode GE 1) then ov_config_for_realtime    ;adjust config for realtime (if not in test case)
    ;---------- Configuration Consistency Checks -----------------------------------------------
    n_types =!ov_config.n_types
    if (min(!ov_config.atype[0:n_types-1]) LT 0) OR (max(!ov_config.atype[0:n_types-1]) GT 3) $  
        then print_exit, 'ov_config: atype value out of range'    ; check if values are reasonable
       
    if (min(!ov_config.jtype[0:n_types-1]) LT 1) OR (max(!ov_config.jtype[0:n_types-1]) GT 6) $   
        then print_exit, 'ov_config: jtype value out of range' 
 
    curr_julian_date = SYSTIME( /JULIAN)
    julian_date_1960=julday(1, 1, 1960, 1, 1, 1)       
    if ((!ov_config.epoch_start GT !ov_config.epoch_end) $
        OR (!ov_config.epoch_start LT julian_date_1960) $
        OR (!ov_config.epoch_end GT curr_julian_date))  $
     then print_exit, 'ov_config: error with epoch dates'

end

;---------------------------------------------------------------------------
;---------------------------------------------------------------------------
pro ov_create_config_struct            ;create structure with dummy variables
    compile_opt idl2 ; DEFINT32 strictarr
;    p = path_sep()   ; the path separator for this OS
    dummy_array=intarr(20)-1
    ov_config = { $         ;--------- define config structure with dummy values -----------------------                                       
        do_test: -1,          $   ; 0: don't run tests,  
                                  ; 1: text file tests for energy fluxes, 2: plot file tests for energy fluxes
                                  ; 3: text file tests for number fluxes, 4: plot file tests for number fluxes
        mode: -1,             $   ; 0: historical, 1: real-time, 2: SWPC real-time                               
        epoch_start: julday(1,1,2000,1,1,0), $      ;julday(mon, day, yr, hour, min, sec)
        epoch_end: julday(1,1,2000,1,1,0), $       ;julday(mon, day, yr, hour, min, sec)
        epoch_delta_minute: -1, $
        do_psfile: -1,       $    ; 0: don't create plot file(s), 1: create plot file(s)
        do_txtfile: -1,       $   ; 0: don't create data text file(s), 1: create text file(s)
        do_interpolation: -1,$    ; 0: no interpolation, 1: do simple interpolation in midnight region
        n_types:-1,          $    ; number of types of output files: should match number of
                                  ; elements in atype0 and jtype0 arrays
        atype: dummy_array, $     ;auroral types: 0=diff, 1=mono, 2=wave, 3=ions
        jtype: dummy_array, $     ;output types: 1=electron energy flux, 2=ion energy flux,
                                  ;              3=e- number flux,       4=ion number flux, 
                                  ;              5=e- average energy,    6=ion average energy
        show_model_stats:-1, $    ; 0: don't print model stats on draw_je plots (e.g. number of  
                                  ;   satellite days and empirical range, 1: print stats on plots)
                                  
        do_database:0, $          ;Used for real time system. 1=save model run information in  
                                  ;   MySQL database. 0=don't use database.
        proxy_host: '', $         ;For use of real time system behind a firewall.  Otherwise leave blank.                     
        proxy_port: '', $         ;For use of real time system behind a firewall.  Otherwise leave blank.                           

        dir_col: ' ', $           ; color file directory
        dir_log: ' ', $           ; log file directory
        dir_omni2: ' ', $         ; omni2 data directory
        dir_output: ' ', $        ; output directory
        dir_premodel: ' ', $      ; premodel directory
        dir_testdata: ' ' $      ; test_data directory                                         
    }
    defsysv, '!ov_config', ov_config

end
;---------------------------------------------------------------------------
pro ov_config_for_realtime     ;Set some !ov_config parameters if running real time
  compile_opt idl2 ; DEFINT32 strictarr

  case !ov_config.mode of    
  1: begin
       !ov_config.do_psfile = 1
       !ov_config.do_txtfile = 1
       !ov_config.show_model_stats = 0  ; don't imprint model stats on plot output
     end
  
  2: begin                              ;SWPC case
     end
  endcase
end    
;----------------------------------------------------------------------------
pro ov_config_for_test
    compile_opt idl2 ; DEFINT32 strictarr
;Set !ov_config parameters if doing a test case, i.e., !ov_config.do_test GE 1.
;DO NOT ADJUST THE PARAMETERS IN THIS ROUTINE.

    !ov_config.epoch_delta_minute=60.0   ;set test parameters 
    !ov_config.do_interpolation=0
    !ov_config.n_types=4             ; number of types of output files: should match number of 
                                         ;   elements in atype0 and jtype0 arrays 
    !ov_config.atype=[0, 1, 2, 3]    ; auroral types: 0=diff, 1=mono, 2=wave, 3=ions

         

    !ov_config.jtype=[1, 1, 1, 2]    ; output types:  1=electron energy flux, 2=ion energy flux, 
                                         ;   3=e- number flux, 4=ion number flux, 5=e- average energy
                                         ;   6=ion average energy
   if ((!ov_config.do_test EQ 3) OR (!ov_config.do_test EQ 4)) then $
       !ov_config.jtype=[3, 3, 3, 4]
    
    ; Set test parameters for test file tests (cases 1 and 3)
    if ((!ov_config.do_test EQ 1) OR (!ov_config.do_test EQ 3)) then begin
        !ov_config.epoch_start=julday( 10, 24, 2011,  14,  0,  0 ) ; julday(mon, day, yr, hour, min, sec)
        !ov_config.epoch_end=julday( 10, 25, 2011, 3, 59,  0 )
        !ov_config.do_psfile=0                                     ; don't create plot files
        !ov_config.do_txtfile=1                                    ; create data text files
    endif
    
    ; Set test parameters for test plot files (test cases 2 and 4)
    if ((!ov_config.do_test EQ 2) OR (!ov_config.do_test EQ 4)) then begin                     
        !ov_config.epoch_start=julday( 10, 25, 2011,  1,  0,  0 ) ; julday(mon, day, yr, hour, min, sec)
        !ov_config.epoch_end=julday( 10, 25, 2011, 1, 0,  1 )
        !ov_config.do_psfile=1                                    ; create plot files
        !ov_config.do_txtfile=0                                   ; don't create data text files
    endif
end    
        
;---------------------------------------------------------------------------
pro print_exit, str, err_num=err_num                             ;print string, close files, and exit
    compile_opt idl2 ; DEFINT32 strictarr
    
    caldat,systime(1, /julian), month0, day0, year0, hour0, minute0    ;fill with numeric variables
    time_stamp=string(year0, month0, day0, hour0, minute0, $                  ;ISO 8601 e.g., 2012-07-02T06:39Z
      format='(i4, "-", i2.2, "-", i2.2, "T", i2.2, ":", i2.2, "Z")')
    if (!ov_config.mode EQ 2) then prefix=timestamp+':ERROR:' $
    else prefix='**** '
    
    print   
    print, prefix, str                                                      
    if ( n_elements( err_num ) ne 0 ) then $                     ;print error number and message if provided
      print, prefix, err_num, ' ', !ERROR_STATE.MSG
    print
    
    close, /all 
    stop
end

