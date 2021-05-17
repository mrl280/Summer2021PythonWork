@'e:\IDL\stats_hok\Eregion\programs_2013\kat_package\plotlib.pro'
@'e:\IDL\peter\work_2018\INV_cly_comparison\polygone_plot.pro

pro scatter_10_12mhz

ang=(findgen(13)/6)*3.14
usersym,0.4*cos(ang),0.4*sin(ang);,/fill


year=2016

for mm=3,3 do begin
month=mm
ff0=+string(format='(i02)',month)+'_10_12_diff_RKN'
;ff0=+string(format='(i02)',month)+'_10_12_diff_RKN_summer'
;ff0=+string(format='(i02)',month)+'_10_12_diff_RKN_equinox'
;ff0=+string(format='(i02)',month)+'_10_12_diff_RKN_equinox_plus_october'
;ff0=+string(format='(i02)',month)+'_10_12_diff_RKN_15_21_UT_equinox'
;ff0=+string(format='(i02)',month)+'_10_12_diff_RKN_15_21_ut_summer'
;ff0=+string(format='(i02)',month)+'_10_12_diff_RKN_all_data' ;; entire database for picked observations
ff0=+string(format='(i02)',month)+'_10_12_diff_RKN_15_21_ut_all_data'; entire dta base for all observations

;ff0='06_06_10_12_diff_RKN_14_23_ut' ;; extended period for March 06

filename0='e:\idl\michael\10_12_mhz_COMPARISON\rkn\2016\'+string(format='(i02)',month)+'\'+string(ff0)

filename=filename0+'.txt


file_ps='e:\idl\michael\10_12_mhz_COMPARISON\rkn\2016\03\scatter_'+string(ff0)+'.ps'

  SET_PLOT,'PS';
   DEVICE,/port,xsize=5.5,ysize=10,xoffset=0.5, yoffset=0.1,/inches,FILENAME=file_ps, $
   FONT_SIZE=17,/COLOR

  LOADCT,12 ; load the correct color table

  !P.multi=[0,2,4]

;!P.charthick=3

file0=filename;'e:\idl\peter\5_min_12MHz_data\2016'+string(format='(I02)', month_select)+string(format='(I02)', day_select)+'.txt'
queue0=FILE_LINES(file0)
raw0=FLTARR(14,queue0-1)
;print,'Number of 12MHz data points=',queue0
fh0=101
header0=' '
OPENR,fh0,file0

readf,fh0,header0
READF,fh0,raw0
CLOSE,fh0
month1=raw0[0,*]
day1 =raw0[1,*]
UT1 =raw0[3,*]
vel_10MHz1=raw0[4,*]
vel_12MHz1=raw0[8,*]
num_points_10_1=raw0[2,*]
num_points_12_1=raw0[6,*]
gate= raw0[10,*]

UT_min=15
UT_max=21


for nn=0,3 do begin

;gg_min=gate(0)+nn*3
;gg_max=gate(0)+nn*3+2

gg_min=gate(0)+3+nn*6
gg_max=gate(0)+3+nn*6+5

test_for_gates=where(gate ge gg_min and gate le gg_max and UT1 ge UT_min and UT1 le UT_max and abs(vel_10MHz1) gt 10 and abs(vel_12MHz1) gt 10 and num_points_10_1 gt 2 and num_points_12_1 gt 2, count)
print,count

plot,findgen(24), -24+2*findgen(24),psym=4,symsize=1.0, YRANGE=[-1200,1200],YSTYLE=1,xRANGE=[-1200, 500],xSTYLE=1,xtitle='LOS Velocity (10 MHz) [m/s]',ytitle='LOS Velocity (12 MHz) [m/s]',$
/nodata, xminor=5,font=0;, title=ff0;,/isotropic;, title='2016'+string(format='(I02)', month_select)+string(format='(I02)', day_select) ;col=100,
;if count gt 0 then plots, vel_10MHz1(test_for_gates),vel_12MHz1(test_for_gates),psym=8;,symsize=0.2
if count gt 0 then plots, vel_10MHz1(test_for_gates),vel_12MHz1(test_for_gates),psym=3;,symsize=0.2

if count gt 0 then xyouts, -1000,1000,'n='+string(format='(i4)',n_elements(vel_10MHz1(test_for_gates))),font=0,charsize=0.6
if count gt 0 then xyouts, -1000,800,'gg: '+string(format='(i02)',gg_min)+'-'+string(format='(i02)',gg_max),font=0,charsize=0.5
;xyouts, -200,1000,'UT:'+string(format='(i02)',UT_min)+'-'+string(format='(i02)',UT_max),font=0,charsize=0.6

if nn eq 0 then xyouts, 300,-1000,'a',font=0,charsize=1.0
if nn eq 1 then xyouts, 300,-1000,'b',font=0,charsize=1.0
if nn eq 2 then xyouts, 300,-1000,'c',font=0,charsize=1.0
if nn eq 3 then xyouts, 300,-1000,'d',font=0,charsize=1.0



oplot, [-1500,1500],[-1500,1500],linestyle=0,col=100
oplot, [00,00],[-1500,1500],linestyle=2
oplot, [-1500,1500],[00,00],linestyle=2

oplot, [-1100,-100],[-400,-400],linestyle=0,col=200
oplot, [-400,-400], [-1100,-100],linestyle=0,col=200

;oplot, [-1100,-400],[-200,-200],linestyle=0,col=200
endfor
endfor

print, ' DONE'

if !d.name eq 'PS' or !d.name eq 'ps' then device,/close


END
