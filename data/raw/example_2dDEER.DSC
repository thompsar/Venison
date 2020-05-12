#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	CPLX
XTYP	IDX
YTYP	IDX
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
IIFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	754
XMIN	-68.000000
XWID	6024.000000
YPTS	32
YMIN	-0.000000
YWID	124.000000
*
*	Documentational Text:
*
TITL	'example_2dDEER'
IRNAM	'Intensity'
IINAM	'Intensity'
XNAM	'Time'
YNAM	'Time'
IRUNI	''
IIUNI	''
XUNI	'ns'
YUNI	'ns'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    andrew
DATE    12/13/17
TIME    08:45:22
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    PLS
OXS1    TADC
AXS1    B0VL
AXS2    NONE
AXS3    
MWPW    0.01256
A1CT    1.1904
A1SW    0.005
MWFQ    3.34377e+10
AVGS    3
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     cwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
Power              12.56 mW
PowerAtten         6 dB
QValue             300

.DVC     fieldCtrl, 1.0

AllegroMode        False
CenterField        11904.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     12140.000 G
SweepDirection     Up
SweepWidth         50.0 G

.DVC     freqCounter, 1.0

FrequencyMon       33.437700 GHz
QMonitBridge       Off

.DVC     ftBridge, 1.0

Attenuation        60.00 dB
ELDORAtt           0 dB
FrequencyA         33.507300 GHz
VideoBW            200 MHz
VideoGain          45 dB

.DVC     ftEpr, 1.0

AutoTimeOut        True
AveragesPerScan    1000
ELDORFreqStart     32.997718 GHz
ELDORFreqWidth     0.763949 GHz
FTAcqModeSlct      Run from PulseSPEL
PCycleAllowed      True
PCycleOn           False
PPExtTrg           False
PPExtTrgSlope      Rising
PlsSPELEXPSlct     2D-4P ELDOR d1avg D2O
PlsSPELGlbTxt      \
;\n\
;       ELDOR package\n\
;\n\
;     PulseSPEL general variables, definitions & conventions\n\
\n\
begin defs\n\
\n\
; Variables\n\
; Comments:\n\
\n\
p0   = 16      ; 90 pulse length\n\
p1   = 32      ; 180 pulse length\n\
p2   = 32      ; inversion pulse\n\
p3   = 32      ; ELDOR pulse width\n\
\n\
pg   = 67      ; Integrator Gate Width\n\
d0   = 572       ; Detection Trigger offset time\n\
\n\
d1   = 128     ; tau\n\
d2   = 6000    ; tau 2 (set equal to d1 for SIFTER and DQ-EPR)\n\
d3   = 60      ; T\n\
d4   = 200     ; time to start before axis=0 (added to tau2 in DQ and in SIFTER)\n\
d5   = 0       ;\n\
\n\
d9   = 0       ; DAF (Delay After Flash) delay pulse sequence until\n\
               ; after external event (eg. laser pulse)\n\
d30  = 8       ; 1st time scale increment\n\
d31  = 4       ; 2nd time scale increment\n\
\n\
a    = 1000    ; (RESERVED) number of transient averages (TRANS)\n\
               ; if greater than 1024, a must be a multiple of 1024\n\
h    = 1000    ; (CONVENTION) number of shots per point (INTG and SPT)\n\
n    = 8       ; (CONVENTION) number of sweeps to accumulate\n\
\n\
srt  = 1960 * srtu ; (RESERVED) SRT - Shot Repetition Time (srtu=1.02 us)\n\
\n\
; VARIABLES for HS_ENDOR\n\
\n\
rf1  = 10 MHz  ; frequency of RF-Source 1 in kHz\n\
df21 = 10000   ; RF sweep width in kHz\n\
\n\
df11 = 40      ; frequency increment in kHz of rf1 (calculated in sequence)\n\
\n\
p10  = 10000   ; RF pulse length\n\
d10  = 10000   ; MW delay during RF pulse\n\
d12  = 1000    ; delay before RF pulse\n\
d13  = 3000    ; delay after RF pulse\n\
\n\
; end of HS_ENDOR variables\n\
\n\
b    = 1       ;\n\
c    = 1       ;\n\
r    = 1       ;\n\
s    = 1       ;\n\
t    = 1       ;\n\
\n\
w    = 800     ; (CONVENTION) Split point of sweeps (points)\n\
               ; i.e. pulse separation where TWT gate pulse can be split.\n\
dx   = 0       ; (RESERVED) 1st time scale starting value\n\
dy   = 0       ; (RESERVED) 2nd time scale starting value\n\
\n\
eldor = 33507.33 MHz \n\
end defs\n\
\n\
; Note:\n\
;       Variables Marked 'RESERVED' = variable may be used for other times,\n\
;                                       however only this variable can change\n\
;                                        the indicated parameter.\n\
;                                       eg. only 'a' will change the number of\n\
;                                           transients averaged in SpecJet\n\
\n\
\n
PlsSPELLISTSlct    2-step
PlsSPELPhPrgEx     Normal
PlsSPELPrg         
PlsSPELPrgTxt      \
;\n\
;\n\
;4-Pulse ELDOR\n\
;\n\
;       M. Pannier, S. Veit, A. Godt, G. Jeschke, and H. W. Spiess\n\
;               J. Magn Res. 2000, v. 142, p. 331\n\
;\n\
;\n\
begin defs              ; dimension of data array [sx,sy]\n\
 dim s[512]             ;       for 2 pulse Field Sweep set-up\n\
 dim1 s[512]            ;       for 2 pulse Field Sweep\n\
 dim2 s[64]             ;       for Optimization of ELDOR pulse\n\
 dim3 s[512]            ;       for 4-pulse ELDOR setup\n\
 dim4 s[330]            ;       for 1D 4-Pulse ELDOR\n\
 dim5 s[500, 32]        ;       for 2D 4-Pulse ELDOR\n\
 dim6 s[754, 32]        ;       for 2D 4-Pulse ELDOR d1 averaging for D2O\n\
 dim7 s[500, 32]        ;       for 2D 4-Pulse ELDOR d1 averaging for H2O\n\
end defs\n\
;\n\
;\n\
begin lists "none"\n\
 ph1 +x\n\
 ph2 +x\n\
 ph3 +x                 ; not used\n\
 ph4 +x\n\
 asg1 +a\n\
 bsg1 +b\n\
end lists\n\
\n\
begin lists1 "2-step"\n\
 ph1 +x -x\n\
 ph2 +x\n\
 ph3 +x                 ; not used\n\
 ph4 +x\n\
 asg1 +a -a\n\
 bsg1 +b -b\n\
end lists1\n\
\n\
begin lists2 "8-step"\n\
; This phase cycle removes all echoes except the refocused echo\n\
 ph1 +x -x\n\
 ph2 +x +x -x -x +y +y -y -y\n\
 ph3 +x                 ; not used\n\
 ph4 +x\n\
 asg1 +a -a +a -a -a +a -a +a\n\
 bsg1 +b -b +b -b -b +b -b +b\n\
end lists2\n\
\n\
\n\
;\n\
; Standing Echo for 2-Pulse Set-Up\n\
;\n\
begin exp "2P ESE Setup" [TRANS QUAD]\n\
        ; QUAD detection of Transient\n\
\n\
    shot i=1 to a               ; accumulation loop\n\
      d9                        ; DAF (Delay After Flash)\n\
      p3 [*]                    ; Inversion Pulse (not used)\n\
      d3                        ; T\n\
      p0 [ph1]                  ;\n\
      d1                        ; tau1\n\
      p1 [ph2]                  ;\n\
      d1                        ; tau1\n\
      d0                        ; constant acquisition delay\n\
      dig [sg1]                 ; acquisition\n\
    next i                      ; end of accumulation loop\n\
\n\
end exp\n\
\n\
;\n\
; 2-pulse Field Sweep\n\
;\n\
begin exp1 "2P ESE Field Sweep" [INTG QUAD lists lists1]\n\
        ; QUAD detection of Transient\n\
\n\
bsweep x=1 to sx\n\
    shot i=1 to h               ; accumulation loop\n\
      d9                        ; DAF (Delay After Flash)\n\
      p3 [*]                    ; Inversion Pulse (not used)\n\
      d3                        ; T\n\
      p0 [ph1]                  ;\n\
      d1                        ; tau1\n\
      p1 [ph2]                  ;\n\
      d1                        ; tau1\n\
      d0                        ; constant acquisition delay\n\
      acq [sg1]                 ; acquisition\n\
    next i                      ; end of accumulation loop\n\
next x\n\
end exp1\n\
\n\
;\n\
; ELDOR Optimization\n\
;\n\
begin exp2 "ELDOR pulse Optimization" [SPT QUAD lists lists1]\n\
        ; QUAD detection of Transient\n\
p4=p3                           ; Store initial inversion pulse length\n\
\n\
 for k=1 to n                   ; averaging loop\n\
  totscans(n)                   ; output number of scans to do\n\
\n\
  p3=p4                         ; starting inversion pulse length and assignment of x-axis\n\
  dx=p4\n\
  \n\
  sweep x=1 to sx\n\
    shot i=1 to h               ; accumulation loop\n\
      d9                        ; DAF (Delay After Flash)\n\
      p3 [ELDOR]                ; Inversion Pulse (ELDOR = Observe Frequency)\n\
      dx                        ; increment x-axis\n\
      d3                        ; delay after inversion pulse\n\
      p0 [ph1]                  ; 2nd pulse\n\
      d1                        ; tau1\n\
      p1 [ph2]                  ; 3rd pulse\n\
      d1                        ; tau1\n\
      d0                        ; constant acquisition delay\n\
      acq [sg1]                 ; acquisition\n\
    next i                      ; end of accumulation loop\n\
    p3=p3+d30\n\
    dx=dx+d30\n\
  next x\n\
\n\
  scansdone(k)\n\
 next k\n\
\n\
end exp2\n\
\n\
;\n\
; 4-Pulse ELDOR Set-Up\n\
;\n\
begin exp3 "4P ELDOR Setup" [TRANS QUAD]\n\
        ; QUAD detection of Transient\n\
\n\
d7=d1+d2                        ; calculate separation\n\
d7=d7-d3\n\
\n\
    shot i=1 to a               ; accumulation loop\n\
      d9                        ; DAF (Delay After Flash)\n\
      p0 [ph1]                  ; 1st pulse\n\
      d1                        ; tau1\n\
      p1 [ph2]                  ; 2nd pulse\n\
      d3                        ; T + d30 * sx < d1 + d2\n\
      p3 [ELDOR]                ; ELDOR pulse\n\
      d7                        ; tau2+tau1-T\n\
      p1 [ph4]                  ; 3rd pulse\n\
      d2                        ; tau2\n\
      d0                        ; constant acquisition delay\n\
      dig [sg1]                 ; acquisition\n\
    next i                      ; end of accumulation loop\n\
\n\
end exp3\n\
\n\
;\n\
; 4-Pulse ELDOR\n\
;\n\
begin exp4 "1D-4P ELDOR"  [INTG QUAD]\n\
        ; QUAD detection with Integration\n\
\n\
d8=d3                           ; store T\n\
 for k=1 to n                   ; averaging loop\n\
   totscans (n)                 ; output of total number of scans\n\
\n\
   d3=d8                        ; recall T for averaging loop\n\
   d7=d1+d2                     ; calculate separation\n\
   d7=d7-d3\n\
   dx=d3-d1                     ; and define x-axis\n\
\n\
   sweep x=1 to sx              ; T sweep loop\n\
    shot i=1 to h               ; accumulation loop\n\
      d9                        ; DAF (Delay After Flash)\n\
      p0 [ph1]                  ; 1st pulse\n\
      d1                        ; tau1\n\
      p1 [ph2]                  ; 2nd pulse\n\
      d3                        ; T + d30 * sx-1 < d1 + d2\n\
      p3 [ELDOR]                ; ELDOR pulse\n\
      d7                        ; tau1+tau2-T\n\
      p1 [ph4]                  ; 3rd pulse\n\
      d2                        ; tau2\n\
      d0                        ; constant acquisition delay\n\
      acq [sg1]\n\
    next i                      ; end of accumulation loop\n\
\n\
    d3=d3+d30                   ; increment T\n\
    d7=d7-d30                   ; decrement time between pulse 3 and 4\n\
    dx=dx+d30                   ; increment x axis\n\
   next x                       ; end of T sweep loop\n\
\n\
  scansdone(k)                  ; output number of scans done\n\
 next k\n\
\n\
end exp4\n\
\n\
;\n\
; 4-Pulse ELDOR\n\
;\n\
begin exp5 "2D-4P ELDOR"  [INTG QUAD]\n\
        ; QUAD detection with Integration\n\
        ; Y axis is a dummy axis so that averaging occurs in slices\n\
        ; If the cooling water gets to warm and the bridge shuts down\n\
        ; Simply extracted the good slices and perform a sum projection\n\
\n\
d8=d3                           ; store T\n\
for y = 1 to sy\n\
 for k=1 to n                   ; averaging loop\n\
   totscans (n)                 ; output of total number of scans\n\
\n\
   d3=d8                        ; recall T for averaging loop\n\
   d7=d1+d2                     ; calculate separation\n\
   d7=d7-d3\n\
   dx=d3-d1                     ; and define x-axis\n\
\n\
   sweep x=1 to sx              ; T sweep loop\n\
    shot i=1 to h               ; accumulation loop\n\
      d9                        ; DAF (Delay After Flash)\n\
      p0 [ph1]                  ; 1st pulse\n\
      d1                        ; tau1\n\
      p1 [ph2]                  ; 2nd pulse\n\
      d3                        ; T + d30 * sx-1 < d1 + d2\n\
      p3 [ELDOR]                ; ELDOR pulse\n\
      d7                        ; tau1+tau2-T\n\
      p1 [ph4]                  ; 3rd pulse\n\
      d2                        ; tau2\n\
      d0                        ; constant acquisition delay\n\
      acq [sg1]\n\
    next i                      ; end of accumulation loop\n\
\n\
    d3=d3+d30                   ; increment T\n\
    d7=d7-d30                   ; decrement time between pulse 3 and 4\n\
    dx=dx+d30                   ; increment x axis\n\
   next x                       ; end of T sweep loop\n\
\n\
  scansdone(k)                  ; output number of scans done\n\
 next k\n\
dy = dy + d31\n\
next y\n\
\n\
end exp5\n\
\n\
;\n\
; 4-Pulse ELDOR\n\
;\n\
begin exp6 "2D-4P ELDOR d1avg D2O"  [INTG QUAD]\n\
; QUAD detection with Integration\n\
; Y axis is a dummy axis so that averaging occurs in slices\n\
; If the cooling water gets to warm and the bridge shuts down\n\
; Simply extracted the good slices and perform a sum projection\n\
\n\
d8=d3                           ; store T\n\
d11=d1               ;store d1\n\
for y = 1 to sy\n\
\n\
for k=1 to 8                    ; d1 averaging loop\n\
totscans (n)                    ; output of total number of scans\n\
\n\
d12=k*16\n\
d12=d12-16\n\
d3=d8+d12                       ; recall T for averaging loop\n\
d1=d11+d12          ; recall d1\n\
d7=d1+d2                        ; calculate separation\n\
d7=d7-d3\n\
dx=d3-d1                        ; and define x-axis\n\
\n\
sweep x=1 to sx         ; T sweep loop\n\
shot i=1 to h           ; accumulation loop\n\
d9                      ; DAF (Delay After Flash)\n\
p0 [ph1]                        ; 1st pulse\n\
d1                      ; tau1\n\
p1 [ph2]                        ; 2nd pulse\n\
d3                      ; T + d30 * sx-1 < d1 + d2\n\
p3 [ELDOR]              ; ELDOR pulse\n\
d7                      ; tau1+tau2-T\n\
p1 [ph4]                        ; 3rd pulse\n\
d2                      ; tau2\n\
d0                      ; constant acquisition delay\n\
acq [sg1]\n\
next i                  ; end of accumulation loop\n\
\n\
d3=d3+d30                       ; increment T\n\
d7=d7-d30                       ; decrement time between pulse 3 and 4\n\
dx=dx+d30                       ; increment x axis\n\
next x                  ; end of T sweep loop\n\
\n\
scansdone(k)                    ; output number of scans done\n\
next k\n\
dy = dy + d31\n\
next y\n\
\n\
end exp6\n\
\n\
\n\
\n\
\n\
;\n\
; 4-Pulse ELDOR\n\
;\n\
begin exp7 "2D-4P ELDOR d1avg H2O"  [INTG QUAD]\n\
; QUAD detection with Integration\n\
; Y axis is a dummy axis so that averaging occurs in slices\n\
; If the cooling water gets to warm and the bridge shuts down\n\
; Simply extracted the good slices and perform a sum projection\n\
;Now with added d1 averaging\n\
\n\
d8=d3                           ; store T\n\
d11=d1               ;store d1\n\
for y = 1 to sy\n\
\n\
for k=1 to 8                    ; d1 averaging loop\n\
totscans (n)                    ; output of total number of scans\n\
\n\
d12=k*8\n\
d12=d12-8\n\
d3=d8+d12                       ; recall T for averaging loop\n\
d1=d11+d12          ; recall d1\n\
d7=d1+d2                        ; calculate separation\n\
d7=d7-d3\n\
dx=d3-d1                        ; and define x-axis\n\
\n\
sweep x=1 to sx         ; T sweep loop\n\
shot i=1 to h           ; accumulation loop\n\
d9                      ; DAF (Delay After Flash)\n\
p0 [ph1]                        ; 1st pulse\n\
d1                      ; tau1\n\
p1 [ph2]                        ; 2nd pulse\n\
d3                      ; T + d30 * sx-1 < d1 + d2\n\
p3 [ELDOR]              ; ELDOR pulse\n\
d7                      ; tau1+tau2-T\n\
p1 [ph4]                        ; 3rd pulse\n\
d2                      ; tau2\n\
d0                      ; constant acquisition delay\n\
acq [sg1]\n\
next i                  ; end of accumulation loop\n\
\n\
d3=d3+d30                       ; increment T\n\
d7=d7-d30                       ; decrement time between pulse 3 and 4\n\
dx=dx+d30                       ; increment x axis\n\
next x                  ; end of T sweep loop\n\
\n\
scansdone(k)                    ; output number of scans done\n\
next k\n\
dy = dy + d31\n\
next y\n\
\n\
end exp7\n\
\n\
\n\
\n
Psd1               {2;33,4;0} 0,0,484,660,6836,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,40,40,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd10              {2;33,4;0} 0,484,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd11              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd12              {2;33,4;0} 0,484,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd13              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd14              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd15              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd16              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd17              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd18              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd19              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd2               {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd20              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd21              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd22              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd23              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd24              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd25              {2;33,4;0} 0,484,768,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,44,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd26              {2;33,4;0} 0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6852,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd27              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd28              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd29              {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd3               {2;33,4;0} 0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6852,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd4               {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6808,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd5               {2;33,4;0} 0,4,12928,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,67,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd6               {2;33,4;0} 0,484,660,6836,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd7               {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd8               {2;33,4;0} 0,484,484,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd9               {2;33,4;0} 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
QuadDetect         True
RF1Prg             {2;33,3;0.000[MHz]} 0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
RF2Prg             {2;33,3;0.000[MHz]} 0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
ReplaceMode        Off
ShotRepTime        1999.20 us
ShotsPLoop         1000
SmoothAllowed      False
SmoothPoints       1
SptProgress        33 %
StochMode          False
SweepsPExp         8
TriggerTimeOut     1668 s
XAxisQuant         Time
XSpecRes           754
YAxisQuant         Magnetic Field
YSpecRes           32

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         3
NbScansDone        3
NbScansToDo        8
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       1

.DVC     routeTrg, 1.0


.DVC     sigChanSmall, 1.0


.DVC     transRec, 1.0


*
************************************************************
