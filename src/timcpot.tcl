#!/usr/bin/env wish8.6
#!/usr/bin/env wish
#
# Usage: timcpot.tcl data.dir(/evtid/evtid.staid.channel.SAC) out.dir(/evtid.pickfile) sta.loc.txt(name stla stlo stel)
#
# Author: zhangfengxue336@163.com
#
# sac, ttimes, evtlocate, hilbertsac, ps2raster 


package require Thread
package require Ttrace

proc main {argc argv} {
    global bp1 bp2 evtinter recstalim
    global stwin ltwin stltsw stltlim
    global kswin kssw klim slim
    global aicwin aicsw aiclim
    global shiftmark shift1 shift2
    global rawdatsw polysw insw
    global model
    global phasetype
    global datdir outdir stafile
    global sc
    global stafx stafy
    global wavefx wavefy wavezoomfx wavezoomfy
    global z1 z2
    global dxx dyy dzz
    global evtdephand evtdepmax
    global pickerrmax stalim
    global channelsw rotatesw
    global vibnum vibdist
    global colorband colorn
    global waveformsw envelopesw instfreqsw
    global readthreadnum calculatethreadnum
    global singlewaveheight

    if {$argc<3} {
	#puts " Version: R1.0  Author: zhangfengxue336@163.com";
	puts " Usage: timcpot.tcl data.dir(/evtid/evtid.staid.channel.SAC) out.dir(/evtid.pickfile) sta.loc.txt(name stla stlo stel)"
	exit;
    }
    #set colorn 5
    #set colorband(0) "#A60000"; set colorband(1) "#E10000"; set colorband(2) "#FF4B00";
    #set colorband(3) "#FF8400"; set colorband(4) "#FFC600"; set colorband(5) "#FFFA98";

    set colorn 5
    set colorband(0) "#FFA500"; set colorband(1) "#FF7F00"; set colorband(2) "#FF7256";
    set colorband(3) "#FF6347"; set colorband(4) "#FF4500"; set colorband(5) "#FF0000";


    # default values
    set vibnum 5;      set vibdist 100;
    set bp1 5;         set bp2 5;
    set  stwin 2.0;    set ltwin 2.0;    set stltsw 1;    set stltlim 100.0;
    set  kswin 2.0;                      set   kssw 0;    set    klim   2.0;  set slim 2.0;
    set aicwin 2.0;                      set  aicsw 0;    set  aiclim   2.0;
    set shiftmark o;   set shift1 -10.0; set shift2 180.0;
    set evtinter 10.0; set recstalim 10;

    set z1 0.0;        set z2 20.0;
    set dxx 0.5; set dyy 0.5; set dzz 1.0
    set evtdepmax 40

    set readthreadnum 3;  set calculatethreadnum 3

    frame .menu -borderwidth 1.5 -relief groove
    frame .menuline1 -borderwidth 1.0 -relief groove
    frame .menuline2 -borderwidth 1.0 -relief groove
    frame .menuline3 -borderwidth 1.0 -relief groove
    frame .leftcolumn -borderwidth 1.5 -relief groove
    frame .righttop -borderwidth 1.5 -relief groove
    frame .rightbottom -borderwidth 1.5 -relief groove
    ############################ 
    frame .stltframe -borderwidth 1.5 -relief raised
    label .lstwin -text "stwin:"
    text  .tstwin -width 4 -height 1.1 -bg white -fg blue
          .tstwin insert 1.0 $stwin
    label .lltwin -text "ltwin:"
    text  .tltwin -width 4 -height 1.1 -bg white -fg blue
          .tltwin insert 1.0 $ltwin
    label .lstltsw -text "stltsw:"
    text  .tstltsw -width 4 -height 1.1 -bg white -fg blue
          .tstltsw insert 1.0 $stltsw
    label .lstltlim -text "stltlim:"
    text  .tstltlim -width 4 -height 1.1 -bg white -fg blue
          .tstltlim insert 1.0 $stltlim
    pack .lstwin .tstwin .lltwin .tltwin .lstltsw .tstltsw .lstltlim .tstltlim -side left -in .stltframe
    pack .stltframe -side left -in .menuline1

    frame .ksframe -borderwidth 1.5 -relief raised
    label .lkswin -text "kswin:"
    text  .tkswin -width 4 -height 1.1 -bg white -fg blue
          .tkswin insert 1.0 $kswin
    label .lkssw -text "kssw:"
    text  .tkssw -width 4 -height 1.1 -bg white -fg blue
          .tkssw insert 1.0 $kssw
    label .lklim -text "klim:"
    text  .tklim -width 4 -height 1.1 -bg white -fg blue
          .tklim insert 1.0 $klim
    label .lslim -text "slim:"
    text  .tslim -width 4 -height 1.1 -bg white -fg blue
          .tslim insert 1.0 $slim
    pack .lkswin .tkswin .lkssw .tkssw .lklim .tklim .lslim .tslim -side left -in .ksframe
    pack .ksframe -side left -in .menuline1

    frame .aicframe -borderwidth 1.5 -relief raised
    label .laicwin -text "aicwin:"
    text  .taicwin -width 4 -height 1.1 -bg white -fg blue
          .taicwin insert 1.0 $aicwin
    label .laicsw -text "aicsw:"
    text  .taicsw -width 4 -height 1.1 -bg white -fg blue
          .taicsw insert 1.0 $aicsw
    label .laiclim -text "aiclim:"
    text  .taiclim -width 4 -height 1.1 -bg white -fg blue
          .taiclim insert 1.0 $aiclim
    pack .laicwin .taicwin .laicsw .taicsw .laiclim .taiclim -side left -in .aicframe
    pack .aicframe -side left -in .menuline1

    frame .evtparaframe -borderwidth 1.5 -relief raised
    label .levtinter -text "evtinter:"
    text  .tevtinter -width 4 -height 1.2 -bg white -fg red
          .tevtinter insert 1.0 $evtinter
    label .lrecstalim -text "stalim:"
    text  .trecstalim -width 4 -height 1.2 -bg white -fg red
          .trecstalim insert 1.0 $recstalim
    pack .levtinter .tevtinter .lrecstalim .trecstalim -side left -in .evtparaframe
    pack .evtparaframe -side left -in .menuline1

    frame .framethread -borderwidth 1.5 -relief groove
    label .lthreadread -text "Parallel.Read:"
    text  .tthreadread -width 3 -height 1.2 -bg white -fg red
          .tthreadread insert 1.0 $readthreadnum
    label .lthreadcalc -text "Parallel.Calc:"
    text  .tthreadcalc -width 3 -height 1.2 -bg white -fg red
          .tthreadcalc insert 1.0 $calculatethreadnum
    pack  .lthreadread .tthreadread .lthreadcalc .tthreadcalc -side left -in .framethread
    pack  .framethread -side left -in .menuline2

    frame .vibframe -borderwidth 1.5 -relief groove
    label .lvibnum -text "Circle.Amt:" -width 8
    text  .tvibnum -width 3 -height 1.2 -bg white -fg red
          .tvibnum insert 1.0 $vibnum
    label .lvibdist -text "Dist.Step:" -width 7
    text  .tvibdist -width 5 -height 1.2 -bg white -fg red
          .tvibdist insert 1.0 $vibdist
    #label .bp1 -text $bp1 -width 4
    #label .bp  -text "-"  -width 1
    #label .bp2 -text $bp2 -width 4
    pack .lvibnum .tvibnum .lvibdist .tvibdist -side left -in .vibframe
    pack .vibframe -side left -in .menuline2

    frame .evtwaveframe -borderwidth 1.5 -relief groove
    label .lshiftmark -text "shift mark: $shiftmark" -relief raised -bg orange -width 13
    label .lshift1 -text "shift1:"
    text  .tshift1 -width 5 -height 1.2 -bg white -fg red
          .tshift1 insert 1.0 $shift1
    label .lshift2 -text "shift2:"
    text  .tshift2 -width 5 -height 1.2 -bg white -fg red
          .tshift2 insert 1.0 $shift2
    pack .lshiftmark .lshift1 .tshift1 .lshift2 .tshift2 -side left -in .evtwaveframe
    pack .evtwaveframe -side left -in .menuline2

    menu .menushiftmark -tearoff 0
         .menushiftmark add command -label "set to b" -command {.lshiftmark configure -text "shift mark: b"; set shiftmark b}
         .menushiftmark add command -label "set to e" -command {.lshiftmark configure -text "shift mark: e"; set shiftmark e}
         .menushiftmark add command -label "set to o" -command {.lshiftmark configure -text "shift mark: o"; set shiftmark o}
         .menushiftmark add command -label "set to a" -command {.lshiftmark configure -text "shift mark: a"; set shiftmark a}
         .menushiftmark add command -label "set to f" -command {.lshiftmark configure -text "shift mark: f"; set shiftmark f}
         .menushiftmark add command -label "set to t0" -command {.lshiftmark configure -text "shift mark: t0"; set shiftmark t0}
         .menushiftmark add command -label "set to t1" -command {.lshiftmark configure -text "shift mark: t1"; set shiftmark t1}
         .menushiftmark add command -label "set to t2" -command {.lshiftmark configure -text "shift mark: t2"; set shiftmark t2}
         .menushiftmark add command -label "set to t3" -command {.lshiftmark configure -text "shift mark: t3"; set shiftmark t3}
         .menushiftmark add command -label "set to t4" -command {.lshiftmark configure -text "shift mark: t4"; set shiftmark t4}
         .menushiftmark add command -label "set to t5" -command {.lshiftmark configure -text "shift mark: t5"; set shiftmark t5}
         .menushiftmark add command -label "set to t6" -command {.lshiftmark configure -text "shift mark: t6"; set shiftmark t6}
         .menushiftmark add command -label "set to t7" -command {.lshiftmark configure -text "shift mark: t7"; set shiftmark t7}
         .menushiftmark add command -label "set to t8" -command {.lshiftmark configure -text "shift mark: t8"; set shiftmark t8}
         .menushiftmark add command -label "set to t9" -command {.lshiftmark configure -text "shift mark: t9"; set shiftmark t9}
    #frame .framezoom -borderwidth 1.5 -relief groove
    #button .scUp -text "Up" -command {scaleTrace 2} -height 1 -width 1 -bg orange
    #button .scDn -text "Dn" -command {scaleTrace 0.5} -height 1 -width 1 -bg orange
    #pack .scDn .scUp -side left -in .framezoom
    #pack .framezoom -side left -in .menuline2

    #frame .framewave -borderwidth 1.5 -relief groove -bg orange
    #set rawdatsw 0
    #radiobutton .rawdat0 -variable rawdatsw -text "Rn" -value 0 -height 1 -width 2
    #radiobutton .rawdat1 -variable rawdatsw -text "Ry" -value 1 -height 1 -width 2
    #set polysw 0;set insw 1;
    #radiobutton .poly0 -variable polysw -command {polymove} -text "Sn" -value 0 -height 1 -width 2
    #radiobutton .poly1 -variable polysw -command {polymove} -text "Sy" -value 1 -height 1 -width 2
    #pack .rawdat0 .rawdat1 .poly0 .poly1 -side left -in .framewave
    #pack .framewave -side left -in .menuline2

    set model "iasp91"
    frame .framemodel -borderwidth 1.5 -relief groove
    radiobutton .ak135  -variable model -text "ak135"  -value "ak135"
    radiobutton .iasp91 -variable model -text "iasp91" -value "iasp91"
    pack .ak135 .iasp91 -side left -in .framemodel
    pack .framemodel -side left -in .menuline2



    label .lshowallpicks -text "check picks" -relief raised -bg orange -width 11
    pack  .lshowallpicks -side left -in .menuline2
    label .labelzoomcurve -justify center -text "check curve" -relief raised -width 11 -bg orange
    pack  .labelzoomcurve -side left -in .menuline2
    ############################ 
    frame .framesta -borderwidth 1.5 -relief groove
    label .lsw -text "sw:"
    text  .tsw  -width 4 -height 1.2 -bg lightgray -fg red
    set sc 3
    label .lsc -text "Sta.Circle:"
    text  .tsc  -width 2 -height 1.2 -bg white -fg red
	  .tsc insert 1.0 $sc
    pack .lsc .tsc -in .framesta -side left
    pack .framesta -side left -in .menuline3

    label .llocationpara -text "locate para:" -relief raised -bg orange -width 10
    pack .llocationpara -side left -in .menuline3
    menu .menulocationpara -tearoff 0
	 .menulocationpara add command -label "0.50 x 0.50 x 5.0" -command {setlocationpara 0.50 0.50 5.0}
	 .menulocationpara add command -label "0.50 x 0.50 x 2.0" -command {setlocationpara 0.50 0.50 2.0}
	 .menulocationpara add command -label "0.50 x 0.50 x 1.0" -command {setlocationpara 0.50 0.50 1.0}
	 .menulocationpara add command -label "0.20 x 0.20 x 1.0" -command {setlocationpara 0.20 0.20 1.0}
	 #.menulocationpara add command -label "0.10 x 0.10 x 1.0" -command {setlocationpara 0.10 0.10 1.0}
	 #.menulocationpara add command -label "0.05 x 0.05 x 0.5" -command {setlocationpara 0.05 0.05 0.5}
	 #.menulocationpara add command -label "0.02 x 0.02 x 0.2" -command {setlocationpara 0.02 0.02 0.2}
	 #.menulocationpara add command -label "0.01 x 0.01 x 0.1" -command {setlocationpara 0.01 0.01 0.1}
    frame .framexy -borderwidth 1.5 -relief groove
    text .tsx1 -width 8 -height 1.2 -bg lightgray -fg red
    text .tdxx -width 5 -height 1.2 -bg white -fg red
	 .tdxx insert 1.0 $dxx
    text .tsx2 -width 8 -height 1.2 -bg lightgray -fg red

    text .tsy1 -width 8 -height 1.2 -bg lightgray -fg red
    text .tdyy -width 5 -height 1.2 -bg white -fg red
	 .tdyy insert 1.0 $dyy
    text .tsy2 -width 8 -height 1.2 -bg lightgray -fg red

    text .tsz1 -width 4 -height 1.2 -bg white -fg red
	 .tsz1 insert 1.0 $z1
    text .tdzz -width 4 -height 1.2 -bg white -fg red
	 .tdzz insert 1.0 $dzz
    text .tsz2 -width 4 -height 1.2 -bg white -fg red
	 .tsz2 insert 1.0 $z2


    pack .tsx1 .tdxx .tsx2 .tsy1 .tdyy .tsy2 .tsz1 .tdzz .tsz2 -in .framexy -side left
    pack .framexy -side left -in .menuline3

    frame .frameevlocate -borderwidth 1 -relief groove
    set pickerrmax 40;
    label .lPKerr -text "ErrLimit:"
    text .tPKerr -width 3 -height 1.2 -bg white -fg red
         .tPKerr insert 1.0 $pickerrmax
    set stalim 4;
    label .lstalim -text "StaLimit:";
    text .tstalim -width 2 -height 1.2 -bg white -fg red
         .tstalim insert 1.0 $stalim
    pack .lPKerr .tPKerr .lstalim .tstalim -side left -in .frameevlocate
    pack .frameevlocate -side left -in .menuline3

    #menu .menuwaveheight -tearoff 0
	# .menuwaveheight add command -label "150" -command {set singlewaveheight 150; .twaveH delete 1.0 end; .twaveH insert 1.0 $singlewaveheight;}
	# .menuwaveheight add command -label "100" -command {set singlewaveheight 100; .twaveH delete 1.0 end; .twaveH insert 1.0 $singlewaveheight;}
	# .menuwaveheight add command -label "50"  -command {set singlewaveheight 50;  .twaveH delete 1.0 end; .twaveH insert 1.0 $singlewaveheight;}
    #set singlewaveheight 100;
    #frame .framewaveH -borderwidth 1 -relief groove
    #label .labelwaveH -justify center -text "Wave H:" -relief raised -width 7 -bg orange
    #text .twaveH -width 3 -height 1.2 -bg white -fg red
    #     .twaveH insert 1.0 $singlewaveheight;
    #pack .labelwaveH .twaveH -side left -in .framewaveH
    #pack .framewaveH -side left -in .menuline3

    #### Information Display
    label .information -text "Display Screen" -bg white -fg red -width 70
    pack .information -side left -in .menuline3
    #### Information Display
    ############################
    ############################ 
    pack  .menuline2 .menuline3 -side top -in .menu -fill x

    pack .menu -side top -fill x;    # main control

    frame .leftcolumnsubtop -borderwidth 1.0 -relief groove
    set stafx 315; set stafy 320
    frame .framestafig -borderwidth 1.0 -relief groove

    scrollbar .stafigurey -command ".stafigure yview"
    pack .stafigurey -side right -in .framestafig -fill y

    frame .frameevdp -borderwidth 1.0 -relief groove
    label .labelevdpmax -text "MXDP" -fg red
    text .textevdp -width 4 -height 1.2 -bg white -fg red
	 .textevdp insert 1.0 $evtdepmax
    label .labelevdp -text "EVDP"
    scale .scaleevtdep -from 0 -to $evtdepmax -length 280 -width 12 -variable evtdephand -command {set evtdephand}
    pack .labelevdpmax .textevdp .labelevdp .scaleevtdep -in .frameevdp -side top
    pack .frameevdp -side right -in .framestafig

    frame .framelendge -borderwidth 1.0 -relief groove
    canvas .lendge -background white -width 114 -height 20
    .lendge create oval 2 7 8 13 -fill red
    .lendge create text 10 10 -text "Picked" -anchor w
    .lendge create oval 58 7 64 13 -fill green
    .lendge create text 66 10 -text "Unpick" -anchor w
    label .labelfree  -width 5 -relief raise -fg black   -text "+:free"
    label .labelsac   -width 5 -relief raise -fg magenta -text "+:file"
    label .labelrloc  -width 7 -relief raise -fg blue    -text "+:locate"
    label .labelfinal -width 6 -relief raise -fg darkorange  -text "+:final"
    pack .lendge .labelfree .labelsac .labelrloc .labelfinal -in .framelendge -side left
    pack .framelendge -in .framestafig -side bottom -fill x

    canvas .stafigure -background white -width $stafx -height $stafy -closeenough 1 \
                      -xscrollcommand ".stafigurex set" -yscrollcommand ".stafigurey set"
    scrollbar .stafigurex -command ".stafigure xview" -orient horizontal

    pack .stafigurex .stafigure -side bottom -in .framestafig -fill x

    pack .framestafig -in .leftcolumnsubtop
    pack .leftcolumnsubtop -in .leftcolumn -side top

    ########################
    frame .leftcolumnsubmiddle -borderwidth 1.0 -relief groove
    frame .framepara -borderwidth 0.0 -relief groove
    frame .framepara1 -borderwidth 1.0 -relief groove

    frame .framefilter -borderwidth 1.0 -relief groove
    set tt "Notes\n";
    set tt "${tt}bp1=bp2: unfiltered"
    label .lfilter -justify left -fg red -text $tt
    pack .lfilter -in .framefilter -side top

    frame .frameband -borderwidth 1.0 -relief groove
    label .lbp1 -text "bp1:" -width 3
    text  .tbp1 -width 4 -height 1.2 -bg white -fg red
          .tbp1 insert 1.0 $bp1
    label .lbp2 -text "bp2:" -width 3
    text  .tbp2 -width 4 -height 1.2 -bg white -fg red
          .tbp2 insert 1.0 $bp2
    pack .lbp1 .tbp1 .lbp2 .tbp2 -in .frameband -side left
    pack .frameband -in .framefilter -side left
    pack .framefilter -in .framepara1 -side left

    frame .frameprecmd -borderwidth 1.0 -relief groove
    set tt "SAC commands are here:\n"
    set tt "${tt}one command one line"
    label .lpre -justify left -fg black -text $tt
    pack .lpre -in .frameprecmd -side top 
    set preline "rmean\nrtrend\ndecimate 2";
    text .tpre -bg white -fg black -width 24 -height 4 \
	       -yscrollcommand ".tprey set"
         .tpre insert 1.0 $preline
    scrollbar .tprey -command ".tpre yview"
    pack .tpre .tprey -in .frameprecmd -side left  -fill y

    pack .frameprecmd -in .framepara1 -side top
    pack .framepara1 -in .framepara -side top -fill x

    ########################
    frame .framepara2 -borderwidth 1.0 -relief groove
    frame .framechannel  -borderwidth 1.0 -relief groove
    frame .framechannel1 -borderwidth 1.0 -relief groove
    frame .framechannel2 -borderwidth 1.0 -relief groove
    set tt "Channel: BH{ZEN}.SAC"
    label .ltext -justify left -fg black -text $tt
    set channelsw "VH"
    radiobutton .cvert      -variable channelsw -text "V"    -value "V"
    radiobutton .chori      -variable channelsw -text "H"    -value "H"
    radiobutton .cverthori  -variable channelsw -text "V&H"  -value "VH"
    pack .ltext .cvert .chori .cverthori -in .framechannel1 -side left
    set tt "Rotate Horizontal channel"
    label .lrotate -justify left -fg black -text $tt
    set rotatesw "YES"
    radiobutton .rotateyes -variable rotatesw -text "YES" -value "YES"
    radiobutton .rotateno  -variable rotatesw -text "NO"  -value "NO"
    pack .lrotate .rotateyes .rotateno -in .framechannel2 -side left

    pack .framechannel1 .framechannel2 -in .framechannel -side top -fill x
    #pack .framechannel -in .leftcolumn -side top -fill x
    pack .framechannel -in .framepara2 -side top -fill x
    pack .framepara2 -in .framepara -side top -fill x
    pack .framepara -in .leftcolumnsubmiddle -fill x
    pack .leftcolumnsubmiddle -in .leftcolumn -side top -fill x
    ########################
    frame .frameevtid -borderwidth 1.0 -relief groove
    label .evtnote -justify center -fg black -text "Event Directory List"
    listbox .evtid -width 28 -exportselection 0 -selectforeground red \
                     -xscrollcommand ".evtidx set" -yscrollcommand ".evtidy set"
    scrollbar .evtidy -command ".evtid yview"
    scrollbar .evtidx -command ".evtid xview" -orient horizontal
    pack .evtnote .evtidx -side top -in .frameevtid -fill x
    pack .evtidy .evtid -side right -in .frameevtid -fill y
    pack .frameevtid -in .leftcolumn -fill y -side left

    #frame .frameevtlist -borderwidth 1.0 -relief groove -bg orange
    #listbox .evtlist -width 15 -exportselection 0 -selectforeground red \
    #                -xscrollcommand ".evtlistx set" -yscrollcommand ".evtlisty set"
    #scrollbar .evtlisty -command ".evtlist yview"
    #scrollbar .evtlistx -command ".evtlist xview" -orient horizontal
    #pack .evtlistx -in .frameevtlist -side top -fill x
    #pack .evtlisty .evtlist -in .frameevtlist -side right  -fill y
    #pack .frameevtlist -in .leftcolumn -fill y -side left

    frame .framefilelist -borderwidth 1.0 -relief groove
    label .fileindicator -justify center -fg black -text "Output File List"
    pack .fileindicator -in .framefilelist -side top -fill x

    frame .framechanid -borderwidth 1.0 -relief groove
    label .colorbk1 -justify left -bg lightgray -width 3 -text "No"  -relief raised
    label .colorbk2 -justify left -bg yellow -width 3 -text "Yes" -relief raised
    pack  .colorbk1 .colorbk2 -in .framechanid -side left -fill x
    label .chanZ -justify center -fg black -text "ZZ"
    label .chanE -justify center -fg black -text "ET"
    label .chanN -justify center -fg black -text "NR"
    pack .chanN .chanE .chanZ -in .framechanid -side right -fill x
    pack .framechanid -in .framefilelist -side top -fill x

    frame .framename -borderwidth 1.0 -relief groove
    label .filePKTRLOC  -justify left -fg black -text "PKTRLOC"
    label .filePKTFINAL -justify left -fg black -text "PKTFINAL"
    label .fileEVTRLOC  -justify left -fg black -text "EVTRLOC"
    label .fileEVTFINAL -justify left -fg black -text "EVTFINAL"
    label .fileEVTSACF  -justify left -fg black -text "EVTSACF"
    label .fileSTATION  -justify left -fg black -text "STATION"
    label .fileWAVE     -justify left -fg black -text "WAVE"
    label .filePKTCMP   -justify left -fg black -text "PKTCMP"
    pack .filePKTRLOC .filePKTFINAL .fileEVTRLOC .fileEVTFINAL .fileEVTSACF .fileSTATION .fileWAVE .filePKTCMP -in .framename -side top
    pack .framename -in .framefilelist -side left -fill y

    set labelbg lightgray
    set labelwd 2
    frame .frameindicatorz -borderwidth 1.0 -relief groove
    label .labelBHZPKTRLOC  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHZPKTFINAL -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHZEVTRLOC  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHZEVTFINAL -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHZEVTSACF  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHZSTATION  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHZWAVE     -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHZPKTCMP   -justify left -bg $labelbg -width $labelwd -relief raised
    pack .labelBHZPKTRLOC .labelBHZPKTFINAL .labelBHZEVTRLOC .labelBHZEVTFINAL .labelBHZEVTSACF .labelBHZSTATION .labelBHZWAVE .labelBHZPKTCMP -in .frameindicatorz -side top
    pack .frameindicatorz -in .framefilelist -side left -fill y

    frame .frameindicatore -borderwidth 1.0 -relief groove
    label .labelBHEPKTRLOC  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHEPKTFINAL -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHEEVTRLOC  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHEEVTFINAL -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHEEVTSACF  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHESTATION  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHEWAVE     -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHEPKTCMP   -justify left -bg $labelbg -width $labelwd -relief raised
    pack .labelBHEPKTRLOC .labelBHEPKTFINAL .labelBHEEVTRLOC .labelBHEEVTFINAL .labelBHEEVTSACF .labelBHESTATION .labelBHEWAVE .labelBHEPKTCMP -in .frameindicatore -side top
    pack .frameindicatore -in .framefilelist -side left -fill y

    frame .frameindicatorn -borderwidth 1.0 -relief groove
    label .labelBHNPKTRLOC  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHNPKTFINAL -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHNEVTRLOC  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHNEVTFINAL -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHNEVTSACF  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHNSTATION  -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHNWAVE     -justify left -bg $labelbg -width $labelwd -relief raised
    label .labelBHNPKTCMP   -justify left -bg $labelbg -width $labelwd -relief raised
    pack .labelBHNPKTRLOC .labelBHNPKTFINAL .labelBHNEVTRLOC .labelBHNEVTFINAL .labelBHNEVTSACF .labelBHNSTATION .labelBHNWAVE .labelBHNPKTCMP -in .frameindicatorn -side top
    pack .frameindicatorn -in .framefilelist -side left -fill y

    pack .framefilelist -in .leftcolumn -side top 


    pack .leftcolumn -side left -fill y;    # main control

    set wavefx 324; set wavefy 380
    frame .framewaveall -borderwidth 1.0 -relief groove
    frame .framewavez -borderwidth 1.0 -relief groove
    frame .framewavee -borderwidth 1.0 -relief groove
    frame .framewaven -borderwidth 1.0 -relief groove

    #######################  UD
    frame .framewavez1 -borderwidth 1.0 -relief groove
    label .zchannel -justify left -text "Z:\nZ:" -fg red
    button .loadz    -height 1 -width 2 -bg orange -text "Load"   -command {showini;loadchannel BHZ}
    button .saveptz  -height 1 -width 3 -bg orange -text "SavePT" -command {outpicktime BHZ}
    button .locatez  -height 1 -width 3 -bg orange -text "Locate" -command {outpicktime BHZ;relocate BHZ}
    button .zoomoutz -height 1 -width 1 -bg orange -text "A-"     -command {waveallscalezoomy BHZ 0.5}
    button .zoominz  -height 1 -width 1 -bg orange -text "A+"     -command {waveallscalezoomy BHZ 1.5}
    button .delz     -height 1 -width 1 -bg orange -text "Del"    -command {deletesub BHZ}
    button .figoutz  -height 1 -width 3 -bg orange -text "OutAll" -command {outall BHZ}
    pack .zchannel .loadz .saveptz .locatez .zoomoutz .zoominz .delz .figoutz -in .framewavez1 -side left

    frame .framewavez2 -borderwidth 1.0 -relief groove
    #button .zoomout -text "zoom-" -command {scaletraceinallz 0.5} -height 1 -width 5 -bg orange
    #button .zoomin  -text "zoom+" -command {scaletraceinallz 2.0} -height 1 -width 5 -bg orange
    label .zinfofile -justify left -text "Location-File" -width 41 -fg magenta
    label .zinfonew1 -justify left -text "Location-Rloc" -fg blue
    label .zinfonew2 -justify left -text "Location-Final" -fg darkorange
    pack  .zinfofile .zinfonew1 .zinfonew2 -in .framewavez2 -side top -fill x

    canvas .waveallz -background white -width $wavefx -height $wavefy -closeenough 1 \
                    -xscrollcommand ".waveallzx set" -yscrollcommand ".waveally set"
    scrollbar .waveally -command ".waveallz yview"
    scrollbar .waveallzx -command ".waveallz xview" -orient horizontal
    pack .framewavez1 .framewavez2 .waveallz -in .framewavez -side top -fill x

    ####################### EW
    frame .framewavee1 -borderwidth 1.0 -relief groove
    label .echannel -justify left -text "E:\nT:" -fg red
    button .loade    -height 1 -width 2 -bg orange -text "Load"   -command {showini;loadchannel BHE}
    button .savepte  -height 1 -width 3 -bg orange -text "SavePT" -command {outpicktime BHE}
    button .locatee  -height 1 -width 3 -bg orange -text "Locate" -command {
	if {[tk_messageBox -type yesno -message "Relocation through S-wave data.\nPlease ensure data having enough accuracy!\nContinue?"]=="yes"} {outpicktime BHE;relocate BHE}
    }
    button .zoomoute -height 1 -width 1 -bg orange -text "A-"     -command {waveallscalezoomy BHE 0.5}
    button .zoomine  -height 1 -width 1 -bg orange -text "A+"     -command {waveallscalezoomy BHE 1.5}
    button .dele     -height 1 -width 1 -bg orange -text "Del"    -command {deletesub BHE}
    button .figoute  -height 1 -width 3 -bg orange -text "OutAll" -command {outall BHE}
    pack .echannel .loade .savepte .locatee .zoomoute .zoomine .dele .figoute -in .framewavee1 -side left

    frame .framewavee2 -borderwidth 1.0 -relief groove
    label .einfofile -justify left -text "Location-File" -width 41 -fg magenta
    label .einfonew1 -justify left -text "Location-Rloc" -fg blue
    label .einfonew2 -justify left -text "Location-Final" -fg darkorange
    pack .einfofile .einfonew1 .einfonew2 -in .framewavee2 -side top -fill x

    canvas .wavealle -background white -width $wavefx -height $wavefy -closeenough 1
    pack .framewavee1 .framewavee2 .wavealle -in .framewavee -side top -fill x

    ####################### NS
    frame .framewaven1 -borderwidth 1.0 -relief groove
    label .nchannel -justify left -text "N:\nR:" -fg red
    button .loadn    -height 1 -width 2 -bg orange -text "Load"   -command {showini;loadchannel BHN}
    button .saveptn  -height 1 -width 3 -bg orange -text "SavePT" -command {outpicktime BHN}
    button .locaten  -height 1 -width 3 -bg orange -text "Locate" -command {
	if {[tk_messageBox -type yesno -message "Relocation through S-wave data.\nPlease ensure data having enough accuracy!\nContinue?"]=="yes"} {outpicktime BHN;relocate BHN}
    }
    button .zoomoutn -height 1 -width 1 -bg orange -text "A-"     -command {waveallscalezoomy BHN 0.5}
    button .zoominn  -height 1 -width 1 -bg orange -text "A+"     -command {waveallscalezoomy BHN 1.5}
    button .deln     -height 1 -width 1 -bg orange -text "Del"    -command {deletesub BHN}
    button .figoutn  -height 1 -width 3 -bg orange -text "OutAll" -command {outall BHN}
    pack .nchannel .loadn .saveptn .locaten .zoomoutn .zoominn .deln .figoutn -in .framewaven1 -side left

    frame .framewaven2 -borderwidth 1.0 -relief groove
    label .ninfofile -justify left -text "Location-File" -width 41 -fg magenta
    label .ninfonew1 -justify left -text "Location-Rloc" -fg blue
    label .ninfonew2 -justify left -text "Location-Final" -fg darkorange
    pack .ninfofile .ninfonew1 .ninfonew2 -in .framewaven2 -side top -fill x

    canvas .wavealln -background white -width $wavefx -height $wavefy -closeenough 1
    pack .framewaven1 .framewaven2 .wavealln -in .framewaven -side top -fill x
    #######################

    pack .waveally .framewaven .framewavee .framewavez -in .framewaveall -side right -fill y
    pack .framewaveall -in .righttop -side left; # main control-sub

    pack .righttop -side top -fill x;    # main control

    frame .frameshowhide -borderwidth 1.0 -relief groove
    #label .labelshowwaveall -justify center -text "Click here to show wave profile sorted by dist" -relief raised -width 40 -bg orange
    #label .labelhidewaveall -justify center -text "Click here to hide wave profile sorted by dist" -relief raised -width 40 -bg orange
    label .labelshowwaveall -justify center -text "Shrink measuring workspace" -relief raised -width 25 -bg orange
    label .labelhidewaveall -justify center -text "Enlarge measuring workspace" -relief raised -width 25 -bg orange
    pack .labelhidewaveall -in .frameshowhide -side right

    label .labeltext -justify center -text "Waveform information:" -width 18
    set waveformsw 1
    checkbutton .checkwaveform -text "Waveform" -variable waveformsw -command {waveformswsub}
    set envelopesw 0
    checkbutton .checkwaveenvelope -text "Envelope" -variable envelopesw -command {envelopeswsub}
    set instfreqsw 0
    checkbutton .checkwaveinstfreq -text "Instfreq" -variable instfreqsw -command {instfreqswsub}
    pack .labeltext .checkwaveform .checkwaveenvelope .checkwaveinstfreq -in .frameshowhide -side left

    menu .menuwaveheight -tearoff 0
	 .menuwaveheight add command -label "150" -command {set singlewaveheight 150; .twaveH delete 1.0 end; .twaveH insert 1.0 $singlewaveheight; updatewavezoom}
	 .menuwaveheight add command -label "100" -command {set singlewaveheight 100; .twaveH delete 1.0 end; .twaveH insert 1.0 $singlewaveheight; updatewavezoom}
	 .menuwaveheight add command -label "50"  -command {set singlewaveheight 50;  .twaveH delete 1.0 end; .twaveH insert 1.0 $singlewaveheight; updatewavezoom}
	 .menuwaveheight add command -label "30"  -command {set singlewaveheight 30;  .twaveH delete 1.0 end; .twaveH insert 1.0 $singlewaveheight; updatewavezoom}
    set singlewaveheight 100;
    #frame .framewaveH -borderwidth 1 -relief groove
    label .labelwaveH -justify center -text "Trace H:" -relief raised -width 8 -bg orange
    text .twaveH -width 3 -height 1.2 -bg white -fg red
         .twaveH insert 1.0 $singlewaveheight;
    #pack .labelwaveH .twaveH -side left -in .framewaveH
    pack .labelwaveH .twaveH -in .frameshowhide -side left


    label .labelzoomdn -justify center -text "A-" -relief raised -width 4 -bg orange
    label .labelzoomup -justify center -text "A+" -relief raised -width 4 -bg orange
    pack .labelzoomdn .labelzoomup -in .frameshowhide -side left

    set phasetype "P"
    radiobutton .phaseP -variable phasetype -text "P" -value "P"
    radiobutton .phaseS -variable phasetype -text "S" -value "S"
    pack .phaseP .phaseS -in .frameshowhide -side left

    label .labelapb -justify center -text "AutoPicker" -relief raised -width 9 -bg orange
    label .labelapp -justify center -text ">" -relief raised -width 1 -bg orange
    pack .labelapb .labelapp -in .frameshowhide -side left


    pack .frameshowhide -side top -fill x; # main control-sub

    set wavezoomfx 1000; set wavezoomfy 1000; # if wavezoomfx is modified then some parameter should be modified, such as
    frame .framewavezoom -borderwidth 1.0 -relief groove
    canvas .wavezoom -background white -width $wavezoomfx -height $wavezoomfy -closeenough 1 \
                     -xscrollcommand ".wavezoomx set"  -yscrollcommand ".wavezoomy set"
    scrollbar .wavezoomy -command "wavezoomverticalmove"
    #scrollbar .wavezoomx -command "wavezoomstainfomove;.wavezoom xview" -orient horizontal
    scrollbar .wavezoomx -command "wavezoomhorizontalmove" -orient horizontal
    pack .wavezoomx -in .framewavezoom -side bottom -fill x
    pack .wavezoomy .wavezoom -in .framewavezoom -side right -fill y

    pack .framewavezoom -in .rightbottom -side left; # main control-sub

    frame .framewavezoomreplot -borderwidth 1.0 -relief groove
    canvas .wavezoomreplot -background white -width $wavezoomfx -height $wavezoomfy -closeenough 1
    pack .wavezoomreplot -in .framewavezoomreplot -side right -fill y


    pack .rightbottom -side top -fill x;    # main control

    #############
    bindunit
    bindsel
    bindplot
    #############
    set line $argv
    regsub -all {[[:blank:]]+} $line " " line
    set datdir  [lindex [split $line " "] 0]
    set outdir  [lindex [split $line " "] 1]
    set stafile [lindex [split $line " "] 2]
    # Add eventid to .eventid
    set datdirlist [glob -nocomplain $datdir/*]
    set datdirlist [lsort $datdirlist]
    foreach tem $datdirlist {
 	set evtid [lindex [split $tem "/"] end]
	.evtid insert end $evtid
    }


    readstafile $stafile
    setstafigpara_plotsta_evtfree 1;

    #stationinitial
}
#################################################
proc setlocationpara {xval yval zval} {
    global dxx dyy dzz;

    set dxx $xval; .tdxx delete 1.0 end; .tdxx insert 1.0 $dxx;
    set dyy $yval; .tdyy delete 1.0 end; .tdyy insert 1.0 $dyy;
    set dzz $zval; .tdzz delete 1.0 end; .tdzz insert 1.0 $dzz;

}
#################################################
proc fileindicator {channel} {
    global datdir outdir evtidcur

    set PKTRLOC  "$outdir/$evtidcur.$channel.PICK.FORLOC";
    set PKTFINAL "$outdir/$evtidcur.$channel.PICK.FINAL";
    set EVTRLOC  "$outdir/$evtidcur.$channel.EVENT.RLOC";
    set EVTFINAL "$outdir/$evtidcur.$channel.EVENT.FINAL";
    set EVTSACF  "$outdir/$evtidcur.$channel.EVENT.SACF";
    set STATION  "$outdir/$evtidcur.$channel.STATION.PS";
    set WAVE     "$outdir/$evtidcur.$channel.WAVE.PS";
    set PKTCMP   "$outdir/$evtidcur.$channel.COMPARISON.PS";
    if {[file exists $PKTRLOC] ==1} {.label${channel}PKTRLOC  configure -bg yellow;} else {.label${channel}PKTRLOC  configure -bg lightgray;}
    if {[file exists $PKTFINAL]==1} {.label${channel}PKTFINAL configure -bg yellow;} else {.label${channel}PKTFINAL configure -bg lightgray;}
    if {[file exists $EVTRLOC] ==1} {.label${channel}EVTRLOC  configure -bg yellow;} else {.label${channel}EVTRLOC  configure -bg lightgray;}
    if {[file exists $EVTFINAL]==1} {.label${channel}EVTFINAL configure -bg yellow;} else {.label${channel}EVTFINAL configure -bg lightgray;}
    if {[file exists $EVTSACF] ==1} {.label${channel}EVTSACF  configure -bg yellow;} else {.label${channel}EVTSACF  configure -bg lightgray;}
    if {[file exists $STATION] ==1} {.label${channel}STATION  configure -bg yellow;} else {.label${channel}STATION  configure -bg lightgray;}
    if {[file exists $WAVE]    ==1} {.label${channel}WAVE     configure -bg yellow;} else {.label${channel}WAVE     configure -bg lightgray;}
    if {[file exists $PKTCMP]  ==1} {.label${channel}PKTCMP   configure -bg yellow;} else {.label${channel}PKTCMP   configure -bg lightgray;}
}
#################################################
proc deletesub {channel} {
    global datdir outdir evtidcur

    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}
    set ids [$canvasname find withtag wave]
    if {$ids==""} {return}

    set filelist [glob -nocomplain $outdir/$evtidcur.$channel.*]
    foreach element $filelist {
	#puts $element
	exec rm $element
    }
    fileindicator $channel

    if {$channel=="BHZ"} {set infofile .zinfofile; set infonew1 .zinfonew1; set infonew2 .zinfonew2;}
    if {$channel=="BHE"} {set infofile .einfofile; set infonew1 .einfonew1; set infonew2 .einfonew2;}
    if {$channel=="BHN"} {set infofile .ninfofile; set infonew1 .ninfonew1; set infonew2 .ninfonew2;}
    $infonew1 configure -text ""; $infonew2 configure -text ""
}
#################################################
proc outall {channel} {
    global datdir outdir evtidcur
    global year jday hour min sec msec mag dt
    global sfnum sfname sflat sflon sfele sfpick sfpt sfphase sfsnr
    global sfgcarc sfdist sfbaz sfaz
    global evtlon evtlat evtdep evtdate evttime evtomark
    global evtlonsac evtlatsac evtdepsac evtdatesac evttimesac evtomarksac
    global evtomarkvalsw
    global wavezoomfx wavezoomfy

    if {$channel=="BHZ"} {set canvasname .waveallz; set infonew2 .zinfonew2;}
    if {$channel=="BHE"} {set canvasname .wavealle; set infonew2 .einfonew2;}
    if {$channel=="BHN"} {set canvasname .wavealln; set infonew2 .ninfonew2;}
    set ids [$canvasname find withtag wave]
    if {$ids==""} {return}

    # output event information finally
    if {$evtomarkvalsw==0} {
	set chooseval [tk_messageBox -type yesno -message "Event position is manually marked, shoot time is not identified. Continue to write out all results?"]
	if {$chooseval=="no"} {return;}
    }
    set eventfinalfile "$outdir/$evtidcur.$channel.EVENT.FINAL"
    set fout [open $eventfinalfile w]
	puts $fout [format "%s %s %s %10.4f %10.4f %6.1f" $evtidcur $evtdate $evttime $evtlon $evtlat $evtdep]
    close $fout
    # Show message in label
    set evtlontext [format "%9.4f" $evtlon]
    set evtlattext [format "%8.4f" $evtlat]
    set evtdeptext [format "%5.1f" $evtdep]
    set infotext "$evtdate $evttime $evtlontext $evtlattext $evtdeptext"
    $infonew2 configure -text "$infotext"

    # Draw circle in station figure & mark event cross to this point
    drawgcarc $evtlon $evtlat orange 4
    .stafigure delete evt.4
    markevtpos $evtlon $evtlat 4
    set swinwaveallz [.waveallz find withtag wave]
    set swinwavealle [.wavealle find withtag wave]
    set swinwavealln [.wavealln find withtag wave]
    set swinwavezoom [.wavezoom find withtag wave]
    # reorder waveform profile
    # The next five "if then" sentences may be not needed, as those have been achieved in loadchannel subroutine. Just leave them here for check unforeseen problems. 
    if {$swinwaveallz!="" || $swinwavealle!="" || $swinwavealln!=""} {
	caldistindeg $evtlon $evtlat
	sortbydist
	synthtraveltime
    }
    if {$swinwaveallz!=""} {
	waveallini BHZ
	reorderwaveinall BHZ
	.waveallz delete synthtmark
	if {$evtomarkvalsw==1} {plotsynthtmarkinwaveall BHZ}
    }
    if {$swinwavealle!=""} {
	waveallini BHE
	reorderwaveinall BHE
	.wavealle delete synthtmark
	if {$evtomarkvalsw==1} {plotsynthtmarkinwaveall BHE}
    }
    if {$swinwavealln!=""} {
	waveallini BHN
	reorderwaveinall BHN
	.wavealln delete synthtmark
	if {$evtomarkvalsw==1} {plotsynthtmarkinwaveall BHN}
    }
    if {$swinwavezoom!=""} {
	reorderwaveinzoom
	.wavezoom delete synthtmark
	if {$evtomarkvalsw==1} {plotsynthtmarkinwavezoom .wavezoom}
    }

    # output figure
    set wavepsfile "$outdir/$evtidcur.$channel.WAVE.PS"
    set stationpsfile "$outdir/$evtidcur.$channel.STATION.PS"
    $canvasname postscript -file $wavepsfile
    .stafigure  postscript -file $stationpsfile

    # Check if ps2raster is available
    set ps2rasteravailable 0;
    set cmd "command -v ps2raster";
    try {
	exec sh -c $cmd
    } trap CHILDSTATUS {msg opt} {
	if {[lindex [dict get $opt -errorcode] 2] == 1} {
	    set ps2rasteravailable 0;
	}
    }

    if {$ps2rasteravailable==1} {
	exec ps2raster -A -E600 $wavepsfile
	exec ps2raster -A -E600 $stationpsfile
    }

    set swexist [winfo exists .topshowallpicks]
    if {$swexist==1} {
	if { [.topshowallpicks.picksfigure find withtag line]!=""} {
	    set comparisonchartpsfile "$outdir/$evtidcur.$channel.COMPARISON.PS"
	    .topshowallpicks.picksfigure postscript -file $comparisonchartpsfile
	    if {$ps2rasteravailable==1} {exec ps2raster -A -E600 $comparisonchartpsfile;}
	}
    }
    set swexist [winfo exists .topzoomcurve]
    if {$swexist==1} {
	if { [.topzoomcurve.curvefigure find withtag line]!=""} {
	    zoomcurvefigureinipre $channel
	    set zoomcurvepsfile "$outdir/$evtidcur.$channel.CURVE.PS"
	    .topzoomcurve.curvefigure postscript -file $zoomcurvepsfile
	    if {$ps2rasteravailable==1} {exec ps2raster -A -E600 $zoomcurvepsfile;}
	}
    }

    if {$swinwavezoom!=""} {
	set xratio1 [lindex [.wavezoom xview] 0]
	set xratio2 [lindex [.wavezoom xview] 1]
	set yratio1 [lindex [.wavezoom yview] 0]
	set yratio2 [lindex [.wavezoom yview] 1]
	set xval    [expr $wavezoomfx*$xratio1]
	#set widval  [expr $wavezoomfx*($xratio2-$xratio1)]
	set widval  [expr $wavezoomfx*(1.0-$xratio1)]
	set yval    [expr $wavezoomfy*$yratio1]
	#set heival  [expr $wavezoomfy*($yratio2-$yratio1)]
	#set heival  [expr $wavezoomfy*(1.0-$yratio1)]

	set wavezoompsfile "$outdir/$evtidcur.$channel.WAVEZOOM.PS"
	set height [expr $wavezoomfy/72*2.53995]
	#.wavezoom postscript -file $wavezoompsfile -x $xval -y $yval -width $widval -height $heival -pagey ${heivalpoi}c -pageanchor nw -fontmap "-*-Times-Bold-r-Normal--*-*-*"
	#.wavezoom postscript -file $wavezoompsfile -x $xval -y 0.0 -width $widval -height $wavezoomfy -pagey ${height}c -pageanchor nw;
	if {$ps2rasteravailable==1} {exec ps2raster -A -E600 $wavezoompsfile;}
    }

    # output picktime with event information
    set secs [expr $sec+$msec/1000.0]
    #puts "taddtcl $year $jday $hour $min $secs $evtomark";
    set line [taddtcl $year $jday $hour $min $secs $evtomark]
    set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
    set evtyear [lindex $line 0]; set evtjday [lindex $line 1]; set evthour [lindex $line 2];
    set evtmin  [lindex $line 3]; set secline [lindex $line 4];
    set seclist [split $secline "\."]
    set evtsec  [lindex $seclist 0];
    set val2    [scan $evtsec "%d"];
    set evtmsec [expr int(($secline-$val2)*1000.0)]
    set evtmsec [format "%03d" $evtmsec];

    # Bug for format
    set evtyearval [scan $evtyear "%d"];
    set evtjdayval [scan $evtjday "%d"];
    set evthourval [scan $evthour "%d"];
    set evtminval  [scan $evtmin  "%d"];
    set evtsecval  [scan $evtsec  "%d"];
    set evtmsecval [scan $evtmsec "%d"];

    # Calculate SNR
    signalnoiseratio $channel
    set pickfinalfile "$outdir/$evtidcur.$channel.PICK.FINAL"
    set fout [open $pickfinalfile w]
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	puts $fout [format "%8s %10.4f %10.4f %7.4f %1d %10.3f %1s %10.3f %04d %03d %02d %02d %02d %03d \
                            %.3f %04d %03d %02d %02d %02d %03d \
                            %10.4f %10.4f %6.1f %4.1f \
                            %s %7.4f %10.4f %10.4f %10.4f" \
                            $sfname($i) $sflat($sfname($i)) $sflon($sfname($i)) $sfele($sfname($i)) \
                            $sfpick($sfname($i).$channel) $sfpt($sfname($i).$channel) $sfphase($sfname($i).$channel) $sfsnr($sfname($i).$channel) \
                            $year $jday $hour $min $sec $msec \
                            $evtomark $evtyearval $evtjdayval $evthourval $evtminval $evtsecval $evtmsecval \
                            $evtlat $evtlon $evtdep $mag \
                            $evtidcur $dt $sfgcarc($sfname($i)) $sfbaz($sfname($i)) $sfaz($sfname($i))]
    }
    puts $fout "> $evtidcur"
    close $fout

    # output event information from SACfile
    set eventsacfile "$outdir/$evtidcur.$channel.EVENT.SACF"
    set fout [open $eventsacfile w]
	puts $fout [format "%s %s %s %10.4f %10.4f %6.1f" $evtidcur $evtdatesac $evttimesac $evtlonsac $evtlatsac $evtdepsac]
    close $fout

    fileindicator $channel
}
#################################################
proc signalnoiseratio {channel} {
    global sfnum sfdist sfpick sfpt sfphase
    global sfname
    global sfsnr
    #global chanidcur
    global sacdatwinb sacdatwine sacdat dt
    global sacdatfileb sacdatfilee

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($sfname($i).$channel) == 1} {
	    set twinl [expr $sfpt($sfname($i).$channel)-2.0]
	    set twinr [expr $sfpt($sfname($i).$channel)+2.0]
	    if {$twinl<$sacdatfileb($sfname($i).$channel)} {set twinl $sacdatfileb($sfname($i).$channel)}
	    if {$twinr>$sacdatfilee($sfname($i).$channel)} {set twinr $sacdatfilee($sfname($i).$channel)}
	    set twinhalf [expr ($twinr-$twinl)/2.0]
	    # Adjust twinl and twinr
	    #set twinl [expr $sfpt($sfname($i).$channel)-$twinhalf]
	    #set twinr [expr $sfpt($sfname($i).$channel)+$twinhalf]
	    # Locate to data-index point
	    set datbegint [expr int(($twinl-$sacdatfileb($sfname($i).$channel))/$dt)+1]
	    set datendint [expr int(($twinr-$sacdatfileb($sfname($i).$channel))/$dt)+1]
	    set datmidint [expr int(($sfpt($sfname($i).$channel)-$sacdatfileb($sfname($i).$channel))/$dt)+1]

	    # First part
	    set win1 0.0
	    for {set n $datbegint} {$n<$datmidint} {incr n 1} {
		set nn [expr $n-1]
		set sgval [lindex $sacdat($sfname($i).$channel) $nn]
		set win1 [expr $win1+$sgval*$sgval]
	    }
	    # Second part
	    set win2 0.0
	    for {set n $datmidint} {$n<=$datendint} {incr n 1} {
		set nn [expr $n-1]
		set sgval [lindex $sacdat($sfname($i).$channel) $nn]
		set win2 [expr $win2+$sgval*$sgval]
	    }

	    # Calculate SNR
	    if {$win1!=0.0} {
		set sfsnr($sfname($i).$channel) [expr $win2/$win1]
	    } else {
		set sfsnr($sfname($i).$channel) -12345.0
	    }
	} else {
	    set sfsnr($sfname($i).$channel) -12345.0
	}
    }
}
#################################################
proc scaletraceinallz {mul} {
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale

    set ids [.waveallz find withtag wave]
    if {$ids==""} {return}
    .waveallz scale wave $waveallbd 0 $mul 1
    #.waveallz configure -scrollregion [list 0 0 $wavefx $wavefy]
}
#################################################
proc waveallscalezoomy {channel mul} {
    global sfyaxisposinwaveall
    global sfnum dtname sfdist

    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}

    set ids [$canvasname find withtag wave]
    if {$ids==""} {return}

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	$canvasname scale wave.$dtname($i) 0 $sfyaxisposinwaveall($dtname($i).$channel) 1 $mul
    }
}
proc wavezoomscalezoomy {mul} {
    global sfyaxisposinwavezoom
    global chanidcur
    global sfnum dtname sfdist

    set ids [.wavezoom find withtag wave]
    if {$ids==""} {return}

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	.wavezoom scale         wave.$dtname($i) 0 $sfyaxisposinwavezoom($dtname($i).$chanidcur) 1 $mul
	.wavezoom scale waveenvelope.$dtname($i) 0 $sfyaxisposinwavezoom($dtname($i).$chanidcur) 1 $mul
	.wavezoom scale waveinstfreq.$dtname($i) 0 $sfyaxisposinwavezoom($dtname($i).$chanidcur) 1 $mul
    }
}
#################################################
proc setstafigpara_plotsta_evtfree {evtfreesw} {
    global stafx stafy
    global sc swidth sheight bd
    global xscale yscale x1 x2 y1 y2
    global stanum
    global stname stla stlo
    global cursorcross

    .tsx1 delete 1.0 end;.tsx1 insert 1.0 $x1;
    .tsx2 delete 1.0 end;.tsx2 insert 1.0 $x2;
    .tsy1 delete 1.0 end;.tsy1 insert 1.0 $y1;
    .tsy2 delete 1.0 end;.tsy2 insert 1.0 $y2;

    set bd 15;
    set swidth [expr $stafx-$bd-$bd]
    .tsw delete 1.0 end; .tsw insert 1.0 $swidth
    set sheight [expr ($y2-$y1)/($x2-$x1)*$swidth]
    #set sheight $swidth
    #if {$sheight<[expr 320-$bd-$bd]} {set sheight [expr 320-$bd-$bd]}
    #set tem [expr $sheight/$swidth];puts $tem
    set stafy [expr $sheight+$bd+$bd];
    #if {$stafy<320} {set stafy 320;}
    .stafigure configure -width $stafx -height $stafy
    .stafigure configure -scrollregion [list 0 0 $stafx $stafy]

    set xscale [expr  $swidth/($x2-$x1)];   # Define the number of each degree, dpi/deg
    set yscale [expr $sheight/($y2-$y1)];   # Define the number of each degree, dpi/deg

    #set id [lindex [.stafigure find withtag evt] 0]
    .stafigure delete line postext station gcarc

    .stafigure create line [expr $stafx-$bd-$swidth]                $bd  [expr $stafx-$bd]          $bd                -width 1 -tags line -fill blue
    .stafigure create line [expr $stafx-$bd-$swidth] [expr $bd+$sheight] [expr $stafx-$bd]         [expr $bd+$sheight] -width 1 -tags line -fill blue
    .stafigure create line [expr $stafx-$bd-$swidth]                $bd  [expr $stafx-$bd-$swidth] [expr $bd+$sheight] -width 1 -tags line -fill blue
    .stafigure create line [expr $stafx-$bd]                        $bd  [expr $stafx-$bd]         [expr $bd+$sheight] -width 1 -tags line -fill blue
    set tx1 [format "%.2f" $x1];set tx2 [format "%.2f" $x2];
    set ty1 [format "%.2f" $y1];set ty2 [format "%.2f" $y2];
    .stafigure create text [expr $stafx-$bd*1.6-$swidth]  $bd                 -text "($tx1° $ty2°)"  -tags postext -anchor sw -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    .stafigure create text [expr $stafx-$bd/2.0]         [expr $bd+$sheight]  -text "($tx2° $ty1°)"  -tags postext -anchor ne -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    .stafigure create text [expr $stafx-$bd*1.6-$swidth] [expr $bd+$sheight]  -text "($tx1° $ty1°)"  -tags postext -anchor nw -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    .stafigure create text [expr $stafx-$bd/2.0]          $bd                 -text "($tx2° $ty2°)"  -tags postext -anchor se -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    set distline [distaztcl $y1 $x1 $y1 $x2]
    set distwidth [lindex $distline 0]; set distwidth [expr $distwidth*6371.0*acos(-1.0)/180.0]; set distwidth1 [format "%.2f" $distwidth]
    set distline [distaztcl $y2 $x1 $y2 $x2]
    set distwidth [lindex $distline 0]; set distwidth [expr $distwidth*6371.0*acos(-1.0)/180.0]; set distwidth2 [format "%.2f" $distwidth]
    set distline [distaztcl $y1 $x1 $y2 $x1]
    set distheight [lindex $distline 0]; set distheight [expr $distheight*6371.0*acos(-1.0)/180.0]; set distheight [format "%.2f" $distheight]
    set xc [expr $stafx/2.0]; set yc [expr $stafy/2.0]
    .stafigure create text $xc                $bd                 -text "$distwidth2 km"  -tags postext -anchor s           -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    .stafigure create text $xc                [expr $bd+$sheight] -text "$distwidth1 km"  -tags postext -anchor n           -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    .stafigure create text $bd                $yc                 -text "$distheight km" -tags postext -anchor s -angle 90 -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    .stafigure create text [expr $bd+$swidth] $yc                 -text "$distheight km" -tags postext -anchor n -angle 90 -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;


    for {set i 1} {$i<=$stanum} {incr i} {
	set line [stafigurexy2poixy $stlo($i) $stla($i) $x1 $y1 $xscale $yscale $bd $sheight]
	#set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
	set curx [lindex $line 0]; set cury [lindex $line 1];
	set kstnm $stname($i)
	.stafigure creat oval [expr $curx-$sc] [expr $cury-$sc] [expr $curx+$sc] [expr $cury+$sc] \
		                  -tags [list station station.$kstnm] -outline gray
	#set kstnm [string range $kstnm 1 3]
	#.stafigure creat text $curx $cury -text $kstnm -tags [list kstnm kstnm.$kstnm] -fill blue -anchor c;
    }

    .stafigure create line 0 0 $stafx 0 -width 1 -tags [list linex1 linexy line] -fill green -state $cursorcross
    .stafigure create line 0 0 $stafx 0 -width 1 -tags [list linex2 linexy line] -fill green -state $cursorcross
    .stafigure create line 0 0 0 $stafy -width 1 -tags [list liney1 linexy line] -fill green -state $cursorcross
    .stafigure create line 0 0 0 $stafy -width 1 -tags [list liney2 linexy line] -fill green -state $cursorcross

    .stafigure delete evt evterr
    if {$evtfreesw==1} {
	set lonc [expr ($x1+$x2)/2.0]; set latc [expr ($y1+$y2)/2.0]
	markevtpos $lonc $latc 0
    }

}
#################################################
proc loadstapara {} {
    global stafill stawidth staoutline taglist

    set ids [.stafigure find withtag station]
    set taglist {}
    foreach element $ids {
	set tag [lindex [.stafigure gettags $element] 1]
	lappend taglist $tag
	set stafill($tag) [.stafigure itemcget $element -fill]
	set stawidth($tag) [.stafigure itemcget $element -width]
	set staoutline($tag) [.stafigure itemcget $element -outline]
    }

}
#################################################
proc resetupstapara {} {
    global stafill stawidth staoutline taglist
    foreach element $taglist {
	.stafigure itemconfigure $element -fill $stafill($element)
	.stafigure itemconfigure $element -width $stawidth($element)
	.stafigure itemconfigure $element -outline $staoutline($element)
    }
}
#################################################
proc loadevtpara {} {
    global evtnum temevtlon temevtlat 
    global gcarcexist gcarccrossnum gcarccolor
    global evtlonget evtlatget

    array unset temevtlon; array unset temevtlat
    set ids [.stafigure find withtag evt]
    if {$ids != ""} {
	set evtnum -10
	foreach element $ids {
	    set tag [lindex [.stafigure gettags $element] 1];
	    set num [lindex [split $tag "."] 1];
	    if {$evtnum<$num} {set evtnum $num}
	    getevtpos $element; set temevtlon($num) $evtlonget; set temevtlat($num) $evtlatget
	}
    }
    # Load information for gcarc
    set gcarcexist 0
    set ids [.stafigure find withtag gcarc]
    if {$ids != ""} {
	set gcarcexist 1
	foreach element $ids {
	    set gcarccrossnum [lindex [.stafigure gettags $element] 1];
	    set gcarccolor [.stafigure itemcget $element -fill]
	    break;
	}
    }
}
#################################################
proc resetupevtpara {} {
    global evtnum temevtlon temevtlat
    global gcarcexist gcarccrossnum gcarccolor

    foreach nn [array names temevtlon] {
	if {$nn==0} {
	    markevtpos $temevtlon($nn) $temevtlat($nn) 0
	} else {
	    markevtpos $temevtlon($nn) $temevtlat($nn) $nn
	}
    }
    # Reset gcarc
    if {$gcarcexist==1} {
	drawgcarc $temevtlon($gcarccrossnum) $temevtlat($gcarccrossnum) $gcarccolor $gcarccrossnum
    }
}
#################################################
proc wavezoomhorizontalmove {args} {
    global wavezoomfx
    if {[lindex $args 0]=="scroll"} {.wavezoom xview [lindex $args 0] [lindex $args 1] [lindex $args 2]}
    if {[lindex $args 0]=="moveto"} {.wavezoom xview [lindex $args 0] [lindex $args 1]}

    set ids [.wavezoom find withtag statext]
    if {$ids==""} {return}
    set tt [.wavezoom coords [lindex $ids 0]]
    set xx [lindex $tt 0]
    set xfractor [lindex [.wavezoom xview] 0]
    set delx [expr $wavezoomfx*$xfractor-$xx]
    .wavezoom move statext $delx 0
    plotscaleinwavezoom .wavezoom
}
proc wavezoomverticalmove {args} {
    if {[lindex $args 0]=="scroll"} {.wavezoom yview [lindex $args 0] [lindex $args 1] [lindex $args 2]}
    if {[lindex $args 0]=="moveto"} {.wavezoom yview [lindex $args 0] [lindex $args 1]}

    # Error control
    set ids [.wavezoom find withtag wave]; if {$ids==""} {return}
    plotscaleinwavezoom .wavezoom
}
#################################################
proc loadwavezoomstainfopicktmark {} {
    global stamarkposy stamarktext stamarkcolor stamarkfont stamarktagmatch stamarktagid statagidlist
    global tmarkposx tmarkposy1 tmarkposy2 tmarkwidth tmarktagmatch tmarktagid tmarktagidlist
    global wavezoomx1 wavezoomxscale wavezoombd

    set ids [.wavezoom find withtag statext]
    set statagidlist {}
    foreach element $ids {
	set tags [.wavezoom gettags $element]
	set tag [lindex $tags 1]
	lappend statagidlist $tag
	set stamarktagmatch($tag) [lindex $tags 0]
	set stamarktagid($tag)    [lindex $tags 1]
	set stamarkposy($tag) [lindex [.wavezoom coords $element] 1]
	set stamarktext($tag) [.wavezoom itemcget $element -text]
	set stamarkcolor($tag) [.wavezoom itemcget $element -fill]
	set stamarkfont($tag) [.wavezoom itemcget $element -font]
    }

    set ids [.wavezoom find withtag tmark]
    set tmarktagidlist {}
    foreach element $ids {
	set tags [.wavezoom gettags $element]
	set tag [lindex $tags 1]
	lappend tmarktagidlist $tag
	set tmarktagmatch($tag) [lindex $tags 0]
	set tmarktagid($tag)    [lindex $tags 1]
	set tt [.wavezoom coords $element];
	set tmarkposy1($tag) [lindex $tt 1]; set tmarkposy2($tag) [lindex $tt 3]
	set curx [lindex $tt 0]
	set tmarkposx($tag) [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
	set tmarkwidth($tag) [.wavezoom itemcget $element -width]
    }
}
proc resetupwavezoomstainfopicktmark {canvasname} {
    global stamarkposy stamarktext stamarkcolor stamarkfont stamarktagmatch stamarktagid statagidlist
    global wavezoomfx
    global tmarkposx tmarkposy1 tmarkposy2 tmarkwidth tmarktagmatch tmarktagid tmarktagidlist
    global wavezoomx1 wavezoomxscale wavezoombd

    $canvasname delete statext
    set xratio [lindex [.wavezoom xview] 0]
    set poix [expr $xratio*$wavezoomfx]
    foreach tag $statagidlist {
	$canvasname create text $poix $stamarkposy($tag) -text $stamarktext($tag) -tags [list $stamarktagmatch($tag) $stamarktagid($tag)] -fill $stamarkcolor($tag) -font $stamarkfont($tag) -anchor w
    }

    $canvasname delete tmark
    foreach tag $tmarktagidlist {
	set curx [wavezoomx2poix $tmarkposx($tag) $wavezoomx1 $wavezoomxscale $wavezoombd]
	$canvasname create line $curx $tmarkposy1($tag) $curx $tmarkposy2($tag) -width $tmarkwidth($tag) -tags [list $tmarktagmatch($tag) $tmarktagid($tag)] -fill red
    }
}
#################################################


#################################################
proc getevtpos {tag} {
    global evtlonget evtlatget
    global sc swidth sheight bd
    global xscale yscale x1 x2 y1 y2


    #set id [lindex [.stafigure find withtag tag] 0]
    #if {$id != ""} {
	set line [.stafigure coords $tag];
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line
	set evtpoix [lindex [split $line " "] 0];
	set evtpoiy [lindex [split $line " "] 1];
	set line [stafigurepoixy2xy $evtpoix $evtpoiy $x1 $y2 $xscale $yscale $bd];
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line
	set evtlonget [lindex [split $line " "] 0];
	set evtlatget [lindex [split $line " "] 1];

	set tagtem [lindex [.stafigure gettags $tag] 1];
	set num [lindex [split $tagtem "."] 1];
	return $num
    #}
}
#################################################
proc markevtpos {lonpos latpos number} {
    #global evtlon evtlat
    global sc swidth sheight bd
    global xscale yscale x1 x2 y1 y2
    global isrelocate
    global evtherrsd evtverrsd evtherrrms evtverrrms

    #if {$delsw==1} {
	#.stafigure delete evt
    #}

    set line [stafigurexy2poixy $lonpos $latpos $x1 $y1 $xscale $yscale $bd $sheight];set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line
    set evtpoix [lindex [split $line " "] 0];
    set evtpoiy [lindex [split $line " "] 1];
    .stafigure delete evt.$number
    if {$number==0} {
	# location for free
	.stafigure create text $evtpoix $evtpoiy -text "+" -tags [list evt evt.$number] -anchor c -font bold -fill black;
    } elseif {$number==1} {
	# location from SACfile
	.stafigure create text $evtpoix $evtpoiy -text "+" -tags [list evt evt.$number] -anchor c -font bold -fill magenta;
    } elseif {$number==4} {
	# location from FINALfile
	.stafigure create text $evtpoix $evtpoiy -text "+" -tags [list evt evt.$number] -anchor c -font bold -fill darkorange;
    } elseif {$number==2} {
	# location from EVTRLOCfile (SD)
	.stafigure create text $evtpoix $evtpoiy -text "+" -tags [list evt evt.$number] -anchor c -font bold -fill blue;
	set line [stafigurexy2poixy [expr $lonpos-$evtherrsd] $latpos $x1 $y1 $xscale $yscale $bd $sheight];set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line
	set poidiff [lindex [split $line " "] 0]; set poidiff [expr abs($poidiff-$evtpoix)]
	.stafigure delete evterr.sd
	.stafigure create oval [expr $evtpoix-$poidiff] [expr $evtpoiy-$poidiff] [expr $evtpoix+$poidiff] [expr $evtpoiy+$poidiff] -tags [list evterr evterr.sd] -outline gray -fill gray -stipple gray12
    } elseif {$number==3} {
	# location from EVTRLOCfile (RMS)
	.stafigure create text $evtpoix $evtpoiy -text "+" -tags [list evt evt.$number] -anchor c -font bold -fill blue;
	set line [stafigurexy2poixy [expr $lonpos-$evtherrrms] $latpos $x1 $y1 $xscale $yscale $bd $sheight];set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line
	set poidiff [lindex [split $line " "] 0]; set poidiff [expr abs($poidiff-$evtpoix)]
	.stafigure delete evterr.rms
	.stafigure create oval [expr $evtpoix-$poidiff] [expr $evtpoiy-$poidiff] [expr $evtpoix+$poidiff] [expr $evtpoiy+$poidiff] -tags [list evterr evterr.rms] -outline gray -fill gray -stipple gray12
    }
}
#################################################

#################################################
proc readstafile stafile {
    global x1 x2 y1 y2
    global stanum
    global stname stla stlo

    set y1 [exec awk "BEGIN{tem= 180} {if(\$2<tem)tem=\$2} END{print tem}" $stafile]
    set y2 [exec awk "BEGIN{tem=-180} {if(\$2>tem)tem=\$2} END{print tem}" $stafile]
    set x1 [exec awk "BEGIN{tem= 180} {if(\$3<tem)tem=\$3} END{print tem}" $stafile]
    set x2 [exec awk "BEGIN{tem=-180} {if(\$3>tem)tem=\$3} END{print tem}" $stafile]

    set fsort [open $stafile]
	set stanum 0;
	while {[gets $fsort line]>=0} {
	    set stanum [expr $stanum+1]
	    set line [string trim $line]
	    regsub -all {[[:blank:]]+} $line " " line
	    set stname($stanum) [lindex [split $line " "] 0]
	    set stla($stanum) [lindex [split $line " "] 1]
	    set stlo($stanum) [lindex [split $line " "] 2]
	    if {[string first "." $stname($stanum)]>=0} {
		tk_messageBox -type ok -message "Station name cannot contain dot, please check and confirm. $stname($stanum)";
		exit;
	    }
	}
    close $fsort
}

#################################################

#################################################
proc bindplot {} {
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global evtlon evtlat evtdep evtdate evttime evtomark
    global evtlonget evtlatget
    global evtdatesac evttimesac evtomarksac evtlonsac evtlatsac evtdepsac
    global evtdatesd  evttimesd  evtomarksd  evtlonsd  evtlatsd  evtdepsd evtlonrms evtlatrms evtdeprms
    global evtdatefin evttimefin evtomarkfin evtlonfin evtlatfin evtdepfin
    global evtdephand evtdepmax
    global sc swidth sheight bd
    global xscale yscale x1 x2 y1 y2
    global sfnum dtname sfdist sfgcarc
    global sfpt sfpick sfphase sfmansw
    global phasetype
    global zoomdragcontrol zoompress2control allpress2control
    global model
    global chanidcur
    global tmarkcurvesw
    global evtomarkvalsw
    global evtcrosspress
    global evterrpress

    set zoomdragcontrol 1
    set zoompress2control 1
    set allpress2control 1
    set evtomarkvalsw 1
    set evtcrosspress 0
    set evterrpress 0

    .stafigure bind evt.0 <1> {
	set curX [%W canvasx %x]
	set curY [%W canvasy %y]
    }
    .stafigure bind evt.0 <B1-Motion> {
	set temX [%W canvasx %x];    set temY [%W canvasy %y]
	set delX [expr $temX-$curX]; set delY [expr $temY-$curY]
	.stafigure move current $delX $delY
	set curX $temX
	set curY $temY
    }
    .stafigure bind evt <ButtonPress-3> {
	# disable measure distance
	set evtcrosspress 1

	# get & save evtpos in evtlon evtlat, evtdep is set up according to numval
	# num=0, free; num=1, evtdepsac; num=2, evtdepsd; num=3, evtdeprms
	set num [getevtpos current];
	# set evtdep and vibcolor
	if {$num==0} {set evtlon $evtlonget; set evtlat $evtlatget; set evtdep $evtdephand; set vibcolor black}
	if {$num==1} {set evtlon $evtlonsac; set evtlat $evtlatsac; set evtdep $evtdepsac;  set vibcolor magenta}
	if {$num==2} {set evtlon $evtlonsd;  set evtlat $evtlatsd;  set evtdep $evtdepsd;   set vibcolor blue}
	if {$num==3} {set evtlon $evtlonrms; set evtlat $evtlatrms; set evtdep $evtdeprms;  set vibcolor blue}
	if {$num==4} {set evtlon $evtlonfin; set evtlat $evtlatfin; set evtdep $evtdepfin;  set vibcolor orange}
	if {$evtdep>$evtdepmax} {
	    set evtdepmax $evtdep
	    .textevdp insert 1.0 $evtdepmax
	    .scaleevtdep configure -to $evtdepmax
	}
	set evtdephand $evtdep
	drawgcarc $evtlon $evtlat $vibcolor $num
	#markevtpos #if this work again, evtpos will have a small change due to accuracy

	#if {[.wavezoom find withtag wave]==""} {return}
	set swinwaveallz [.waveallz find withtag wave]
	set swinwavealle [.wavealle find withtag wave]
	set swinwavealln [.wavealln find withtag wave]
	set swinwavezoom [.wavezoom find withtag wave]

	if {$swinwaveallz!="" || $swinwavealle!="" || $swinwavealln!=""} {
	    caldistindeg $evtlon $evtlat
	    sortbydist
	    synthtraveltime
	}

	set evtomarkvalsw 1;
	if {$num==0} {set evtdate "0000/00/00"; set evttime "00:00:00.000"; set evtomark "-12345.0";  set evtomarkvalsw 0;}
	if {$num==1} {set evtdate $evtdatesac;  set evttime $evttimesac;    set evtomark $evtomarksac}
	if {$num==2} {set evtdate $evtdatesd;   set evttime $evttimesd;     set evtomark $evtomarksd}
	if {$num==3} {set evtdate $evtdatesd;   set evttime $evttimesd;     set evtomark $evtomarksd}
	if {$num==4} {set evtdate $evtdatefin;  set evttime $evttimefin;    set evtomark $evtomarkfin}

	if {$chanidcur=="BHZ"} {set infonew1 .zinfonew1;}
	if {$chanidcur=="BHE"} {set infonew1 .einfonew1;}
	if {$chanidcur=="BHN"} {set infonew1 .ninfonew1;}
	if {$num==2 || $num==3} {
	    set evtlonout [formatsub94 $evtlon]; set evtlatout [formatsub84 $evtlat]; set evtdepout [formatsub51 $evtdep];
	    $infonew1 configure -text "$evtdate $evttime $evtlonout $evtlatout $evtdepout"
	}

	## if $swinwavezoom!="", then chanidcur should have value
	#if {$swinwavezoom!=""} {
	    ##### set evtomark
	    #set seln -1
	    #for {set i 1} {$i<=$sfnum} {incr i 1} {
		#if {$sfpick($dtname($i).$chanidcur)==1} {set seln $i;break;}
	    #}
	    #if {$seln==-1} {
		#tk_messageBox -type ok -message "Warning, can not plot synthetic time, because there is no picked time in the current select channel."
	    #} else {
		#set line [exec modeltime.pl $model all $evtdep $sfgcarc($dtname($seln)) | tr "\n" " "]
		#set line [modeltimetcl $model all $evtdep $sfgcarc($dtname($seln))]
		#set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
		#set ttlist [split $line " "]
		#set onsettime [lindex $ttlist 1]
		#set evtomark [expr $sfpt($dtname($seln).$chanidcur)-$onsettime]
		######
		#set evtomarkvalsw 1;
	    #}
	#}

	if {$swinwaveallz!=""} {
	    waveallini BHZ
	    reorderwaveinall BHZ
	    .waveallz delete synthtmark
	    if {$evtomarkvalsw==1} {plotsynthtmarkinwaveall BHZ}
	}
	if {$swinwavealle!=""} {
	    waveallini BHE
	    reorderwaveinall BHE
	    .wavealle delete synthtmark
	    if {$evtomarkvalsw==1} {plotsynthtmarkinwaveall BHE}
	}
	if {$swinwavealln!=""} {
	    waveallini BHN
	    reorderwaveinall BHN
	    .wavealln delete synthtmark
	    if {$evtomarkvalsw==1} {plotsynthtmarkinwaveall BHN}
	}
	if {$swinwavezoom!=""} {
	    reorderwaveinzoom
	    .wavezoom delete synthtmark
	    if {$evtomarkvalsw==1} {plotsynthtmarkinwavezoom .wavezoom}
	}
	if {$chanidcur != ""} {checkpicksub $chanidcur;zoomcurvefigureinipre $chanidcur}
	#waveallzini
	#reorderwaveinallz
	#plotsynthtmarkinwaveall $channel
	#if {$swinwavezoom==""} {return;}
	#reorderwaveinzoom
	#plotsynthtmarkinwavezoom .wavezoom
	##plotwaveall
	##plottmarkinwaveall
    }
    .stafigure bind evterr <ButtonPress-3> {
	set evterrpress 1
	.stafigure itemconfigure current -state hidden
    }
    .stafigure bind evterr <ButtonPress-2> {
	.stafigure itemconfigure current -state hidden
    }
    .stafigure bind gcarc <2> {
	.stafigure delete gcarc
    }
    .stafigure bind station <2> {
	set id   [.stafigure find withtag current]
	set tag [lindex [.stafigure gettags $id] 1]
	set kstnm [lindex [split $tag "."] 1]
	commonselect $kstnm 1
    }
    .wavezoom bind synthtmark <Motion> {
	.wavezoom delete movelabel
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set tag [lindex [.wavezoom gettags current] 1]
	.wavezoom create text $curx $cury -text $tag -tag movelabel -fill blue -anchor s
    }
    .wavezoomreplot bind synthtmark <Motion> {
	.wavezoomreplot delete movelabel
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set tag [lindex [.wavezoomreplot gettags current] 1]
	.wavezoomreplot create text $curx $cury -text $tag -tag movelabel -fill blue -anchor s
    }
    .waveallz bind synthtmark <Motion> {
	.waveallz delete movelabel
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set tag [lindex [.waveallz gettags current] 1]
	.waveallz create text $curx $cury -text $tag -tag movelabel -fill blue -anchor s
    }
    .wavealle bind synthtmark <Motion> {
	.wavealle delete movelabel
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set tag [lindex [.wavealle gettags current] 1]
	.wavealle create text $curx $cury -text $tag -tag movelabel -fill blue -anchor s
    }
    .wavealln bind synthtmark <Motion> {
	.wavealln delete movelabel
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set tag [lindex [.wavealln gettags current] 1]
	.wavealln create text $curx $cury -text $tag -tag movelabel -fill blue -anchor s
    }
    .wavezoom       bind movelabel <Motion> {.wavezoom       delete movelabel}
    .wavezoomreplot bind movelabel <Motion> {.wavezoomreplot delete movelabel}
    .waveallz       bind movelabel <Motion> {.waveallz       delete movelabel}
    .wavealle       bind movelabel <Motion> {.wavealle       delete movelabel}
    .wavealln       bind movelabel <Motion> {.wavealln       delete movelabel}
    #################################################
    .wavezoomreplot bind tmark <ButtonPress-1> {
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
    }
    .wavezoomreplot bind tmark <B1-Motion> {
	set temx [%W canvasx %x]
	set temy [%W canvasy %y]
	set delx [expr $temx-$curx]
	set dely [expr $temy-$cury]
	.wavezoomreplot move current $delx 0
	set curx $temx;
	set cury $temy;

	# shift tmark in wavezoom continue
	set tag [lindex [.wavezoomreplot gettags current] 1]
	.wavezoom move $tag $delx 0

	# Prepare for shift tmark in waveall continue
	set tag [lindex [.wavezoomreplot gettags current] 1]
	set kstnm [lindex [split $tag "."] 1]
	set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
	set sfpt($kstnm.$chanidcur) $tt

	# shift tmark in waveall continue
	if {$tmarkcurvesw($chanidcur)==1} {
	    if {$chanidcur=="BHZ"} {set canvasname .waveallz}
	    if {$chanidcur=="BHE"} {set canvasname .wavealle}
	    if {$chanidcur=="BHN"} {set canvasname .wavealln}
	    set curxold [lindex [$canvasname coords tmark.$kstnm] 0]
	    set tt [waveallxy2poixy $sfpt($kstnm.$chanidcur) 0 $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
	    set curxnew [lindex $tt 0]
	    set delx [expr $curxnew-$curxold]
	    $canvasname move tmark.$kstnm $delx 0
	} else {
	    plotcurve $chanidcur
	}
	colorstation $chanidcur
	plotcurdots $chanidcur
	getpickcurrentevtdatrange $chanidcur
	plotzoomcurve $chanidcur
    }
    .wavezoomreplot bind tmark <ButtonRelease-1> {
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];

	# Set sfpt lastly
	set tag [lindex [.wavezoomreplot gettags current] 1]
	set kstnm [lindex [split $tag "."] 1]
	set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
	set sfpt($kstnm.$chanidcur) $tt
	set sfphase($kstnm.$chanidcur) $phasetype
	set sfmansw($kstnm.$chanidcur) 1

	# shift tmark in waveall lastly
	if {$tmarkcurvesw($chanidcur)==1} {
	    if {$chanidcur=="BHZ"} {set canvasname .waveallz}
	    if {$chanidcur=="BHE"} {set canvasname .wavealle}
	    if {$chanidcur=="BHN"} {set canvasname .wavealln}
	    set curxold [lindex [$canvasname coords tmark.$kstnm] 0]
	    set tt [waveallxy2poixy $sfpt($kstnm.$chanidcur) 0 $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
	    set curxnew [lindex $tt 0]
	    set delx [expr $curxnew-$curxold]
	    $canvasname move tmark.$kstnm $delx 0
	} else {
	    plotcurve $chanidcur
	}
	colorstation $chanidcur
	plotcurdots $chanidcur
	getpickcurrentevtdatrange $chanidcur
	plotzoomcurve $chanidcur
    }
    .wavezoomreplot bind tmark <ButtonPress-2> {set zoompress2control 0;}
    .wavezoomreplot bind tmark <ButtonRelease-2> {
	set zoompress2control 1
	if {$chanidcur=="BHZ"} {set canvasname .waveallz}
	if {$chanidcur=="BHE"} {set canvasname .wavealle}
	if {$chanidcur=="BHN"} {set canvasname .wavealln}
	set tag [lindex [.wavezoomreplot gettags current] 1]
	.wavezoomreplot delete $tag
	.wavezoom delete $tag
	$canvasname delete $tag
	set kstnm [lindex [split $tag "."] 1]
	set sfpick($kstnm.$chanidcur) 0
	set sfmansw($kstnm.$chanidcur) 0
	set sfphase($kstnm.$chanidcur) "U"
	.stafigure      itemconfigure station.$kstnm -fill green
	.wavezoom       itemconfigure statext.$kstnm -fill green3
	.wavezoomreplot itemconfigure statext.$kstnm -fill green3
	plotcurve $chanidcur
	colorstation $chanidcur
	plotcurdots $chanidcur
	getpickcurrentevtdatrange $chanidcur
	plotzoomcurve $chanidcur
    }
    .wavezoomreplot bind statext <ButtonPress-2> {set zoompress2control 0;}
    .wavezoomreplot bind statext <ButtonRelease-2> {
	set zoompress2control 1
	set id   [.wavezoomreplot find withtag current]
	set tag [lindex [.wavezoomreplot gettags $id] 1]
	set kstnm [lindex [split $tag "."] 1]
	commonselect $kstnm 0
    }
    .wavezoomreplot bind wave <ButtonPress-2> {set zoompress2control 0;}
    .wavezoomreplot bind wave <ButtonRelease-2> {
	set zoompress2control 1
	set id   [.wavezoomreplot find withtag current]
	set tag [lindex [.wavezoomreplot gettags $id] 1]
	set kstnm [lindex [split $tag "."] 1]
	commonselect $kstnm 0
    }
    ################################################
    .wavezoom bind tmark <ButtonPress-1> {
	set zoomdragcontrol 0
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
    }
    .wavezoom bind tmark <B1-Motion> {
	set temx [%W canvasx %x]
	set temy [%W canvasy %y]
	set delx [expr $temx-$curx]
	set dely [expr $temy-$cury]
	.wavezoom move current $delx 0
	set curx $temx;
	set cury $temy;

	# Prepare for shift tmark in waveall continue
	set tag [lindex [.wavezoom gettags current] 1]
	set kstnm [lindex [split $tag "."] 1]
	set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
	set sfpt($kstnm.$chanidcur) $tt

	# shift tmark in waveall continue
	if {$tmarkcurvesw($chanidcur)==1} {
	    if {$chanidcur=="BHZ"} {set canvasname .waveallz}
	    if {$chanidcur=="BHE"} {set canvasname .wavealle}
	    if {$chanidcur=="BHN"} {set canvasname .wavealln}
	    set curxold [lindex [$canvasname coords tmark.$kstnm] 0]
	    set tt [waveallxy2poixy $sfpt($kstnm.$chanidcur) 0 $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
	    set curxnew [lindex $tt 0]
	    set delx [expr $curxnew-$curxold]
	    $canvasname move tmark.$kstnm $delx 0
	} else {
	    plotcurve $chanidcur
	}
	colorstation $chanidcur
	plotcurdots $chanidcur
	getpickcurrentevtdatrange $chanidcur
	plotzoomcurve $chanidcur
	#.wavezoom move   tmark.$dtname($i) 0 $dely
	#set stamarkposy($tag) [lindex [.wavezoom coords $element] 1]
	#set stamarktext($tag) [.wavezoom itemcget $element -text]
	#.wavezoom itemconfigure wave.$tag -fill magenta
	#set tmarkposx($tag) [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
    }
    .wavezoom bind tmark <ButtonRelease-1> {
	set zoomdragcontrol 1
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];

	# Set sfpt lastly
	set tag [lindex [.wavezoom gettags current] 1]
	set kstnm [lindex [split $tag "."] 1]
	set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
	set sfpt($kstnm.$chanidcur) $tt
	set sfphase($kstnm.$chanidcur) $phasetype
	set sfmansw($kstnm.$chanidcur) 1

	# shift tmark in waveall lastly
	if {$tmarkcurvesw($chanidcur)==1} {
	    if {$chanidcur=="BHZ"} {set canvasname .waveallz}
	    if {$chanidcur=="BHE"} {set canvasname .wavealle}
	    if {$chanidcur=="BHN"} {set canvasname .wavealln}
	    set curxold [lindex [$canvasname coords tmark.$kstnm] 0]
	    set tt [waveallxy2poixy $sfpt($kstnm.$chanidcur) 0 $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
	    set curxnew [lindex $tt 0]
	    set delx [expr $curxnew-$curxold]
	    $canvasname move tmark.$kstnm $delx 0
	} else {
	    plotcurve $chanidcur
	}
	colorstation $chanidcur
	plotcurdots $chanidcur
	getpickcurrentevtdatrange $chanidcur
	plotzoomcurve $chanidcur
    }

    .wavezoom bind tmark <ButtonPress-2> {set zoompress2control 0;}
    .wavezoom bind tmark <ButtonRelease-2> {
	set zoompress2control 1
	if {$chanidcur=="BHZ"} {set canvasname .waveallz}
	if {$chanidcur=="BHE"} {set canvasname .wavealle}
	if {$chanidcur=="BHN"} {set canvasname .wavealln}
	set tag [lindex [.wavezoom gettags current] 1]
	.wavezoom delete $tag
	$canvasname delete $tag
	set kstnm [lindex [split $tag "."] 1]
	set sfpick($kstnm.$chanidcur) 0
	set sfmansw($kstnm.$chanidcur) 0
	set sfphase($kstnm.$chanidcur) "U"
	.stafigure itemconfigure station.$kstnm -fill green
	.wavezoom  itemconfigure statext.$kstnm -fill green3
	plotcurve $chanidcur
	colorstation $chanidcur
	plotcurdots $chanidcur
	getpickcurrentevtdatrange $chanidcur
	plotzoomcurve $chanidcur
    }
    .wavezoom bind statext <ButtonPress-2> {set zoompress2control 0;}
    .wavezoom bind statext <ButtonRelease-2> {
	set zoompress2control 1
	set id   [.wavezoom find withtag current]
	set tag [lindex [.wavezoom gettags $id] 1]
	set kstnm [lindex [split $tag "."] 1]
	commonselect $kstnm 0
    }
    .wavezoom bind wave <ButtonPress-2> {set zoompress2control 0;}
    .wavezoom bind wave <ButtonRelease-2> {
	set zoompress2control 1
	set id   [.wavezoom find withtag current]
	set tag [lindex [.wavezoom gettags $id] 1]
	set kstnm [lindex [split $tag "."] 1]
	commonselect $kstnm 0
    }
    .waveallz bind wave <ButtonPress-2> {set allpress2control 0;}
    .waveallz bind wave <ButtonRelease-2> {
	set allpress2control 1
	set id   [.waveallz find withtag current]
	set tag [lindex [.waveallz gettags $id] 1]
	set kstnm [lindex [split $tag "."] 1]
	commonselect $kstnm 1
    }
    .wavealle bind wave <ButtonPress-2> {set allpress2control 0;}
    .wavealle bind wave <ButtonRelease-2> {
	set allpress2control 1
	set id   [.wavealle find withtag current]
	set tag [lindex [.wavealle gettags $id] 1]
	set kstnm [lindex [split $tag "."] 1]
	commonselect $kstnm 1
    }
    .wavealln bind wave <ButtonPress-2> {set allpress2control 0;}
    .wavealln bind wave <ButtonRelease-2> {
	set allpress2control 1
	set id   [.wavealln find withtag current]
	set tag [lindex [.wavealln gettags $id] 1]
	set kstnm [lindex [split $tag "."] 1]
	commonselect $kstnm 1
    }
}
#################################################
proc commonselect {tag movesw} {
    # movesw==1 wavezoom move to yview 
    global sfnum dtname sfdist
    global backcolor

    set swinwaveallz [.waveallz find withtag wave]
    set swinwavealle [.wavealle find withtag wave]
    set swinwavealln [.wavealln find withtag wave]
    set swinwavezoom [.wavezoom find withtag wave]
    set swinwavezoomreplot [.wavezoomreplot find withtag wave]
    set swinstafigure [.stafigure find withtag station]
    if {$swinwaveallz==""} {return}

    set selnum -1	
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {[string compare $dtname($i) $tag]==0} {set selnum $i;break;}
    }
    if {$selnum==-1} {return}

    if {$swinwaveallz!=""} {
	switch [.waveallz itemcget wave.$tag -fill] {
	    gray    {.waveallz itemconfigure wave.$tag -fill magenta}
	    magenta {.waveallz itemconfigure wave.$tag -fill gray}
	}
    }
    if {$swinwavealle!=""} {
	switch [.wavealle itemcget wave.$tag -fill] {
	    gray    {.wavealle itemconfigure wave.$tag -fill magenta}
	    magenta {.wavealle itemconfigure wave.$tag -fill gray}
	}
    }
    if {$swinwavealln!=""} {
	switch [.wavealln itemcget wave.$tag -fill] {
	    gray    {.wavealln itemconfigure wave.$tag -fill magenta}
	    magenta {.wavealln itemconfigure wave.$tag -fill gray}
	}
    }
    if {$swinwavezoom!=""} {
	# Move wavezoom yview
	if {$movesw==1} {
	    set yratio [expr (($sfnum+1.0)-$selnum-2.0)/$sfnum]
	    .wavezoom yview moveto $yratio
	    plotscaleinwavezoom .wavezoom
	}
	# change color periodly
	switch [.wavezoom itemcget wave.$tag -fill] {
	    gray    {.wavezoom itemconfigure wave.$tag -fill magenta}
	    magenta {.wavezoom itemconfigure wave.$tag -fill gray}
	}
	#switch [.wavezoom itemcget statext.$tag -fill] {
	#    green    {.wavezoom itemconfigure statext.$tag -fill magenta;   set backcolor green}
	#    red      {.wavezoom itemconfigure statext.$tag -fill magenta;   set backcolor red}
	#    magenta  {.wavezoom itemconfigure statext.$tag -fill $backcolor;}
	#}
	switch [.wavezoom itemcget statext.$tag -font] {
	    TkDefaultFont    {.wavezoom itemconfigure statext.$tag -font bold}
	    bold    {.wavezoom itemconfigure statext.$tag -font TkDefaultFont}
	}
    }
    if {$swinstafigure!=""} {
	switch [.stafigure itemcget station.$tag -outline] {
	    gray    {.stafigure itemconfigure station.$tag -width 1.5 -outline purple}
	    purple {.stafigure itemconfigure station.$tag -width 1 -outline gray}
	}
    }
    if {$swinwavezoomreplot!=""} {
	switch [.wavezoomreplot itemcget wave.$tag -fill] {
	    gray    {.wavezoomreplot itemconfigure wave.$tag -fill magenta}
	    magenta {.wavezoomreplot itemconfigure wave.$tag -fill gray}
	}
	switch [.wavezoomreplot itemcget statext.$tag -font] {
	    TkDefaultFont    {.wavezoomreplot itemconfigure statext.$tag -font bold}
	    bold    {.wavezoomreplot itemconfigure statext.$tag -font TkDefaultFont}
	}
    }


}
#################################################
proc reorderwaveinall {channel} {
    global sfyaxisposinwaveall sfyaxisposinwavezoom
    global sfnum dtname sfdist
    global sacdatwinb sacdatwine sacdat dt
    global waveallbd wavezoombd
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfname sflat sflon sfele sfpick sfpt

    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	# Find the y-position
	set line [waveallxy2poixy 0 $sfdist($i) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
	set tmarky [lindex $line 1]
	#puts $dtname($i),$sfyaxisposinwaveall($dtname($i)),$tmarky

	# Calculate distance & move, include wave and tmark
	set dely [expr $tmarky-$sfyaxisposinwaveall($dtname($i).$channel)]
	$canvasname move  wave.$dtname($i) 0 $dely
	$canvasname move tmark.$dtname($i) 0 $dely

	# Save for the next
	set sfyaxisposinwaveall($dtname($i).$channel) $tmarky
    }
    plottmarkrangeinwaveall $channel
    plotcurve $channel
}
#################################################
proc reorderwaveinzoom {} {
    global sfyaxisposinwaveall sfyaxisposinwavezoom
    global sfnum dtname sfdist sfaz
    global sacdatwinb sacdatwine sacdat dt
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfname sflat sflon sfele sfpick sfpt
    global chanidcur

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	# Find the y-position
	set tmarky [expr $singlewaveheight*($sfnum+1-$i)]

	# Calculate distance & move, include wave, text, and tmark
	set dely [expr $tmarky-$sfyaxisposinwavezoom($dtname($i).$chanidcur)]
	.wavezoom move    wave.$dtname($i) 0 $dely
	.wavezoom move statext.$dtname($i) 0 $dely
	.wavezoom move   tmark.$dtname($i) 0 $dely
	.wavezoom move waveenvelope.$dtname($i) 0 $dely
	.wavezoom move waveinstfreq.$dtname($i) 0 $dely

	# Reload distance text
	set kstnm [lindex [.wavezoom itemcget statext.$dtname($i) -text] 0]
	set disttext [format "D:%6.2f" $sfdist($i)]
	set aztext [format "AZ:%6.2f" $sfaz($dtname($i))]
	set tt "$kstnm $disttext $aztext\n#$i"
	.wavezoom itemconfigure statext.$dtname($i) -text $tt
	
	# Save for the next
	set sfyaxisposinwavezoom($dtname($i).$chanidcur) $tmarky
    }
    plottmarkrangeinwavezoom;
}
#################################################
proc drawgcarc {lon lat color num} {
    global vibnum vibdist

    .stafigure delete gcarc
    set circlewidth 1
    #if {$color=="orange"} {set circlewidth 2}
    for {set i 1} {$i<=$vibnum} {incr i 1} {
	#if {} # fill circles with different colors
	set points [gcarcdat $lon $lat [expr ($vibdist/10.0*0.0899)*$i]]; #10km=0.0899
	.stafigure create line $points -width $circlewidth -tags [list gcarc $num] -fill $color
    }
}
#################################################
proc getpickcurrentevtdatrange {channel} {
    global pickcurlimt1 pickcurlimt2 pickcurlimd1 pickcurlimd2
    global sfnum dtname sfdist
    global sfname sflat sflon sfele sfpick sfpt
    global evtomark

    set pickcurlimt1 1.e+20; set pickcurlimt2 -1.e+20;
    set pickcurlimd1 1.e+20; set pickcurlimd2 -1.e+20;

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($dtname($i).$channel)==1} {
	    set timeval [expr $sfpt($dtname($i).$channel)-$evtomark]
	    if {$timeval<=$pickcurlimt1} {set pickcurlimt1 $timeval}
	    if {$timeval>=$pickcurlimt2} {set pickcurlimt2 $timeval}
	    if {$sfdist($i)<=$pickcurlimd1} {set pickcurlimd1 $sfdist($i)}
	    if {$sfdist($i)>=$pickcurlimd2} {set pickcurlimd2 $sfdist($i)}
	}
    }
}

proc zoomcurvefigureini {channel pickcurt1 pickcurt2 pickcurd1 pickcurd2} {
    global pickcurfigrx1 pickcurfigrx2 pickcurfigry1 pickcurfigry2
    global wavefx wavefy
    global pickcurfigbd pickcurfigfx pickcurfigfy
    global pickcurfigx1 pickcurfigx2 pickcurfigy1 pickcurfigy2 pickcurfigxscale pickcurfigyscale

    set swexist [winfo exists .topzoomcurve]; if {$swexist==0} {return}
    .topzoomcurve.curvefigure delete line curve curvedot
    if {$pickcurt1>=$pickcurt2 || $pickcurd1>=$pickcurd2} {return;}
    # recorde dat boundary
    set pickcurfigrx1 $pickcurt1; set pickcurfigrx2 $pickcurt2;
    set pickcurfigry1 $pickcurd1; set pickcurfigry2 $pickcurd2;

    set pickcurfigbd 25.0
    set pickcurfigfx $wavefx; set pickcurfigfy $wavefy

    set shiftdx [expr ($pickcurt2-$pickcurt1)/20.0]; #set shiftdx 0.0;
    set shiftdy [expr ($pickcurd2-$pickcurd1)/20.0]; #set shiftdy 0.0;
    set pickcurfigx1 [expr $pickcurt1-$shiftdx]; set pickcurfigx2 [expr $pickcurt2+$shiftdx];
    set pickcurfigy1 [expr $pickcurd1-$shiftdy]; set pickcurfigy2 [expr $pickcurd2+$shiftdy];
    set pickcurfigxscale [expr ($pickcurfigfx-$pickcurfigbd*2.0)/($pickcurfigx2-$pickcurfigx1)];
    set pickcurfigyscale [expr ($pickcurfigfy-$pickcurfigbd*2.0)/($pickcurfigy2-$pickcurfigy1)];

    set line [waveallxy2poixy $pickcurfigx1 $pickcurfigy1 $pickcurfigx1 $pickcurfigy1 $pickcurfigxscale $pickcurfigyscale $pickcurfigbd $pickcurfigfy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx1 [lindex [split $line " "] 0]; set ly1 [lindex [split $line " "] 1]
    set line [waveallxy2poixy $pickcurfigx2 $pickcurfigy1 $pickcurfigx1 $pickcurfigy1 $pickcurfigxscale $pickcurfigyscale $pickcurfigbd $pickcurfigfy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx2 [lindex [split $line " "] 0]; set ly2 [lindex [split $line " "] 1]
    .topzoomcurve.curvefigure create line $lx1 $ly1 $lx2 $ly2 -width 2 -tags [list line] -fill blue
    set lxmarkx [expr ($lx1+$lx2)/2.0];

    set line [waveallxy2poixy $pickcurfigx1 $pickcurfigy2 $pickcurfigx1 $pickcurfigy1 $pickcurfigxscale $pickcurfigyscale $pickcurfigbd $pickcurfigfy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx2 [lindex [split $line " "] 0]; set ly2 [lindex [split $line " "] 1]
    .topzoomcurve.curvefigure create line $lx1 $ly1 $lx2 $ly2 -width 2 -tags [list line] -fill blue
    set lymarky [expr ($ly1+$ly2)/2.0];
    if {$channel=="BHZ"} {set titletext "UD or Z";}
    if {$channel=="BHE"} {set titletext "EW or T";}
    if {$channel=="BHN"} {set titletext "NS or R";}
    .topzoomcurve.curvefigure create text $lxmarkx $pickcurfigbd           -text $titletext -tags [list line] -fill red            -anchor s -font -*-times-bold-r-normal-*-*-150-*-*-*-*-*-*;
    .topzoomcurve.curvefigure create text $lxmarkx [expr $pickcurfigfy-13] -text "Time(s)"  -tags [list line] -fill blue           -anchor n -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    .topzoomcurve.curvefigure create text 13 $lymarky                      -text "Dist(km)" -tags [list line] -fill blue -angle 90 -anchor s -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;

    set intervalx [expr round(($pickcurfigx2-$pickcurfigx1)/5.0/10.0)]; # split to 5 parts and get the nearest 10-based number
    set intervalx [expr $intervalx*10.0]
    set intervalmarkx [expr ($pickcurfigx2-$pickcurfigx1)/5.0-$intervalx]
    for {set xx [expr int(ceil($pickcurfigx1))]} {$xx<=[expr int($pickcurfigx2)]} {incr xx 1} {
	# in the case of small-limited
        if {$intervalx<=0.0} {
	    if {$intervalmarkx>=2.5} {set intervalx 2.0;}
	    if {$intervalmarkx<=2.5} {set intervalx 1.0;}
	}
	if {fmod($xx,$intervalx)==0.0} {
	    set line [waveallxy2poixy $xx $pickcurfigy1 $pickcurfigx1 $pickcurfigy1 $pickcurfigxscale $pickcurfigyscale $pickcurfigbd $pickcurfigfy]
	    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
	    set lx [lindex [split $line " "] 0]; set ly [lindex [split $line " "] 1]
	    .topzoomcurve.curvefigure creat oval [expr $lx-2] [expr $ly-2] [expr $lx+2] [expr $ly+2] -tags [list line] -fill red -outline red
	    .topzoomcurve.curvefigure creat text $lx [expr $ly+2] -text $xx -tags [list line] -fill blue -anchor n -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;;
	}
    }

    set intervaly [expr round(($pickcurfigy2-$pickcurfigy1)/5.0/10.0)]; # split to 5 parts and get the nearest 10-based number
    set intervaly [expr $intervaly*10.0]
    set intervalmarky [expr ($pickcurfigy2-$pickcurfigy1)/5.0-$intervaly]
    for {set yy [expr int(ceil($pickcurfigy1))]} {$yy<=[expr int($pickcurfigy2)]} {incr yy 1} {
	# in the case of small-limited
        if {$intervaly<=0.0} {
	    if {$intervalmarky>=2.5} {set intervaly 2.0;}
	    if {$intervalmarky<=2.5} {set intervaly 1.0;}
	}
	if {fmod($yy,$intervaly)==0.0} {
	    set line [waveallxy2poixy $pickcurfigx1 $yy $pickcurfigx1 $pickcurfigy1 $pickcurfigxscale $pickcurfigyscale $pickcurfigbd $pickcurfigfy]
	    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
	    set lx [lindex [split $line " "] 0]; set ly [lindex [split $line " "] 1]
	    .topzoomcurve.curvefigure creat oval [expr $lx-2] [expr $ly-2] [expr $lx+2] [expr $ly+2] -tags [list line] -fill red -outline red
	    .topzoomcurve.curvefigure creat text [expr $lx-2] $ly -text $yy -tags [list line] -fill blue -angle 90 -anchor s -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*; # -font fixed; "-*-times-medium-r-normal-*-18-*-*-*-*-*-iso8859-1";
	}
    }

}
proc plotzoomcurve {channel} {
    global sfnum dtname sfdist
    global sfname sflat sflon sfele sfpick sfpt
    global pickcurfigbd pickcurfigfx pickcurfigfy
    global pickcurfigx1 pickcurfigx2 pickcurfigy1 pickcurfigy2 pickcurfigxscale pickcurfigyscale
    global evtomark
    global curvepoix curvepoiy

    set swexist [winfo exists .topzoomcurve]
    if {$swexist==0} {return}
    if {$evtomark==-12345.0} {return}

    if {[.topzoomcurve.curvefigure find withtag line]==""} {return}
    .topzoomcurve.curvefigure delete curve curvedot
    set aa {}
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($dtname($i).$channel)==1} {
	    set tt [expr $sfpt($dtname($i).$channel)-$evtomark]
	    set tt [waveallxy2poixy $tt $sfdist($i) $pickcurfigx1 $pickcurfigy1 $pickcurfigxscale $pickcurfigyscale $pickcurfigbd $pickcurfigfy]
	    set tmarkx [lindex $tt 0]; set curvepoix($dtname($i).$channel) $tmarkx
	    set tmarky [lindex $tt 1]; set curvepoiy($dtname($i).$channel) $tmarky
	    .topzoomcurve.curvefigure create oval [expr $tmarkx-2] [expr $tmarky-2] [expr $tmarkx+2] [expr $tmarky+2] -tags [list curvedot $dtname($i)] -fill red -outline black
	    lappend aa $tmarkx $tmarky
	}
    }
    if {[llength $aa] <4} {return}
    .topzoomcurve.curvefigure create line $aa -width 2 -tags [list curve] -fill red
}
proc zoomcurvesub {poix poiy scale} {
    global pickcurlimt1 pickcurlimt2 pickcurlimd1 pickcurlimd2
    global pickcurfigbd pickcurfigfx pickcurfigfy
    global pickcurfigx1 pickcurfigx2 pickcurfigy1 pickcurfigy2 pickcurfigxscale pickcurfigyscale
    global chanidcur

    if {[.topzoomcurve.curvefigure find withtag line]==""} {return}

    set tt [waveallpoixy2xy $poix $poiy $pickcurfigx1 $pickcurfigy2 $pickcurfigxscale $pickcurfigyscale $pickcurfigbd]
    set lon [lindex $tt 0];
    set lat [lindex $tt 1];

    # Zoomin/out Range
    set xx1 [expr $lon-($lon-$pickcurfigx1)*$scale]; if {$xx1<$pickcurlimt1} {set xx1 $pickcurlimt1}
    set xx2 [expr $lon+($pickcurfigx2-$lon)*$scale]; if {$xx2>$pickcurlimt2} {set xx2 $pickcurlimt2}
    set yy1 [expr $lat-($lat-$pickcurfigy1)*$scale]; if {$yy1<$pickcurlimd1} {set yy1 $pickcurlimd1}
    set yy2 [expr $lat+($pickcurfigy2-$lat)*$scale]; if {$yy2>$pickcurlimd2} {set yy2 $pickcurlimd2}

    # Replot Range;
    zoomcurvefigureini $chanidcur $xx1 $xx2 $yy1 $yy2
    plotzoomcurve $chanidcur
}
proc zoomcurvefigureinipre {channel} {
    global pickcurlimt1 pickcurlimt2 pickcurlimd1 pickcurlimd2

    set swexist [winfo exists .topzoomcurve]
    if {$swexist==1} {
	getpickcurrentevtdatrange $channel
	zoomcurvefigureini $channel $pickcurlimt1 $pickcurlimt2 $pickcurlimd1 $pickcurlimd2
	plotzoomcurve $channel
    }
}
#################################################
proc loadbindfortopzoomcurve {} {
    global pickcurfigbd pickcurfigfx pickcurfigfy
    global pickcurfigx1 pickcurfigx2 pickcurfigy1 pickcurfigy2 pickcurfigxscale pickcurfigyscale
    global pickcurlimt1 pickcurlimt2 pickcurlimd1 pickcurlimd2
    global pickcurfigrx1 pickcurfigrx2 pickcurfigry1 pickcurfigry2
    global chanidcur
    global curvepoix curvepoiy
    global sfnum dtname sfdist
    global sfname sflat sflon sfele sfpick sfpt

    bind .topzoomcurve.curvefigure <Motion> {
	if {[.topzoomcurve.curvefigure find withtag line]==""} {return}
	set curx [%W canvasx %x]
	set cury [%W canvasy %y]
	set tt [waveallpoixy2xy $curx $cury $pickcurfigx1 $pickcurfigy2 $pickcurfigxscale $pickcurfigyscale $pickcurfigbd]
	set t1 [lindex $tt 0]; set t1 [formatsub3 $t1]
	set t2 [lindex $tt 1]; set t2 [formatsub3 $t2]
	.information configure -text "Poix:$curx Poiy:$cury Time-Dist:$t1 $t2"
    }
    bind .topzoomcurve.curvefigure <Button-4> {
	set curx [%W canvasx %x]
	set cury [%W canvasy %y]
	zoomcurvesub $curx $cury 0.7
    }
    bind .topzoomcurve.curvefigure <Button-5> {
	set curx [%W canvasx %x]
	set cury [%W canvasy %y]
	zoomcurvesub $curx $cury 1.2
    }
    bind .topzoomcurve.curvefigure <ButtonPress-1> {
	set curX %x
	set curY %y
	# The same as above two lines
	#set curX [%W canvasx %x]
	#set curY [%W canvasy %y]
	%W config -cursor hand1
    }
    bind .topzoomcurve.curvefigure <B1-Motion> {
	set delX [expr %x-$curX]
	set delY [expr %y-$curY]
	set curX %x
	set curY %y
	# The same as above two lines
	#set curX [%W canvasx %x]
	#set curY [%W canvasy %y]

	if {[.topzoomcurve.curvefigure find withtag line]==""} {return}
	# Move range
	set londel [expr $delX/$pickcurfigxscale]
	set latdel [expr $delY/$pickcurfigyscale]
	# The next four lines: if comment, it will appear as finger drag in android system (move & zoom). if uncomment, only move
	#if {[expr $pickcurfigrx1-$londel]<$pickcurlimt1 || [expr $pickcurfigrx1-$londel]>$pickcurlimt2} {return}
	#if {[expr $pickcurfigrx2-$londel]<$pickcurlimt1 || [expr $pickcurfigrx2-$londel]>$pickcurlimt2} {return}
	#if {[expr $pickcurfigry1+$latdel]<$pickcurlimd1 || [expr $pickcurfigry1+$latdel]>$pickcurlimd2} {return}
	#if {[expr $pickcurfigry2+$latdel]<$pickcurlimd1 || [expr $pickcurfigry2+$latdel]>$pickcurlimd2} {return}
	set xx1 [expr $pickcurfigrx1-$londel]; if {$xx1<$pickcurlimt1} {set xx1 $pickcurlimt1}
	set xx2 [expr $pickcurfigrx2-$londel]; if {$xx2>$pickcurlimt2} {set xx2 $pickcurlimt2}
	set yy1 [expr $pickcurfigry1+$latdel]; if {$yy1<$pickcurlimd1} {set yy1 $pickcurlimd1}
	set yy2 [expr $pickcurfigry2+$latdel]; if {$yy2>$pickcurlimd2} {set yy2 $pickcurlimd2}

	# Replot Range;
	zoomcurvefigureini $chanidcur $xx1 $xx2 $yy1 $yy2
	plotzoomcurve $chanidcur
    }
    bind .topzoomcurve.curvefigure <ButtonRelease-1> {
	%W config -cursor arrow
    }
    bind .topzoomcurve.curvefigure <Double-Button-1> {
	if {[.topzoomcurve.curvefigure find withtag line]==""} {return}
	set curX %x
	set curY %y
	set distmin 100;
	for {set i 1} {$i<=$sfnum} {incr i 1} {
	    if {$sfpick($dtname($i).$chanidcur)==1} {
		set delx [expr $curvepoix($dtname($i).$chanidcur)-$curX]
		set dely [expr $curvepoiy($dtname($i).$chanidcur)-$curY]
		set disttem [expr sqrt($delx*$delx+$dely*$dely)]
		if {$disttem<$distmin} {
		    set distmin $disttem
		    set seln $i
		}
	    }
	}

	if {$distmin<5} {
	    commonselect $dtname($seln) 1
	}
    }
    .topzoomcurve.curvefigure bind curvedot <2> {
	set tagtem [lindex [.topzoomcurve.curvefigure gettags current] 1]
	commonselect $tagtem 1
    }

}
proc loadbindfortopshowallpicks {} {
    global pickfigbd pickfigfx pickfigfy
    global pickfigx1 pickfigx2 pickfigy1 pickfigy2 pickfigxscale pickfigyscale
    global chanidcur
    global pickspoix pickspoiy
    global sfnum dtname sfdist
    global sfname sflat sflon sfele sfpick sfpt

    bind .topshowallpicks.picksfigure <Motion> {
	if {[.topshowallpicks.picksfigure find withtag line]==""} {return}
	set curx [%W canvasx %x]
	set cury [%W canvasy %y]
	set tt [waveallpoixy2xy $curx $cury $pickfigx1 $pickfigy2 $pickfigxscale $pickfigyscale $pickfigbd]
	set t1 [lindex $tt 0]; set t1 [formatsub3 $t1]
	set t2 [lindex $tt 1]; set t2 [formatsub3 $t2]
	.information configure -text "Poix:$curx Poiy:$cury Dist-Time:$t1 $t2"
    }
    bind .topshowallpicks.picksfigure <Double-Button-1> {
	if {[.topshowallpicks.picksfigure find withtag line]==""} {return}
	set curX %x
	set curY %y
	set distmin 100;
	for {set i 1} {$i<=$sfnum} {incr i 1} {
	    if {$sfpick($dtname($i).$chanidcur)==1} {
		set delx [expr $pickspoix($dtname($i).$chanidcur)-$curX]
		set dely [expr $pickspoiy($dtname($i).$chanidcur)-$curY]
		set disttem [expr sqrt($delx*$delx+$dely*$dely)]
		if {$disttem<$distmin} {
		    set distmin $disttem
		    set seln $i
		}
	    }
	}

	if {$distmin<5} {
	    commonselect $dtname($seln) 1
	}
    }
    .topshowallpicks.picksfigure bind curdots <2> {
	set tagtem [lindex [.topshowallpicks.picksfigure gettags current] 1]
	commonselect $tagtem 1
    }
}

#################################################
proc getpickdatandrange {channel} {
    global outdir
    global pickt1 pickt2 pickd1 pickd2
    global pickdat
    global distlim
    global evtlat evtlon evtdep

    set pickfilelist [glob -nocomplain $outdir/*.$channel.PICK.FINAL]
    set lines "";
    foreach element $pickfilelist {
	if {$distlim != "ALL"} {
	    set posline [exec awk "{if(NR==1){if(\$7\"\" == \"U\" || \$7\"\" == \"S\" || \$7\"\" == \"P\"){print \$22,\$23}else{print \$20,\$21}}}" $element]
	    #set posline [exec awk "{if(NR==1){print \$20,\$21}}" $element]
	    set argslist [split $posline " "]
	    set lattem [lindex $argslist 0]; set lontem [lindex $argslist 1];
	    set distline [distaztcl $evtlat $evtlon $lattem $lontem]
	    set distline [string trim $distline];regsub -all {[[:blank:]]+} $distline " " distline
	    set dist [lindex $distline 0]
	    set dist [expr $dist*6371.0*acos(-1.0)/180.0]
	    if {$dist>$distlim} {continue;}
	}
	set line [exec awk "{if(\$5==1){if(\$7\"\" == \"U\" || \$7\"\" == \"S\" || \$7\"\" == \"P\"){print \$6-\$15,\$28}else{print \$6-\$13,\$26}}}" $element]
	#set line [exec awk "{if(\$5==1){print \$6-\$13,\$26}}" $element]
	set lines "$lines\n$line"
    }

    set pickdat ""
    set pickt1  1.e+20;
    set pickt2 -1.e+20;
    set pickd1  1.e+20;
    set pickd2 -1.e+20;
    set lines [split $lines "\n"]
    set linenum [llength $lines]
    for {set i 0} {$i<$linenum} {incr i 1} {
	set line [lindex $lines $i];
	if {$line==""} {continue;}
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	set linelist [split $line " "];
	set t [lindex $linelist 0];
	set dist [lindex $linelist 1];
	set dist [expr $dist*6371.0*acos(-1.0)/180.0]
	if {$t<$pickt1} {set pickt1 $t}
	if {$t>$pickt2} {set pickt2 $t}
	if {$dist<$pickd1} {set pickd1 $dist}
	if {$dist>$pickd2} {set pickd2 $dist}
	lappend pickdat $t $dist
    }

}

proc plotallpicks {channel} {
    global pickt1 pickt2 pickd1 pickd2
    global pickdat
    global pickfigbd pickfigfx pickfigfy wavefx wavefy
    global pickfigx1 pickfigx2 pickfigy1 pickfigy2 pickfigxscale pickfigyscale

    getpickdatandrange $channel
    .topshowallpicks.picksfigure delete line dots curdots
    if {$pickt1>=$pickt2 || $pickd1>=$pickd2} {return;}
    set pickfigbd 25.0
    set pickfigfx $wavefx; set pickfigfy $wavefy

    set shiftdx [expr ($pickd2-$pickd1)/20.0];
    set shiftdy [expr ($pickt2-$pickt1)/20.0];
    set pickfigx1 [expr $pickd1-$shiftdx]; set pickfigx2 [expr $pickd2+$shiftdx];
    set pickfigy1 [expr $pickt1-$shiftdy]; set pickfigy2 [expr $pickt2+$shiftdy];
    set pickfigxscale [expr ($pickfigfx-$pickfigbd*2.0)/($pickfigx2-$pickfigx1)];
    set pickfigyscale [expr ($pickfigfy-$pickfigbd*2.0)/($pickfigy2-$pickfigy1)];

    set line [waveallxy2poixy $pickfigx1 $pickfigy1 $pickfigx1 $pickfigy1 $pickfigxscale $pickfigyscale $pickfigbd $pickfigfy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx1 [lindex [split $line " "] 0]; set ly1 [lindex [split $line " "] 1]
    set line [waveallxy2poixy $pickfigx2 $pickfigy1 $pickfigx1 $pickfigy1 $pickfigxscale $pickfigyscale $pickfigbd $pickfigfy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx2 [lindex [split $line " "] 0]; set ly2 [lindex [split $line " "] 1]
    .topshowallpicks.picksfigure create line $lx1 $ly1 $lx2 $ly2 -width 2 -tags [list line] -fill blue
    set lxmarkx [expr ($lx1+$lx2)/2.0];

    set line [waveallxy2poixy $pickfigx1 $pickfigy2 $pickfigx1 $pickfigy1 $pickfigxscale $pickfigyscale $pickfigbd $pickfigfy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx2 [lindex [split $line " "] 0]; set ly2 [lindex [split $line " "] 1]
    .topshowallpicks.picksfigure create line $lx1 $ly1 $lx2 $ly2 -width 2 -tags [list line] -fill blue
    set lymarky [expr ($ly1+$ly2)/2.0];
    if {$channel=="BHZ"} {set titletext "UD or Z";}
    if {$channel=="BHE"} {set titletext "EW or T";}
    if {$channel=="BHN"} {set titletext "NS or R";}
    .topshowallpicks.picksfigure create text $lxmarkx $pickfigbd           -text $titletext -tags [list line] -fill red            -anchor s -font -*-times-bold-r-normal-*-*-150-*-*-*-*-*-*;
    .topshowallpicks.picksfigure create text $lxmarkx [expr $pickfigfy-13] -text "Dist(km)" -tags [list line] -fill blue           -anchor n -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    .topshowallpicks.picksfigure create text 13 $lymarky                   -text "Time(s)"  -tags [list line] -fill blue -angle 90 -anchor s -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;

    set intervalx [expr round(($pickfigx2-$pickfigx1)/5.0/10.0)]; # split to 5 parts and get the nearest 10-based number
    set intervalx [expr $intervalx*10.0]
    set intervalmarkx [expr ($pickfigx2-$pickfigx1)/5.0-$intervalx]
    for {set xx [expr int(ceil($pickfigx1))]} {$xx<=[expr int($pickfigx2)]} {incr xx 1} {
	# in the case of small-limited
        if {$intervalx<=0.0} {
	    if {$intervalmarkx>=2.5} {set intervalx 2.0;}
	    if {$intervalmarkx<=2.5} {set intervalx 1.0;}
	}
	if {fmod($xx,$intervalx)==0.0} {
	    set line [waveallxy2poixy $xx $pickfigy1 $pickfigx1 $pickfigy1 $pickfigxscale $pickfigyscale $pickfigbd $pickfigfy]
	    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
	    set lx [lindex [split $line " "] 0]; set ly [lindex [split $line " "] 1]
	    .topshowallpicks.picksfigure creat oval [expr $lx-2] [expr $ly-2] [expr $lx+2] [expr $ly+2] -tags [list line] -fill red -outline red
	    .topshowallpicks.picksfigure creat text $lx [expr $ly+2] -text $xx -tags [list line] -fill blue -anchor n -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;;
	}
    }

    set intervaly [expr round(($pickfigy2-$pickfigy1)/5.0/10.0)]; # split to 5 parts and get the nearest 10-based number
    set intervaly [expr $intervaly*10.0]
    set intervalmarky [expr ($pickfigy2-$pickfigy1)/5.0-$intervaly]
    for {set yy [expr int(ceil($pickfigy1))]} {$yy<=[expr int($pickfigy2)]} {incr yy 1} {
	# in the case of small-limited
        if {$intervaly<=0.0} {
	    if {$intervalmarky>=2.5} {set intervaly 2.0;}
	    if {$intervalmarky<=2.5} {set intervaly 1.0;}
	}
	if {fmod($yy,$intervaly)==0.0} {
	    set line [waveallxy2poixy $pickfigx1 $yy $pickfigx1 $pickfigy1 $pickfigxscale $pickfigyscale $pickfigbd $pickfigfy]
	    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
	    set lx [lindex [split $line " "] 0]; set ly [lindex [split $line " "] 1]
	    .topshowallpicks.picksfigure creat oval [expr $lx-2] [expr $ly-2] [expr $lx+2] [expr $ly+2] -tags [list line] -fill red -outline red
	    .topshowallpicks.picksfigure creat text [expr $lx-2] $ly -text $yy -tags [list line] -fill blue -angle 90 -anchor s -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*; # -font fixed; "-*-times-medium-r-normal-*-18-*-*-*-*-*-iso8859-1";
	}
    }
    #### plot pick dat
    set datnum [llength $pickdat]
    for {set i 0} {$i<$datnum} {incr i 2} {
	set tt [lindex $pickdat $i]; set dist [lindex $pickdat [expr $i+1]]
	set line [waveallxy2poixy $dist $tt $pickfigx1 $pickfigy1 $pickfigxscale $pickfigyscale $pickfigbd $pickfigfy]
	set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
	set lx [lindex [split $line " "] 0]; set ly [lindex [split $line " "] 1]
	.topshowallpicks.picksfigure creat oval [expr $lx-0.5] [expr $ly-0.5] [expr $lx+0.5] [expr $ly+0.5] -tags [list dots] -fill black
    }
}
proc checkpicksub {channel} {
    global plotdotsavailablesw

    set swexist [winfo exists .topshowallpicks]
    if {$swexist==1} {
	plotallpicks $channel
	set plotdotsavailablesw 1
	plotcurdots $channel
    }
}
#################################################
proc waveallini {channel} {
    global dtname sfnum
    global sfgcarc sfdist
    global sacdatwinb sacdatwine
    global shift1 shift2
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global vibnum vibdist


    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}

    $canvasname delete line
    set waveallbd 25.0

    set waveallx1($channel) $sacdatwinb($channel); set waveallx2($channel) $sacdatwine($channel)
    set shiftdy [expr ($sfdist($sfnum)-$sfdist(1))/20.0];
    # In the case of all station have a equal dist
    if {$shiftdy==0.0} {set shiftdy 10.0}
    set waveally1($channel) [expr $sfdist(1)-$shiftdy];
    set waveally2($channel) [expr $sfdist($sfnum)+$shiftdy]
    set waveallxscale($channel) [expr ($wavefx-$waveallbd*2.0)/($waveallx2($channel)-$waveallx1($channel))]
    set waveallyscale($channel) [expr ($wavefy-$waveallbd*2.0)/($waveally2($channel)-$waveally1($channel))]

    set line [waveallxy2poixy $waveallx1($channel) $waveally1($channel) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx1 [lindex [split $line " "] 0]; set ly1 [lindex [split $line " "] 1]
    set line [waveallxy2poixy $waveallx2($channel) $waveally1($channel) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx2 [lindex [split $line " "] 0]; set ly2 [lindex [split $line " "] 1]
    $canvasname create line $lx1 $ly1 $lx2 $ly2 -width 2 -tags [list line] -fill blue
    set lxmarkx [expr ($lx1+$lx2)/2.0];

    set line [waveallxy2poixy $waveallx1($channel) $waveally2($channel) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx2 [lindex [split $line " "] 0]; set ly2 [lindex [split $line " "] 1]
    $canvasname create line $lx1 $ly1 $lx2 $ly2 -width 2 -tags [list line] -fill blue
    set lymarky [expr ($ly1+$ly2)/2.0];
    if {$channel=="BHZ"} {set titletext "UD or Z";}
    if {$channel=="BHE"} {set titletext "EW or T";}
    if {$channel=="BHN"} {set titletext "NS or R";}
    $canvasname create text $lxmarkx $waveallbd           -text $titletext -tags [list line] -fill red            -anchor s -font -*-times-bold-r-normal-*-*-150-*-*-*-*-*-*;
    $canvasname create text $lxmarkx [expr $wavefy-13]    -text "Time(s)"  -tags [list line] -fill blue           -anchor n -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
    $canvasname create text 13 $lymarky                   -text "Dist(km)" -tags [list line] -fill blue -angle 90 -anchor s -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;

    set intervalx [expr round(($waveallx2($channel)-$waveallx1($channel))/5.0/10.0)]; # split to 5 parts and get the nearest 10-based number
    set intervalx [expr $intervalx*10.0]
    set intervalmarkx [expr ($waveallx2($channel)-$waveallx1($channel))/5.0-$intervalx]
    for {set xx [expr int(ceil($waveallx1($channel)))]} {$xx<=[expr int($waveallx2($channel))]} {incr xx 1} {
	# in the case of small-limited
        if {$intervalx<=0.0} {
	    if {$intervalmarkx>=2.5} {set intervalx 2.0;}
	    if {$intervalmarkx<=2.5} {set intervalx 1.0;}
	}
	if {fmod($xx,$intervalx)==0.0} {
	    set line [waveallxy2poixy $xx $waveally1($channel) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
	    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
	    set lx [lindex [split $line " "] 0]; set ly [lindex [split $line " "] 1]
	    $canvasname creat oval [expr $lx-2] [expr $ly-2] [expr $lx+2] [expr $ly+2] -tags [list line] -fill red -outline red
	    $canvasname creat text $lx [expr $ly+2] -text $xx -tags [list line] -fill blue -anchor n -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*;
	}
    }

    set intervaly [expr round(($waveally2($channel)-$waveally1($channel))/5.0/10.0)]; # split to 5 parts and get the nearest 10-based number
    set intervaly [expr $intervaly*10.0]
    set intervalmarky [expr ($waveally2($channel)-$waveally1($channel))/5.0-$intervaly]
    for {set yy [expr int(ceil($waveally1($channel)))]} {$yy<=[expr int($waveally2($channel))]} {incr yy 1} {
	# in the case of small-limited
        if {$intervaly<=0.0} {
	    if {$intervalmarky>=2.5} {set intervaly 2.0;}
	    if {$intervalmarky<=2.5} {set intervaly 1.0;}
	}
	if {fmod($yy,$intervaly)==0.0} {
	    set line [waveallxy2poixy $waveallx1($channel) $yy $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
	    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
	    set lx [lindex [split $line " "] 0]; set ly [lindex [split $line " "] 1]
	    $canvasname creat oval [expr $lx-2] [expr $ly-2] [expr $lx+2] [expr $ly+2] -tags [list line] -fill red -outline red
	    $canvasname creat text [expr $lx-2] $ly -text $yy -tags [list line] -fill blue -angle 90 -anchor s -font -*-TimesRoman-bold-r-normal-*-*-90-*-*-*-*-*-*; # -font fixed; "-*-times-medium-r-normal-*-18-*-*-*-*-*-iso8859-1";
	}
    }

    #set line [expr ($sfdist($sfnum)-$sfdist(1))/$vibnum]
    
}
#################################################
proc wavezoomini {} {
    global dtname sfnum
    global sfgcarc sfdist
    global shift1 shift2
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale

    .wavezoom delete wave statext tmark tmarkrange synthtmark ruler line movelabel
    .wavezoom delete waveenvelope waveinstfreq
    .wavezoom configure -scrollregion [list 0 0 1000 100]
    set wavezoomfx 1000.0
    #set singlewaveheight 100;
    set wavezoombd 0
    #set wavezoomx1 [expr $sacdatbval+$shift1]; set wavezoomx2 [expr $sacdatbval+$shift2]
    #set wavezoomy1 0; set wavezoomy2 [expr $singlewaveheight*($sfnum+1)]
    #set wavezoomxscale [expr ($wavezoomfx-$wavezoombd*2.0)/($wavezoomx2-$wavezoomx1)]
    #set wavezoomyscale 1.0; #unused
}

#################################################
proc waveallxy2poixy {wx wy wx1 wy1 xscale yscale bd wheight} {
    set poix [expr ($wx-$wx1)*$xscale+$bd];
    set poiy [expr $wheight-(($wy-$wy1)*$yscale+$bd)];
    return "$poix $poiy";
}
proc waveallpoixy2xy {poix poiy wx1 wy2 xscale yscale bd} {
    set x [expr $wx1+($poix-$bd)/$xscale]; #set x [format "%8.2f" $x]
    set y [expr $wy2-($poiy-$bd)/$yscale]; #set y [format "%6.2f" $y]
    return "$x $y"
}
#################################################
proc wavezoomx2poix {wx wx1 xscale bd} {
    set poix [expr ($wx-$wx1)*$xscale+$bd];
    return "$poix";
}
#################################################
proc wavezoompoix2x {poix wx1 xscale bd} {
    set x [expr $wx1+($poix-$bd)/$xscale]; #set x [format "%10.2f" $x]
    return "$x";
}
#################################################
proc formatsub1 {a} {set a [format "%.1f" $a]; return $a}
proc formatsub2 {a} {set a [format "%.2f" $a]; return $a}
proc formatsub3 {a} {set a [format "%.3f" $a]; return $a}
proc formatsub4 {a} {set a [format "%.4f" $a]; return $a}
proc formatsub94 {a} {set a [format "%9.4f" $a]; return $a}
proc formatsub84 {a} {set a [format "%8.4f" $a]; return $a}
proc formatsub51 {a} {set a [format "%5.1f" $a]; return $a}
#################################################
proc loadchannel {channel} {
    global loadpresw loadstabz loadstabe loadstabn
    global datdir outdir evtidcur
    global sacfilenum
    global sfnum
    global evtlon evtlat evtdep evtdate evttime evtomark
    global evtdatesac evttimesac evtomarksac evtlonsac evtlatsac evtdepsac
    global evtdatesd  evttimesd  evtomarksd  evtlonsd  evtlatsd  evtdepsd evtlonrms evtlatrms evtdeprms
    global evtdatefin evttimefin evtomarkfin evtlonfin evtlatfin evtdepfin
    global evtdepmax evtdephand
    global chanidcur
    global isrelocate
    global channelsw
    global evtomarkvalsw

    if {$loadpresw==0} {return;}
    if {$channel=="BHZ" && $channelsw=="H"} {tk_messageBox -type ok -message "This channel is not available";return;}
    if {$channel=="BHE" && $channelsw=="V"} {tk_messageBox -type ok -message "This channel is not available";return;}
    if {$channel=="BHN" && $channelsw=="V"} {tk_messageBox -type ok -message "This channel is not available";return;}
    if {$channel=="BHZ"} {set infofile .zinfofile; set infonew1 .zinfonew1; set infonew2 .zinfonew2; set loadstabz [expr $loadstabz+1]; set loadstasw $loadstabz;}
    if {$channel=="BHE"} {set infofile .einfofile; set infonew1 .einfonew1; set infonew2 .einfonew2; set loadstabe [expr $loadstabe+1]; set loadstasw $loadstabe;}
    if {$channel=="BHN"} {set infofile .ninfofile; set infonew1 .ninfonew1; set infonew2 .ninfonew2; set loadstabn [expr $loadstabn+1]; set loadstasw $loadstabn;}
    $infofile configure -text ""; $infonew1 configure -text ""; $infonew2 configure -text ""
    set chanidcur $channel

    set ctrlsw [getsacinfo $channel]; if {$ctrlsw==1} {return}
    setstafigpara_plotsta_evtfree 1;
    ###puts [exec date]
    ###puts "loadstation"
    #puts [clock milliseconds]
    loadstation $channel $loadstasw
    #puts [clock milliseconds]
    ###puts [exec date]
    if {$sfnum != $sacfilenum} {
	set message "Station amount in exported results and those of SACfiles data do not match with each other.\n$sfnum vs $sacfilenum.\nMaybe some data has been added or removed.\nDo you wish to renew?";
	set chooseval [tk_messageBox -type yesno -message $message];
	if {$chooseval=="no"} {return;}
	# If yes, then renew
	renewstation $channel
	#tk_messageBox -type ok -message "Station number ($sfnum) does not match SACfile ($sacfilenum), please check and confirm"
	#return
    }
    #################################
    # loadsw=0, from SACfile; loadsw=1, from EVENT.RLOCfile; loadsw=2 from EVENT.FINALfile
    # loadevtdatetime from SACfile
    loadevtdatetime $channel 0
    set evtdate $evtdatesac
    set evttime $evttimesac
    set evtomark $evtomarksac
    set infotext "$evtdate $evttime"

    # loadevtpos from SACfile
    set errsw [loadevtpos $channel 0]; if {$errsw==1} {return}
    set evtlon $evtlonsac
    set evtlat $evtlatsac
    set evtdep $evtdepsac
    .stafigure delete evt.1
    markevtpos $evtlon $evtlat 1

    set evtlontext [format "%9.4f" $evtlon]
    set evtlattext [format "%8.4f" $evtlat]
    set evtdeptext [format "%5.1f" $evtdep]
    set infotext "$infotext $evtlontext $evtlattext $evtdeptext"
    $infofile configure -text "$infotext"

    set evtinfosw "SACfile"
    ##############################################
    #load eventinfomartion from EVENT.RLOCfile

    set eventfile "$outdir/$evtidcur.$channel.EVENT.RLOC";
    if {[file exists $eventfile]==1} {
	set isrelocate 1
	# load evtdate evttime evtomark, save in evtdatesd  evttimesd  evtomarksd
	loadevtdatetime $channel 1
	set evtdate $evtdatesd
	set evttime $evttimesd
	set evtomark $evtomarksd

	# load evtlon evtlat evtdep, save in ......
	loadevtpos $channel 1
	set evtlon $evtlonsd
	set evtlat $evtlatsd
	set evtdep $evtdepsd
	.stafigure delete evt.2
	markevtpos $evtlonsd $evtlatsd 2
	.stafigure delete evt.3
	markevtpos $evtlonrms $evtlatrms 3

	# Show message in label
	set evtlontext [format "%9.4f" $evtlonsd]
	set evtlattext [format "%8.4f" $evtlatsd]
	set evtdeptext [format "%5.1f" $evtdepsd]
        set infotext "$evtdate $evttime $evtlontext $evtlattext $evtdeptext"
	$infonew1 configure -text "$infotext"
	set evtlontext [format "%9.4f" $evtlonrms]
	set evtlattext [format "%8.4f" $evtlatrms]
	set evtdeptext [format "%5.1f" $evtdeprms]
        set infotext "$evtdate $evttime $evtlontext $evtlattext $evtdeptext"
	#$infonew2 configure -text "$infotext"
	set evtinfosw "EVTRLOCfile"
    } else {
	set isrelocate 0
    }
    ##############################################
    #load eventinfomartion from EVENT.FINALfile

    set eventfile "$outdir/$evtidcur.$channel.EVENT.FINAL";
    if {[file exists $eventfile]==1} {
	# load evtdate evttime evtomark, save in evtdatefin  evttimefin  evtomarkfin
	loadevtdatetime $channel 2
	set evtdate $evtdatefin
	set evttime $evttimefin
	set evtomark $evtomarkfin

	# load evtlon evtlat evtdep, save in ......
	loadevtpos $channel 2
	set evtlon $evtlonfin
	set evtlat $evtlatfin
	set evtdep $evtdepfin
	.stafigure delete evt.4
	markevtpos $evtlon $evtlat 4

	# Show message in label
	set evtlontext [format "%9.4f" $evtlon]
	set evtlattext [format "%8.4f" $evtlat]
	set evtdeptext [format "%5.1f" $evtdep]
        set infotext "$evtdate $evttime $evtlontext $evtlattext $evtdeptext"
	$infonew2 configure -text "$infotext"
	set evtinfosw "EVTFINALfile"
    }

    if {$evtdep>$evtdepmax} {
	set evtdepmax $evtdep
	.textevdp insert 1.0 $evtdepmax
	.scaleevtdep configure -to $evtdepmax
    }
    ## set evdp to evtdepscale 
    set evtdephand $evtdep
    if {$evtinfosw=="SACfile"} {set vibcolor magenta; set temnum 1}
    if {$evtinfosw=="EVTRLOCfile"} {set vibcolor blue; set temnum 2}
    if {$evtinfosw=="EVTFINALfile"} {set vibcolor orange; set temnum 4}
    drawgcarc $evtlon $evtlat $vibcolor $temnum
    caldistindeg $evtlon $evtlat
    sortbydist
    set evtomarkvalsw 1
    ##########################################
    ###puts [exec date]
    ###puts "readsacdat"
    set errsw [readsacdat $channel]; if {$errsw==1} {return}
    ###puts [exec date]
    waveallini $channel
    wavezoomini
    ###puts [exec date]
    ###puts "plotwaveinall"
    #puts [clock milliseconds]
    plotwaveinall $channel
    #puts [clock milliseconds]

    ###puts [exec date]
    plottmarkinwaveall $channel;
    plotcurve $channel;
    plottmarkrangeinwaveall $channel;
    ###puts [exec date]
    ###puts "synthtraveltime"
    synthtraveltime
    ###puts [exec date]
    plotsynthtmarkinwaveall $channel

    ## reorder other channel waveform profile if appear
    set swinwaveallz [.waveallz find withtag wave]
    set swinwavealle [.wavealle find withtag wave]
    set swinwavealln [.wavealln find withtag wave]

    if {$channel=="BHZ"} {
	if {$swinwavealle!=""} {
	    waveallini BHE
	    reorderwaveinall BHE
	    plotsynthtmarkinwaveall BHE
	}
	if {$swinwavealln!=""} {
	    waveallini BHN
	    reorderwaveinall BHN
	    plotsynthtmarkinwaveall BHN
	}
    }
    if {$channel=="BHE"} {
	if {$swinwaveallz!=""} {
	    waveallini BHZ
	    reorderwaveinall BHZ
	    plotsynthtmarkinwaveall BHZ
	}
	if {$swinwavealln!=""} {
	    waveallini BHN
	    reorderwaveinall BHN
	    plotsynthtmarkinwaveall BHN
	}
    }
    if {$channel=="BHN"} {
	if {$swinwaveallz!=""} {
	    waveallini BHZ
	    reorderwaveinall BHZ
	    plotsynthtmarkinwaveall BHZ
	}
	if {$swinwavealle!=""} {
	    waveallini BHE
	    reorderwaveinall BHE
	    plotsynthtmarkinwaveall BHE
	}
    }
    ## fileindicator
    fileindicator $channel

    # Check picks
    checkpicksub $channel
    # Check current curve
    zoomcurvefigureinipre $channel
}
proc updatewavezoom {} {
    global chanidcur
    global timewinleft timewinright

    if {[.wavezoom find withtag wave]==""} {return}
    waveallsub $timewinleft $timewinright $chanidcur;
}
#############################
proc waveallsub {curxl curxr chanid} {
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale waveallbd
    global sacdatwinb sacdatwine
    global wavezoomfx wavezoombd
    global wavezoomx1 wavezoomx2 wavezoomxscale
    global chanidcur

    set timel $curxl; set timer $curxr;
    if {$curxr<$curxl} {set timel $curxr; set timer $curxl;}
    #### If out of boundary the set to boundary
    set tt [waveallpoixy2xy $timel 0 $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
    set wavezoomx1 [lindex $tt 0]; if {$wavezoomx1<$sacdatwinb($chanid)} {set wavezoomx1 $sacdatwinb($chanid)}
    set tt [waveallpoixy2xy $timer 0 $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
    set wavezoomx2 [lindex $tt 0]; if {$wavezoomx2>$sacdatwine($chanid)} {set wavezoomx2 $sacdatwine($chanid)}

    ##### If select-area is totally located in blank area,return
    if {$wavezoomx1>$sacdatwine($chanid)} {return}
    if {$wavezoomx2<$sacdatwinb($chanid)} {return}

    #wavezoomini; If add this line, re-select from the waveall will move wavezoom yview to 0.0
    set wavezoomfx 1000.0
    set wavezoomxscale [expr ($wavezoomfx-$wavezoombd*2.0)/($wavezoomx2-$wavezoomx1)]
    ###puts [exec date]
    ###puts "plotwavezoom"
    #puts [clock milliseconds]
    plotwavezoom $chanid
    #puts [clock milliseconds]

    ###puts [exec date]
    set chanidcur $chanid
    plottmarkinwavezoom
    plottmarkrangeinwavezoom
    plotscaleinwavezoom .wavezoom
    ###puts [exec date]
    ###puts "plotsynthtmarkinwavezoom"
    plotsynthtmarkinwavezoom .wavezoom
    ###puts [exec date]
    # Check picks
    checkpicksub $chanid
    # Check current curve
    zoomcurvefigureinipre $chanid
}
proc filelabelsub {labelname} {
    global outdir evtidcur

    if {$evtidcur==""} {return;}
    set ll [string length $labelname];
    set lastindex [expr $ll-1];
    set channel [string range $labelname 6 8];
    set lastpart [string range $labelname 9 $lastindex];

    if {$lastpart=="PKTRLOC"}  {set FILE "$outdir/$evtidcur.$channel.PICK.FORLOC";}
    if {$lastpart=="PKTFINAL"} {set FILE "$outdir/$evtidcur.$channel.PICK.FINAL";}
    if {$lastpart=="EVTRLOC"}  {set FILE "$outdir/$evtidcur.$channel.EVENT.RLOC";}
    if {$lastpart=="EVTFINAL"} {set FILE "$outdir/$evtidcur.$channel.EVENT.FINAL";}
    if {$lastpart=="EVTSACF"}  {set FILE "$outdir/$evtidcur.$channel.EVENT.SACF";}
    if {$lastpart=="STATION"}  {set FILE "$outdir/$evtidcur.$channel.STATION.PS";}
    if {$lastpart=="WAVE"}     {set FILE "$outdir/$evtidcur.$channel.WAVE.PS";}
    if {$lastpart=="PKTCMP"}   {set FILE "$outdir/$evtidcur.$channel.COMPARISON.PS";}

    if {[file exists $FILE]==1} {exec rm $FILE; $labelname configure -bg lightgray}

    fileindicator BHZ;fileindicator BHE;fileindicator BHN
}
#################################################
proc bindsel {} {
    global bp1 bp2
    global datdir evtdir
    global evtidcur evthourdir outdir
    global stafx stafy
    global sc swidth sheight bd
    global xscale yscale x1 x2 y1 y2
    global evtlon evtlat evtdep evtdate evttime
    global sfnum sacfilenum dtname sfdist sfmansw
    global sacdatwinb sacdatwine sacdat dt
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global zoomdragcontrol zoompress2control allpress2control
    global evtdephand
    global channelsw rotatesw
    global loadpresw loadstabz loadstabe loadstabn
    global chanidcur
    global synthtmarkinwavezoomsw synthtmarkinwaveallsw
    global tmarkcurvesw
    global distmeasure evtcrosspress evterrpress
    global cursorcross
    global readthreadnum calculatethreadnum
    global plotdotsavailablesw
    global keyxnum
    global zoomupdnavailable
    global pickcurlimt1 pickcurlimt2 pickcurlimd1 pickcurlimd2
    global stwin ltwin stltlim
    global timewinleft timewinright
    global forcesw forcewin
    global methodsw
    global distlim

    set distlim 100
    set methodsw "STALTA"
    set forcesw "Yes"; set forcewin 5.0
    set chanidcur ""
    set evtidcur ""
    set distmeasure 0; set evtcrosspress 0
    set cursorcross hidden
    set loadpresw 0
    set synthtmarkinwavezoomsw 0
    set synthtmarkinwaveallsw(BHZ) 1;
    set synthtmarkinwaveallsw(BHE) 1;
    set synthtmarkinwaveallsw(BHN) 1;
    set tmarkcurvesw(BHZ) 0; set tmarkcurvesw(BHE) 0; set tmarkcurvesw(BHN) 0;
    set plotdotsavailablesw 0

    set keyxnum 0

    set zoomupdnavailable 1
    #bind all <Enter> {focus %W}

    bind .labelBHZPKTRLOC  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHEPKTRLOC  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHNPKTRLOC  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHZPKTFINAL <Double-Button-1> {filelabelsub %W;}
    bind .labelBHEPKTFINAL <Double-Button-1> {filelabelsub %W;}
    bind .labelBHNPKTFINAL <Double-Button-1> {filelabelsub %W;}
    bind .labelBHZEVTRLOC  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHEEVTRLOC  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHNEVTRLOC  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHZEVTFINAL <Double-Button-1> {filelabelsub %W;}
    bind .labelBHEEVTFINAL <Double-Button-1> {filelabelsub %W;}
    bind .labelBHNEVTFINAL <Double-Button-1> {filelabelsub %W;}
    bind .labelBHZEVTSACF  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHEEVTSACF  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHNEVTSACF  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHZSTATION  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHESTATION  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHNSTATION  <Double-Button-1> {filelabelsub %W;}
    bind .labelBHZWAVE     <Double-Button-1> {filelabelsub %W;}
    bind .labelBHEWAVE     <Double-Button-1> {filelabelsub %W;}
    bind .labelBHNWAVE     <Double-Button-1> {filelabelsub %W;}
    bind .labelBHZPKTCMP   <Double-Button-1> {filelabelsub %W;}
    bind .labelBHEPKTCMP   <Double-Button-1> {filelabelsub %W;}
    bind .labelBHNPKTCMP   <Double-Button-1> {filelabelsub %W;}


    bind .lshiftmark <ButtonRelease-1> {
	tk_popup .menushiftmark %X %Y
    }

    bind .llocationpara <ButtonRelease-1> {
	tk_popup .menulocationpara %X %Y
    }
    bind .labelwaveH <ButtonRelease-1> {
	tk_popup .menuwaveheight %X %Y
    }

    bind .lshowallpicks <ButtonRelease-1> {
	.lshowallpicks configure -bg orange
	set plotdotsavailablesw 0
	set swexist [winfo exists .topshowallpicks]
	if {$swexist==0} {
	    set topname [toplevel .topshowallpicks];
	    wm title $topname "Consistency Check"
	    frame $topname.frametop
	    label $topname.distlable  -justify center -text "Radius:"
	    radiobutton $topname.dist1 -text "100km"  -value "100"   -variable distlim  -command {if {$chanidcur!=""} {checkpicksub $chanidcur}}
	    radiobutton $topname.dist2 -text "200km"  -value "200"   -variable distlim  -command {if {$chanidcur!=""} {checkpicksub $chanidcur}}
	    radiobutton $topname.dist3 -text "500km"  -value "500"   -variable distlim  -command {if {$chanidcur!=""} {checkpicksub $chanidcur}}
	    radiobutton $topname.dist4 -text "ALL"    -value "ALL"   -variable distlim  -command {if {$chanidcur!=""} {checkpicksub $chanidcur}}
	    #radiobutton $topname.dist5 -text "" -value "input" -variable distlim
	    #text $topname.disttext -width 3 -height 1.2 -bg white -fg red
	    pack $topname.distlable $topname.dist1 $topname.dist2 $topname.dist3 $topname.dist4 -in $topname.frametop -side left
	    canvas $topname.picksfigure -background white -width $wavefx -height $wavefy -closeenough 1
	    pack $topname.frametop $topname.picksfigure -in $topname -side top
	    loadbindfortopshowallpicks
	}
	if {$chanidcur==""} {return}
	plotallpicks $chanidcur
	set plotdotsavailablesw 1
	plotcurdots $chanidcur
    }
    bind .labelzoomcurve <ButtonRelease-1> {
	.labelzoomcurve configure -bg orange
	set swexist [winfo exists .topzoomcurve]
	if {$swexist==0} {
	    toplevel .topzoomcurve
	    wm title .topzoomcurve "Traveltime curve"
	    canvas .topzoomcurve.curvefigure -background white -width $wavefx -height $wavefy -closeenough 1
	    pack .topzoomcurve.curvefigure -in .topzoomcurve -side top
	    loadbindfortopzoomcurve
	}
	if {$chanidcur==""} {return}
	getpickcurrentevtdatrange $chanidcur
	zoomcurvefigureini $chanidcur $pickcurlimt1 $pickcurlimt2 $pickcurlimd1 $pickcurlimd2
	plotzoomcurve $chanidcur
    }
    bind .lshowallpicks <ButtonPress-1>   {
	.lshowallpicks configure -bg lightgrey
    }
    bind .labelzoomcurve <ButtonPress-1>   {
	.labelzoomcurve configure -bg lightgrey
    }
    bind .labelapb <ButtonPress-1>   {
	.labelapb configure -bg lightgrey
    }
    bind .labelapb <ButtonRelease-1>   {
	.labelapb configure -bg orange
	if {[.wavezoom find withtag wave]==""} {return}
	loadandpreparedatandcallpick
    }
    bind .labelapp <ButtonPress-1>  {
	.labelapp configure -bg lightgrey
    }
    bind .labelapp <ButtonRelease-1> {
	.labelapp configure -bg orange
	set app [toplevel .menuapp]
	wm title $app "AutoPicker Parameter"

	frame $app.frame0
	frame $app.frame1
	frame $app.framespace1
	frame $app.frame2
	frame $app.frame3
	frame $app.frame4
	frame $app.framespace2
	frame $app.frame5
	pack $app.frame0 $app.frame1 $app.framespace1 $app.frame2 $app.frame3 $app.frame4 $app.framespace2 $app.frame5 -side top -fill x

	label $app.marktext -text "Force picks to follow predictions?"
	radiobutton $app.forceyes -variable forcesw -text "Yes" -value "Yes"
	radiobutton $app.forceno  -variable forcesw -text "No"  -value "No"
	pack $app.marktext $app.forceyes $app.forceno -side left -in $app.frame0

	if {$forcesw=="Yes"} {set colortem "white";     set statustem "normal"}
	if {$forcesw=="No"}  {set colortem "lightgray"; set statustem "disable"}
	label $app.timeshiftlabel -text "If force picks to follow predictions,\nthen set a shift for time window:"
	text $app.timeshifttext -width 8 -height 1.3 -bg $colortem -fg blue
	     $app.timeshifttext insert 1.0 $forcewin
	     $app.timeshifttext configure  -state $statustem
	label $app.temtext1 -text "secs"
	pack $app.timeshifttext $app.timeshiftlabel -side right -in $app.frame1

	label $app.textspace1 -text "     "
	pack $app.textspace1 -side right -in $app.framespace1

	label $app.method -text "Method for AutoPicker:"
	radiobutton $app.stalta -variable methodsw -text "STA/LTA" -value "STALTA"
	radiobutton $app.aic    -variable methodsw -text "AIC"     -value "AIC"
	pack $app.method $app.stalta $app.aic -side left -in $app.frame2


	if {$methodsw=="STALTA"} {set colortem "white";     set statustem "normal"}
	if {$methodsw=="AIC"}    {set colortem "lightgray"; set statustem "disable"}
	label $app.llta -text "   LTA-win:"
	text $app.tlta -width 8 -height 1.1 -bg $colortem -fg blue
	     $app.tlta insert 1.0 $ltwin
	     $app.tlta configure  -state $statustem
	label $app.temtext3 -text "secs"
	pack $app.tlta $app.llta -side right -in $app.frame3

	label $app.lsta -text "STA-win:"
	text $app.tsta -width 8 -height 1.1 -bg $colortem -fg blue
	     $app.tsta insert 1.0 $stwin
	     $app.tsta configure  -state $statustem
	label $app.temtext2 -text "secs  "
	pack $app.tsta $app.lsta -side right -in $app.frame3



	label $app.lrat -text "lta-to-sta-ratio-limit:"
	text $app.trat -width 8 -height 1.1 -bg white -fg blue
	     $app.trat insert 1.0 $stltlim
	#pack $app.trat $app.lrat -side right -in $app.frame4

	label $app.textspace2 -text "     "
	pack $app.textspace2 -side right -in $app.framespace2

	frame $app.framebutton
	button $app.ok     -height 1 -width 6 -bg orange -text "OK"     -command {
		set forcewin [$app.timeshifttext get 1.0 end]; set forcewin [string trim $forcewin];
		set stwin    [$app.tsta get 1.0 end];          set stwin    [string trim $stwin];
		set ltwin    [$app.tlta get 1.0 end];          set ltwin    [string trim $ltwin];
		set stltlim  [$app.trat get 1.0 end];          set stltlim  [string trim $stltlim];
		destroy $app
	}
	button $app.cancel -height 1 -width 6 -bg orange -text "Cancel" -command {destroy $app}
	pack $app.cancel $app.ok -side right -in $app.framebutton
	pack $app.framebutton -side top -in $app.frame5


	focus $app
	grab $app
	wm transient $app .

	#raise $app
	#tkwait window $app
	bind $app.forceyes <Button-1> {$app.timeshifttext configure -bg white     -state normal}
	bind $app.forceno  <Button-1> {$app.timeshifttext configure -bg lightgray -state disable}
	bind $app.stalta   <Button-1> {
	    $app.tsta configure -bg white -state normal;
	    $app.tlta configure -bg white -state normal;
	}
	bind $app.aic      <Button-1> {
	    $app.tsta configure -bg lightgray  -state disable;
	    $app.tlta configure -bg lightgray  -state disable;
	}
    }

    bind .scaleevtdep <Button-4> {set evtdephand [expr $evtdephand-1];}
    bind .scaleevtdep <Button-5> {set evtdephand [expr $evtdephand+1];}

    bind .evtid <ButtonRelease-1> {
	set evtidcur [.evtid get [.evtid curselection]]
	exec rm -rf ./swap
	exec mkdir ./swap
	#if {[file exists ./swap]==0} {
	    #exec mkdir swap
	#} else {
	    #exec rm -rf ./swap; This line cannot work, unknown reason
	#}

	## read band pass filter parameter
	set bp1 [.tbp1 get 1.0 end]; set bp1 [string trim $bp1];
	set bp2 [.tbp2 get 1.0 end]; set bp2 [string trim $bp2];


	## read SAC commands
	set tt [.tpre get 1.0 end]
	set linenum [exec echo $tt | wc]
	set linenum [string trim $linenum];regsub -all {[[:blank:]]+} $linenum " " linenum
	set linenum [lindex $linenum 0];
	set linenum [expr $linenum-1]
	set precmds {}
	for {set i 1} {$i<=$linenum} {incr i 1} {
	    set tt [.tpre get $i.0 $i.end]
	    set precmds "$precmds\n$tt";
	}
	set filelist [glob -nocomplain $datdir/$evtidcur/*.BHZ.SAC]

	# set thread num & create thread, then set them to wait
	for {set i 0} {$i<$readthreadnum} {incr i 1} {
	    set workid($i) [thread::create {
		package require Ttrace;
		proc preparesac {acmds} {
		    #exec echo "$acmds" | sac
		    set fout [open "| sac >/dev/null" w]
		    puts $fout "$acmds";
		    close $fout
		}
		set sacio [open "| sac >/dev/null" w]
		thread::wait;
	    }]
	}
	# obtain main thread id
	set mainid [thread::id]
	###puts [exec date]
	###puts "preparedat"
	if {$channelsw=="V" || $channelsw=="VH"} {
	    set filenum 0
	    foreach BHZelement $filelist {
		# Check if three channel
		set BHNelement $BHZelement
		set BHEelement $BHZelement
		set BHNelement [string replace $BHNelement end-4 end-4 "N"]
		set BHEelement [string replace $BHEelement end-4 end-4 "E"]
		###### 0-Not exist; 1-exist
		if {[file exists $BHNelement]==1 && [file exists $BHEelement]==1} {
		    set rcmds "r $BHZelement";
		    if {$bp1!=$bp2} {set bpcmds "bp co $bp1 $bp2 p 2\n";}
		    if {$bp1==$bp2} {set bpcmds "";} 
		    set wcmds "w dir swap over";
		    set envelopecmds "envelope\nw dir swap append .ENVELOPE";
		    set acmds "$rcmds\n$precmds\n$bpcmds\n$wcmds\n$envelopecmds\n";
		    #exec echo "$acmds" | sac
		    incr filenum 1;
		    set imark [expr int(fmod($filenum,$readthreadnum))]
		    #thread::send -async $workid($imark) "preparesac {$acmds\nq}"
		    thread::send -async $workid($imark) "puts \$sacio {$acmds}"
		}
	    }
	}
	if {$channelsw=="H" || $channelsw=="VH"} {
	    set filenum 0
	    foreach element $filelist {
		set BHNelement $element
		set BHEelement $element
		set BHNelement [string replace $BHNelement end-4 end-4 "N"]
		set BHEelement [string replace $BHEelement end-4 end-4 "E"]
		###### 0-Not exist; 1-exist
		if {[file exists $BHNelement]==1 && [file exists $BHEelement]==1} {
		    set rcmds "r $BHNelement $BHEelement";
		    if {$rotatesw=="YES"} {
			# Some datfile doesn't have same npts, however b is different, so using CUTERR to fill zero
			set line [saclsttcl b e npts f $BHNelement]; regsub -all {[[:blank:]]+} $line " " line;
			set begn [lindex $line 1]; set endn [lindex $line 2]; set nptsn [lindex $line 3];
			set line [saclsttcl b e npts f $BHEelement]; regsub -all {[[:blank:]]+} $line " " line;
			set bege [lindex $line 1]; set ende [lindex $line 2]; set nptse [lindex $line 3];
			set beg $begn; if {$bege<$beg} {set beg $bege;}
			set end $endn; if {$ende>$end} {set end $ende;}
			set cutcmds "cuterr fillz\n cut b $beg $end\n"
			set rotatecmd "rotate to gcp\n";
		    }
		    if {$rotatesw=="NO"}  {
			set cutcmds ""
			set rotatecmd "";
		    }
		    if {$bp1!=$bp2} {set bpcmds "bp co $bp1 $bp2 p 2\n";}
		    if {$bp1==$bp2} {set bpcmds "";} 
		    set wcmds "w dir swap over";
		    set envelopecmds "envelope\nw dir swap append .ENVELOPE";
		    set acmds "$cutcmds\n$rcmds\n$precmds\n$rotatecmd\n$bpcmds\n$wcmds\n$envelopecmds\n";
		    #exec echo "$acmds" | sac
		    incr filenum 1;
		    set imark [expr int(fmod($filenum,$readthreadnum))]
		    #thread::send -async $workid($imark) "preparesac {$acmds\nq}"
		    thread::send -async $workid($imark) "puts \$sacio {$acmds}"
		}
	    }
	}
	# send stop command to each thread
	for {set i 0} {$i<$readthreadnum} {incr i 1} {
	    #thread::send -async $workid($i) "flush \$sacio"
	    thread::send -async $workid($i) "puts \$sacio q\n"
	    thread::send -async $workid($i) "thread::release $workid($i)"
	}
	while {[llength [thread::names]] > 1} {after 200}
	###puts [exec date]
	set loadpresw 1; set loadstabz 0; set loadstabe 0; set loadstabn 0;
        fileindicator BHZ;
        fileindicator BHE;
        fileindicator BHN;
	setstafigpara_plotsta_evtfree 1
	showini
	labelini
	.waveallz delete all
	.wavealle delete all
	.wavealln delete all
	.wavezoom delete all
	###puts [exec date]
	###puts "hilbert"
	####### run hilbertsac for file in swap dir

	# if not use this, hilbert will not conduct on all expected file properly, when filenum is about great than 240*3
	if {$channelsw=="H" || $channelsw=="V"} {set muln 1;}
	if {$channelsw=="VH"} {set muln 3;}
	set filelist [glob -nocomplain swap/$evtidcur.*.SAC]
	while {[llength $filelist]< [expr $filenum*$muln]} {
	    after 500
	    set filelist [glob -nocomplain swap/$evtidcur.*.SAC]; 
	}
	# set thread num & create thread, then set them to wait
	for {set i 0} {$i<$readthreadnum} {incr i 1} {
	    set workid($i) [thread::create {package require Ttrace; thread::wait;}]
	    set fileargv($i) {}
	}
	# obtain main thread id
	set mainid [thread::id]

	set i 0
	foreach element $filelist {
	    #exec hilbertsac $element $element 1.1
	    incr i 1;
	    set imark [expr int(fmod($i,$readthreadnum))]
	    lappend  fileargv($imark) $element
	    #thread::send -async $workid($imark) "exec hilbertsac $element swap 1.1"
	    #puts $i,$imark,$workid($imark)
	}
	# hilbert, and send stop command to each thread
	for {set i 0} {$i<$readthreadnum} {incr i 1} {
	    thread::send -async $workid($i) "exec hilbertsac $fileargv($i) swap 1.1"
	    thread::send -async $workid($i) "thread::release $workid($i)"
	}
	while {[llength [thread::names]] > 1} {after 200}
	###puts [exec date]

	resetstafigview
	### if necessary, delete current dot's signs in topshow
	set swexist [winfo exists .topshowallpicks.picksfigure]
	if {$swexist==1} {.topshowallpicks.picksfigure delete curdots}
	set swexist [winfo exists .topzoomcurve.curvefigure]
	if {$swexist==1} {.topzoomcurve.curvefigure delete curve curvedot}
    }

    bind .stafigure <Control-Button-4> {%W yview scroll -2 units}
    bind .stafigure <Control-Button-5> {%W yview scroll +2 units}
    bind .stafigure <Motion> {
	set curx [%W canvasx %x];   # zero-point from top-left canvas corner
	set cury [%W canvasy %y];
	set curx [expr round($curx)]
	set cury [expr round($cury)]
	#set curx %x;               # zero-point from top-left displayarea corner
	#set cury %y
        ############################################
	#set lat [expr $y2-($cury-$bd)/$yscale];
	#set lon [expr $x1+($curx-$bd)/$xscale];
	##set lat [format "%6.2f" $lat]; # This two line cannot workable, Don't know reason
	##set lon [format "%6.2f" $lon]; # This two line cannot workable, Don't know reason
        ############################################
	set line [stafigurepoixy2xy $curx $cury $x1 $y2 $xscale $yscale $bd];
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line
	set lon [lindex [split $line " "] 0]; set lon [formatsub4 $lon]
	set lat [lindex [split $line " "] 1]; set lat [formatsub4 $lat]
	.information configure -text "Poix:$curx Poiy:$cury Pos:$lon $lat"
	.stafigure coords linex1 0 $cury [expr $curx-2] $cury
	.stafigure coords linex2 [expr $curx+2] $cury $stafx $cury
	.stafigure coords liney1 $curx 0 $curx [expr $cury-2]
	.stafigure coords liney2 $curx [expr $cury+2] $curx $stafy

	# Show event error
	set evterrpress 0
	.stafigure itemconfigure evterr -state normal

	set evtcrosspress 0
	if {$distmeasure==0} {return}
	if {$distmeasure==1} {
	    .stafigure itemconfigure linexy -state hidden
	    #.stafigure itemconfigure linex2 -state hidden
	    #.stafigure itemconfigure liney1 -state hidden
	    #.stafigure itemconfigure liney2 -state hidden
	}
	# for measuring distance
	set distline [distaztcl $lat0 $lon0 $lat $lon]
	set distline [string trim $distline];regsub -all {[[:blank:]]+} $distline " " distline
	set dist [lindex $distline 0]
	set dist [expr $dist*acos(-1.0)*6371.0/180.0]
	#set dist [format "%.3f" $dist]; # This line cannot workable, Don't know reason
	if {[info exist dist0]==1} {set dist [expr $dist+$dist0]}
	set dist [formatsub3 $dist]
	.stafigure delete tagdistance
	.stafigure create text $curx $cury -text "$dist km" -tags [list tagdistance] -anchor s -fill red
	.stafigure create line $curx0 $cury0 $curx $cury -width 1 -tags [list tagdistance] -fill red 
    }

    #bind .stafigure <1> {
	#.stafigure coords evt.0 $curx $cury
    #}
    bind .stafigure <Enter> {focus %W}
    bind .stafigure <Key-x> {
	# This process needs that the widget should firstly get focus
	set statussw [.stafigure itemcget linex1 -state]
	if {$statussw==""} {
	    .stafigure itemconfigure linexy -state hidden
	    .wavezoom itemconfigure linexy -state hidden
	    set cursorcross hidden
	} elseif {$statussw=="normal"} {
	    .stafigure itemconfigure linexy -state hidden
	    .wavezoom itemconfigure linexy -state hidden
	    set cursorcross hidden
	} elseif {$statussw=="hidden"} {
	    .stafigure itemconfigure linexy -state normal
	    .wavezoom itemconfigure linexy -state normal
	    set cursorcross normal
	}
    }
    bind .wavezoom <Key-x> {
	set statussw [.stafigure itemcget linex1 -state]
	if {$statussw==""} {
	    .stafigure itemconfigure linexy -state hidden
	    .wavezoom itemconfigure linexy -state hidden
	    set cursorcross hidden
	} elseif {$statussw=="normal"} {
	    .stafigure itemconfigure linexy -state hidden
	    .wavezoom itemconfigure linexy -state hidden
	    set cursorcross hidden
	} elseif {$statussw=="hidden"} {
	    .stafigure itemconfigure linexy -state normal
	    .wavezoom itemconfigure linexy -state normal
	    set cursorcross normal
	}
    }
    bind .stafigure <ButtonPress-1> {
	# The following line for drag
	set curX %x
	set curY %y
	%W config -cursor hand1

	set swmarkevt 1
    }
    bind .stafigure <B1-Motion> {

	set delX [expr %x-$curX]
	set delY [expr %y-$curY]
	set curX %x
	set curY %y

	loadstapara
	loadevtpara

	set londel [expr $delX/$xscale]
	set latdel [expr $delY/$yscale]
	if {[expr $x1-$londel]<-180.0 || [expr $x1-$londel]>180.0} {return}
	if {[expr $x2-$londel]<-180.0 || [expr $x2-$londel]>180.0} {return}
	if {[expr $y1+$latdel]< -90.0 || [expr $y1+$latdel]> 90.0} {return}
	if {[expr $y2+$latdel]< -90.0 || [expr $y2+$latdel]> 90.0} {return}
	set x1 [expr $x1-$londel];if {$x1<-180.0} {set x1 -180.0;}
	set x2 [expr $x2-$londel];if {$x2> 180.0} {set x2  180.0;}
	set y1 [expr $y1+$latdel];if {$y1< -90.0} {set y1  -90.0;}
	set y2 [expr $y2+$latdel];if {$y2>  90.0} {set y2   90.0;}
	#.stafigure move gcarc $delX $delY

	setstafigpara_plotsta_evtfree 0
	resetupstapara
	resetupevtpara


	set swmarkevt 0
    }
    bind .stafigure <ButtonRelease-1> {
	
	if {$swmarkevt==1} {.stafigure coords evt.0 $curx $cury}
	%W config -cursor arrow
    }
    bind .stafigure <Button-4> {
	### load & save evtparameters
	loadevtpara

	### load & save station parameters
	loadstapara

	set curx [%W canvasx %x]
	set cury [%W canvasy %y]
	zoom $curx $cury 0.7
	setstafigpara_plotsta_evtfree 0;

	### re-setup evtparameters
	resetupevtpara

	### re-setup station parameters
	resetupstapara
    }
    bind .stafigure <Button-5> {
	### load & save evtparameters
	loadevtpara

	### load & save station parameters
	loadstapara

	set curx [%W canvasx %x]
	set cury [%W canvasy %y]
	zoom $curx $cury 1.2
	setstafigpara_plotsta_evtfree 0;

	### re-setup evtparameters
	resetupevtpara

	### re-setup station parameters
	resetupstapara
    }
    bind .stafigure <Key-space> {
	resetstafigview
    }
    bind .stafigure <Double-Button-3> {
	#set ndist [expr $ndist+1]
	set curx1 [%W canvasx %x]
	set cury1 [%W canvasy %y]
	
	.stafigure delete tagdistance tagdistancepre
	set distmeasure 0
	# End measurement & reset variables
	if {[info exist dist]==1} {unset dist};
	if {[info exist dist0]==1} {unset dist0};
	if {[info exist curx0]==1} {unset curx0};
	if {[info exist cury0]==1} {unset cury0};
	#unset curx0; unset cury0
    }
    bind .stafigure <3> {
	# The following line for drag through B3
	#set curX %x
	#set curY %y
	#%W config -cursor hand1

	if {$evtcrosspress==1 || $evterrpress==1} {return;}
	# The following line for measuring distance
	###### plot polygonal line segment - fixed
	set curx [%W canvasx %x]
	set cury [%W canvasy %y]
	if {[info exist curx0]==1 && [info exist cury0]==1} {
	    .stafigure create line $curx0 $cury0 $curx $cury -width 1 -tags [list tagdistancepre] -fill red
	}
	#if {$evtcrosspress==1} {return;}
	###### for first point of new polygonal line segment - moving
	set curx0 [%W canvasx %x]
	set cury0 [%W canvasy %y]
	set line [stafigurepoixy2xy $curx0 $cury0 $x1 $y2 $xscale $yscale $bd];
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line
	set lon0 [lindex [split $line " "] 0];
	set lat0 [lindex [split $line " "] 1];
	if {[info exist dist]==1} {set dist0 $dist;}
	set distmeasure 1
    }
    bind .stafigure <Control-Button-3> {

	#if {$distmeasure!=1} {return;}
	#set curx [%W canvasx %x]
	#set cury [%W canvasy %y]
	#set line [stafigurepoixy2xy $curx $cury $x1 $y2 $xscale $yscale $bd];
	#set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line
	#set lon [lindex [split $line " "] 0];
	#set lat [lindex [split $line " "] 1];
	#.stafigure create line $curx0 $cury0 $curx $cury -width 1 -tags [list tagdistancepre] -fill black
	#set dist0 $dist
	#set curx0 $curx
	#set cury0 $cury
	#set lon0 $lon
	#set lat0 $lat

    }
    bind .stafigure <B3-Motion> {
	# This part has been replaced by <B1-Motion>, So Return;
	return;

	# if drag then delete tagdistance & set distmeasure to 0
	.stafigure delete tagdistance
	set distmeasure 0

	set delX [expr %x-$curX]
	set delY [expr %y-$curY]
	set curX %x
	set curY %y

	loadstapara
	loadevtpara

	set londel [expr $delX/$xscale]
	set latdel [expr $delY/$yscale]
	if {[expr $x1-$londel]<-180.0 || [expr $x1-$londel]>180.0} {return}
	if {[expr $x2-$londel]<-180.0 || [expr $x2-$londel]>180.0} {return}
	if {[expr $y1+$latdel]< -90.0 || [expr $y1+$latdel]> 90.0} {return}
	if {[expr $y2+$latdel]< -90.0 || [expr $y2+$latdel]> 90.0} {return}
	set x1 [expr $x1-$londel];if {$x1<-180.0} {set x1 -180.0;}
	set x2 [expr $x2-$londel];if {$x2> 180.0} {set x2  180.0;}
	set y1 [expr $y1+$latdel];if {$y1< -90.0} {set y1  -90.0;}
	set y2 [expr $y2+$latdel];if {$y2>  90.0} {set y2   90.0;}
	#.stafigure move gcarc $delX $delY

	setstafigpara_plotsta_evtfree 0
	resetupstapara
	resetupevtpara
    }
    #bind .stafigure <ButtonRelease-3> {
	#%W config -cursor arrow
    #}
    ####################
    bind .waveallz <Enter> {focus %W}
    bind .wavealle <Enter> {focus %W}
    bind .wavealln <Enter> {focus %W}
    bind .waveallz <Key-space> {waveallshowhide .waveallz}
    bind .wavealle <Key-space> {waveallshowhide .wavealle}
    bind .wavealln <Key-space> {waveallshowhide .wavealln}

    #################### Right click to view the nearest wave
    bind .waveallz <3> {
	if {[.wavezoom find withtag wave]==""} {return}
	set cury  [%W canvasy %y]
	set chanid BHZ
	if {$chanid!=$chanidcur} {return}
	set tt [waveallpoixy2xy 0 $cury $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
	set ydist [lindex $tt 1]
	set distmin 100000.0
	for {set i 1} {$i<=$sfnum} {incr i 1} {
	    if {[expr abs($sfdist($i)-$ydist)]<=$distmin} {set distmin [expr abs($sfdist($i)-$ydist)]; set selnum $i}
	}
	set yratio [expr (($sfnum+1.0)-$selnum-2.0)/$sfnum]
	.wavezoom yview moveto $yratio
	plotscaleinwavezoom .wavezoom
    }
    bind .wavealle <3> {
	if {[.wavezoom find withtag wave]==""} {return}
	set cury  [%W canvasy %y]
	set chanid BHE
	if {$chanid!=$chanidcur} {return}
	set tt [waveallpoixy2xy 0 $cury $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
	set ydist [lindex $tt 1]
	set distmin 100000.0
	for {set i 1} {$i<=$sfnum} {incr i 1} {
	    if {[expr abs($sfdist($i)-$ydist)]<=$distmin} {set distmin [expr abs($sfdist($i)-$ydist)]; set selnum $i}
	}
	set yratio [expr (($sfnum+1.0)-$selnum-2.0)/$sfnum]
	.wavezoom yview moveto $yratio
	plotscaleinwavezoom .wavezoom
    }
    bind .wavealln <3> {
	if {[.wavezoom find withtag wave]==""} {return}
	set cury  [%W canvasy %y]
	set chanid BHN
	if {$chanid!=$chanidcur} {return}
	set tt [waveallpoixy2xy 0 $cury $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
	set ydist [lindex $tt 1]
	set distmin 100000.0
	for {set i 1} {$i<=$sfnum} {incr i 1} {
	    if {[expr abs($sfdist($i)-$ydist)]<=$distmin} {set distmin [expr abs($sfdist($i)-$ydist)]; set selnum $i}
	}
	set yratio [expr (($sfnum+1.0)-$selnum-2.0)/$sfnum]
	.wavezoom yview moveto $yratio
	plotscaleinwavezoom .wavezoom
    }
    #################### Right click to view the nearest wave
    #################### Double click to switch tmark or curve
    bind .waveallz <Double-Button-1> {
	if {[.waveallz find withtag wave]==""} {return}
	if {[.waveallz find withtag tmark]==""} {
	    .waveallz delete curve
	    set tmarkcurvesw(BHZ) 1
	    plottmarkinwaveall BHZ
	} else {
	    .waveallz delete tmark
	    set tmarkcurvesw(BHZ) 0
	    plotcurve BHZ
	}
    }
    bind .wavealle <Double-Button-1> {
	if {[.wavealle find withtag wave]==""} {return}
	if {[.wavealle find withtag tmark]==""} {
	    .wavealle delete curve
	    set tmarkcurvesw(BHE) 1
	    plottmarkinwaveall BHE
	} else {
	    .wavealle delete tmark
	    set tmarkcurvesw(BHE) 0
	    plotcurve BHE
	}
    }
    bind .wavealln <Double-Button-1> {
	if {[.wavealln find withtag wave]==""} {return}
	if {[.wavealln find withtag tmark]==""} {
	    .wavealln delete curve
	    set tmarkcurvesw(BHN) 1
	    plottmarkinwaveall BHN
	} else {
	    .wavealln delete tmark
	    set tmarkcurvesw(BHN) 0
	    plotcurve BHN
	}
    }
    #################### Double click to switch tmark or curve
    #################### Press x to show wavezoom
    bind .waveallz <Key-x> {
	if {[.waveallz find withtag wave]==""} {return}
	.waveallz delete box boxtem; .wavealle delete box boxtem; .wavealln delete box boxtem
	if {$keyxnum==0} {
	    set keyxnum 1
	    set curxl [%W canvasx %x];
	    .waveallz create line [list $curxl 0 $curxl $wavefy] -width 1 -fill lightgreen -tags boxtem
	} elseif {$keyxnum==1} {
	    set keyxnum 0
	    set curxr [%W canvasx %x];
	    if {$curxl==$curxr} {return;}
	    if {[expr abs($curxr-$curxl)<10]} {tk_messageBox -type ok -message "window is too short";return}
	    .waveallz create rect $curxl 0                         $curxr $waveallbd                -fill lightgreen -outline lightgreen -stipple gray12 -tag box
	    .waveallz create rect $curxl $waveallbd                $curxr [expr $wavefy-$waveallbd]                  -outline lightgreen -stipple gray12 -tag [list box boxtem]
	    .waveallz create rect $curxl [expr $wavefy-$waveallbd] $curxr $wavefy                   -fill lightgreen -outline lightgreen -stipple gray12 -tag box
	    set chanid BHZ
	    set timewinleft $curxl; set timewinright $curxr;
	    waveallsub $curxl $curxr $chanid;
	}
    }
    bind .wavealle <Key-x> {
	if {[.wavealle find withtag wave]==""} {return}
	.waveallz delete box boxtem; .wavealle delete box boxtem; .wavealln delete box boxtem
	if {$keyxnum==0} {
	    set keyxnum 1
	    set curxl [%W canvasx %x];
	    .wavealle create line [list $curxl 0 $curxl $wavefy] -width 1 -fill lightgreen -tags boxtem
	} elseif {$keyxnum==1} {
	    set keyxnum 0
	    set curxr [%W canvasx %x];
	    if {[expr abs($curxr-$curxl)<10]} {tk_messageBox -type ok -message "window is too short";return}
	    .wavealle create rect $curxl 0                         $curxr $waveallbd                -fill lightgreen -outline lightgreen -stipple gray12 -tag box
	    .wavealle create rect $curxl $waveallbd                $curxr [expr $wavefy-$waveallbd]                  -outline lightgreen -stipple gray12 -tag [list box boxtem]
	    .wavealle create rect $curxl [expr $wavefy-$waveallbd] $curxr $wavefy                   -fill lightgreen -outline lightgreen -stipple gray12 -tag box
	    set chanid BHE
	    set timewinleft $curxl; set timewinright $curxr;
	    waveallsub $curxl $curxr $chanid;
	}
    }
    bind .wavealln <Key-x> {
	if {[.wavealln find withtag wave]==""} {return}
	.waveallz delete box boxtem; .wavealle delete box boxtem; .wavealln delete box boxtem
	if {$keyxnum==0} {
	    set keyxnum 1
	    set curxl [%W canvasx %x];
	    .wavealln create line [list $curxl 0 $curxl $wavefy] -width 1 -fill lightgreen -tags boxtem
	} elseif {$keyxnum==1} {
	    set keyxnum 0
	    set curxr [%W canvasx %x];
	    if {[expr abs($curxr-$curxl)<10]} {tk_messageBox -type ok -message "window is too short";return}
	    .wavealln create rect $curxl 0                         $curxr $waveallbd                -fill lightgreen -outline lightgreen -stipple gray12 -tag box
	    .wavealln create rect $curxl $waveallbd                $curxr [expr $wavefy-$waveallbd]                  -outline lightgreen -stipple gray12 -tag [list box boxtem]
	    .wavealln create rect $curxl [expr $wavefy-$waveallbd] $curxr $wavefy                   -fill lightgreen -outline lightgreen -stipple gray12 -tag box
	    set chanid BHN
	    set timewinleft $curxl; set timewinright $curxr;
	    waveallsub $curxl $curxr $chanid;
	}
    }
    #################### Press x to show wavezoom
    #################### Drag to show wavezoom
    bind .waveallz <ButtonPress-1> {
	set curxl [%W canvasx %x];
	set ifdrag 0
    }
    bind .waveallz <B1-Motion> {
	if {[.waveallz find withtag wave]==""} {return}
	set ifdrag 1
	set curxr [%W canvasx %x];
	.waveallz delete box boxtem; .wavealle delete box boxtem; .wavealln delete box boxtem
	.waveallz create rect $curxl 0                         $curxr $waveallbd                -fill lightgreen -outline lightgreen -stipple gray12 -tag box
	.waveallz create rect $curxl $waveallbd                $curxr [expr $wavefy-$waveallbd]                  -outline lightgreen -stipple gray12 -tag [list box boxtem]
	.waveallz create rect $curxl [expr $wavefy-$waveallbd] $curxr $wavefy                   -fill lightgreen -outline lightgreen -stipple gray12 -tag box
    }
    bind .waveallz <ButtonRelease-1> {
	if {$ifdrag==0} {return}
	# if window is too small then return
	if {[expr abs($curxr-$curxl)<10]} {tk_messageBox -type ok -message "window is too short";return}
	set chanid BHZ
	if {$ifdrag==1} {
	    set timewinleft $curxl; set timewinright $curxr;
	    waveallsub $curxl $curxr $chanid;
	    #set timel $curxl; set timer $curxr;
	    #if {$curxr<$curxl} {set timel $curxr; set timer $curxl;}
	    #set tt [waveallpoixy2xy $timel 0 $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
	    #set wavezoomx1 [lindex $tt 0]; if {$wavezoomx1<$sacdatwinb($chanid)} {set wavezoomx1 $sacdatwinb($chanid)}
	    #set tt [waveallpoixy2xy $timer 0 $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
	    #set wavezoomx2 [lindex $tt 0]; if {$wavezoomx2>$sacdatwine($chanid)} {set wavezoomx2 $sacdatwine($chanid)}

	    ###### If select-area is totally located in blank area,return
	    #if {$wavezoomx1>$sacdatwine($chanid)} {return}
	    #if {$wavezoomx2<$sacdatwinb($chanid)} {return}

	    ##wavezoomini; If add this line, re-select from the waveall will move wavezoom yview to 0.0
	    #set wavezoomfx 1000.0
	    #set wavezoomxscale [expr ($wavezoomfx-$wavezoombd*2.0)/($wavezoomx2-$wavezoomx1)]
	    #plotwavezoom $chanid
	    #set chanidcur $chanid
	    #plottmarkinwavezoom
	    #plottmarkrangeinwavezoom
	    #plotscaleinwavezoom .wavezoom
	    #plotsynthtmarkinwavezoom .wavezoom
	    ## Check picks
	    #checkpicksub $chanid
	}
    }
    #########
    bind .wavealle <ButtonPress-1> {
	set curxl [%W canvasx %x];
	set ifdrag 0
    }
    bind .wavealle <B1-Motion> {
	if {[.wavealle find withtag wave]==""} {return}
	set ifdrag 1
	set curxr [%W canvasx %x];
	.waveallz delete box boxtem; .wavealle delete box boxtem; .wavealln delete box boxtem
	.wavealle create rect $curxl 0                         $curxr $waveallbd                -fill lightgreen -outline lightgreen -stipple gray12 -tag box
	.wavealle create rect $curxl $waveallbd                $curxr [expr $wavefy-$waveallbd]                  -outline lightgreen -stipple gray12 -tag [list box boxtem]
	.wavealle create rect $curxl [expr $wavefy-$waveallbd] $curxr $wavefy                   -fill lightgreen -outline lightgreen -stipple gray12 -tag box
    }
    bind .wavealle <ButtonRelease-1> {
	if {$ifdrag==0} {return}
	# if window is too small then return
	if {[expr abs($curxr-$curxl)<10]} {tk_messageBox -type ok -message "window is too short";return}
	set chanid BHE
	if {$ifdrag==1} {
	    set timewinleft $curxl; set timewinright $curxr;
	    waveallsub $curxl $curxr $chanid;
	}
    }

    #########
    bind .wavealln <ButtonPress-1> {
	set curxl [%W canvasx %x];
	set ifdrag 0
    }
    bind .wavealln <B1-Motion> {
	if {[.wavealln find withtag wave]==""} {return}
	set ifdrag 1
	set curxr [%W canvasx %x];
	.waveallz delete box boxtem; .wavealle delete box boxtem; .wavealln delete box boxtem
	.wavealln create rect $curxl 0                         $curxr $waveallbd                -fill lightgreen -outline lightgreen -stipple gray12 -tag box
	.wavealln create rect $curxl $waveallbd                $curxr [expr $wavefy-$waveallbd]                  -outline lightgreen -stipple gray12 -tag [list box boxtem]
	.wavealln create rect $curxl [expr $wavefy-$waveallbd] $curxr $wavefy                   -fill lightgreen -outline lightgreen -stipple gray12 -tag box
    }
    bind .wavealln <ButtonRelease-1> {
	if {$ifdrag==0} {return}
	# if window is too small then return
	if {[expr abs($curxr-$curxl)<10]} {tk_messageBox -type ok -message "window is too short";return}
	set chanid BHN
	if {$ifdrag==1} {
	    set timewinleft $curxl; set timewinright $curxr;
	    waveallsub $curxl $curxr $chanid;
	}
    }
    #################### Drag to show wavezoom
    #bind .wavealln <Enter> {focus %W}
    #bind .wavealln <Key-space> {
	#pack forget .waveallz
	#pack forget .wavealle
	#pack forget .wavealln
	#pack .labelshowwavealln -in .framewaven -side top -fill x
    #}
    bind .labelshowwaveall <ButtonRelease-1> {
	showini
	#pack .waveallz -in .framewavez -side top -fill x
	#pack .wavealle -in .framewavee -side top -fill x
	#pack .wavealln -in .framewaven -side top -fill x
	#pack forget .labelshowwaveall
	#pack .labelhidewaveall -in .frameshowhide -side right -fill x
	#pack .framestafig -in .leftcolumnsubtop
	#pack .framepara -in .leftcolumnsubmiddle -fill x
    }
    bind .labelhidewaveall <ButtonRelease-1> {
	pack forget .waveallz
	pack forget .wavealle
	pack forget .wavealln
	pack forget .labelhidewaveall
	pack .labelshowwaveall -in .frameshowhide -side right -fill x
	pack forget .framepara; .leftcolumnsubmiddle configure -height 1.0

	if {$chanidcur==""   } {pack .waveallz -in .leftcolumnsubmiddle -side top}
	if {$chanidcur=="BHZ"} {pack .waveallz -in .leftcolumnsubmiddle -side top}
	if {$chanidcur=="BHE"} {pack .wavealle -in .leftcolumnsubmiddle -side top}
	if {$chanidcur=="BHN"} {pack .wavealln -in .leftcolumnsubmiddle -side top}
	#pack forget .framestafig
	#if {$chanidcur==""   } {pack .waveallz -in .leftcolumnsubtop -side top}
	#if {$chanidcur=="BHZ"} {pack .waveallz -in .leftcolumnsubtop -side top}
	#if {$chanidcur=="BHE"} {pack .wavealle -in .leftcolumnsubtop -side top}
	#if {$chanidcur=="BHN"} {pack .wavealln -in .leftcolumnsubtop -side top}

    }

    ####################
    bind .labelzoomup <ButtonPress-1>   {
	if {$zoomupdnavailable==0} {return;}
	.labelzoomup configure -bg lightgrey
    }
    bind .labelzoomdn <ButtonPress-1>   {
	if {$zoomupdnavailable==0} {return;}
	.labelzoomdn configure -bg lightgrey
    }
    bind .labelzoomup <ButtonRelease-1> {
	if {$zoomupdnavailable==0} {return;}
	.labelzoomup configure -bg orange; wavezoomscalezoomy 2.0
    }
    bind .labelzoomdn <ButtonRelease-1> {
	if {$zoomupdnavailable==0} {return;}
	.labelzoomdn configure -bg orange; wavezoomscalezoomy 0.5
    }

    #################### 
    bind .wavezoom <Button-4> {
	if {[.wavezoom find withtag wave]==""} {return}
	%W yview scroll -2 units;
	plotscaleinwavezoom .wavezoom
    }
    bind .wavezoom <Button-5> {
	if {[.wavezoom find withtag wave]==""} {return}
	%W yview scroll +2 units;
	plotscaleinwavezoom .wavezoom
    }
    bind .wavezoom <Double-Button-1> {
	if {[.wavezoom find withtag wave]==""} {return}
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set seln [expr round($cury/$singlewaveheight)]
	set seln [expr $sfnum+1-$seln]
	if {$seln<1 || $seln>$sfnum} {return}
	addpick .wavezoom $curx $seln red
	set sfmansw($dtname($seln).$chanidcur) 1

#	if {$chanidcur=="BHZ"} {set canvasname .waveallz}
#	if {$chanidcur=="BHE"} {set canvasname .wavealle}
#	if {$chanidcur=="BHN"} {set canvasname .wavealln}
#	if {[.wavezoom find withtag tmark.$dtname($seln)]!=""} {
#	    # Move tmark to new position
#	    set curx0 [lindex [.wavezoom coords tmark.$dtname($seln)] 0]
#	    set delx [expr $curx-$curx0]
#	    .wavezoom move tmark.$dtname($seln) $delx 0
#
#	    # recorder picktime
#	    set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
#	    set sfpt($dtname($seln).$chanidcur) $tt
#
#	    # shift tmark in waveall, or adjust curve
#	    if {$tmarkcurvesw($chanidcur)==1} {
#		set curxold [lindex [$canvasname coords tmark.$dtname($seln)] 0]
#		set tt [waveallxy2poixy $sfpt($dtname($seln).$chanidcur) $sfdist($seln) $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
#		set curxnew [lindex $tt 0]
#		set delx [expr $curxnew-$curxold]
#		$canvasname move tmark.$dtname($seln) $delx 0
#	    } else {
#		plotcurve $chanidcur
#	    }
#
#	} else {
#	    # Create a tmark
#	    set tmarkx $curx;
#	    #set tmarky $sfyaxisposinwavezoom($dtname($seln))
#	    set tmarky [expr $singlewaveheight*($sfnum+1-$seln)]
#	    set tmarkhalflength [expr $singlewaveheight*0.3]
#	    .wavezoom create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 2 -tags [list tmark tmark.$dtname($seln)] -fill red
#
#	    # recorder picktime and change sfpick
#	    set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
#	    set sfpt($dtname($seln).$chanidcur) $tt
#	    set sfpick($dtname($seln).$chanidcur) 1
#
#	    # Add a new tmark in waveall, or adjust curve
#	    if {$tmarkcurvesw($chanidcur)==1} {
#		set tt [waveallxy2poixy $sfpt($dtname($seln).$chanidcur) $sfdist($seln) $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
#		set tmarkx [lindex $tt 0]
#		set tmarky [lindex $tt 1]
#		#set tmarky $sfyaxisposinwaveall($dtname($seln))
#		#set tmarkhalflength [expr $waveallyscale*0.75]
#		set tmarkhalflength 7
#		$canvasname create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 1 -tags [list tmark tmark.$dtname($seln)] -fill red
#	    } else {
#		plotcurve $chanidcur
#	    }
#	    .wavezoom  itemconfigure statext.$dtname($seln) -fill red
#	    #.stafigure itemconfigure station.$dtname($seln) -fill red
#	}
#	colorstation $chanidcur
#	plotcurdots $chanidcur
#	getpickcurrentevtdatrange $chanidcur
#	plotzoomcurve $chanidcur
    }
    bind .wavezoomreplot <Double-Button-1> {
	if {[.wavezoomreplot find withtag wave]==""} {return}
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set seln [expr round($cury/$singlewaveheight)]
	set seln [expr $sfnum+1-$seln]
	if {$seln<1 || $seln>$sfnum} {return}
	addpick .wavezoomreplot $curx $seln red
	set sfmansw($dtname($seln).$chanidcur) 1

#	if {$chanidcur=="BHZ"} {set canvasname .waveallz}
#	if {$chanidcur=="BHE"} {set canvasname .wavealle}
#	if {$chanidcur=="BHN"} {set canvasname .wavealln}
#	if {[.wavezoomreplot find withtag tmark.$dtname($seln)]!=""} {
#	    # Move tmark to new position
#	    set curx0 [lindex [.wavezoomreplot coords tmark.$dtname($seln)] 0]
#	    set delx [expr $curx-$curx0]
#	    .wavezoomreplot move tmark.$dtname($seln) $delx 0
#	    .wavezoom       move tmark.$dtname($seln) $delx 0
#
#	    # recorder picktime
#	    set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
#	    set sfpt($dtname($seln).$chanidcur) $tt
#
#	    # shift tmark in waveall, or adjust curve
#	    if {$tmarkcurvesw($chanidcur)==1} {
#		set curxold [lindex [$canvasname coords tmark.$dtname($seln)] 0]
#		set tt [waveallxy2poixy $sfpt($dtname($seln).$chanidcur) $sfdist($seln) $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
#		set curxnew [lindex $tt 0]
#		set delx [expr $curxnew-$curxold]
#		$canvasname move tmark.$dtname($seln) $delx 0
#	    } else {
#		plotcurve $chanidcur
#	    }
#
#	} else {
#	    # Create a tmark
#	    set tmarkx $curx;
#	    #set tmarky $sfyaxisposinwavezoom($dtname($seln))
#	    set tmarky [expr $singlewaveheight*($sfnum+1-$seln)]
#	    set tmarkhalflength [expr $singlewaveheight*0.3]
#	    .wavezoomreplot create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 2 -tags [list tmark tmark.$dtname($seln)] -fill red
#	    .wavezoom       create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 2 -tags [list tmark tmark.$dtname($seln)] -fill red
#
#	    # recorder picktime and change sfpick
#	    set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
#	    set sfpt($dtname($seln).$chanidcur) $tt
#	    set sfpick($dtname($seln).$chanidcur) 1
#
#	    # Add a new tmark in waveall, or adjust curve
#	    if {$tmarkcurvesw($chanidcur)==1} {
#		set tt [waveallxy2poixy $sfpt($dtname($seln).$chanidcur) $sfdist($seln) $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
#		set tmarkx [lindex $tt 0]
#		set tmarky [lindex $tt 1]
#		#set tmarky $sfyaxisposinwaveall($dtname($seln))
#		#set tmarkhalflength [expr $waveallyscale*0.75]
#		set tmarkhalflength 7
#		$canvasname create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 1 -tags [list tmark tmark.$dtname($seln)] -fill red
#	    } else {
#		plotcurve $chanidcur
#	    }
#	    .wavezoomreplot itemconfigure statext.$dtname($seln) -fill red
#	    .wavezoom       itemconfigure statext.$dtname($seln) -fill red
#	    #.stafigure      itemconfigure station.$dtname($seln) -fill red
#	}
#	colorstation $chanidcur
#	plotcurdots $chanidcur
#	getpickcurrentevtdatrange $chanidcur
#	plotzoomcurve $chanidcur
    }
    bind .wavezoom <ButtonPress-1> {
	set curxl [%W canvasx %x];
	set ifdrag 0
    }
    bind .wavezoom <B1-Motion> {
	if {[.wavezoom find withtag wave]==""} {return}
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
	set tt [formatsub3 $tt]
	.information configure -text "Poix:$curx Poiy:$cury Time:$tt"
	
	if {$zoomdragcontrol==0} {return}
	set ifdrag 1
	set curxr [%W canvasx %x];
	.wavezoom delete box
	.wavezoom create rect $curxl 0 $curxr $wavezoomfy -fill gray -outline white -stipple gray12 -tag box
    }
    bind .wavezoom <ButtonRelease-1> {
	###### if not drag then return
	if {$ifdrag==0} {return}
	###### if drag delete box
	.wavezoom delete box
	###### if not found wave then return
	if {[.wavezoom find withtag wave]==""} {return}
	set curxr [%W canvasx %x];
	###### if window is too small then return
	if {[expr abs($curxr-$curxl)<10]} {return}

	loadwavezoomstainfopicktmark
	# sort curxl curxr
	if {$curxr<$curxl} {set tem   $curxr; set curxr $curxl; set curxl $tem; }

	set scalemul [expr 1000.0/($curxr-$curxl)]
	set xratio [expr $curxl/$wavezoomfx]
	set wavezoomfx [expr $wavezoomfx*$scalemul]
	set wavezoomxscale [expr ($wavezoomfx-$wavezoombd*2.0)/($wavezoomx2-$wavezoomx1)]
	.wavezoom configure -scrollregion [list 0 0 $wavezoomfx $wavezoomfy]
	.wavezoom scale wave 0 0 $scalemul 1
	.wavezoom xview moveto [expr $xratio]
	if {[.wavezoom find withtag waveenvelope]!=""} {.wavezoom scale waveenvelope 0 0 $scalemul 1}
	if {[.wavezoom find withtag waveinstfreq]!=""} {.wavezoom scale waveinstfreq 0 0 $scalemul 1}
	resetupwavezoomstainfopicktmark  .wavezoom
	plottmarkrangeinwavezoom
	plotscaleinwavezoom .wavezoom
	plotsynthtmarkinwavezoom .wavezoom
    }
    bind .wavezoom <Control-Button-4> {
	if {[.wavezoom find withtag wave]==""} {return}
	loadwavezoomstainfopicktmark
	set scalemul 2.0
	set curx [%W canvasx %x];
	set curscreenx %x
	set xratio [expr $curx/$wavezoomfx]

	set wavezoomfx [expr $wavezoomfx*$scalemul]
	set wavezoomxscale [expr ($wavezoomfx-$wavezoombd*2.0)/($wavezoomx2-$wavezoomx1)]
	.wavezoom configure -scrollregion [list 0 0 $wavezoomfx $wavezoomfy]
	.wavezoom scale wave 0 0 $scalemul 1
	.wavezoom xview moveto [expr ($xratio*$wavezoomfx-$curscreenx)/$wavezoomfx]
	if {[.wavezoom find withtag waveenvelope]!=""} {.wavezoom scale waveenvelope 0 0 $scalemul 1}
	if {[.wavezoom find withtag waveinstfreq]!=""} {.wavezoom scale waveinstfreq 0 0 $scalemul 1}
	resetupwavezoomstainfopicktmark  .wavezoom
	plottmarkrangeinwavezoom
	plotscaleinwavezoom .wavezoom
	plotsynthtmarkinwavezoom .wavezoom
    }
    bind .wavezoom <Control-Button-5> {
	if {[.wavezoom find withtag wave]==""} {return}
	loadwavezoomstainfopicktmark
	set scalemul 0.5
	set curx [%W canvasx %x];
	set curscreenx %x
	set xratio [expr $curx/$wavezoomfx]

	if {$wavezoomfx<=1000} {set scalemul 1.0;}
	if {[expr $wavezoomfx*$scalemul]<1000.0} {
	    set scalemul [expr 1000.0/$wavezoomfx]
	}
	set wavezoomfx [expr $wavezoomfx*$scalemul]
	set wavezoomxscale [expr ($wavezoomfx-$wavezoombd*2.0)/($wavezoomx2-$wavezoomx1)]
	.wavezoom configure -scrollregion [list 0 0 $wavezoomfx $wavezoomfy]
	.wavezoom scale wave 0 0 $scalemul 1
	.wavezoom xview moveto [expr ($xratio*$wavezoomfx-$curscreenx)/$wavezoomfx]
	if {[.wavezoom find withtag waveenvelope]!=""} {.wavezoom scale waveenvelope 0 0 $scalemul 1}
	if {[.wavezoom find withtag waveinstfreq]!=""} {.wavezoom scale waveinstfreq 0 0 $scalemul 1}
	resetupwavezoomstainfopicktmark  .wavezoom
	plottmarkrangeinwavezoom
	plotscaleinwavezoom .wavezoom
	plotsynthtmarkinwavezoom .wavezoom
    }
    bind .wavezoom <3> {
	if {[.wavezoom find withtag wave]==""} {return}
	loadwavezoomstainfopicktmark
	set scalemul 0.5
	set curx [%W canvasx %x];
	set curscreenx %x
	set xratio [expr $curx/$wavezoomfx]

	if {$wavezoomfx<=1000} {set scalemul 1.0;}
	if {[expr $wavezoomfx*$scalemul]<1000.0} {
	    set scalemul [expr 1000.0/$wavezoomfx]
	}
	set wavezoomfx [expr $wavezoomfx*$scalemul]
	set wavezoomxscale [expr ($wavezoomfx-$wavezoombd*2.0)/($wavezoomx2-$wavezoomx1)]
	.wavezoom configure -scrollregion [list 0 0 $wavezoomfx $wavezoomfy]
	.wavezoom scale wave 0 0 $scalemul 1
	set xratio [expr ($xratio*$wavezoomfx-$curscreenx)/$wavezoomfx]
	.wavezoom xview moveto $xratio
	if {[.wavezoom find withtag waveenvelope]!=""} {.wavezoom scale waveenvelope 0 0 $scalemul 1}
	if {[.wavezoom find withtag waveinstfreq]!=""} {.wavezoom scale waveinstfreq 0 0 $scalemul 1}
	resetupwavezoomstainfopicktmark .wavezoom
	plottmarkrangeinwavezoom
	plotscaleinwavezoom .wavezoom
	plotsynthtmarkinwavezoom .wavezoom
    }
    bind .wavezoom <Enter> {focus %W}
    bind .wavezoom <Key-space> {
	if {[.wavezoom find withtag wave]==""} {return}
	replot
	.labelzoomup configure -bg lightgray
	.labelzoomdn configure -bg lightgray
	set zoomupdnavailable 0
    }
    bind .wavezoomreplot <Enter> {focus %W}
    bind .wavezoomreplot <Key-space> {
	pack forget .framewavezoomreplot
	pack .framewavezoom -in .rightbottom -side left
	.labelzoomup configure -bg orange
	.labelzoomdn configure -bg orange
	set zoomupdnavailable 1
    }


    #bind .wavezoom <Enter> {%W config -cursor tcross}
    bind .wavezoom <Motion> {
	if {[.wavezoom find withtag wave]==""} {return}
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd];
	set tt [formatsub3 $tt]
	.information configure -text "Poix:$curx Poiy:$cury Time:$tt"
	.wavezoom coords linex1 0 $cury [expr $curx-2] $cury
	.wavezoom coords linex2 [expr $curx+2] $cury $wavezoomfx $cury
	.wavezoom coords liney1 $curx 0 $curx [expr $cury-2]
	.wavezoom coords liney2 $curx [expr $cury+2] $curx $wavezoomfy
    }
    bind .wavezoomreplot <Motion> {
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
	set tt [formatsub3 $tt]
	.information configure -text "Poix:$curx Poiy:$cury Time:$tt"
    }
    bind .wavezoom <2> {
	if {[.wavezoom find withtag wave]==""} {return}
	if {$zoompress2control==0} {return}
	if {[.wavezoom find withtag synthtmark]==""} {
	    set synthtmarkinwavezoomsw 1
	    plotsynthtmarkinwavezoom .wavezoom
	} else {
	    set synthtmarkinwavezoomsw 0
	    .wavezoom delete synthtmark
	}
    }
    bind .wavezoomreplot <2> {
	if {[.wavezoomreplot find withtag wave]==""} {return}
	if {$zoompress2control==0} {return}
	if {[.wavezoomreplot find withtag synthtmark]==""} {
	    set synthtmarkinwavezoomsw 1
	    plotsynthtmarkinwavezoom .wavezoomreplot
	    plotsynthtmarkinwavezoom .wavezoom
	} else {
	    set synthtmarkinwavezoomsw 0
	    .wavezoomreplot delete synthtmark
	    .wavezoom       delete synthtmark
	}
    }
    ################################ Middle click to turn on/off synthetic time
    bind .waveallz <2> {
	if {[.waveallz find withtag wave]==""} {return}
	if {$allpress2control==0} {return}
	if {[.waveallz find withtag synthtmark]==""} {
	    set synthtmarkinwaveallsw(BHZ) 1
	    plotsynthtmarkinwaveall BHZ
	} else {
	    set synthtmarkinwaveallsw(BHZ) 0
	    .waveallz delete synthtmark
	}
    }
    bind .wavealle <2> {
	if {[.wavealle find withtag wave]==""} {return}
	if {$allpress2control==0} {return}
	if {[.wavealle find withtag synthtmark]==""} {
	    set synthtmarkinwaveallsw(BHE) 1
	    plotsynthtmarkinwaveall BHE
	} else {
	    set synthtmarkinwaveallsw(BHE) 0
	    .wavealle delete synthtmark
	}
    }
    bind .wavealln <2> {
	if {[.wavealln find withtag wave]==""} {return}
	if {$allpress2control==0} {return}
	if {[.wavealln find withtag synthtmark]==""} {
	    set synthtmarkinwaveallsw(BHN) 1
	    plotsynthtmarkinwaveall BHN
	} else {
	    set synthtmarkinwaveallsw(BHN) 0
	    .wavealln delete synthtmark
	}
    }
    ################################ Middle click to turn on/off synthetic time
    ################################ Moving and showing x-y-coordinate value
    bind .waveallz <Motion> {
	if {[.waveallz find withtag line]==""} {return}
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set chanid BHZ
	set tt [waveallpoixy2xy $curx $cury $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
	set t1 [lindex $tt 0]; set t1 [formatsub3 $t1]
	set t2 [lindex $tt 1]; set t2 [formatsub3 $t2]
	.information configure -text "Poix:$curx Poiy:$cury Time-Dist:$t1 $t2"
    }
    bind .wavealle <Motion> {
	if {[.wavealle find withtag line]==""} {return}
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set chanid BHE
	set tt [waveallpoixy2xy $curx $cury $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
	set t1 [lindex $tt 0]; set t1 [formatsub3 $t1]
	set t2 [lindex $tt 1]; set t2 [formatsub3 $t2]
	.information configure -text "Poix:$curx Poiy:$cury Time-Dist:$t1 $t2"
    }
    bind .wavealln <Motion> {
	if {[.wavealln find withtag line]==""} {return}
	set curx [%W canvasx %x];
	set cury [%W canvasy %y];
	set chanid BHN
	set tt [waveallpoixy2xy $curx $cury $waveallx1($chanid) $waveally2($chanid) $waveallxscale($chanid) $waveallyscale($chanid) $waveallbd]
	set t1 [lindex $tt 0]; set t1 [formatsub3 $t1]
	set t2 [lindex $tt 1]; set t2 [formatsub3 $t2]
	.information configure -text "Poix:$curx Poiy:$cury Time-Dist:$t1 $t2"
    }
    ################################ Moving and showing x-y-coordinate value
}
############################################################
proc resetstafigview {} {
    global x1 y1 x2 y2
    global stafile
    ### load & save evtparameters
    loadevtpara

    ### load & save station parameters
    loadstapara

    set y1 [exec awk "BEGIN{tem= 180} {if(\$2<tem)tem=\$2} END{print tem}" $stafile]
    set y2 [exec awk "BEGIN{tem=-180} {if(\$2>tem)tem=\$2} END{print tem}" $stafile]
    set x1 [exec awk "BEGIN{tem= 180} {if(\$3<tem)tem=\$3} END{print tem}" $stafile]
    set x2 [exec awk "BEGIN{tem=-180} {if(\$3>tem)tem=\$3} END{print tem}" $stafile]
    setstafigpara_plotsta_evtfree 0;

    ### re-setup evtparameters
    resetupevtpara

    ### re-setup station parameters
    resetupstapara

}
############################################################
proc addpick {wavezoomname curx seln tmarkcolor} {
    global sfnum dtname sfdist sfpick sfpt sfphase
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global chanidcur
    global tmarkcurvesw
    global phasetype
    global pickcurlimt1 pickcurlimt2 pickcurlimd1 pickcurlimd2

    if {$chanidcur=="BHZ"} {set waveallname .waveallz}
    if {$chanidcur=="BHE"} {set waveallname .wavealle}
    if {$chanidcur=="BHN"} {set waveallname .wavealln}

    if {[$wavezoomname find withtag tmark.$dtname($seln)]!=""} {
	# Move tmark to new position
	set curx0 [lindex [$wavezoomname coords tmark.$dtname($seln)] 0]
	set delx [expr $curx-$curx0]
	$wavezoomname move tmark.$dtname($seln) $delx 0
	if {$wavezoomname == ".wavezoomreplot"} {
	    .wavezoom move tmark.$dtname($seln) $delx 0
	}

	# recorder picktime
	set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
	set sfpt($dtname($seln).$chanidcur) $tt
	set sfphase($dtname($seln).$chanidcur) $phasetype

	# shift tmark in waveall, or adjust curve
	if {$tmarkcurvesw($chanidcur)==1} {
	    set curxold [lindex [$waveallname coords tmark.$dtname($seln)] 0]
	    set tt [waveallxy2poixy $sfpt($dtname($seln).$chanidcur) $sfdist($seln) $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
	    set curxnew [lindex $tt 0]
	    set delx [expr $curxnew-$curxold]
	    $waveallname move tmark.$dtname($seln) $delx 0
	} else {
	    plotcurve $chanidcur
	}

    } else {
	# Create a tmark
	set tmarkx $curx;
	#set tmarky $sfyaxisposinwavezoom($dtname($seln))
	set tmarky [expr $singlewaveheight*($sfnum+1-$seln)]
	set tmarkhalflength [expr $singlewaveheight*0.3]
	$wavezoomname create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 2 -tags [list tmark tmark.$dtname($seln)] -fill $tmarkcolor
	if {$wavezoomname == ".wavezoomreplot"} {
	    .wavezoom create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 2 -tags [list tmark tmark.$dtname($seln)] -fill $tmarkcolor
	}

	# recorde picktime and change sfpick
	set tt [wavezoompoix2x $curx $wavezoomx1 $wavezoomxscale $wavezoombd]
	set sfpt($dtname($seln).$chanidcur) $tt
	set sfphase($dtname($seln).$chanidcur) $phasetype
	set sfpick($dtname($seln).$chanidcur) 1

	# Add a new tmark in waveall, or adjust curve
	if {$tmarkcurvesw($chanidcur)==1} {
	    set tt [waveallxy2poixy $sfpt($dtname($seln).$chanidcur) $sfdist($seln) $waveallx1($chanidcur) $waveally1($chanidcur) $waveallxscale($chanidcur) $waveallyscale($chanidcur) $waveallbd $wavefy]
	    set tmarkx [lindex $tt 0]
	    set tmarky [lindex $tt 1]
	    #set tmarky $sfyaxisposinwaveall($dtname($seln))
	    #set tmarkhalflength [expr $waveallyscale*0.75]
	    set tmarkhalflength 7
	    $waveallname create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 1 -tags [list tmark tmark.$dtname($seln)] -fill red
	} else {
	    plotcurve $chanidcur
	}
	$wavezoomname  itemconfigure statext.$dtname($seln) -fill red
	if {$wavezoomname == ".wavezoomreplot"} {
	    .wavezoom  itemconfigure statext.$dtname($seln) -fill red
	}

    }
    colorstation $chanidcur
    plotcurdots $chanidcur
    getpickcurrentevtdatrange $chanidcur
    #set temi 0;
    #for {set i 1} {$i<=$sfnum} {incr i 1} {set temi [expr $temi+$sfpick($dtname($i).$chanidcur)]}
    #if {$temi>=2} {zoomcurvefigureini $chanidcur $pickcurlimt1 $pickcurlimt2 $pickcurlimd1 $pickcurlimd2}
    zoomcurvefigureini $chanidcur $pickcurlimt1 $pickcurlimt2 $pickcurlimd1 $pickcurlimd2
    plotzoomcurve $chanidcur

}
############################################################
proc loadandpreparedatandcallpick {} {
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallbd wavezoombd singlewaveheight
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfnum dtname sfpick
    global zoomupdnavailable

    global sacdatwinb sacdatwine sacdat dt
    global sacdatfileb sacdatfilee
    global chanidcur
    global timelist
    global evtomark
    global phasetype sfphase sfmansw

    global stwin ltwin stltlim
    global forcesw forcewin methodsw



    set tt [.wavezoom xview]
    set xfractor1 [lindex $tt 0]; set xfractor2 [lindex $tt 1];
    set xpoil [expr round($wavezoomfx*$xfractor1)]
    set xpoir [expr round($wavezoomfx*$xfractor2)]
    set xl [wavezoompoix2x $xpoil $wavezoomx1 $wavezoomxscale $wavezoombd]
    set xr [wavezoompoix2x $xpoir $wavezoomx1 $wavezoomxscale $wavezoombd]

    set tt [.wavezoom yview]
    set yfractor1 [lindex $tt 0]; set yfractor2 [lindex $tt 1];
    set ypoit [expr round($wavezoomfy*$yfractor1)]
    set ypoib [expr round($wavezoomfy*$yfractor2)]
    set yt [expr int($ypoit/$singlewaveheight)+1];if {$yt<1} {set yt 1;}
    set yb [expr int($ypoib/$singlewaveheight)+0];if {$yb>$sfnum} {set yb $sfnum;}
    set y1 [expr $sfnum+1-$yb]
    set y2 [expr $sfnum+1-$yt]

    for {set i $y1} {$i<=$y2} {incr i 1} {
	#if {$sfpick($dtname($i).$chanidcur)==1 && $sfphase($dtname($i).$chanidcur)==$phasetype} {continue}
	if {$sfmansw($dtname($i).$chanidcur)==1 && $sfphase($dtname($i).$chanidcur)==$phasetype} {continue}

	# Adjust the begin time for each file
	if {$xl<$sacdatfileb($dtname($i).$chanidcur)} {
	    set datbegtime $sacdatfileb($dtname($i).$chanidcur)
	} elseif {$xl>=$sacdatfileb($dtname($i).$chanidcur) && $xl<$sacdatfilee($dtname($i).$chanidcur)} {
	    set datbegtime $xl
	} elseif {$xl>=$sacdatfilee($dtname($i).$chanidcur)} {
	    continue;
	}

	if {$xr<=$sacdatfileb($dtname($i).$chanidcur)} {
	    continue;
	} elseif {$xr>$sacdatfileb($dtname($i).$chanidcur) && $xr<=$sacdatfilee($dtname($i).$chanidcur)} {
	    set datendtime $xr
	} elseif {$xr>$sacdatfilee($dtname($i).$chanidcur)} {
	    set datendtime $sacdatfilee($dtname($i).$chanidcur)
	}

	# Optimize the time window (According to predictions If forcesw==Yes)
	if {$forcesw=="Yes"} {
	set jend [llength $timelist($dtname($i))];
	set swtp 0; set swts 0;
	for {set j 0} {$j<$jend} {incr j 3} {
	    set tmarklabel [lindex $timelist($dtname($i)) $j];
	    if {[string index $tmarklabel 0]=="P" && $swtp==0} {
		set n [expr $j+1]; set tp [lindex $timelist($dtname($i)) $n];
		set swtp 1;
	    }

	    if {[string index $tmarklabel 0]=="S" && $swts==0} {
		set n [expr $j+1]; set ts [lindex $timelist($dtname($i)) $n];
		set swts 1;
	    }
	    if {$swtp==1 && $swts==1} {break;}
	}	
	#set tpleft [expr $tp/(1.0+0.15)];      set tpright [expr $tp/(1.0-0.15)];
	#set tsleft [expr $ts/(1.0+0.15)];      set tsright [expr $ts/(1.0-0.15)];
	set tpleft [expr $tp-$forcewin];      set tpright [expr $tp+$forcewin];
	set tsleft [expr $ts-$forcewin];      set tsright [expr $ts+$forcewin];
	set tpwinbeg [expr $evtomark+$tpleft]; set tpwinend [expr $evtomark+$tpright];
	set tswinbeg [expr $evtomark+$tsleft]; set tswinend [expr $evtomark+$tsright];
	if {$tpwinend>$tswinbeg} {
	    set tpwinend [expr ($tpwinend+$tswinbeg)/2.0];
	    set tswinbeg [expr ($tpwinend+$tswinbeg)/2.0];
	}
	if {$phasetype=="P"} {
	    if {$datbegtime>=$tpwinbeg} {set datbegtime $datbegtime}
	    if {$datbegtime< $tpwinbeg} {set datbegtime $tpwinbeg}
	    if {$datendtime<=$tpwinbeg} {continue;}

	    if {$datbegtime>=$tpwinend} {continue;}
	    if {$datendtime> $tpwinend} {set datendtime $tpwinend}
	    if {$datendtime<=$tpwinend} {set datendtime $datendtime}
	} elseif {$phasetype=="S"} {
	    if {$datbegtime>=$tswinbeg} {set datbegtime $datbegtime}
	    if {$datbegtime< $tswinbeg} {set datbegtime $tswinbeg}
	    if {$datendtime<=$tswinbeg} {continue;}

	    if {$datbegtime>=$tswinend} {continue;}
	    if {$datendtime> $tswinend} {set datendtime $tswinend}
	    if {$datendtime<=$tswinend} {set datendtime $datendtime}
	}

	}

	# Obtain data-index range (for time)
	set datbegint [expr int(($datbegtime-$sacdatfileb($dtname($i).$chanidcur))/$dt)+1]
	set datendint [expr int(($datendtime-$sacdatfileb($dtname($i).$chanidcur))/$dt)+1]

	# Copy windat
	set sacdattem [lrange $sacdat($dtname($i).$chanidcur) [expr $datbegint-1] [expr $datendint-1]]

	# Prepare data
	set ampy1 {}; set ampy2 {};
	set npts [expr $datendint-$datbegint+1]
	for {set l 0} {$l<$npts} {incr l 1} {
	    set sgval [lindex $sacdattem $l]
	    lappend ampy1 $sgval
	    lappend ampy2 [expr $sgval*$sgval]
	}

	if {$methodsw=="STALTA"} {
	    ####  BEG STA/LTA
	    # Refine time window if necessary
	    set stwin0 $stwin;
	    set ltwin0 $ltwin;
	    if {[expr $datendtime-$datbegtime]<=[expr $stwin+$ltwin]} {
		set stwin0 [expr ($datendtime-$datbegtime)/2.5]
		set ltwin0 [expr ($datendtime-$datbegtime)/2.5]
	    }

	    set stwinint [expr int($stwin0/$dt)]
	    set ltwinint [expr int($ltwin0/$dt)]
	    set tt [stalta $ampy1 $ampy2 $npts $stwinint $ltwinint]
	    if {$tt==""} {continue;}
	    set slvalue [lindex $tt 0]
	    set selint [lindex $tt 1]
	    set seltime [expr $datbegtime+$selint*$dt]
	    ####  END STA/LTA
	} elseif {$methodsw=="AIC"} {
	    ####  BEG AIC
	    #set aicwin0 $stwin
	    #if {[expr $datendtime-$datbegtime]<$stwin} {
		#set aicwin0 [expr $datendtime-$datbegtime]
	    #}
	    #set aicwinint [expr int($aicwin0/$dt)]
	    set tt [aicsub $ampy1 $ampy2 $npts]
	    if {$tt==""} {continue;}
	    set aicvalue [lindex $tt 0]
	    set selint [lindex $tt 1]
	    set seltime [expr $datbegtime+$selint*$dt]
	    ####  END AIC
	}

	#### make a mark
	set curx [wavezoomx2poix $seltime $wavezoomx1 $wavezoomxscale $wavezoombd]
	if {$phasetype == "P"} {set tmarkcolor "red"}
	if {$phasetype == "S"} {set tmarkcolor "red"}
	if {$zoomupdnavailable==1} {
	    addpick .wavezoom $curx $i $tmarkcolor
	} else {
	    addpick .wavezoomreplot $curx $i $tmarkcolor
	}



    }
}
#################################################
proc aicsub {ampy1 ampy2 npts} {

    if {$npts<=4} {return "";}
    set aicvalue 1.0e20; set selint -1
    set bsm1 0.0; set bsm2 0.0;
    set esm1 0.0; set esm2 0.0;
    for {set i 1} {$i<=1} {incr i 1} {
	set n [expr $i-1];
	set sgval1 [lindex $ampy1 $n]
	set sgval2 [lindex $ampy2 $n]

	set bsm1 [expr $bsm1+$sgval1]
	set bsm2 [expr $bsm2+$sgval2]
    }
    for {set i 2} {$i<=$npts} {incr i 1} {
	set n [expr $i-1];
	set sgval1 [lindex $ampy1 $n]
	set sgval2 [lindex $ampy2 $n]

	set esm1 [expr $esm1+$sgval1]
	set esm2 [expr $esm2+$sgval2]
    }

    set checka [expr $bsm1+$esm1]
    set checkb [expr $bsm2+$esm2]
    set bxave [expr $bsm1/1.0]
    set btem [expr ($bsm2-1.0*$bxave*$bxave)/1.0]
    set exave [expr $esm1/($npts-1)]
    set etem [expr ($esm2-($npts-1)*$exave*$exave)/($npts-1-1)]
    if {$btem==0.0} {
	set btem 0.0;
    } elseif {$btem>0.0} {
	set btem [expr log($btem)]
    } else {
	puts "Message from AIC: Unknown Warning1";
    }
    if {$etem==0.0} {
	set etem 0.0;
    } elseif {$etem>0.0} {
	set etem [expr log($etem)]
    } else {
	puts "Message from AIC: Unknown Warning2";
    }
    set tem [expr 1.0*$btem+($npts-1-1)*$etem]
    #set aicvalue $tem

    for {set i 2} {$i<=[expr $npts-1]} {incr i 1} {
	set n [expr $i-1];
	set sgval1 [lindex $ampy1 $n]
	set sgval2 [lindex $ampy2 $n]
	set bsm1 [expr $bsm1+$sgval1]
	set bsm2 [expr $bsm2+$sgval2]
	set esm1 [expr $esm1-$sgval1]
	set esm2 [expr $esm2-$sgval2]
	#if {$checka!=[expr $bsm1+$esm1] || $checkb!=[expr $bsm2+$esm2]} {
	    #set tem1 [expr $bsm1+$esm1-$checka]
	    #set tem2 [expr $bsm2+$esm2-$checkb]
	    #puts "Message from AIC: Unknown Warning3 $tem1 $tem2";
	#}
	set bxave [expr $bsm1/$i]
	set btem [expr ($bsm2-$i*$bxave*$bxave)/($i-1)]
	set exave [expr $esm1/($npts-$i)]
	set divtem [expr $npts-$i-1]
	if {$divtem==0} {set divtem 1}
	set etem [expr ($esm2-($npts-$i)*$exave*$exave)/$divtem]

	if {$btem==0.0} {
	    set btem 0.0;
	} elseif {$btem>0.0} {
	    set btem [expr log($btem)]
	} else {
	    #puts "Message from AIC: Unknown Warning4, $btem";
	    continue;
	}
	if {$etem==0.0} {
	    set etem 0.0;
	} elseif {$etem>0.0} {
	    set etem [expr log($etem)]
	} else {
	    #puts "Message from AIC: Unknown Warning5, $etem";
	    continue;
	}
	set tem [expr $i*$btem+($npts-$i-1)*$etem]
	if {$tem<$aicvalue} {
	   set aicvalue $tem;
	   set selint $n
	}
    }

    if {$selint==-1} {
	return "";
    } else {
	return "$aicvalue $selint"
    }

}
#################################################
proc stalta {ampy1 ampy2 npts stwinint ltwinint} {

    if {$npts<=4} {return "";}
    set slvalue 0.0; set selint -1
    set win1 0.0
    for {set i 0} {$i<$stwinint} {incr i 1} {
	set sgval [lindex $ampy2 $i]
	set win1 [expr $win1+$sgval]
    }
    set win2 0.0
    for {set i $stwinint} {$i< [expr $stwinint+$ltwinint-1]} {incr i 1} {
	set sgval [lindex $ampy2 $i]
	set win2 [expr $win2+$sgval]
    }


    if {$win1!=0.0} {
	set ratio [expr ($win2/$ltwinint)/($win1/$stwinint)]
	set slvalue $ratio
	set selint $stwinint
    }

    set rangeint [expr $npts-$stwinint-$ltwinint]
    for {set i 1} {$i<=$rangeint} {incr i 1} {
	set sgval1 [lindex $ampy2 [expr $i-1]]
	set sgval2 [lindex $ampy2 [expr $i-1+$stwinint]]
	set sgval3 [lindex $ampy2 [expr $i-1+$stwinint+$ltwinint]]
	set win1 [expr $win1-$sgval1]
	set win1 [expr $win1+$sgval2]
	set win2 [expr $win2-$sgval2]
	set win2 [expr $win2+$sgval3]

	if {$win1==0.0} {continue}
	set ratio [expr ($win2/$ltwinint)/($win1/$stwinint)]
	if {$slvalue<=$ratio} {
	    set slvalue $ratio
	    set selint [expr $i+$stwinint]
	}

    }
    if {$selint==-1} {
	return "";
    } else {
	return "$slvalue $selint"
    }
}
#################################################
proc waveallshowhide {canvasname} {

    if {[$canvasname find withtag wave]==""} {return}
    set statussw [$canvasname itemcget wave -state]
    if {$statussw==""} {
	$canvasname itemconfigure wave -state hidden
    } elseif {$statussw=="normal"} {
	$canvasname itemconfigure wave -state hidden
    } elseif {$statussw=="hidden"} {
	$canvasname itemconfigure wave -state normal
    }
}
#################################################
proc showini {} {
    pack .waveallz -in .framewavez -side top -fill x
    pack .wavealle -in .framewavee -side top -fill x
    pack .wavealln -in .framewaven -side top -fill x
    pack forget .labelshowwaveall
    pack .labelhidewaveall -in .frameshowhide -side right -fill x
    #pack .framestafig -in .leftcolumnsubtop
    pack .framepara -in .leftcolumnsubmiddle -fill x
}
proc labelini {} {
    .zinfofile configure -text "Location-File"
    .zinfonew1 configure -text "Location-Rloc"
    .zinfonew2 configure -text "Location-Final";
    .einfofile configure -text "Location-File"
    .einfonew1 configure -text "Location-Rloc"
    .einfonew2 configure -text "Location-Final";
    .ninfofile configure -text "Location-File"
    .ninfonew1 configure -text "Location-Rloc"
    .ninfonew2 configure -text "Location-Final";

}
#################################################
proc plotscaleinwavezoom {canvasname} {
    global wavefx wavefy wavezoomfx wavezoomfy
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale

    $canvasname delete ruler
    set xratio2 [lindex [.wavezoom xview] 1]
    set yratio1 [lindex [.wavezoom yview] 0]
    set xright [expr $xratio2*$wavezoomfx]
    set ytop  [expr $yratio1*$wavezoomfy]

    set xx1 [expr $xright-30-$wavezoomxscale]
    set xx2 [expr $xright-30]
    set xx0 [expr ($xx1+$xx2)/2.0]
    set yy  [expr $ytop+20]
    $canvasname create line [list $xx1 $yy $xx2 $yy] -width 2 -tags [list ruler] -fill red
    $canvasname create text $xx0 $yy -text "1 sec" -tags [list ruler] -anchor n -fill red
}
#################################################
proc plottmarkinwavezoom {} {
    global sfnum dtname sfdist
    global sacdatwinb sacdatwine sacdat dt
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfnum sfname sflat sflon sfele sfpick sfpt
    global chanidcur

    .wavezoom delete tmark
    set tmarkhalflength [expr $singlewaveheight*0.3]
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($dtname($i).$chanidcur)==1} {
	    #Find the x-y-axis pos, tmarky is also equal to sfyaxisposinwavezoom
	    set tmarkx [wavezoomx2poix $sfpt($dtname($i).$chanidcur) $wavezoomx1 $wavezoomxscale $wavezoombd]
	    set tmarky [expr $singlewaveheight*($sfnum+1-$i)]
	    .wavezoom create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 2 -tags [list tmark tmark.$dtname($i)] -fill red
	}
    }

}
#################################################
proc plottmarkrangeinwavezoom {} {
    global sfnum dtname sfdist
    global sacdatwinb sacdatwine sacdat dt
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfnum sfname sflat sflon sfele sfpick sfpt
    global chanidcur

    return;
    .wavezoom delete tmarkrange
    set tmarkxmin  100000000.0; set tmarkxminpoix 0
    set tmarkxmax -100000000.0; set tmarkxmaxpoix 0
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($dtname($i).$chanidcur)==1} {
	    if {$sfpt($dtname($i).$chanidcur)<=$tmarkxmin} {set tmarkxmin $sfpt($dtname($i).$chanidcur);}
	    if {$sfpt($dtname($i).$chanidcur)>=$tmarkxmax} {set tmarkxmax $sfpt($dtname($i).$chanidcur);}
	}
    }
    set tmarkxminpoix [wavezoomx2poix [expr $tmarkxmin-0.5] $wavezoomx1 $wavezoomxscale $wavezoombd]
    set tmarkxmaxpoix [wavezoomx2poix [expr $tmarkxmax+0.5] $wavezoomx1 $wavezoomxscale $wavezoombd]
    set tmarkxminpoix [expr $tmarkxminpoix-0]
    set tmarkxmaxpoix [expr $tmarkxmaxpoix+0]
    .wavezoom create line $tmarkxminpoix 0 $tmarkxminpoix $wavezoomfy -width 1 -tags [list tmarkrange tmark.tmarkmin] -fill red
    .wavezoom create line $tmarkxmaxpoix 0 $tmarkxmaxpoix $wavezoomfy -width 1 -tags [list tmarkrange tmark.tmarkmax] -fill red
}
#################################################
proc plotsynthtmarkinwavezoom {canvasname} {
    global sfnum dtname sfdist sfgcarc
    global sacdatwinb sacdatwine sacdat dt
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfname sflat sflon sfele sfpick sfpt
    global model timelist
    global evtomark
    global synthtmarkinwavezoomsw

    $canvasname delete synthtmark
    if {$synthtmarkinwavezoomsw==0} {return;}
    set tmarkhalflength [expr $singlewaveheight*0.20]
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set iend [llength $timelist($dtname($i))]
	for {set j 0} {$j<$iend} {incr j 3} {
	    set tmarklabel [lindex $timelist($dtname($i)) $j];
	    set n [expr $j+1]; set tt [lindex $timelist($dtname($i)) $n];
	    set tt [expr $tt+$evtomark];
	    if {$tt>=$wavezoomx1 && $tt<=$wavezoomx2} {
		set tmarkx [wavezoomx2poix $tt $wavezoomx1 $wavezoomxscale $wavezoombd]
		set tmarky [expr $singlewaveheight*($sfnum+1-$i)]
		$canvasname create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 1.5 \
                                      -tags [list synthtmark $tmarklabel] -fill cyan
	    }
	}
    }
}
#################################################
proc plotsynthtmarkinwaveall {channel} {
    global sfnum dtname sfdist sfgcarc
    global sacdatwinb sacdatwine sacdat dt
    global waveallbd wavezoombd
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfname sflat sflon sfele sfpick sfpt
    global model timelist
    global evtomark
    global synthtmarkinwaveallsw

    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}
    $canvasname delete synthtmark
    if {$synthtmarkinwaveallsw($channel)==0} {return;}
    #set tmarkhalflength [expr $waveallyscale($channel)*0.55]
    set tmarkhalflength 5
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set iend [llength $timelist($dtname($i))]
	for {set j 0} {$j<$iend} {incr j 3} {
	    set tmarklabel [lindex $timelist($dtname($i)) $j]
	    set n [expr $j+1]; set tt [lindex $timelist($dtname($i)) $n]
	    set tt [expr $tt+$evtomark];
	    if {$tt>=$waveallx1($channel) && $tt<=$waveallx2($channel)} {
		set line [waveallxy2poixy $tt $sfdist($i) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
		set tmarkx [lindex $line 0]; set tmarky [lindex $line 1];
		$canvasname create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 1 \
                                      -tags [list synthtmark $tmarklabel] -fill cyan
	    }
	}
    }

}

#################################################
proc plottmarkinwaveall {channel} {
    global sfnum dtname sfdist
    global sacdatwinb sacdatwine sacdat dt
    global waveallbd wavezoombd
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfname sflat sflon sfele sfpick sfpt
    global tmarkcurvesw
    # tmarkcurvesw=1, plot tmark; tmarkcurvesw=0, plot curve; 

    if {$tmarkcurvesw($channel)==0} {return}
    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}
    $canvasname delete tmark
    #set tmarkhalflength [expr $waveallyscale*0.75]
    set tmarkhalflength 7
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($dtname($i).$channel)==1} {
	    #Find the x-y-axis pos, tmarky is also equal to sfyaxisposinwaveall
	    set line [waveallxy2poixy $sfpt($dtname($i).$channel) $sfdist($i) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
	    set tmarkx [lindex $line 0]; set tmarky [lindex $line 1];
	    $canvasname create line $tmarkx [expr $tmarky-$tmarkhalflength] $tmarkx [expr $tmarky+$tmarkhalflength] -width 1 -tags [list tmark tmark.$dtname($i)] -fill red
	}
    }

}
#################################################
proc plottmarkrangeinwaveall {channel} {
    global sfnum dtname sfdist
    global sacdatwinb sacdatwine sacdat dt
    global waveallbd wavezoombd
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfname sflat sflon sfele sfpick sfpt

    return;
    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}

    $canvasname delete tmarkrange
    set tmarkxmin  100000000.0; set tmarkxminpoix 0
    set tmarkxmax -100000000.0; set tmarkxmaxpoix 0
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($dtname($i).$channel)==1} {
	    if {$sfpt($dtname($i).$channel)<=$tmarkxmin} {set tmarkxmin $sfpt($dtname($i).$channel);}
	    if {$sfpt($dtname($i).$channel)>=$tmarkxmax} {set tmarkxmax $sfpt($dtname($i).$channel);}
	}
    }
    set tt [waveallxy2poixy [expr $tmarkxmin-0.5] 0 $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
    set tmarkxminpoix [lindex $tt 0]
    set tt [waveallxy2poixy [expr $tmarkxmax+0.5] 0 $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
    set tmarkxmaxpoix [lindex $tt 0]
    set tmarkxminpoix [expr $tmarkxminpoix-0]
    set tmarkxmaxpoix [expr $tmarkxmaxpoix+0]
    $canvasname create line $tmarkxminpoix $waveallbd $tmarkxminpoix [expr $wavefy-$waveallbd] -width 1 -tags [list tmarkrange tmark.tmarkmin] -fill red
    $canvasname create line $tmarkxmaxpoix $waveallbd $tmarkxmaxpoix [expr $wavefy-$waveallbd] -width 1 -tags [list tmarkrange tmark.tmarkmax] -fill red
}
#################################################
proc plotcurve {channel} {
    global sfnum dtname sfdist
    global waveallbd wavezoombd
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfname sflat sflon sfele sfpick sfpt
    global tmarkcurvesw
    # tmarkcurvesw=1, plot tmark; tmarkcurvesw=0, plot curve; 

    if {$tmarkcurvesw($channel)==1} {return}
    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}
    $canvasname delete curve
    set aa {}
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($dtname($i).$channel)==1} {
	    set tt [waveallxy2poixy $sfpt($dtname($i).$channel) $sfdist($i) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
	    set tmarkx [lindex $tt 0]
	    set tmarky [lindex $tt 1]
	    lappend aa $tmarkx $tmarky
	}
    }
    if {[llength $aa] <4} {return}
    $canvasname create line $aa -width 1 -tags [list curve] -fill red
}
#################################################
proc plotcurdots {channel} {
    global plotdotsavailablesw
    global sfnum dtname sfdist
    global sfname sflat sflon sfele sfpick sfpt
    global pickfigbd pickfigfx pickfigfy
    global pickfigx1 pickfigx2 pickfigy1 pickfigy2 pickfigxscale pickfigyscale
    global evtomark
    global pickspoix pickspoiy
    
    set swexist [winfo exists .topshowallpicks]
    if {$plotdotsavailablesw==0 || $swexist==0} {return}
    if {$evtomark==-12345.0} {return}

    .topshowallpicks.picksfigure delete curdots
    if {[.topshowallpicks.picksfigure find withtag line]==""} {return}
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($dtname($i).$channel)==1} {
	    set tt [expr $sfpt($dtname($i).$channel)-$evtomark]
	    set line [waveallxy2poixy $sfdist($i) $tt $pickfigx1 $pickfigy1 $pickfigxscale $pickfigyscale $pickfigbd $pickfigfy]
	    set dist [lindex $line 0]; set pickspoix($dtname($i).$channel) $dist
	    set tt [lindex $line 1];   set pickspoiy($dtname($i).$channel) $tt
	    .topshowallpicks.picksfigure creat oval [expr $dist-2.5] [expr $tt-2.5] [expr $dist+2.5] [expr $tt+2.5] -tags [list curdots $dtname($i)] -fill magenta -outline magenta3
	}
    }
}
#################################################
proc plotwaveinall {channel} {
    global sfnum dtname sfdist
    global sfyaxisposinwaveall
    global sacdatwinb sacdatwine sacdat dt
    global sacdatfileb sacdatfilee
    global waveallbd wavezoombd
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global readthreadnum calculatethreadnum

    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}

    # set thread num & create thread, then set them to wait
    for {set i 0} {$i<$calculatethreadnum} {incr i 1} {
	set workid($i) [thread::create {package require Ttrace; thread::wait;}]
    }
    # obtain main thread id
    set mainid [thread::id]

    $canvasname delete wave box
    for {set i 1} {$i<=$sfnum} {incr i 1} {

	#plotwaveinallsub $channel $i
	# Find the first point x-val & y-val, save in lx,ly
	set line [waveallxy2poixy $sacdatfileb($dtname($i).$channel) $sfdist($i) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
	set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
	set lx [lindex [split $line " "] 0]; set ly [lindex [split $line " "] 1]

	# Save vertical position in sfyaxispos
	set sfyaxisposinwaveall($dtname($i).$channel) $ly

	###### workable with tk8.5.13
	#plotwavesub $canvasname $sacdat($dtname($i).$channel) $lx $ly \
	#            $sacdatfileb($dtname($i).$channel) $waveallxscale($channel) $dt 15.0 wave $dtname($i) gray

	###### parallel, workable with tk8.6.12
	#set id [thread::create {package require Ttrace;thread::wait;}]
	#set mainid [thread::id]
	#thread::send -async $id "plotwavesubthread $mainid $canvasname {$sacdat($dtname($i).$channel)} $lx $ly \
        #                         $sacdatfileb($dtname($i).$channel) $waveallxscale($channel) $dt 15.0 wave $dtname($i) gray"

	set imark [expr int(fmod($i,$calculatethreadnum))]
	thread::send -async $workid($imark) "plotwavesubthread $mainid $canvasname {$sacdat($dtname($i).$channel)} $lx $ly \
                                             $sacdatfileb($dtname($i).$channel) $waveallxscale($channel) $dt 15.0 wave $dtname($i) gray"
    }
    # send stop command to each thread
    for {set i 0} {$i<$calculatethreadnum} {incr i 1} {
	thread::send -async $workid($i) "thread::release $workid($i)"
    }
    while {[llength [thread::names]] > 1} {after 200}

	# Set amplitude scale
	#set amx [minmaxabs $sacdat($dtname($i).$channel)]
	#if {$amx>0} {
	#    #set poipam [expr $waveallyscale($channel)/$amx]
	#    set poipam [expr 15.0/$amx]
	#} else {
	#    set poipam 0
	#}

	#set poipam [expr $poipam*1.0]

	#set n 0
	#set t $sacdatfileb($dtname($i).$channel)
	#set aa {}
	#foreach point $sacdat($dtname($i).$channel) {
	    #set xx($n) [expr $lx+$waveallxscale($channel)*($t-$sacdatfileb($dtname($i).$channel))]
	    #set yy($n) [expr $ly-$poipam*$point]
	    #lappend aa $xx($n) $yy($n)
	    #incr n 1
	    #set t [expr $t+$dt]
	#}
        #if {[llength $aa] <4} {continue;}
	#$canvasname create line $aa -width 1 -tags [list wave wave.$dtname($i)] -fill gray


	########### this parallel is not workable
	#    set n 6
	#    for {set i 1} {$i<=$n} {incr i 1} {set var$i "";}
	#    for {set i 1} {$i<=$sfnum} {incr i 1} {
	#	set imark [expr int(fmod($i,$n))]
	#	if {$imark==0} {set imark $n}
	#	lappend var$imark $i
	#    }

	#    foreach i $var1 j $var2 k $var3 l $var4 m $var5 n $var6 {
	#	if {$i!=""} {plotwaveinallsub $channel $i}
	#	if {$j!=""} {plotwaveinallsub $channel $j}
	#	if {$k!=""} {plotwaveinallsub $channel $k}
	#	if {$l!=""} {plotwaveinallsub $channel $l}
	#	if {$m!=""} {plotwaveinallsub $channel $m}
	#	if {$n!=""} {plotwaveinallsub $channel $n}
	#    }
}
########################################################
ttrace::eval {
    proc plotwavesubthread {mainid canvas dat firpoix firpoiy begtime xscale dt traceheight tagmatch kstnm color} {

	set amx -1.e+20
	foreach element $dat {
	    set absval [expr abs($element)]
	    if {$absval > $amx} {set amx $absval}
	}

	#set amx [minmaxabsthread $dat]
	if {$amx>0} {
	    set poipam [expr $traceheight/$amx]
	} else {
	    set poipam 0
	}
	set poipam [expr $poipam*1.0]

	set n 0
	set t $begtime
	set aa {}
	foreach point $dat {
	    set xx($n) [expr $firpoix+$xscale*($t-$begtime)]
	    set yy($n) [expr $firpoiy-$poipam*$point]
	    lappend aa $xx($n) $yy($n)
	    incr n 1
	    set t [expr $t+$dt]
	}
	if {[llength $aa]<4} {return;}
	thread::send -async $mainid "canvasplotsub $canvas {$aa} $tagmatch $kstnm $color"
	#thread::exit
    }
    ############################################################
    proc canvasplotsub {canvas aa tagmatch kstnm color} {
	$canvas create line $aa -width 0.5 -tags [list $tagmatch $tagmatch.$kstnm] -fill $color
    }
    ############################################################
    proc threadrelease {workid} {
	thread::release $workid
    }
}
#################################################
proc plotwaveinallsub {channel ii} {
    global sfnum dtname sfdist
    global sfyaxisposinwaveall
    global sacdatwinb sacdatwine sacdat dt
    global sacdatfileb sacdatfilee
    global waveallbd wavezoombd
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale

    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}

    # Find the first point x-val & y-val, save in lx,ly
    set line [waveallxy2poixy $sacdatfileb($dtname($ii).$channel) $sfdist($ii) $waveallx1($channel) $waveally1($channel) $waveallxscale($channel) $waveallyscale($channel) $waveallbd $wavefy]
    set line [string trim $line]; regsub -all {[[:blank:]]+} $line " " line
    set lx [lindex [split $line " "] 0]; set ly [lindex [split $line " "] 1]

    # Save vertical position in sfyaxispos
    set sfyaxisposinwaveall($dtname($ii).$channel) $ly

    
    plotwavesub $canvasname $sacdat($dtname($ii).$channel) $lx $ly \
                $sacdatfileb($dtname($ii).$channel) $waveallxscale($channel) $dt 15.0 wave $dtname($ii) gray
}
#################################################
proc replot {} {
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallbd wavezoombd singlewaveheight
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global sfnum dtname
    global wavecolor

    pack forget .framewavezoom
    pack .framewavezoomreplot -in .rightbottom -side left
    .wavezoomreplot configure -scrollregion [list 0 0 $wavezoomfx $wavezoomfy]

    set tt [.wavezoom xview]
    set xfractor1 [lindex $tt 0];set xfractor2 [lindex $tt 1];
    set xpoil [expr round($wavezoomfx*$xfractor1)]
    set xpoir [expr round($wavezoomfx*$xfractor2)]
    set xl [wavezoompoix2x $xpoil $wavezoomx1 $wavezoomxscale $wavezoombd]
    set xr [wavezoompoix2x $xpoir $wavezoomx1 $wavezoomxscale $wavezoombd]

    set tt [.wavezoom yview]
    set yfractor1 [lindex $tt 0];set yfractor2 [lindex $tt 1];
    set ypoit [expr round($wavezoomfy*$yfractor1)]
    set ypoib [expr round($wavezoomfy*$yfractor2)]
    set yt [expr int($ypoit/$singlewaveheight)+1];if {$yt<1} {set yt 1;}
    set yb [expr int($ypoib/$singlewaveheight)+1];if {$yb>$sfnum} {set yb $sfnum;}
    set y1 [expr $sfnum+1-$yb]
    set y2 [expr $sfnum+1-$yt]
    for {set i $y1} {$i<=$y2} {incr i 1} {
	set wavecolor($i) [.wavezoom itemcget wave.$dtname($i) -fill]
    }

    .wavezoomreplot xview moveto $xfractor1
    .wavezoomreplot yview moveto $yfractor1
    .wavezoomreplot delete notes
    .wavezoomreplot creat text [expr $xpoil+500] [expr $ypoit+1] -text "Press Key-space to return" -tags [list notes] -fill blue -font bold -anchor n;

    replotwavewin $xl $xr $y1 $y2
    loadwavezoomstainfopicktmark
    resetupwavezoomstainfopicktmark .wavezoomreplot
    plotscaleinwavezoom .wavezoomreplot
    plotsynthtmarkinwavezoom .wavezoomreplot
}
#################################################
proc replotwavewin {xl xr y1 y2} {
    global sfnum dtname sfdist sfaz
    global sfpt sfpick
    global sfyaxisposinwavezoom
    global sacdatwinb sacdatwine sacdat dt
    global sacdatenvelope sacdatinstfreq
    global sacdatfileb sacdatfilee
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global wavecolor
    global chanidcur
    global waveformsw envelopesw instfreqsw

    .wavezoomreplot delete wave statext waveenvelope waveinstfreq

    #set wavezoomx1int [expr int(($xl-$sacdatwinb)/$dt)+1]
    #set wavezoomx2int [expr int(($xr-$sacdatwinb)/$dt)+1]

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$i>=$y1 && $i<=$y2} {
	    # Adjust the begin time for each file
	    if {$xl<$sacdatfileb($dtname($i).$chanidcur)} {
		set datbegtime $sacdatfileb($dtname($i).$chanidcur)
	    } elseif {$xl>=$sacdatfileb($dtname($i).$chanidcur) && $xl<$sacdatfilee($dtname($i).$chanidcur)} {
		set datbegtime $xl
	    } elseif {$xl>=$sacdatfilee($dtname($i).$chanidcur)} {
		continue;
	    }

	    if {$xr<=$sacdatfileb($dtname($i).$chanidcur)} {
		continue;
	    } elseif {$xr>$sacdatfileb($dtname($i).$chanidcur) && $xr<=$sacdatfilee($dtname($i).$chanidcur)} {
		set datendtime $xr;
	    } elseif {$xr>$sacdatfilee($dtname($i).$chanidcur)} {
		set datendtime $sacdatfilee($dtname($i).$chanidcur)
	    }
	    set datbegint [expr int(($datbegtime-$sacdatfileb($dtname($i).$chanidcur))/$dt)+1]
	    set datendint [expr int(($datendtime-$sacdatfileb($dtname($i).$chanidcur))/$dt)+1]

	    # Find the first point x-val & y-val, save in lx,ly
	    set lx [wavezoomx2poix $datbegtime $wavezoomx1 $wavezoomxscale $wavezoombd]
	    set ly [expr $singlewaveheight*($sfnum+1-$i)]
	    # Copy windat
	    set sacdattem [lrange $sacdat($dtname($i).$chanidcur) [expr $datbegint-1] [expr $datendint-1]]

	    # Set amplitude scale
	    #set amx [minmaxabs $sacdattem]
	    #if {$amx>0} {
		#set poipam [expr $singlewaveheight/$amx]
	    #} else {
		#set poipam 0
	    #}
	    #set poipam [expr $poipam*1.0]

	    #set n 0
	    #set t $datbegtime
	    #set aa {}
	    #foreach point $sacdattem {
		#set xx($n) [expr $lx+$wavezoomxscale*($t-$datbegtime)]
		#set yy($n) [expr $ly-$poipam*$point]
		#lappend aa $xx($n) $yy($n)
		#incr n 1
		#set t [expr $t+$dt]
	    #}
	    #.wavezoomreplot create line $aa -width 1 -tags [list wave wave.$dtname($i)] -fill $wavecolor($i)
	    plotwavesub .wavezoomreplot $sacdattem $lx $ly \
                           $datbegtime $wavezoomxscale $dt $singlewaveheight wave $dtname($i) $wavecolor($i)

	    ####### plot envelope ########
	    set sacdattem [lrange $sacdatenvelope($dtname($i).$chanidcur) [expr $datbegint-1] [expr $datendint-1]]
	    plotwavesub .wavezoomreplot $sacdattem $lx $ly \
                           $datbegtime $wavezoomxscale $dt $singlewaveheight waveenvelope $dtname($i) orange
	    ####### plot instfreq ########
	    set sacdattem [lrange $sacdatinstfreq($dtname($i).$chanidcur) [expr $datbegint-1] [expr $datendint-1]]
	    plotwavesub .wavezoomreplot $sacdattem $lx $ly \
                           $datbegtime $wavezoomxscale $dt $singlewaveheight waveinstfreq $dtname($i) peru
	}
    }

    if {$waveformsw==0} {.wavezoomreplot itemconfigure wave         -state hidden}
    if {$envelopesw==0} {.wavezoomreplot itemconfigure waveenvelope -state hidden}
    if {$instfreqsw==0} {.wavezoomreplot itemconfigure waveinstfreq -state hidden}

}
#################################################
proc plotwavezoom {channel} {
    global sfnum dtname sfdist sfaz
    global sfpt sfpick
    global sfyaxisposinwavezoom
    global sacdatwinb sacdatwine sacdat dt
    global sacdatenvelope sacdatinstfreq
    global sacdatfileb sacdatfilee
    global waveallbd wavezoombd singlewaveheight
    global wavefx wavefy wavezoomfx wavezoomfy
    global waveallx1 waveallx2 waveally1 waveally2 waveallxscale waveallyscale
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global wavezoomdatbegtime wavezoomdatendtime wavezoombegint wavezoomendint
    global waveformsw envelopesw instfreqsw
    global readthreadnum calculatethreadnum
    global cursorcross

    .wavezoom delete wave statext line waveenvelope waveinstfreq

    set wavezoomfy [expr $singlewaveheight*($sfnum+1)]
    .wavezoom configure -scrollregion [list 0 0 $wavezoomfx $wavezoomfy]

    #set wavezoomx1int [expr int(($wavezoomx1-$sacdatwinb)/$dt)+1]
    #set wavezoomx2int [expr int(($wavezoomx2-$sacdatwinb)/$dt)+1]

    # set thread num & create thread, then set them to wait
    for {set i 0} {$i<$calculatethreadnum} {incr i 1} {
	set workid($i) [thread::create {package require Ttrace; thread::wait;}]
    }
    # obtain main thread id
    set mainid [thread::id]

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	# mark name & dist & az
	set disttext [format "D:%6.2f" $sfdist($i)]
	set aztext [format "AZ:%6.2f" $sfaz($dtname($i))]
	if {$sfpick($dtname($i).$channel)==1} {set color red}
	if {$sfpick($dtname($i).$channel)==0} {set color green3}
	set ly [expr $singlewaveheight*($sfnum+1-$i)]
	.wavezoom create text 0 $ly -text "$dtname($i) $disttext $aztext\n#$i" -tags [list statext statext.$dtname($i)] -fill $color -anchor w;

	# Adjust the begin time for each file
	#if {$wavezoomx1>=$sacdatfilee($dtname($i).$channel)} {continue}
	#if {$wavezoomx2<=$sacdatfileb($dtname($i).$channel)} {continue}
	if {$wavezoomx1<$sacdatfileb($dtname($i).$channel)} {
	    set wavezoomdatbegtime($dtname($i).$channel) $sacdatfileb($dtname($i).$channel)
	} elseif {$wavezoomx1>=$sacdatfileb($dtname($i).$channel) && $wavezoomx1<$sacdatfilee($dtname($i).$channel)} {
	    set wavezoomdatbegtime($dtname($i).$channel) $wavezoomx1
	} elseif {$wavezoomx1>=$sacdatfilee($dtname($i).$channel)} {
	    #continue; # for those waveform is not in window, then set a wrong begin time, which will be token care by lrange, return a bank element
	    set wavezoomdatbegtime($dtname($i).$channel) $wavezoomx1
	}

	if {$wavezoomx2<=$sacdatfileb($dtname($i).$channel)} {
	    #continue; # for those waveform is not in window, then set a wrong end time, which will be token care by lrange, return a bank element
	    set wavezoomdatendtime($dtname($i).$channel) $wavezoomx2
	} elseif {$wavezoomx2>$sacdatfileb($dtname($i).$channel) && $wavezoomx2<=$sacdatfilee($dtname($i).$channel)} {
	    set wavezoomdatendtime($dtname($i).$channel) $wavezoomx2
	} elseif {$wavezoomx2>$sacdatfilee($dtname($i).$channel)} {
	    set wavezoomdatendtime($dtname($i).$channel) $sacdatfilee($dtname($i).$channel)
	}
	set wavezoombegint($dtname($i).$channel) [expr int(($wavezoomdatbegtime($dtname($i).$channel)-$sacdatfileb($dtname($i).$channel))/$dt)+1]
	set wavezoomendint($dtname($i).$channel) [expr int(($wavezoomdatendtime($dtname($i).$channel)-$sacdatfileb($dtname($i).$channel))/$dt)+1]

	# Find the first point x-val & y-val, save in wavezoomfirpoix,wavezoomfirpoiy
	set wavezoomfirpoix [wavezoomx2poix $wavezoomdatbegtime($dtname($i).$channel) $wavezoomx1 $wavezoomxscale $wavezoombd]
	set wavezoomfirpoiy [expr $singlewaveheight*($sfnum+1-$i)]

	# Save vertical position in sfyaxispos
	set sfyaxisposinwavezoom($dtname($i).$channel) $wavezoomfirpoiy

	####### plot waveform ########
	#copy windat
	#set sacdattem [lrange $sacdat($dtname($i).$channel) [expr $wavezoombegint($dtname($i).$channel)-1] [expr $wavezoomendint($dtname($i).$channel)-1]]
	#plotwavesub .wavezoom $sacdattem $wavezoomfirpoix $wavezoomfirpoiy \
        #            $wavezoomdatbegtime($dtname($i).$channel) $wavezoomxscale $dt $singlewaveheight wave $dtname($i) gray

	####### plot envelope ########
	#set sacdattem [lrange $sacdatenvelope($dtname($i).$channel) [expr $wavezoombegint($dtname($i).$channel)-1] [expr $wavezoomendint($dtname($i).$channel)-1]]
	#plotwavesub .wavezoom $sacdattem $wavezoomfirpoix $wavezoomfirpoiy \
        #            $wavezoomdatbegtime($dtname($i).$channel) $wavezoomxscale $dt $singlewaveheight waveenvelope $dtname($i) orange
	#if {$envelopesw==0} {.wavezoom itemconfigure waveenvelope -state hidden}

	####### plot instfreq ########
	#set sacdattem [lrange $sacdatinstfreq($dtname($i).$channel) [expr $wavezoombegint($dtname($i).$channel)-1] [expr $wavezoomendint($dtname($i).$channel)-1]]
	#plotwavesub .wavezoom $sacdattem $wavezoomfirpoix $wavezoomfirpoiy \
        #            $wavezoomdatbegtime($dtname($i).$channel) $wavezoomxscale $dt $singlewaveheight waveinstfreq $dtname($i) peru
	#if {$instfreqsw==0} {.wavezoom itemconfigure waveinstfreq -state hidden}

	set imark [expr int(fmod($i,$calculatethreadnum))]
	####### plot waveform ########
	#copy windat
	set sacdattem [lrange $sacdat($dtname($i).$channel) [expr $wavezoombegint($dtname($i).$channel)-1] [expr $wavezoomendint($dtname($i).$channel)-1]]
	thread::send -async $workid($imark) "plotwavesubthread $mainid .wavezoom {$sacdattem} $wavezoomfirpoix $wavezoomfirpoiy \
                                             $wavezoomdatbegtime($dtname($i).$channel) $wavezoomxscale $dt $singlewaveheight wave $dtname($i) gray"
	####### plot envelope ########
	set sacdattem [lrange $sacdatenvelope($dtname($i).$channel) [expr $wavezoombegint($dtname($i).$channel)-1] [expr $wavezoomendint($dtname($i).$channel)-1]]
	thread::send -async $workid($imark) "plotwavesubthread $mainid .wavezoom {$sacdattem} $wavezoomfirpoix $wavezoomfirpoiy \
                                             $wavezoomdatbegtime($dtname($i).$channel) $wavezoomxscale $dt $singlewaveheight waveenvelope $dtname($i) orange"
	####### plot instfreq ########
	set sacdattem [lrange $sacdatinstfreq($dtname($i).$channel) [expr $wavezoombegint($dtname($i).$channel)-1] [expr $wavezoomendint($dtname($i).$channel)-1]]
	set singlewaveheighttem [expr $singlewaveheight*0.6]
	thread::send -async $workid($imark) "plotwavesubthread $mainid .wavezoom {$sacdattem} $wavezoomfirpoix $wavezoomfirpoiy \
                                             $wavezoomdatbegtime($dtname($i).$channel) $wavezoomxscale $dt $singlewaveheighttem waveinstfreq $dtname($i) peru"

	# Set amplitude scale
	#set amx [minmaxabs $sacdattem]
	#if {$amx>0} {
	#    set poipam [expr $singlewaveheight/$amx]
	#} else {
	#    set poipam 0
	#}
	#set poipam [expr $poipam*1.0]

	#set n 0
	#set t $wavezoomdatbegtime($dtname($i).$channel)
	#set aa {}
	#foreach point $sacdattem {
	#    set xx($n) [expr $wavezoomfirpoix+$wavezoomxscale*($t-$wavezoomdatbegtime($dtname($i).$channel))]
	#    set yy($n) [expr $wavezoomfirpoiy-$poipam*$point]
	#    lappend aa $xx($n) $yy($n)
	#    incr n 1
	#    set t [expr $t+$dt]
	#}
        #if {[llength $aa] <4} {continue;}
	#.wavezoom create line $aa -width 1 -tags [list wave wave.$dtname($i)] -fill gray



    }
    # send stop command to each thread
    for {set i 0} {$i<$calculatethreadnum} {incr i 1} {
	thread::send -async $workid($i) "thread::release $workid($i)"
    }
    while {[llength [thread::names]] > 1} {after 200}

    # show or hide envelope or instfreq
    thread::send -async $mainid "if {$waveformsw==0} {.wavezoom itemconfigure wave         -state hidden}"
    thread::send -async $mainid "if {$envelopesw==0} {.wavezoom itemconfigure waveenvelope -state hidden}"
    thread::send -async $mainid "if {$instfreqsw==0} {.wavezoom itemconfigure waveinstfreq -state hidden}"
    # The next two line work in an unproper queue, It should be work at last
    #if {$envelopesw==0} {.wavezoom itemconfigure waveenvelope -state hidden}
    #if {$instfreqsw==0} {.wavezoom itemconfigure waveinstfreq -state hidden}

    .wavezoom create line 0 0 $wavezoomfx 0 -width 1 -tags [list linex1 linexy line] -fill green -state $cursorcross
    .wavezoom create line 0 0 $wavezoomfx 0 -width 1 -tags [list linex2 linexy line] -fill green -state $cursorcross
    .wavezoom create line 0 0 0 $wavezoomfy -width 1 -tags [list liney1 linexy line] -fill green -state $cursorcross
    .wavezoom create line 0 0 0 $wavezoomfy -width 1 -tags [list liney2 linexy line] -fill green -state $cursorcross
}
#################################################
proc plotwavesub {canvasname datlist firpoix firpoiy firpoitime xscale dt traceheight tagmatchname kstnm color} {

    # Set amplitude scale
    set amx [minmaxabs $datlist]
    if {$amx>0} {
	set poipam [expr $traceheight/$amx]
    } else {
	set poipam 0
    }
    set poipam [expr $poipam*1.0]
    set n 0
    set t $firpoitime
    set aa {}
    foreach point $datlist {
	set xx($n) [expr $firpoix+$xscale*($t-$firpoitime)]
	set yy($n) [expr $firpoiy-$poipam*$point]
	lappend aa $xx($n) $yy($n)
	incr n 1
	set t [expr $t+$dt]
    }
    if {[llength $aa] <4} {return;}
    $canvasname create line $aa -width 1 -tags [list $tagmatchname $tagmatchname.$kstnm] -fill $color
}
#################################################
proc waveformswsub {} {
    global chanidcur
    global waveformsw

    if {[.wavezoom find withtag wave]==""} {return}
    if {$waveformsw==0} {
	.wavezoom itemconfigure wave -state hidden
    } else {
	.wavezoom itemconfigure wave -state normal
    }

    if {[.wavezoomreplot find withtag wave]==""} {return}
    if {$waveformsw==0} {
	.wavezoomreplot itemconfigure wave -state hidden
    } else {
	.wavezoomreplot itemconfigure wave -state normal
    }
}
proc envelopeswsub {} {
    global chanidcur
    global envelopesw
    global sacdatenvelope
    global wavezoomdatbegtime wavezoomdatendtime wavezoombegint wavezoomendint
    global sfnum dtname sfdist sfaz
    global wavezoomx1 wavezoomx2 wavezoomy1 wavezoomy2 wavezoomxscale wavezoomyscale
    global dt 
    global waveallbd wavezoombd singlewaveheight

    if {[.wavezoom find withtag wave]==""} {return}
    if {$envelopesw==0} {
	#.wavezoom delete waveenvelope
	.wavezoom itemconfigure waveenvelope -state hidden
    } else {
	.wavezoom itemconfigure waveenvelope -state normal
	#for {set i 1} {$i<=$sfnum} {incr i 1} {
	#    set sacdattem [lrange $sacdatenvelope($dtname($i).$chanidcur) [expr $wavezoombegint($dtname($i).$chanidcur)-1] [expr $wavezoomendint($dtname($i).$chanidcur)-1]]
	#    set firpoix [wavezoomx2poix $wavezoomdatbegtime($dtname($i).$chanidcur) $wavezoomx1 $wavezoomxscale $wavezoombd]
	#    set firpoiy [expr $singlewaveheight*($sfnum+1-$i)]
	#    plotwavesub .wavezoom $sacdattem $firpoix $firpoiy $wavezoomdatbegtime($dtname($i).$chanidcur) $wavezoomxscale $dt $singlewaveheight waveenvelope $dtname($i) orange
	#}
    }

    if {[.wavezoomreplot find withtag wave]==""} {return}
    if {$envelopesw==0} {
	.wavezoomreplot itemconfigure waveenvelope -state hidden
    } else {
	.wavezoomreplot itemconfigure waveenvelope -state normal
    }
}
proc instfreqswsub {} {
    global chanidcur
    global instfreqsw

    if {[.wavezoom find withtag wave]==""} {return}
    if {$instfreqsw==0} {
	.wavezoom itemconfigure waveinstfreq -state hidden
    } else {
	.wavezoom itemconfigure waveinstfreq -state normal
    }

    if {[.wavezoomreplot find withtag wave]==""} {return}
    if {$instfreqsw==0} {
	.wavezoomreplot itemconfigure waveinstfreq -state hidden
    } else {
	.wavezoomreplot itemconfigure waveinstfreq -state normal
    }
}
#################################################
proc minmax a {
    set min  1.e+20
    set max -1.e+20
    foreach element $a  {
        if {$element > $max} {set max $element}
        if {$element < $min} {set min $element}
    }
    return [expr $max-$min]
}
proc minmaxabs {a} {
    set min  1.e+20
    set max -1.e+20
    foreach element $a  {
	set absval [expr abs($element)]
        if {$absval > $max} {set max $absval}
    }
    return $max
}
#################################################
proc caldistindeg {lon lat} {
    global sfnum sfname sflat sflon sfele sfpick sfpt
    global sfgcarc sfdist sfbaz sfaz
    global readthreadnum calculatethreadnum

    #for {set i 1} {$i<=$sfnum} {incr i 1} {
	#set line [distaztcl $sflat($sfname($i)) $sflon($sfname($i)) $lat $lon]
	#set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	#set linelist [split $line " "];
	#set sfgcarc($sfname($i)) [lindex $linelist 0];
	#set sfbaz($sfname($i)) [lindex $linelist 1];
	#set sfaz($sfname($i)) [lindex $linelist 2];
    #}

    # set thread num & create thread, then set them to wait
    for {set i 0} {$i<$calculatethreadnum} {incr i 1} {
	set workid($i) [thread::create {package require Ttrace; thread::wait;}]
    }
    # obtain main thread id
    set mainid [thread::id]

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set imark [expr int(fmod($i,$calculatethreadnum))]
	thread::send -async $workid($imark) "caldistindegsub $sflat($sfname($i)) $sflon($sfname($i)) $lat $lon $sfname($i)"
    }

    # send stop command to each thread
    for {set i 0} {$i<$calculatethreadnum} {incr i 1} {
	thread::send -async $workid($i) "thread::release $workid($i)"
    }
    while {[llength [thread::names]] > 1} {after 200}

    #copy dat from those (in thread) to those (in main thread)
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set sfgcarc($sfname($i)) [tsv::get sfgcarc $sfname($i)];
	set sfbaz($sfname($i))   [tsv::get sfbaz   $sfname($i)];
	set sfaz($sfname($i))    [tsv::get sfaz    $sfname($i)];
    }
}
ttrace::eval {
    proc caldistindegsub {stla stlo evla evlo kstnm} {
	set line [distaztcl $stla $stlo $evla $evlo]
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	set linelist [split $line " "];
	tsv::set sfgcarc $kstnm [lindex $linelist 0];
	tsv::set sfbaz   $kstnm [lindex $linelist 1];
	tsv::set sfaz    $kstnm [lindex $linelist 2];
    }
}
#################################################
proc sortbydist {} {
    global sfnum sfname sflat sflon sfele sfpick sfpt
    global sfgcarc sfdist sfbaz sfaz
    global dtname

    #copy dtname from sfname
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set dtname($i) $sfname($i)
    }

    # sort by distance
    for {set i 1} {$i<= [expr $sfnum-1]} {incr i 1} {
	for {set j [expr $i+1]} {$j<=$sfnum} {incr j 1} {
	    if {$sfgcarc($dtname($j))<$sfgcarc($dtname($i))} {
		set tem $dtname($j);
		set dtname($j) $dtname($i);
		set dtname($i) $tem;
	    }
	}
    }

    # convert to km
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set sfdist($i) [expr $sfgcarc($dtname($i))*6371.0*acos(-1.0)/180.0]
	#puts $dtname($i),$sfgcarc($dtname($i)),$sfdist($i);
    }
}
#################################################
proc readsacdat {channel} {
    global datdir evtidcur
    global shiftmark shift1 shift2
    global year jday hour min sec msec
    global sfnum dtname
    global sacdatwinb sacdatwine sacdat dt
    global sacdatenvelope sacdatinstfreq
    global sacdatfileb sacdatfilee
    global readthreadnum calculatethreadnum

    # set thread num & create thread, then set them to wait
    for {set i 0} {$i<$readthreadnum} {incr i 1} {
	set workid($i) [thread::create {package require Ttrace; thread::wait;}]
    }
    # obtain main thread id
    set mainid [thread::id]

    set sacdatwinb($channel) 10000000.0
    set sacdatwine($channel) -10000000.0
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set sacfile [glob -nocomplain swap/$evtidcur.$dtname($i).$channel.SAC]
	if {[file exists $sacfile]!=0} {
	    #set line [exec saclst b kztime delta e $shiftmark f $sacfile]
	    set line [saclsttcl kztime delta b e $shiftmark f $sacfile]
	    set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	    set linelist [split $line " "]
	    set t1 [lindex $linelist 1]; set t2 [lindex $linelist 2];
	    set t3 [lindex $linelist 3]; set t4 [lindex $linelist 4];
	    set t5 [lindex $linelist 5]; set t6 [lindex $linelist 6];
	    set t7 [lindex $linelist 7];
	    set sacdatbval [lindex $linelist  8];
	    set sacdateval [lindex $linelist  9];
	    set tmarkval   [lindex $linelist 10];
	}

	set errsw 0
	set t1 [scan $t1 "%d"]; set t2 [scan $t2 "%d"];
	set t3 [scan $t3 "%d"]; set t4 [scan $t4 "%d"];
	set t5 [scan $t5 "%d"]; set t6 [scan $t6 "%d"]; #this line is due to those in getsacinfo, if not scan again, [expr $sec+$msec/1000.0] will show an error message
	if {$t1!=$year || $t2!=$jday || $t3!=$hour} {set errsw 1;}
	if {$t4!=$min  || $t5!=$sec  || $t6!=$msec} {set errsw 1;}
	set t7tem [format "%.3f" $t7]; set dttem [format "%.3f" $dt]; # this line is due to presicion error
	if {$t7tem!=$dttem} {set errsw 1;}
	if {$errsw==1} {
	    tk_messageBox -type ok -message "SACfile $sacfile have different refer-time or data sample, please check and confirm";
	    return 1
	}
	set sacdatfileb($dtname($i).$channel) [expr $tmarkval+$shift1]
	set sacdatfilee($dtname($i).$channel) [expr $tmarkval+$shift2]
	# Find time window
	if {$sacdatfileb($dtname($i).$channel)<=$sacdatwinb($channel)} {set sacdatwinb($channel) $sacdatfileb($dtname($i).$channel)}
	if {$sacdatfilee($dtname($i).$channel)>=$sacdatwine($channel)} {set sacdatwine($channel) $sacdatfilee($dtname($i).$channel)}

	# Adjust dat beg/end-time for each file
	if {$sacdatfileb($dtname($i).$channel)<$sacdatbval} {set sacdatfileb($dtname($i).$channel) $sacdatbval}
	if {$sacdatfilee($dtname($i).$channel)>$sacdateval} {set sacdatfilee($dtname($i).$channel) $sacdateval}


	set imark [expr int(fmod($i,$readthreadnum))]
	## read sacfile dat
	thread::send -async $workid($imark) "readsacdatthreadsub $sacfile $sacdatfileb($dtname($i).$channel) $sacdatfilee($dtname($i).$channel) $dtname($i) sacfile"
	## read envelope dat
	set sacfileenvelope $sacfile.ENVELOPE
	thread::send -async $workid($imark) "readsacdatthreadsub $sacfileenvelope $sacdatfileb($dtname($i).$channel) $sacdatfilee($dtname($i).$channel) $dtname($i) envelope"
	## read instfreq dat
	set sacfileinstfreq $sacfile.instfreq
	thread::send -async $workid($imark) "readsacdatthreadsub $sacfileinstfreq $sacdatfileb($dtname($i).$channel) $sacdatfilee($dtname($i).$channel) $dtname($i) instfreq"

	#continue; # the next code is for reading dat through main thread
	##########################
	## read sacfile dat
	#set mccrsacout [exec mccrsac $sacfile $sacdatfileb($dtname($i).$channel) $sacdatfilee($dtname($i).$channel)]
	#set mccrsacout [readsacdatsub $sacfile $sacdatfileb($dtname($i).$channel) $sacdatfilee($dtname($i).$channel)]
	#set sacdat($dtname($i).$channel) [lreplace $mccrsacout 0 1]

	## read envelope dat
	#set sacfileenvelope $sacfile.ENVELOPE
	#set mccrsacout [exec mccrsac $sacfileenvelope $sacdatfileb($dtname($i).$channel) $sacdatfilee($dtname($i).$channel)]
	#set mccrsacout [readsacdatsub $sacfileenvelope $sacdatfileb($dtname($i).$channel) $sacdatfilee($dtname($i).$channel)]
	#set sacdatenvelope($dtname($i).$channel) [lreplace $mccrsacout 0 1]

	## read instfreq dat
	#set sacfileinstfreq $sacfile.instfreq
	#set mccrsacout [exec mccrsac $sacfileinstfreq $sacdatfileb($dtname($i).$channel) $sacdatfilee($dtname($i).$channel)]
	#set mccrsacout [readsacdatsub $sacfileinstfreq $sacdatfileb($dtname($i).$channel) $sacdatfilee($dtname($i).$channel)]
	#set sacdatinstfreq($dtname($i).$channel) [lreplace $mccrsacout 0 1]

    }
    # send stop command to each thread
    for {set i 0} {$i<$readthreadnum} {incr i 1} {
	thread::send -async $workid($i) "thread::release $workid($i)"
    }
    while {[llength [thread::names]] > 1} {after 200}

    # copy dat from those (in thread) to those (in main thread)
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set sacdat($dtname($i).$channel)         [tsv::get threaddatout $dtname($i).sacfile]
	set sacdatenvelope($dtname($i).$channel) [tsv::get threaddatout $dtname($i).envelope]
	set sacdatinstfreq($dtname($i).$channel) [tsv::get threaddatout $dtname($i).instfreq]
    }
    return 0
}
ttrace::eval {
    proc readsacdatthreadsub {filename begtime endtime kstnm type} {
	#set mccrsacout [exec mccrsac $filename $begtime $endtime]
	set mccrsacout [readsacdatsub $filename $begtime $endtime]
	tsv::set threaddatout $kstnm.$type [lreplace $mccrsacout 0 1]
    }
}
#################################################
ttrace::eval {
proc saclsttcl2 {args} {
    set sacfile [lindex $args end]
    set newargs [lreplace $args end-1 end]

    set f [open $sacfile r]
    fconfigure $f -translation binary
    set aa $sacfile
    foreach element $newargs {
	if {$element == "evla"} {seek $f 140; binary scan [read $f 4] f tt;lappend aa $tt;}
	if {$element == "evlo"} {seek $f 144; binary scan [read $f 4] f tt;lappend aa $tt;}
	if {$element == "evdp"} {seek $f 152; binary scan [read $f 4] f tt;lappend aa $tt;}
	if {$element == "kstnm"} {seek $f 440; binary scan [read $f 8] f tt; set tt [string trim $tt];lappend aa $tt;}
    }
    close $f
    return $aa
}
}
#################################################
ttrace::eval {
proc saclsttcl {args} {
    set sacfile [lindex $args end]
    set newargs [lreplace $args end-1 end]
    
    set f [open $sacfile r]
    fconfigure $f -translation binary
    binary scan [read $f 440] f70i40 varfloat varinteger
    binary scan [read $f 192] a8a16a8a8a8a8a8a8a8a8a8a8a8a8a8a8a8a8a8a8a8a8a8 \
                 kstnm kevnm khole ko ka kt0 kt1 kt2 kt3 kt4 kt5 kt6 kt7 kt8 kt9 kf \
                 kuser0 kuser1 kuser2 kcmpnm knetwk kdatrd kinst
    #binary scan [read $f 24] a8a16  kstnm kevnm
    #binary scan [read $f 24] a8a8a8 khole ko ka
    #binary scan [read $f 24] a8a8a8 kt0 kt1 kt2
    #binary scan [read $f 24] a8a8a8 kt3 kt4 kt5
    #binary scan [read $f 24] a8a8a8 kt6 kt7 kt8
    #binary scan [read $f 24] a8a8a8 kt9 kf  kuser0
    #binary scan [read $f 24] a8a8a8 kuser1 kuser2 kcmpnm
    #binary scan [read $f 24] a8a8a8 knetwk kdatrd kinst
    close $f
    set aa $sacfile
    foreach element $newargs {
	set index [saclsttclsub $element]
	if {$index>=0 && $index<=69} {
	    set tt [lindex $varfloat $index]
	    lappend aa $tt
	}
	if {$index>=70 && $index<=109} {
	    set tt [lindex $varinteger [expr $index-70]]
	    lappend aa $tt
	}
	#if {$index==110} {set kstnm  [exec echo $kstnm  | tr "\n" " "]; set kstnm  [string trim $kstnm];  lappend aa $kstnm}
	# trim should replace "\n", however, it doesn't work in tcl 8.5.13
	if {$index==110} {set kstnm  [string trim $kstnm];  lappend aa $kstnm}
	if {$index==112} {set kevnm  [string trim $kevnm];  lappend aa $kevnm}
	if {$index==116} {set khole  [string trim $khole];  lappend aa $khole}
	if {$index==118} {set ko     [string trim $ko];     lappend aa $ko}
	if {$index==120} {set ka     [string trim $ka];     lappend aa $ka}
	if {$index==122} {set kt0    [string trim $kt0];    lappend aa $kt0}
	if {$index==124} {set kt1    [string trim $kt1];    lappend aa $kt1}
	if {$index==126} {set kt2    [string trim $kt2];    lappend aa $kt2}
	if {$index==128} {set kt3    [string trim $kt3];    lappend aa $kt3}
	if {$index==130} {set kt4    [string trim $kt4];    lappend aa $kt4}
	if {$index==132} {set kt5    [string trim $kt5];    lappend aa $kt5}
	if {$index==134} {set kt6    [string trim $kt6];    lappend aa $kt6}
	if {$index==136} {set kt7    [string trim $kt7];    lappend aa $kt7}
	if {$index==138} {set kt8    [string trim $kt8];    lappend aa $kt8}
	if {$index==140} {set kt9    [string trim $kt9];    lappend aa $kt9}
	if {$index==142} {set kf     [string trim $kf];     lappend aa $kf}
	if {$index==144} {set kuser0 [string trim $kuser0]; lappend aa $kuser0}
	if {$index==146} {set kuser1 [string trim $kuser1]; lappend aa $kuser1}
	if {$index==148} {set kuser2 [string trim $kuser2]; lappend aa $kuser2}
	if {$index==150} {set kcmpnm [string trim $kcmpnm]; lappend aa $kcmpnm}
	if {$index==152} {set knetwk [string trim $knetwk]; lappend aa $knetwk}
	if {$index==154} {set kdatrd [string trim $kdatrd]; lappend aa $kdatrd}
	if {$index==156} {set kinst  [string trim $kinst];  lappend aa $kinst}
	if {$index==200} {
	    set nzyear [lindex $varinteger 0]; set nzyear [format "%04d" $nzyear]
	    set nzjday [lindex $varinteger 1]; set nzjday [format "%03d" $nzjday]
	    set nzhour [lindex $varinteger 2]; set nzhour [format "%02d" $nzhour]
	    set nzmin  [lindex $varinteger 3]; set nzmin  [format "%02d" $nzmin]
	    set nzsec  [lindex $varinteger 4]; set nzsec  [format "%02d" $nzsec]
	    set nzmsec [lindex $varinteger 5]; set nzmsec [format "%03d" $nzmsec]
	    lappend aa $nzyear $nzjday $nzhour $nzmin $nzsec $nzmsec
	}
    }
    return $aa
}
#################################################
proc saclsttclsub {headname} {
    if {[string compare -nocase $headname "delta"] == 0} {return 0}
    if {[string compare -nocase $headname "depmin"] == 0} {return 1}
    if {[string compare -nocase $headname "depmax"] == 0} {return 2}
    if {[string compare -nocase $headname "scale"] == 0} {return 3}
    if {[string compare -nocase $headname "odelta"] == 0} {return 4}
    if {[string compare -nocase $headname "b"] == 0} {return 5}
    if {[string compare -nocase $headname "e"] == 0} {return 6}
    if {[string compare -nocase $headname "o"] == 0} {return 7}
    if {[string compare -nocase $headname "a"] == 0} {return 8}
    if {[string compare -nocase $headname "internal1"] == 0} {return 9}

    if {[string compare -nocase $headname "t0"] == 0} {return 10}
    if {[string compare -nocase $headname "t1"] == 0} {return 11}
    if {[string compare -nocase $headname "t2"] == 0} {return 12}
    if {[string compare -nocase $headname "t3"] == 0} {return 13}
    if {[string compare -nocase $headname "t4"] == 0} {return 14}
    if {[string compare -nocase $headname "t5"] == 0} {return 15}
    if {[string compare -nocase $headname "t6"] == 0} {return 16}
    if {[string compare -nocase $headname "t7"] == 0} {return 17}
    if {[string compare -nocase $headname "t8"] == 0} {return 18}
    if {[string compare -nocase $headname "t9"] == 0} {return 19}

    if {[string compare -nocase $headname "F"] == 0} {return 20}
    if {[string compare -nocase $headname "resp0"] == 0} {return 21}
    if {[string compare -nocase $headname "resp1"] == 0} {return 22}
    if {[string compare -nocase $headname "resp2"] == 0} {return 23}
    if {[string compare -nocase $headname "resp3"] == 0} {return 24}
    if {[string compare -nocase $headname "resp4"] == 0} {return 25}
    if {[string compare -nocase $headname "resp5"] == 0} {return 26}
    if {[string compare -nocase $headname "resp6"] == 0} {return 27}
    if {[string compare -nocase $headname "resp7"] == 0} {return 28}
    if {[string compare -nocase $headname "resp8"] == 0} {return 29}

    if {[string compare -nocase $headname "resp9"] == 0} {return 30}
    if {[string compare -nocase $headname "stla"] == 0} {return 31}
    if {[string compare -nocase $headname "stlo"] == 0} {return 32}
    if {[string compare -nocase $headname "stel"] == 0} {return 33}
    if {[string compare -nocase $headname "stdp"] == 0} {return 34}
    if {[string compare -nocase $headname "evla"] == 0} {return 35}
    if {[string compare -nocase $headname "evlo"] == 0} {return 36}
    if {[string compare -nocase $headname "evel"] == 0} {return 37}
    if {[string compare -nocase $headname "evdp"] == 0} {return 38}
    if {[string compare -nocase $headname "mag"] == 0} {return 39}

    if {[string compare -nocase $headname "user0"] == 0} {return 40}
    if {[string compare -nocase $headname "user1"] == 0} {return 41}
    if {[string compare -nocase $headname "user2"] == 0} {return 42}
    if {[string compare -nocase $headname "user3"] == 0} {return 43}
    if {[string compare -nocase $headname "user4"] == 0} {return 44}
    if {[string compare -nocase $headname "user5"] == 0} {return 45}
    if {[string compare -nocase $headname "user6"] == 0} {return 46}
    if {[string compare -nocase $headname "user7"] == 0} {return 47}
    if {[string compare -nocase $headname "user8"] == 0} {return 48}
    if {[string compare -nocase $headname "user9"] == 0} {return 49}

    if {[string compare -nocase $headname "dist"] == 0} {return 50}
    if {[string compare -nocase $headname "az"] == 0} {return 51}
    if {[string compare -nocase $headname "baz"] == 0} {return 52}
    if {[string compare -nocase $headname "gcarc"] == 0} {return 53}
    if {[string compare -nocase $headname "internal2"] == 0} {return 54}
    if {[string compare -nocase $headname "internal3"] == 0} {return 55}
    if {[string compare -nocase $headname "depmen"] == 0} {return 56}
    if {[string compare -nocase $headname "cmpaz"] == 0} {return 57}
    if {[string compare -nocase $headname "cmpinc"] == 0} {return 58}
    if {[string compare -nocase $headname "XMINIMUM"] == 0} {return 59}

    if {[string compare -nocase $headname "XMAXIMUM"] == 0} {return 60}
    if {[string compare -nocase $headname "YMINIMUM"] == 0} {return 61}
    if {[string compare -nocase $headname "YMAXIMUM"] == 0} {return 62}
    if {[string compare -nocase $headname "ADJTM"] == 0} {return 63}
    if {[string compare -nocase $headname "UNUSED7"] == 0} {return 64}
    if {[string compare -nocase $headname "UNUSED8"] == 0} {return 65}
    if {[string compare -nocase $headname "UNUSED9"] == 0} {return 66}
    if {[string compare -nocase $headname "UNUSED10"] == 0} {return 67}
    if {[string compare -nocase $headname "UNUSED11"] == 0} {return 68}
    if {[string compare -nocase $headname "UNUSED12"] == 0} {return 69}

    if {[string compare -nocase $headname "nzyear"] == 0} {return 70}
    if {[string compare -nocase $headname "nzjday"] == 0} {return 71}
    if {[string compare -nocase $headname "nzhour"] == 0} {return 72}
    if {[string compare -nocase $headname "nzmin"] == 0} {return 73}
    if {[string compare -nocase $headname "nzsec"] == 0} {return 74}
    if {[string compare -nocase $headname "nzmsec"] == 0} {return 75}
    if {[string compare -nocase $headname "nvhdr"] == 0} {return 75}
    if {[string compare -nocase $headname "NORID"] == 0} {return 77}
    if {[string compare -nocase $headname "NEVID"] == 0} {return 78}
    if {[string compare -nocase $headname "NPTS"] == 0} {return 79}

    if {[string compare -nocase $headname "NSPTS"] == 0} {return 80}
    if {[string compare -nocase $headname "NWFID"] == 0} {return 81}
    if {[string compare -nocase $headname "NXSIZE"] == 0} {return 82}
    if {[string compare -nocase $headname "NYSIZE"] == 0} {return 83}
    if {[string compare -nocase $headname "UNUSED15"] == 0} {return 84}
    if {[string compare -nocase $headname "IFTYPE"] == 0} {return 85}
    if {[string compare -nocase $headname "IDEP"] == 0} {return 86}
    if {[string compare -nocase $headname "IZTYPE"] == 0} {return 87}
    if {[string compare -nocase $headname "UNUSED16"] == 0} {return 88}
    if {[string compare -nocase $headname "IINST"] == 0} {return 89}

    if {[string compare -nocase $headname "istreg"] == 0} {return 90}
    if {[string compare -nocase $headname "IEVREG"] == 0} {return 91}
    if {[string compare -nocase $headname "IEVTYP"] == 0} {return 92}
    if {[string compare -nocase $headname "iqual"] == 0} {return 93}
    if {[string compare -nocase $headname "isynth"] == 0} {return 94}
    if {[string compare -nocase $headname "IMAGTYP"] == 0} {return 95}
    if {[string compare -nocase $headname "IMAGSRC"] == 0} {return 96}
    if {[string compare -nocase $headname "UNUSED19"] == 0} {return 97}
    if {[string compare -nocase $headname "UNUSED20"] == 0} {return 98}
    if {[string compare -nocase $headname "UNUSED21"] == 0} {return 99}

    if {[string compare -nocase $headname "unused22"] == 0} {return 100}
    if {[string compare -nocase $headname "unused23"] == 0} {return 101}
    if {[string compare -nocase $headname "unused24"] == 0} {return 102}
    if {[string compare -nocase $headname "unused25"] == 0} {return 103}
    if {[string compare -nocase $headname "unused26"] == 0} {return 104}
    if {[string compare -nocase $headname "leven"] == 0} {return 105}
    if {[string compare -nocase $headname "lpspol"] == 0} {return 106}
    if {[string compare -nocase $headname "lovrok"] == 0} {return 107}
    if {[string compare -nocase $headname "lcalda"] == 0} {return 108}
    if {[string compare -nocase $headname "unused27"] == 0} {return 109}

    if {[string compare -nocase $headname "kstnm"] == 0} {return 110}
    if {[string compare -nocase $headname "kevnm"] == 0} {return 112}
    if {[string compare -nocase $headname "khole"] == 0} {return 116}
    if {[string compare -nocase $headname "ko"] == 0} {return 118}
    if {[string compare -nocase $headname "ka"] == 0} {return 120}
    if {[string compare -nocase $headname "kt0"] == 0} {return 122}
    if {[string compare -nocase $headname "kt1"] == 0} {return 124}
    if {[string compare -nocase $headname "kt2"] == 0} {return 126}
    if {[string compare -nocase $headname "kt3"] == 0} {return 128}
    if {[string compare -nocase $headname "kt4"] == 0} {return 130}

    if {[string compare -nocase $headname "kt5"] == 0} {return 132}
    if {[string compare -nocase $headname "kt6"] == 0} {return 134}
    if {[string compare -nocase $headname "kt7"] == 0} {return 136}
    if {[string compare -nocase $headname "kt8"] == 0} {return 138}
    if {[string compare -nocase $headname "kt9"] == 0} {return 140}
    if {[string compare -nocase $headname "kf"] == 0} {return 142}
    if {[string compare -nocase $headname "kuser0"] == 0} {return 144}
    if {[string compare -nocase $headname "kuser1"] == 0} {return 146}
    if {[string compare -nocase $headname "kuser2"] == 0} {return 148}
    if {[string compare -nocase $headname "kcmpnm"] == 0} {return 150}

    if {[string compare -nocase $headname "knetwk"] == 0} {return 152}
    if {[string compare -nocase $headname "kdatrd"] == 0} {return 154}
    if {[string compare -nocase $headname "kinst"] == 0} {return 156}
    if {[string compare -nocase $headname "kztime"] == 0} {return 200}
}
}
#################################################
ttrace::eval {
#read SAC data between tmin and tmax
proc readsacdatsub {name tmin tmax} {

    set f [open $name r]
    fconfigure $f -translation binary
    binary scan [read $f 440] f70i40 varfloat varinteger

    #binary scan [read $f 280] f70 varfloat
    #binary scan [read $f 160] i40 varinteger
    #binary scan [read $f 24] a8a16  kstnm kevnm
    #binary scan [read $f 24] a8a8a8 khole ko ka
    #binary scan [read $f 24] a8a8a8 kt0 kt1 kt2
    #binary scan [read $f 24] a8a8a8 kt3 kt4 kt5
    #binary scan [read $f 24] a8a8a8 kt6 kt7 kt8
    #binary scan [read $f 24] a8a8a8 kt9 kf  kuser0
    #binary scan [read $f 24] a8a8a8 kuser1 kuser2 kcmpnm
    #binary scan [read $f 24] a8a8a8 knetwk kdatrd kinst

    set dt   [lindex $varfloat   0]
    set b    [lindex $varfloat   5]
    set e    [lindex $varfloat   6]
    set npts [lindex $varinteger 9]

    seek $f 632
    binary scan [read $f [expr $npts*4]] f* vardat
    #if {$tmin<$b || $tmax>$e} {puts "error"}

    set t1num [expr int(($tmin-$b)/$dt)+1]
    set t2num [expr int(($tmax-$b)/$dt)+1]
    set dat [lrange $vardat [expr $t1num-1] [expr $t2num-1]]
    close $f

    set t0 [expr $b+($t1num-1)*$dt]; ####### t0 may have a little shift from tmin
    return "$dt $t0 $dat";
}
}
#################################################
proc loadevtpos {channel sw} {
    # sw==0, load from SACfile; sw==1, load from EVENT.RLOCfile; sw==2, load from EVENT.FINALfile
    # in case of sw==0, return 0, no error, return 1, error

    global datdir outdir evtidcur
    global evtdatesac evttimesac evtomarksac evtlonsac evtlatsac evtdepsac
    global evtdatesd  evttimesd  evtomarksd  evtlonsd  evtlatsd  evtdepsd evtlonrms evtlatrms evtdeprms
    global evtdatefin evttimefin evtomarkfin evtlonfin evtlatfin evtdepfin
    global evtherrsd evtverrsd evtherrrms evtverrrms

    if {$sw==0} {
	set sacfilelist [glob -nocomplain swap/$evtidcur.*.$channel.SAC]
	set sacfilesingle [lindex $sacfilelist 0]
	set line [saclsttcl evlo evla evdp f $sacfilesingle]
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	set linelist [split $line " "];
	set evtlonsac [lindex $linelist 1];
	set evtlatsac [lindex $linelist 2];
	set evtdepsac [lindex $linelist 3];
	if {$evtlonsac==-12345.0 || $evtlatsac==-12345.0 || $evtdepsac==-12345.0} {
	    tk_messageBox -type ok -message "Maybe, event coordination is not defined\nForce all of them to be 0.0\nPlease Note!";
	    set evtlonsac 0.0; set evtlatsac 0.0; set evtdepsac 0.0;
#	    return 1
	}
 	return 0
    } elseif {$sw==1} {
	set eventfile "$outdir/$evtidcur.$channel.EVENT.RLOC";
	#if {[file exists $eventfile]==0} {tk_messageBox -type ok -message "$outdir/$evtidcur.$channel.EVENT.RLOC does not exist"}
	set line [exec grep $evtidcur $eventfile]; regsub -all {[[:blank:]]+} $line " " line
	set evtlonsd [lindex $line 3]
	set evtlatsd [lindex $line 4]
	set evtdepsd [lindex $line 5]
	set evtlonrms [lindex $line 6]
	set evtlatrms [lindex $line 7]
	set evtdeprms [lindex $line 8]
	set evtherrsd [lindex $line 13]
	set evtverrsd [lindex $line 14]
	set evtherrrms [lindex $line 15]
	set evtverrrms [lindex $line 16]
    } elseif {$sw==2} {
	set eventfile "$outdir/$evtidcur.$channel.EVENT.FINAL"
	set line [exec grep $evtidcur $eventfile]; regsub -all {[[:blank:]]+} $line " " line
	set evtlonfin [lindex $line 3]
	set evtlatfin [lindex $line 4]
	set evtdepfin [lindex $line 5]
    }
}
#################################################
proc loadevtdatetime {channel sw} {
    # sw==0, load from SACfile; sw==1, load from EVENT.RLOCfile; sw==2, load from EVENT.FINALfile
    global datdir outdir evtidcur
    global year jday hour min sec msec
    global shiftmark
    global evtdatesac evttimesac evtomarksac evtlonsac evtlatsac evtdepsac
    global evtdatesd  evttimesd  evtomarksd  evtlonsd  evtlatsd  evtdepsd evtlonrms evtlatrms evtdeprms
    global evtdatefin evttimefin evtomarkfin evtlonfin evtlatfin evtdepfin

    if {$sw==0} {
	set sacfilelist [glob -nocomplain swap/$evtidcur.*.$channel.SAC]
	set sacfilesingle [lindex $sacfilelist 0]
	set line [saclsttcl o f $sacfilesingle]
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	set linelist [split $line " "]
	# set evtomarksac
	set evtomarksac [lindex $linelist 1];
	if {$evtomarksac==-12345.0} {
	    #### if omark is undefine in sac file, show a message
	    set message "Maybe, event time is not defined. Force them to be\n 0000/00/00 00:00:00.000";
	    set message "$message\nSynthetic time-mark will be marked at a wrong position.";
	    set message "$message\nIf shift mark is omark, then it will be automatically set to b.\nPlease Note!";
	    tk_messageBox -type ok -message $message;
	    if {$shiftmark=="o"} {.lshiftmark configure -text "shift mark: b"; set shiftmark "b";}
	    #set evtomarksac 0.0;
	    set evtdatesac "0000/00/00";
	    set evttimesac "00:00:00.000";
	} else {
	    set secs [expr $sec+$msec/1000.0]
	    set line [taddtcl $year $jday $hour $min $secs $evtomarksac]
	    set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	    set evtyear [lindex $line 0]; set evtjday [lindex $line 1]; set evthour [lindex $line 2];
	    set evtmin  [lindex $line 3]; set secline [lindex $line 4];
	    set seclist [split $secline "\."]
	    set evtsec  [lindex $seclist 0]; set evtsec [format "%02d" $evtsec];
	    set val2    [scan $evtsec "%d"];
	    set evtmsec [expr int(($secline-$val2)*1000.0)]
	    set evtmsec [format "%03d" $evtmsec];
	    set jday2dateline [jday2datetcl $evtyear $evtjday]
	    set evtmon [lindex $jday2dateline 1]; set evtday [lindex $jday2dateline 2];
	    # set evtdatesac evttimesac
	    set evtdatesac "$evtyear/$evtmon/$evtday"
	    set evttimesac "$evthour:$evtmin:$evtsec.$evtmsec"
	}
    } elseif {$sw==1} {
	set eventfile "$outdir/$evtidcur.$channel.EVENT.RLOC"
	#if {[file exists $eventfile]==0} {tk_messageBox -type ok -message "$outdir/$evtidcur.$channel.EVENT.RLOC does not exist"}
	set line [exec grep $evtidcur $eventfile]; regsub -all {[[:blank:]]+} $line " " line
	set evtdatesd [lindex $line 1];
	set evttimesd [lindex $line 2];
	set line [converttime $evtdatesd $evttimesd]
	set evtyear [lindex $line 0]; set evtjday [lindex $line 1]; set evthour [lindex $line 2];
	set evtmin  [lindex $line 3]; set evtsec  [lindex $line 4]; set evtmsec [lindex $line 5];
	#set evtomarksd [exec t1t2 $year $jday $hour $min $sec $msec $evtyear $evtjday $evthour $evtmin $evtsec $evtmsec]
	set evtomarksd [t1t2tcl $year $jday $hour $min $sec $msec $evtyear $evtjday $evthour $evtmin $evtsec $evtmsec]
    } elseif {$sw==2} {
	set eventfile "$outdir/$evtidcur.$channel.EVENT.FINAL"
	set line [exec grep $evtidcur $eventfile]; regsub -all {[[:blank:]]+} $line " " line
	set evtdatefin [lindex $line 1];
	set evttimefin [lindex $line 2];
	if {$evtdatefin!="0000/00/00"} {
	    set line [converttime $evtdatefin $evttimefin]
	    set evtyear [lindex $line 0]; set evtjday [lindex $line 1]; set evthour [lindex $line 2];
	    set evtmin  [lindex $line 3]; set evtsec  [lindex $line 4]; set evtmsec [lindex $line 5];
	    #set evtomarkfin [exec t1t2 $year $jday $hour $min $sec $msec $evtyear $evtjday $evthour $evtmin $evtsec $evtmsec]
	    set evtomarkfin [t1t2tcl $year $jday $hour $min $sec $msec $evtyear $evtjday $evthour $evtmin $evtsec $evtmsec]
	} else {
	    set evtyear "0000"; set evtjday "000";  set evthour "00";
	    set evtmin  "00";   set evtsec "00";   set evtmsec "000";
	    set evtomarkfin -12345.0
	}
    }

}
#################################################
proc converttime {evtdate evttime} {
    set datelist [split $evtdate "\/"];
    set timelist [split $evttime "\:"];
    set evtyear [lindex $datelist 0]; set evtmon [lindex $datelist 1]; set evtday  [lindex $datelist 2];
    set evthour [lindex $timelist 0]; set evtmin [lindex $timelist 1]; set evtsecs [lindex $timelist 2];
    #set jdayline [exec date2jday $evtyear $evtmon $evtday];
    set jdayline [date2jdaytcl $evtyear $evtmon $evtday];
    set jdayline [string trim $jdayline];regsub -all {[[:blank:]]+} $jdayline " " jdayline
    set evtyear [lindex $jdayline 0]; set evtjday [lindex $jdayline 1]
    set tt [split $evtsecs "\."]
    set evtsec [lindex $tt 0];
    set val2 [scan $evtsec "%d"]; #If not scan, and if evtsec=08, there will be an error, look like invalid octal number
    set evtmsec [expr int(($evtsecs-$val2)*1000.0)];
    set evtmsec [format "%03d" $evtmsec]
    return "$evtyear $evtjday $evthour $evtmin $evtsec $evtmsec"
}
#################################################
proc synthtraveltime {} {
    global sfnum dtname sfgcarc
    global timelist
    global model
    global evtdep
    global readthreadnum calculatethreadnum

    #for {set i 1} {$i<=$sfnum} {incr i 1} {
	#set line [exec modeltime.pl $model all $evtdep $sfgcarc($dtname($i)) | tr "\n" " "]
	#set line [modeltimetcl $model all $evtdep $sfgcarc($dtname($i))]
	#set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	#set timelist($dtname($i)) [split $line " "]
    #}

    # set thread num & create thread, then set them to wait
    for {set i 0} {$i<$readthreadnum} {incr i 1} {
	set workid($i) [thread::create {
		package require Ttrace;
		global rayn raymark;   # This value is shared in each thread
		set rayn 0
		thread::wait;
		}]
	# initial distcmd in each thread
	set distcmd($i) ""
    }
    # obtain main thread id
    set mainid [thread::id]

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set imark [expr int(fmod($i,$readthreadnum))]
	thread::send -async $workid($imark) "recordinfo $dtname($i)"
	set distcmd($imark) "$distcmd($imark)$sfgcarc($dtname($i))\n"
	#thread::send -async $workid($imark) "modeltimetclthreadsub $model all $evtdep $sfgcarc($dtname($i)) $dtname($i)"
    }
    set precmd "$model\nall\n\n$evtdep\n";
    set exitcmd "-99\n-99\n";
    # calculate ttimes, and send stop command to each thread
    for {set i 0} {$i<$readthreadnum} {incr i 1} {
	thread::send -async $workid($i) "getttimes {$precmd} {$distcmd($i)} {$exitcmd}"
	thread::send -async $workid($i) "thread::release $workid($i)"
    }
    while {[llength [thread::names]] > 1} {after 200}

    # copy from synth(in thread) to timelist(in main thread)
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set timelist($dtname($i)) [tsv::get synth $dtname($i)]
    }
}
#################################################
ttrace::eval {
    proc modeltimetclthreadsub {model phase depth distindeg kstnm} {
	#puts $model,$all,$evdep,$gcarc,$kstnm
	set awkcmd "{if(NF==9 && \$2~/\[0-9\]/){print \$3,\$4,\$7}else if(NF==8 && \$1~/\[0-9\]/){print \$2,\$3,\$6}}"
	set lines [exec echo "$model\n$phase\n\n$depth\n$distindeg\n-99\n-99\n" | ttimes | awk $awkcmd | tr "\n" " "]
	set lines [string trim $lines];regsub -all {[[:blank:]]+} $lines " " lines;
	tsv::set synth $kstnm $lines;
    }

    proc recordinfo {kstnm} {
	global rayn raymark
	set rayn [expr $rayn+1]
	set raymark($rayn) $kstnm
    }
    proc getttimes {precmd distcmd exitcmd} {
	global rayn raymark;  # This value is shared in each thread

	set allcmd "$precmd$distcmd$exitcmd"
	set line [exec ttimes << $allcmd]
	set linenew [split $line "\n"]
	set linenum [llength $linenew]
	set marki 0
	for {set i 0} {$i<$linenum} {incr i 1} {
	    set line [lindex $linenew $i];
	    set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	    set fieldnum [llength $line]
	    if {$fieldnum==9 && [regexp {^([0-9]+)$} [lindex $line 1]]==1} {
		set marki [expr $marki+1]
		set timeout($marki) {}
		lappend timeout($marki) [lindex $line 2] [lindex $line 3] [lindex $line 6]
	    } elseif {$fieldnum==8 && [regexp {^([0-9]+)$} [lindex $line 0]]==1} {
		lappend timeout($marki) [lindex $line 1] [lindex $line 2] [lindex $line 5]
	    }	    
	}
	if {$marki!=$rayn} {puts "error from thread for calculating synth"}

	for {set i 1} {$i<=$rayn} {incr i 1} {
	    tsv::set synth $raymark($i) $timeout($i)
	}

    }
}
#################################################
proc getsacinfo {channel} {
    global bp1 bp2
    global datdir evtidcur
    global sacfilenum sacfilename
    global year jday hour min sec msec mag dt

    set sacfilelist [glob -nocomplain swap/$evtidcur.*.$channel.SAC]
    set sacfilenum [llength $sacfilelist]
    if {$sacfilenum==0} {
	tk_messageBox -type ok -message "No available data in swap. If you switch channel parameter, please click evtID first";
	return 1
    }
    for {set i 1} {$i<=$sacfilenum} {incr i 1} {
	set index [expr $i-1]
	set line [lindex $sacfilelist $index];
	set linelist [split $line "."];
	set sacfilename($i) [lindex $linelist end-2];
    }
    set sacfilesingle [lindex $sacfilelist 0]
    #set line [exec saclst kztime delta f $sacfilesingle]
    set line [saclsttcl kztime delta mag f $sacfilesingle]
    set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
    set linelist [split $line " "];
    set year [lindex $linelist 1];  set jday [lindex $linelist 2];  set hour [lindex $linelist 3];
    set min  [lindex $linelist 4];  set sec  [lindex $linelist 5];  set msec [lindex $linelist 6];
    set dt   [lindex $linelist 7];  set mag  [lindex $linelist 8];

    set year [scan $year "%d"];  set jday [scan $jday "%d"];
    set hour [scan $hour "%d"];  set min  [scan $min  "%d"];
    set sec  [scan $sec  "%d"];  set msec [scan $msec "%d"]; # if not scan again, [expr $sec+$msec/1000.0] will show an error message

    set dttext [format "%.3f" $dt];
    set fmax [expr 1.0/$dt/2.0];   set fmax [format "%.3f" $fmax];
    #set bp1f [expr 1.0/$bp1]; set bp2f [expr 1.0/$bp2];
    #puts $bp1f,$bp2f,$fmax;
    if {($bp1>$fmax || $bp2>$fmax) && $bp1!=$bp2} {
	set ctrlsw [tk_messageBox -type yesno -message "Warning:\nData sample rate is $dttext sec, Half frequency is $fmax Hz.\nFilter parameters are $bp1 Hz and $bp2 Hz.\nContinue?"]
	if {$ctrlsw=="no"} {return 1}
    }
    return 0
}
#################################################
proc renewstation {channel} {
    global evtidcur
    global sfnum sfname sflat sflon sfele sfpick sfpt sfphase
    global sacfilenum sacfilename
    global sfnameold

    # make a copy
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set sfnameold($i) $sfname($i);
    }
    # copying and renewing
    for {set i 1} {$i<=$sacfilenum} {incr i 1} {
	set sfname($i) $sacfilename($i);
	set swval [isexist $sfname($i)];
	if {$swval == "no"} {
	    set sacfilesingle "swap/$evtidcur.$sfname($i).$channel.SAC";
	    set line [saclsttcl kstnm stla stlo stel f $sacfilesingle];
	    set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	    set linelist [split $line " "];
	    if {$sfname($i)!=[lindex $linelist 1]} {
		tk_messageBox -type ok -message "The station name are not consistent for $sacfilesingle.\nPlease check!";
		return
	    }
	    set sflat($sfname($i))  [lindex $linelist 2];
	    set sflon($sfname($i))  [lindex $linelist 3];
	    set sfele($sfname($i))  [lindex $linelist 4];
	    set sfpick($sfname($i).$channel) 0
	    set sfmansw($sfname($i).$channel) 0
	    set sfpt($sfname($i).$channel)   -12345.0
	    set sfphase($sfname($i).$channel) "U"
	}
    }
    set sfnum $sacfilenum;

}
proc isexist {staname} {
    global sfnum
    global sfnameold

    # If found, return yes; not found return no;
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$staname == $sfnameold($i)} {
	    return yes;
	}
    }
    return no;
}
#################################################
proc loadstation {channel loadstasw} {
    global datdir outdir evtidcur
    global sfnum sfname sflat sflon sfele sfpick sfpt sfphase sfmansw
    global readthreadnum calculatethreadnum

  if {$loadstasw==1} {

    set pickfile "";
    set pickforloc "$outdir/$evtidcur.$channel.PICK.FORLOC"
    set pickfinal  "$outdir/$evtidcur.$channel.PICK.FINAL"
    if {[file exists $pickforloc]==1} { set pickfile $pickforloc}
    if {[file exists $pickfinal ]==1} { set pickfile $pickfinal }

    if {[file exists $pickfile]==1} {
	set fsort [open $pickfile]
	    set sfnum 0
	    while {[gets $fsort line]>=0} {
		set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
		if {[string range $line 0 0]==">"} {continue;}
		set linelist [split $line " "]
		set sfnum [expr $sfnum+1]
		set sfname($sfnum)          [lindex $linelist 0]
		set sflat($sfname($sfnum))  [lindex $linelist 1]
		set sflon($sfname($sfnum))  [lindex $linelist 2]
		set sfele($sfname($sfnum))  [lindex $linelist 3]
		set sfpick($sfname($sfnum).$channel) [lindex $linelist 4];
		set sfmansw($sfname($sfnum).$channel) $sfpick($sfname($sfnum).$channel)
		set sfpt($sfname($sfnum).$channel)   [lindex $linelist 5]
		set sfphase($sfname($sfnum).$channel) [lindex $linelist 6]
	    }
	close $fsort
    } else {
	set sacfilelist [glob -nocomplain swap/$evtidcur.*.$channel.SAC]

	# set thread num & create thread, then set them to wait
	for {set i 0} {$i<$readthreadnum} {incr i 1} {
	    set workid($i) [thread::create {package require Ttrace; thread::wait;}]
	}
	# obtain main thread id
	set mainid [thread::id]

	set sfnum 0
	foreach element $sacfilelist {
	    set sfnum [expr $sfnum+1]
	    set imark [expr int(fmod($sfnum,$readthreadnum))]
	    thread::send -async $workid($imark) "loadstationsub $element $sfnum"
	}

	# send stop command to each thread
	for {set i 0} {$i<$readthreadnum} {incr i 1} {
	    thread::send -async $workid($i) "thread::release $workid($i)"
	}
	while {[llength [thread::names]] > 1} {after 200}


	# copy dat those in thread to those in main thread
	for {set i 1} {$i<=$sfnum} {incr i 1} {
	    set kstnm [tsv::get sfname $i]
	    set sfname($i) $kstnm
	    set sflat($kstnm) [tsv::get sflat $kstnm]
	    set sflon($kstnm) [tsv::get sflon $kstnm]
	    set sfele($kstnm) [tsv::get sfele $kstnm]
	    set sfpick($kstnm.$channel) 0
	    set sfmansw($kstnm.$channel) 0
	    set sfpt($kstnm.$channel) -12345.0
	    set sfphase($kstnm.$channel) "U"
	}

	#continue
	#set sfnum 0
	#foreach element $sacfilelist {
	#    set sfnum [expr $sfnum+1]
	#    set line [saclsttcl kstnm stla stlo stel f $element]
	#    set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	#    set linelist [split $line " "]
	#    set sfname($sfnum)           [lindex $linelist 1]
	#    set sflat($sfname($sfnum))   [lindex $linelist 2]
	#    set sflon($sfname($sfnum))   [lindex $linelist 3]
	#    set sfele($sfname($sfnum))   [lindex $linelist 4]
	#    set sfpick($sfname($sfnum).$channel)  0
	#    set sfpt($sfname($sfnum).$channel)    -12345.0
	#}


    }

  }
    #for {set i 1} {$i<=$sfnum} {incr i 1} {
	#set id [.stafigure find withtag station.$sfname($i)]
	#if {$id != ""} {
	    #if {$sfpick($sfname($i).$channel)==0} {.stafigure itemconfigure station.$sfname($i) -width 1 -fill green}
	    #if {$sfpick($sfname($i).$channel)==1} {.stafigure itemconfigure station.$sfname($i) -width 1 -fill #a00}
	#}
    #}
    colorstation $channel
}
ttrace::eval {
    proc loadstationsub {sacfile sfnum} {
	set line [saclsttcl kstnm stla stlo stel f $sacfile]
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	set linelist [split $line " "]
	set kstnm [lindex $linelist 1]
	tsv::set sfname $sfnum $kstnm
	tsv::set sflat $kstnm [lindex $linelist 2]
	tsv::set sflon $kstnm [lindex $linelist 3]
	tsv::set sfele $kstnm [lindex $linelist 4]
    }

}
#################################################
proc colorstation {channel} {
    global sfnum sfname sflat sflon sfele sfpick sfpt
    global colorband colorn

    set tmin  1.e+20;
    set tmax -1.e+20;
    for  {set i 1} {$i<=$sfnum} {incr i 1} {
	if {$sfpick($sfname($i).$channel)==1} {
	    if {$sfpt($sfname($i).$channel)>=$tmax} {set tmax $sfpt($sfname($i).$channel)}
	    if {$sfpt($sfname($i).$channel)<=$tmin} {set tmin $sfpt($sfname($i).$channel)}
	}
    }
    set tdt [expr ($tmax-$tmin)/($colorn+1)];

    for {set i 1} {$i<=$sfnum} {incr i 1} {
	set id [.stafigure find withtag station.$sfname($i)]
	if {$id != ""} {
	    if {$sfpick($sfname($i).$channel)==0} {.stafigure itemconfigure station.$sfname($i) -fill green}
	    if {$sfpick($sfname($i).$channel)==1} {
		if {$tdt!=0.0} {
		    set nn [expr int(($sfpt($sfname($i).$channel)-$tmin)/$tdt)]
		} else {
		    set nn $colorn
		}
		if {$nn>$colorn} {set nn $colorn}
		if {$nn<0}       {set nn 0}
		.stafigure itemconfigure station.$sfname($i) -fill $colorband($nn)
	    }
	} else {
	    tk_messageBox -type ok -message "Cannot find station $sfname($i) in station figure, Maybe this station is not included in sta.loc.txt"
	}
    }


}
#################################################
proc relocate {channel} {
    global x1 x2 y1 y2
    global z1 z2
    global dxx dyy dzz
    global datdir outdir evtidcur
    global evtlon evtlat evtdep evtdate evttime evtomark
    global evtomarkvalsw
    global sacdatwinb sacdatwine
    global year jday hour min sec msec
    global pickerrmax stalim
    global model
    global evtdatesd evttimesd evtomarksd evtlonsd evtlatsd evtdepsd evtlonrms evtlatrms evtdeprms
    global evtdepmax evtdephand
    global isrelocate
    global evtherrsd evtverrsd evtherrrms evtverrrms
    global pickcurlimt1 pickcurlimt2 pickcurlimd1 pickcurlimd2

    if {$channel=="BHZ"} {set canvasname .waveallz; set infonew1 .zinfonew1; set infonew2 .zinfonew2;}
    if {$channel=="BHE"} {set canvasname .wavealle; set infonew1 .einfonew1; set infonew2 .einfonew2;}
    if {$channel=="BHN"} {set canvasname .wavealln; set infonew1 .ninfonew1; set infonew2 .ninfonew2;}
    $infonew1 configure -text ""; #$infonew2 configure -text ""

    if {[$canvasname find withtag wave]==""} {return}
    set x1text [format "%7.3f" $x1];    set x2text [format "%7.3f" $x2];
    set y1text [format "%7.3f" $y1];    set y2text [format "%7.3f" $y2];

    .information configure -text "evtlocate $x1text $dxx $x2text $y1text $dyy $y2text $z1 $dzz $z2 PICKfile 0 1 $pickerrmax $stalim 0 0 $model 1"
    set pickfile "$outdir/$evtidcur.$channel.PICK.FORLOC"
    if {[file exists $pickfile]==0} {tk_messageBox -type ok -message "Please save pick time firstly";return;}
    set picknum [exec awk "{sum += \$5}; END {print sum}" $pickfile];
    if {$picknum<3} {tk_messageBox -type ok -message "The amount of picked time is $picknum, less than limited value (3)";return;}
    if {$picknum<$stalim} {tk_messageBox -type ok -message "The amount of picked time is $picknum, less than user-limited value ($stalim).\nPlease decrease StaLimit or increase picks";return;}

    set evtlocateline [exec evtlocate $x1 $dxx $x2 $y1 $dyy $y2 $z1 $dzz $z2 $pickfile 0 1 $pickerrmax $stalim 0 0 $model 1]
    set fout [open "$outdir/$evtidcur.$channel.EVENT.RLOC" w]
	puts $fout $evtlocateline
    close $fout
    set evtdatesd [lindex $evtlocateline 1]
    set evttimesd [lindex $evtlocateline 2]
    ###### set evtdate evttime ###### 
    set evtdate $evtdatesd
    set evttime $evttimesd

    set line [converttime $evtdate $evttime]
    set evtyear [lindex $line 0]; set evtjday [lindex $line 1]; set evthour [lindex $line 2];
    set evtmin  [lindex $line 3]; set evtsec  [lindex $line 4]; set evtmsec [lindex $line 5];
    #set evtomarksd [exec t1t2 $year $jday $hour $min $sec $msec $evtyear $evtjday $evthour $evtmin $evtsec $evtmsec]
    set evtomarksd [t1t2tcl $year $jday $hour $min $sec $msec $evtyear $evtjday $evthour $evtmin $evtsec $evtmsec]
    ###### set evtomark ###### 
    set evtomark $evtomarksd
    set evtomarkvalsw 1;  # in case: the event position is manually marked, evtomarkvalsw=0.

    set isrelocate 1
    set evtlonsd [lindex $evtlocateline 3]
    set evtlatsd [lindex $evtlocateline 4]
    set evtdepsd [lindex $evtlocateline 5]
    set evtherrsd [lindex $evtlocateline 13]
    set evtverrsd [lindex $evtlocateline 14]
    .stafigure delete evt.2
    markevtpos $evtlonsd $evtlatsd 2

    set evtlonrms [lindex $evtlocateline 6]
    set evtlatrms [lindex $evtlocateline 7]
    set evtdeprms [lindex $evtlocateline 8]
    set evtherrrms [lindex $evtlocateline 15]
    set evtverrrms [lindex $evtlocateline 16]
    .stafigure delete evt.3
    markevtpos $evtlonrms $evtlatrms 3

    ###### set evtlon evtlat evtdep ###### 
    set evtlon $evtlonsd;
    set evtlat $evtlatsd;
    set evtdep $evtdepsd
    drawgcarc $evtlon $evtlat blue 2
    if {$evtdep>$evtdepmax} {
	set evtdepmax $evtdep
	.textevdp insert 1.0 $evtdepmax
	.scaleevtdep configure -to $evtdepmax
    }
    set evtdephand $evtdep

    set evtlonsd  [format "%9.4f" $evtlonsd ]; set evtlatsd  [format "%8.4f" $evtlatsd ]; set evtdepsd  [format "%5.1f" $evtdepsd ];
    set evtlonrms [format "%9.4f" $evtlonrms]; set evtlatrms [format "%8.4f" $evtlatrms]; set evtdeprms [format "%5.1f" $evtdeprms];
    $infonew1 configure -text "$evtdate $evttime $evtlonsd $evtlatsd $evtdepsd"
    #$infonew2 configure -text "$evtdate $evttime $evtlonrms $evtlatrms $evtdeprms"

    # reorder the waveform et al.
    set swinwaveallz [.waveallz find withtag wave]
    set swinwavealle [.wavealle find withtag wave]
    set swinwavealln [.wavealln find withtag wave]
    set swinwavezoom [.wavezoom find withtag wave]

    if {$swinwaveallz!="" || $swinwavealle!="" || $swinwavealln!=""} {
	caldistindeg $evtlon $evtlat
	sortbydist
	synthtraveltime	
    }

    if {$swinwaveallz!=""} {
	waveallini BHZ
	reorderwaveinall BHZ
	plotsynthtmarkinwaveall BHZ
    }
    if {$swinwavealle!=""} {
	waveallini BHE
	reorderwaveinall BHE
	plotsynthtmarkinwaveall BHE
    }
    if {$swinwavealln!=""} {
	waveallini BHN
	reorderwaveinall BHN
	plotsynthtmarkinwaveall BHN
    }
    if {$swinwavezoom!=""} {
	reorderwaveinzoom
	plotsynthtmarkinwavezoom .wavezoom
    }

    fileindicator $channel
    plotcurdots $channel

    getpickcurrentevtdatrange $channel
    zoomcurvefigureini $channel $pickcurlimt1 $pickcurlimt2 $pickcurlimd1 $pickcurlimd2
    plotzoomcurve $channel
}
#################################################
proc outpicktime {channel} {
    global datdir outdir evtidcur
    global sfnum sfname sflat sflon sfele sfpick sfpt sfphase sfsnr
    global sfdist
    global year jday hour min sec msec dt
    global evtomark
    #global chanidcur

    if {$channel=="BHZ"} {set canvasname .waveallz}
    if {$channel=="BHE"} {set canvasname .wavealle}
    if {$channel=="BHN"} {set canvasname .wavealln}

    if {[$canvasname find withtag wave]==""} {return}
    #set pickmin 1000000.0
    #for {set i 1} {$i<=$sfnum} {incr i 1} {
	#if {$sfpt($sfname($i).$channel)<=$pickmin} {set pickmin $sfpt($sfname($i).$channel)}
    #}
    #set pickmin [expr $pickmin-1.0]; # unused at present

    # Calculate SNR
    signalnoiseratio $channel
    set fout [open "$outdir/$evtidcur.$channel.PICK.FORLOC" w]
    for {set i 1} {$i<=$sfnum} {incr i 1} {
	puts $fout [format "%8s %10.4f %10.4f %7.4f %1d %10.3f %1s %10.3f %04d %03d %02d %02d %02d %03d" \
                            $sfname($i) $sflat($sfname($i)) $sflon($sfname($i)) $sfele($sfname($i)) \
                            $sfpick($sfname($i).$channel) $sfpt($sfname($i).$channel) $sfphase($sfname($i).$channel) $sfsnr($sfname($i).$channel) \
                            $year $jday $hour $min $sec $msec]
    }
    puts $fout "> $evtidcur"
    close $fout

    fileindicator $channel
}
#################################################
proc findminvalline {fastfile index} {
    set index [expr $index-1]
    set fsort [open $fastfile]
	set minval 10000000.00
	while {[gets $fsort line]>=0} {
	    set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	    set lines [split $line " "]
	    set temval [lindex $lines $index]
	    #set dep [lindex $lines 2]
	    if {$temval<=$minval} {
		set linesel $line
		set minval $temval
	    }
	}
    close $fsort
    return $linesel
}
#################################################
ttrace::eval {
proc distaztcl {stalat stalon evtlat evtlon} {


    if {$stalat==$evtlat && $stalon==$evtlon} {
	set delta 0.0; set baz 0.0; set az 0.0
	set out {}
	lappend out $delta; lappend out $baz; lappend out $az
	return $out
    }
    set pi [expr acos(-1.0)]

    set piby2 [expr $pi/2.0]
    set rad [expr $pi/180.0]
    set sph [expr 1.0/298.257]
    set scolat [expr $piby2-atan((1.0-$sph)*(1.0-$sph)*tan($stalat*$rad))]
    set ecolat [expr $piby2-atan((1.0-$sph)*(1.0-$sph)*tan($evtlat*$rad))]
    set slon [expr $stalon*$rad]
    set elon [expr $evtlon*$rad]

    set a [expr sin($scolat)*cos($slon)]
    set b [expr sin($scolat)*sin($slon)]
    set c [expr cos($scolat)]
    set d [expr sin($slon)]
    set e [expr -cos($slon)]
    set g [expr -$c*$e]
    set h [expr $c*$d]
    set k [expr -sin($scolat)]

    set aa [expr sin($ecolat)*cos($elon)]
    set bb [expr sin($ecolat)*sin($elon)]
    set cc [expr cos($ecolat)]
    set dd [expr sin($elon)]
    set ee [expr -cos($elon)]
    set gg [expr -$cc*$ee]
    set hh [expr $cc*$dd]
    set kk [expr -sin($ecolat)]

    set val [expr $a*$aa+$b*$bb+$c*$cc]
    if {$val<-1.0} {set val -1.0}
    if {$val> 1.0} {set val  1.0}
    set delta [expr acos($val)]
    set delta [expr $delta/$rad]

    set rhs1 [expr ($aa-$d)*($aa-$d)+($bb-$e)*($bb-$e)+$cc*$cc-2.0]
    set rhs2 [expr ($aa-$g)*($aa-$g)+($bb-$h)*($bb-$h)+($cc-$k)*($cc-$k)-2.0]
    set baz  [expr atan2($rhs1,$rhs2)]
    if {$baz<0.0} {set baz [expr $baz+2.0*$pi]}
    set baz [expr $baz/$rad]


    set rhs1 [expr ($a-$dd)*($a-$dd)+($b-$ee)*($b-$ee)+$c*$c-2.0]
    set rhs2 [expr ($a-$gg)*($a-$gg)+($b-$hh)*($b-$hh)+($c-$kk)*($c-$kk)-2.0]
    set az  [expr atan2($rhs1,$rhs2)]
    if {$az<0.0} {set az [expr $az+2.0*$pi]}
    set az [expr $az/$rad]

    set out {}
    lappend out $delta; lappend out $baz; lappend out $az
    return $out
}
}
#################################################

#################################################
proc gcarcdat {lon lat dist} {
    global sc swidth sheight bd
    global xscale yscale x1 x2 y1 y2
    set pi [expr acos(-1.0)]
    set azn [expr $pi*2.0/360.0]
    set lonarc [expr $lon*$pi/180.0]
    set latarc [expr $lat*$pi/180.0]
    set delta  [expr $dist*$pi/180.0]
    set aa {}
    for {set i 1} {$i<=361} {incr i 10} {
	set az [expr $azn*($i-1)]
	set latout [expr asin(sin($latarc)*cos($delta)+cos($latarc)*sin($delta)*cos($az))]
	set tem [expr (cos($delta)-sin($latarc)*sin($latout))/(cos($latarc)*cos($latout))]
	if {$tem<-1.0} {set tem -1.0}
	if {$tem> 1.0} {set tem  1.0}
	set ph     [expr acos($tem)]
	if {$az>=0.0 && $az<=$pi} {
	    set ph [expr abs($ph)]
	} else {
	    set ph [expr -1.0*abs($ph)]
	}
	set lonout [expr $lonarc+$ph]
	set latout [expr $latout*180.0/$pi]
	set lonout [expr $lonout*180.0/$pi]
	if {$lonout> 180.0} {set lonout [expr $lonout-360.0]}
	if {$lonout<-180.0} {set lonout [expr $lonout+360.0]}
	set line [stafigurexy2poixy $lonout $latout $x1 $y1 $xscale $yscale $bd $sheight]
	set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line;
	set poix [lindex [split $line " "] 0];
	set poiy [lindex [split $line " "] 1];
	lappend aa $poix $poiy
    }
    return $aa
}

#################################################
proc zoom {poix poiy scale} {
    global sc swidth sheight bd
    global xscale yscale x1 x2 y1 y2
    
    set line [stafigurepoixy2xy $poix $poiy $x1 $y2 $xscale $yscale $bd];
    set line [string trim $line];regsub -all {[[:blank:]]+} $line " " line
    set lon [lindex [split $line " "] 0];
    set lat [lindex [split $line " "] 1];

    # Avoid the mouse cursor is out of the lon-lat area
    if {$lon<-180.0 || $lon>180.0 || $lat<-90.0 || $lat>90.0} {return;}

    set x1 [expr $lon-($lon-$x1)*$scale];if {$x1<-180.0} {set x1 -180.0;}
    set x2 [expr $lon+($x2-$lon)*$scale];if {$x2> 180.0} {set x2  180.0;}
    set y1 [expr $lat-($lat-$y1)*$scale];if {$y1< -90.0} {set y1  -90.0;}
    set y2 [expr $lat+($y2-$lat)*$scale];if {$y2>  90.0} {set y2   90.0;}
}
#################################################
proc stafigurepoixy2xy {curx cury x1 y2 xscale yscale bd} {
    set lon [expr $x1+($curx-$bd)/$xscale]; #set lon [format "%7.3f" $lon];
    set lat [expr $y2-($cury-$bd)/$yscale]; #set lat [format "%7.3f" $lat];
    return "$lon $lat"
}

proc stafigurexy2poixy {lon lat x1 y1 xscale yscale bd sheight} {
    set curx [expr ($lon-$x1)*$xscale+$bd];
    set cury [expr $sheight-($lat-$y1)*$yscale+$bd];
    return "$curx $cury"
}
#################################################

#################################################
proc uniqevtid a {
    set evtidlist "";
    foreach element $a {
	set evtid [lindex [split $element "/"] end]
	set evtid [lindex [split $evtid "."] 1]
	set is 0
	foreach tem $evtidlist {
	    if {[string compare $tem $evtid] == 0} {
		set is 1;break;
	    }
	}
	if {$is == 0} {lappend evtidlist $evtid;}
    }
    set evtidlist [lsort $evtidlist]
    return $evtidlist
}
#################################################
proc jday2datetcl {year jday} {
    # If not scan, there will be a small error, then m-array will have a wrong index
    set year [scan $year "%d"]
    set jday [scan $jday "%d"]

    set m(1) 31; set m(2) 28; set m(3) 31;  set m(4) 30;  set m(5) 31;  set m(6) 30;
    set m(7) 31; set m(8) 31; set m(9) 30; set m(10) 31; set m(11) 30; set m(12) 31;
    if {[expr fmod($year,4)]==0.0} {set m(2) 29}
    if {[expr fmod($year,4)]==0.0} {
	if {$jday<1 || $jday>366} {tk_messageBox -type ok -message "Jday error:Year:$year, Jday:$jday"}
    } else {
	if {$jday<1 || $jday>365} {tk_messageBox -type ok -message "Jday error:Year:$year, Jday:$jday"}
    }
    set total 0
    for {set i 1} {$i<=12} {incr i 1} {
	set total [expr $total+$m($i)]
	if {$jday<=$total} {
	    set mon $i
	    set day [expr $m($i)-($total-$jday)]
	    break;
	}
    }
    set year [format "%04d" $year]
    set mon  [format "%02d" $mon]
    set day  [format "%02d" $day]
    return "$year $mon $day"
}
#################################################
proc date2jdaytcl {year mon day} {
    # If not scan, there will be a small error, then m-array will have a wrong index
    set year [scan $year "%d"]
    set mon  [scan $mon  "%d"]
    set day  [scan $day  "%d"]

    set m(1) 31; set m(2) 28; set m(3) 31;  set m(4) 30;  set m(5) 31;  set m(6) 30;
    set m(7) 31; set m(8) 31; set m(9) 30; set m(10) 31; set m(11) 30; set m(12) 31;
    if {[expr fmod($year,4)]==0.0} {set m(2) 29}
    if {$mon<1 || $mon>12} {
	tk_messageBox -type ok -message "Month error"
    } else {
	if {$day<1 || $day>$m($mon)} {
	    tk_messageBox -type ok -message "error: Month $mon has $day Days"
	}
    }
    set jday 0
    for {set i 1} {$i<=[expr $mon-1]} {incr i 1} {
	set jday [expr $jday+$m($i)]
    }
    set jday [expr $jday+$day]
    set jday [format "%03d" $jday]
    set year [format "%04d" $year]
    return "$year $jday"
}
#################################################
proc t1t2tcl {year1 jday1 hour1 minu1 sec1 msec1 year2 jday2 hour2 minu2 sec2 msec2} {
    # If not scan, there will be a small error, then m-array will have a wrong index
    set year1 [scan $year1 "%d"]; set jday1 [scan $jday1 "%d"]; set hour1 [scan $hour1 "%d"];
    set minu1 [scan $minu1 "%d"]; set sec1  [scan  $sec1 "%d"]; set msec1 [scan $msec1 "%d"];
    set year2 [scan $year2 "%d"]; set jday2 [scan $jday2 "%d"]; set hour2 [scan $hour2 "%d"];
    set minu2 [scan $minu2 "%d"]; set sec2  [scan  $sec2 "%d"]; set msec2 [scan $msec2 "%d"];

    set sw 0
    if {$year1>$year2} {
	set sw 1;
    } elseif {$year1==$year2} {
	if {$jday1>$jday2} {
	    set sw 1;
	} elseif {$jday1==$jday2} {
	    if {$hour1>$hour2} {
		set sw 1;
	    } elseif {$hour1==$hour2} {
		if {$minu1>$minu2} {
		    set sw 1;
		} elseif {$minu1==$minu2} {
		    if {$sec1>$sec2} {
			set sw 1;
		    } elseif {$sec1==$sec2} {
			if {$msec1>$msec2} {
			    set sw 1;
			}
		    }
		}
	    }
	}
    }
    if {$sw==1} {
	set tem $year1; set year1 $year2; set year2 $tem;
	set tem $jday1; set jday1 $jday2; set jday2 $tem;
	set tem $hour1; set hour1 $hour2; set hour2 $tem;
	set tem $minu1; set minu1 $minu2; set minu2 $tem;
	set tem  $sec1;  set sec1  $sec2;  set sec2 $tem;
	set tem $msec1; set msec1 $msec2; set msec2 $tem;
    }

    set t_day [tytcl $year1 $year2]
    set t_day [expr $t_day-$jday1+$jday2]
    set tem1 [expr $hour1*3600.0+$minu1*60.0+$sec1+$msec1/1000.0];
    set tem2 [expr $hour2*3600.0+$minu2*60.0+$sec2+$msec2/1000.0];
    set t_sec [expr -1.0*$tem1+$tem2]
    set t_sec [expr $t_day*24.0*3600.0+$t_sec]
    if {$sw==1} {set t_sec [expr -1.0*$t_sec]}
    set t_sec [format "%.2f" $t_sec]
    return "$t_sec";
}
#################################################
proc tytcl {year1 year2} {
    set year1 [scan $year1 "%d"];
    set year2 [scan $year2 "%d"];

    set sw 0
    if {$year1>$year2} {set tem $year1; set year1 $year2; set year2 $tem; set sw 1}
    set sum 0
    for {set i $year1} {$i<$year2} {incr i 1} {
	if {[expr fmod($i,4)]==0.0} {
	    set sum [expr $sum+366]
	} else {
	    set sum [expr $sum+365]
	}
    }
    if {$sw==1} {set sum [expr -1.0*$sum]}
    return "$sum"
}
#################################################
proc taddtcl {year jday hour min secs secadd} {

    set year [scan $year "%d"]; set jday [scan $jday "%d"];
    set hour [scan $hour "%d"]; set min  [scan $min  "%d"];
    set secs   [scan $secs   "%f"];
    set secadd [scan $secadd "%f"];

    set secs [expr $secs+$secadd]
    if {$secs<0.0} {
	while {$secs<0.0} {
	    set secs [expr $secs+60.0]; set min [expr $min-1];
	    if {$min <0} {set min  [expr  $min+60]; set hour [expr $hour-1];}
	    if {$hour<0} {set hour [expr $hour+24]; set jday [expr $jday-1];}
	    if {$jday<1} {
		set year [expr $year-1]
		if {fmod($year,4.0)==0.0} {
		    set jday 366;
		} else {
		    set jday 365;
		}
	    }
	}
    } elseif {$secs>=60.0} {
	while {$secs>=60.0} {
	    set secs [expr $secs-60.0]; set min [expr $min+1];
	    if {$min >59} {set min  [expr  $min-60]; set hour [expr $hour+1];}
	    if {$hour>23} {set hour [expr $hour-24]; set jday [expr $jday+1];}
	    if {fmod($year,4.0)==0.0} {set tem 366;} else {set tem 365;}
	    if {$jday>$tem} {set jday [expr $jday-$tem]; set year [expr $year+1];}
	}
    }

    set year [format "%04d" $year]
    set jday [format "%03d" $jday]
    set hour [format "%02d" $hour]
    set min  [format "%02d" $min]
    set secs [format "%.3f" $secs]

    return "$year $jday $hour $min $secs";

}

#################################################
proc modeltimetcl {model phase depth distindeg} {

    set awkcmd "{if(NF==9 && \$2~/\[0-9\]/){print \$3,\$4,\$7}else if(NF==8 && \$1~/\[0-9\]/){print \$2,\$3,\$6}}"
    set lines [exec echo "$model\n$phase\n\n$depth\n$distindeg\n-99\n-99\n" | ttimes | awk $awkcmd | tr "\n" " "]
    return $lines
}
#################################################
proc bindunit {} {
    global stwin ltwin stltsw stltlim
    global kswin kssw klim slim
    global aicwin aicsw aiclim
    global bp1 bp2 evtinter recstalim
    global shift1 shift2
    global sc swidth
    global z1 z2
    global dxx dyy dzz
    global evtdepmax
    global pickerrmax stalim
    global vibnum vibdist
    global readthreadnum calculatethreadnum
    global singlewaveheight

    bind .tstwin     <3> {set stwin [.tstwin get 1.0 end];         set stwin [string trim $stwin];         .information configure -text "set to $stwin";}
    bind .tltwin     <3> {set ltwin [.tltwin get 1.0 end];         set ltwin [string trim $ltwin];         .information configure -text "set to $ltwin";}
    bind .tstltsw    <3> {set stltsw [.tstltsw get 1.0 end];       set stltsw [string trim $stltsw];       .information configure -text "set to $stltsw";}
    bind .tstltlim   <3> {set stltlim [.tstltlim get 1.0 end];     set stltlim [string trim $stltlim];     .information configure -text "set to $stltlim";}
    bind .tkswin     <3> {set kswin [.tkswin get 1.0 end];         set kswin [string trim $kswin];         .information configure -text "set to $kswin";}
    bind .tkssw      <3> {set kssw [.tkssw get 1.0 end];           set kssw [string trim $kssw];           .information configure -text "set to $kssw";}
    bind .tklim      <3> {set klim [.tklim get 1.0 end];           set klim [string trim $klim];           .information configure -text "set to $klim";}
    bind .tslim      <3> {set slim [.tslim get 1.0 end];           set slim [string trim $slim];           .information configure -text "set to $slim";}
    bind .taicwin    <3> {set aicwin [.taicwin get 1.0 end];       set aicwin [string trim $aicwin];       .information configure -text "set to $aicwin";}
    bind .taicsw     <3> {set aicsw [.taicsw get 1.0 end];         set aicsw [string trim $aicsw];         .information configure -text "set to $aicsw";}
    bind .taiclim    <3> {set aiclim [.taiclim get 1.0 end];       set aiclim [string trim $aiclim];       .information configure -text "set to $aiclim";}
    bind .tevtinter  <3> {set evtinter [.tevtinter get 1.0 end];   set evtinter [string trim $evtinter];   .information configure -text "set to $evtinter";}
    bind .trecstalim <3> {set recstalim [.trecstalim get 1.0 end]; set recstalim [string trim $recstalim]; .information configure -text "set to $recstalim";}
    bind .tvibnum    <3> {set vibnum [.tvibnum get 1.0 end];       set vibnum [string trim $vibnum];       .information configure -text "set to $vibnum";}
    bind .tvibdist   <3> {set vibdist [.tvibdist get 1.0 end];     set vibdist [string trim $vibdist];     .information configure -text "set to $vibdist";}

    bind .tbp1       <3> {set bp1 [.tbp1 get 1.0 end];             set bp1 [string trim $bp1];             .information configure -text "set to $bp1";}
    bind .tbp2       <3> {set bp2 [.tbp2 get 1.0 end];             set bp2 [string trim $bp2];             .information configure -text "set to $bp2";}
    bind .tshift1    <3> {set shift1 [.tshift1 get 1.0 end];       set shift1 [string trim $shift1];       .information configure -text "set to $shift1";}
    bind .tshift2    <3> {set shift2 [.tshift2 get 1.0 end];       set shift2 [string trim $shift2];       .information configure -text "set to $shift2";}

    bind .tsc        <3> {set sc [.tsc get 1.0 end];               set sc [string trim $sc];               .information configure -text "set to $sc";}
    bind .tsw        <3> {set swidth [.tsw get 1.0 end];           set swidth [string trim $swidth];       .information configure -text "set to $swidth";}

    bind .tdxx       <3> {set dxx [.tdxx get 1.0 end];             set dxx [string trim $dxx];             .information configure -text "set to $dxx"}
    bind .tdyy       <3> {set dyy [.tdyy get 1.0 end];             set dyy [string trim $dyy];             .information configure -text "set to $dyy"}
    bind .tdzz       <3> {set dzz [.tdzz get 1.0 end];             set dzz [string trim $dzz];             .information configure -text "set to $dzz"}
    bind .tsz1       <3> {set z1 [.tsz1 get 1.0 end];              set z1 [string trim $z1];               .information configure -text "set to $z1"}
    bind .tsz2       <3> {set z2 [.tsz2 get 1.0 end];              set z2 [string trim $z2];               .information configure -text "set to $z2"}
    bind .textevdp   <3> {set evtdepmax [.textevdp get 1.0 end];   set evtdepmax [string trim $evtdepmax]; .information configure -text "set to $evtdepmax"; .scaleevtdep configure -to $evtdepmax;}
    bind .tPKerr     <3> {set pickerrmax [.tPKerr get 1.0 end];    set pickerrmax [string trim $pickerrmax]; .information configure -text "set to $pickerrmax"}
    bind .tstalim    <3> {set stalim [.tstalim get 1.0 end];       set stalim [string trim $stalim];       .information configure -text "set to $stalim"}
    bind .twaveH     <3> {
         set singlewaveheight   [.twaveH get 1.0 end];      set singlewaveheight   [string trim $singlewaveheight];   .information configure -text "set to $singlewaveheight"; updatewavezoom}
    bind .tthreadread <3> {
         set readthreadnum      [.tthreadread get 1.0 end]; set readthreadnum      [string trim $readthreadnum];      .information configure -text "set to $readthreadnum"}
    bind .tthreadcalc <3> {
         set calculatethreadnum [.tthreadcalc get 1.0 end]; set calculatethreadnum [string trim $calculatethreadnum]; .information configure -text "set to $calculatethreadnum"}
}
main $argc $argv
