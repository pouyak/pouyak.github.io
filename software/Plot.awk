#!/usr/bin/awk -f
# Written by Pouya Kheradpour (pouyak@gmail.com)
# Copyright (c) 2005-2015 Massachusetts Institute of Technology
# Released under the MIT license
# https://pouya.kheradpour.com/software/LICENSE.mit

# version 20140120

# INPUT FORMAT:

# X series Y
#      X      Val|Min:Max
#      series name|color:stroke|marker|line|has_legend (color, stroke, marker, line are 1 indexed)
#
#      color/stroke: either name, index, or (r,g,b) in 0->255
#      marker: index (0 off); (c)ircle; (d)iamond; (t)riangle; (u)p triangle; (s)quare; (b)ar [default: c (stat); b (bar)]
#      line: step prefix (s for over and up, S for up and over) index (0 off, 1 solid, dash types 2+); width:dash-array [default: 1 (stat); 0 (bar)]
#
#   scatter plots (Mode=scat)
#      Y      Val1|Err1:Err2|MarkerSize|Text|TextOpts
#      TextOpts not yet implemented
#
#   bar plots (Mode=(s)bar(h))
#      Mode=sbar is stacked (bars overlap each other, not add)
#      Mode=barh is horizontal
#      Y      Val1|Err1:Err2|Val2|Err1:Err2|Text|TextOpts
#      TextOpts
#          (r,g,b) color [default: black]
#          size (floating point) 
#          L/B/Z/C/M/X/A/R (plot left, before bar, zero, plot center, bar center, bar end, after bar, right)
#          S/F (text rotated/not rotated)
#          P for prettynum
#
# %%Newline%% is converted to a newline
# %%Semicolon%% is converted to a semicolon
# %%VBar%% is converted to a |
# %%Tag%%text%%/Tag%% with Super, Sub, Monospace, Underline, Bold, Italic, Color=, Size=
# %%(|X|Y)(Min|Max)%% is replaced with extremes, for Vals or Errs
#
# Val2/Err2 only for bar plots --> Val2 defaults to BarZero

########## TODO
# support for banding/tick lines to go across 
# improve prettybound (esp for log)

########## TEST CODE
# scat plot test:
# (echo 0.2 'a|4|0|3' "0.25|0.2:1"; echo 0.3 'a|4|0|3' 0.9; echo 0.5 b "0.7|0.5:1.02"; echo 0.9 "b||c" 0.2; echo 0.9 "b||d" 0.25;echo 0.9 "b||t" 0.3;echo 0.9 "b||u" 0.35;echo 0.9 "b||s" 0.4;echo "1.04|0.8" b 0.5; echo 0.4 "Acsdfsdf" "0.4||20") | Plot.awk -vXMin=0 -vXMax=1 -vYMin=0 -vYMax=1 -vXTitle="X axis" -vYTitle="Y axis" -vTitle="Main Title" -vLegendWidthDelta=-22 > test.svg

# barh plot test:
# (echo x 3 "0.636364|0.569189:0.703539";echo x 2 "0.371859|0.304539:0.439179";echo x 1 "0.277778|0.215231:0.340325";echo x 0 "0.286432|0.223459:0.349405";echo y 3 "1.09232|0.931304:1.25334";echo y 2 "0.727833|0.586148:0.869518";echo y 1 "0.681211|0.554942:0.80748";echo y 0 "0.721716|0.579578:0.863854";) | Plot.awk -vMode=barh -vMin=0 -vMax=1.0 -vBarZero=0.5 -vXTitle="X axis" -vYTitle="Y axis" -vTitle="Main Title" > test.svg

# box whisker test (5,25,50,75,95%=0.35,0.4,0.6,0.7,0.8; outliers=0.1,0.15):
# (echo x "|red|||0|1" "0.4|0.35|0.7|0.8"; echo x "|red|||0|2" "0.6||0.6"; echo x "|red|c||0" 0.1; echo x "|red|c||0" 0.15) | Plot.awk -vMode=sbar -vMin=0 -vMax=1.0 > test.svg

# roc curve
# cat input.txt | pnstats -o xanNor | awk -F"\t" -vOFS="\t" '
#	BEGIN{print 0, "rand|gray|0|2|0", 0; print 1, "rand|gray|0|2|0", 1}
#
#	{
#		nx = $1 " (AUC=" $2 ")||0|s"
#		n=split($6,A,/;/)
#		print 0, nx, 0
#		for (i=1; i<=n; i++) 
#			print (A[i]-i)/($5-n), nx, i/n
#		print 1, nx, 1
#	}' | Plot.awk -F"\t" -vMin=0 -vMax=1 -vXTitle="1 - Specificity" -vYTitle="Sensitivity" > test.svg


BEGIN{
	# one of scat, {,s}bar{,h}
	if (Mode == "")
		Mode = "scat"

	Mode = tolower(Mode)
}

function floor (x) {return (x == int(x)) ? x : (x > 0 ? int(x) : int(x-1))}
function ceil (x) {return (x == int(x)) ? x : floor(x+1)}

# make min and max reasonably round numbers
# assume ticks divides 10 well (either 5 or 6)
function prettybound(min, max, upper,     x)
{
	if ((max+0) != max || (min+0) != min || min == max)
		return max + (upper ? 1 : -1)
	x = 10^floor(log(max-min)/log(10)-0.5)
	return (upper ? ceil(max/x) : floor(min/x))*x
}

END{
	############# COMMON OPTIONS
	if (LFontSize == "")
		LFontSize = 20

	if (HFontSize == "")
		HFontSize = 25

	if (SFontSize == "")
		SFontSize = 15

	if (Font == "")
		Font = "Calibri"

	# colors produced with http://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
	if (CL == "")
		CL = "(255,0,0) (0,255,0) (255,26,184) (255,211,0) (0,87,0) (131,131,255) (158,79,70) (0,255,193) (0,131,149) (0,0,123) (149,211,79) (246,158,219) (211,17,255) (123,26,105) (246,17,96) (255,193,131) (35,35,8) (140,167,123) (246,131,8) (131,114,0) (114,246,255) (158,193,255) (114,96,123) (158,0,0) (0,79,255) (0,70,149) (211,255,0) (184,79,211) (61,0,26)"


	# 0. No legend
	# 1. Outside, on right
	# 2-4. four corners (TL, TR, BL, BR)
	if (DoLegend == "")
		DoLegend = 1

	if (Quiet == "")
		Quiet = 0

	if (Title == "")
		Title == ""

	# note: (X/Y)Title have positions flipped in {,s}barh
	if (XTitle == "")
		XTitle == ""

	if (YTitle == "")
		YTitle == ""

	# format for values at tick marks
	# if set to #BASE uses BASE^val
	if (AxisFmt == "")
		AxisFmt = "%s"

	if (XAxisFmt == "")
		XAxisFmt = AxisFmt

	if (YAxisFmt == "")
		YAxisFmt = AxisFmt

	# whether series are numbers that should be made pretty
	if (NumericLegend == "")
		NumericLegend = 0

	if (XValTextFlip == "")
		XValTextFlip = 0

	############# OPTIONS FOR ONLY Mode == scat
	if (PSize == "")
		PSize = 2

	if (NumXTick == "")
		NumXTick = 5

	if (NumYTick == "")
		NumYTick = 5

	if (ML == "") # use " " to disable
		ML = "c"

	if (LL == "")
		LL = "1: 1:9,5 1:5,5 1:1,2"

	if (XLog == "")
		XLog = 0
	
	if (YLog == "")
		YLog = 0

	# set to 0 to disable stroke around points
	if (PointStroke == "")
		PointStroke = ""

	# how much of a margin is between the XMin bar and the first tick (must be positive)
	if (XMinTickMargin == "")
		XMinTickMargin = 0

	if (XMaxTickMargin == "")
		XMaxTickMargin = 0

	if (YMinTickMargin == "")
		YMinTickMargin = 0

	if (YMaxTickMargin == "")
		YMaxTickMargin = 0

	############# OPTIONS FOR Mode ~ /bar/
	if (BarWidth == "")
		BarWidth = 35

	if (BarSep == "")
		BarSep = (Mode ~ /sbar/) ? 0 : 50

	if (BarSepLines == "")
		BarSepLines = 0

	if (BarMargin == "")
		BarMargin = (Mode ~ /sbar/) ? 10 : 0

	if (BarZero == "")
		BarZero = 0
	
	if (NumTick != "")
		NumXTick = NumYTick = NumTick

	if (Min != "")
		XMin = YMin = Min

	if (Max != "")
		XMax = YMax = Max

	if (MinTickMargin != "")
		XMinTickMargin = YMinTickMargin = MinTickMargin

	if (MaxTickMargin != "")
		XMaxTickMargin = YMaxTickMargin = MaxTickMargin

	if (Log != "")
		XLog = YLog = Log

	if (Dither == "")
		Dither = 0

	# list of values for which a line should be produced
	if (Grid != "")
		XGrid = YGrid = Grid

	############# X/Y MIN/MAX
	# compute the Min/Max found in the data
	for (i in ValV)
		if ((ValV[i]+0) == ValV[i])
		{
			split(i, A, SUBSEP)
			d=A[2] <= 2 ? A[2] : 2

			if (MinV[d] == "" || ValV[i] < MinV[d])
				MinV[d] = ValV[i]
			if (MaxV[d] == "" || ValV[i] > MaxV[d])
				MaxV[d] = ValV[i]
		}
	for (i in ValE)
		if ((ValE[i]+0) == ValE[i])
		{
			split(i, A, SUBSEP)
			d=A[2] <= 2 ? A[2] : 2

			if (MinV[d] == "" || ValE[i] < MinV[d])
				MinV[d] = ValE[i]
			if (MaxV[d] == "" || ValE[i] > MaxV[d])
				MaxV[d] = ValE[i]
		}

	# use same Min/Max for XY
	if (XYMinMax)
	{
		MinV[1] = MinV[2] = min(MinV[1], MinV[2])
		MaxV[1] = MaxV[2] = min(MaxV[1], MaxV[2])
		XMinTickMargin = YMinTickMargin = max(XMinTickMargin, YMinTickMargin)
		XMaxTickMargin = YMaxTickMargin = max(XMaxTickMargin, YMaxTickMargin)
	}

	# min must be negative of max
	if (SymMinMax)
	{
		for (i=1; i<=2; i++)
		{
			MinV[i] = min(MinV[i], -MaxV[i])
			MaxV[i] = -MinV[i]
		}
		XMinTickMargin = XMaxTickMargin = max(XMinTickMargin, XMaxTickMargin)
		YMinTickMargin = YMaxTickMargin = max(YMinTickMargin, YMaxTickMargin)
	}

	# NOTE --> does not behave well with Log
	if (XMin == "")
		XMin = prettybound(MinV[1] + XMinTickMargin, MaxV[1] - XMaxTickMargin, 0) - XMinTickMargin

	if (XMax == "")
		XMax = prettybound(MinV[1] + XMinTickMargin, MaxV[1] - XMaxTickMargin, 1) + XMaxTickMargin

	if (YMin == "")
		YMin = prettybound(MinV[2] + YMinTickMargin, MaxV[2] - YMaxTickMargin, 0) - YMinTickMargin

	if (YMax == "")
		YMax = prettybound(MinV[2] + YMinTickMargin, MaxV[2] - YMaxTickMargin, 1) + YMaxTickMargin

	# use same Min/Max for XY
	if (XYMinMax)
	{
		XMin = YMin = min(XMin, YMin)
		XMax = YMax = max(XMax, YMax)
	}

	# min must be negative of max
	if (SymMinMax)
	{
		XMin = min(XMin, -XMax)
		XMax = max(-XMin, XMax)
		YMin = min(YMin, -YMax)
		YMax = max(-YMin, YMax)
	}

	############# OPTIONS FOR Mode == scat / bar / sbar
	# if "-" then automatically set to keep aspect ratio same between height/width
	if (BoxHeight == "")
		BoxHeight = 500


	############# OPTIONS FOR Mode == scat / barh / sbarh
	# if "-" then automatically set to keep aspect ratio same between height/width
	if (BoxWidth == "")
		BoxWidth = 500


	#############################################################################
	Cs["blue"] = "(0,0,255)"
	Cs["red"] = "(255,0,0)"
	Cs["white"] = "(255,255,255)"
	Cs["yellow"] = "(255,255,0)"
	Cs["green"] = "(0,255,0)"
	Cs["black"] = "(0,0,0)"
	Cs["orange"] = "(255,127,0)"
	Cs["gray"] = "(175,175,175)"
	Cs["violet"] = Cs["purple"] = "(127,0,255)"

	################# SETUP APPROX WIDTH OF CHARACTERS
	split("! \" # $ % & ' ( ) * + , - . / 0 1 2 3 4 5 6 7 8 9 : ; < = > ? @ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z [ \\ ] ^ _ ` a b c d e f g h i j k l m n o p q r s t u v w x y z { | }", A, / /)
	split("18 26 39 36 55 80 12 20 19 34 37 15 22 15 30 38 36 36 36 38 36 38 37 38 37 15 16 16 36 37 34 65 45 40 41 46 35 34 45 43 14 19 40 33 62 45 49 38 51 41 34 39 45 44 69 40 38 36 21 28 17 36 80 18 33 38 31 36 36 25 36 37 13 14 35 13 58 37 39 38 36 27 28 25 36 35 55 33 35 29 23 21 24", B, / /)
	for (i in A)
		CharWidth[A[i]] = B[i]
	CharWidth[" "] = 12

	################# FILL IN SERIES ATTRIBUTES
	MaxSeriesWidth = 0
	NCL = split(CL, CLA, / /); cC = 1
	NML = split(ML, MLA, / /); cM = 1
	split(LL, LLA, / /);

	for (s=1; s<=NS; s++)
	{
		split(SStr[s], A, /[|]/)

		SName[s] = A[1]

		if (NumericLegend)
			SName[s] = prettynum(SName[s], 0, AxisFmt)

		split(A[2], B, /:/)

		SColor[s] = B[1]
		SStroke[s] = B[2]

		if (SColor[s] == "")
		{
			SColor[s] = CLA[cC]
			if (++cC > NCL)
				cC = 1
		}
		else if (SColor[s] == "0")
			SColor[s] = ""
		else if ((SColor[s]+0) == SColor[s])
			SColor[s] = CLA[SColor[s]]

		if (SColor[s] in Cs)
			SColor[s] = Cs[SColor[s]]

		if (SColor[s] != "")
			SColor[s] = "rgb" SColor[s]

		if (SStroke[s] == "")
			SStroke[s] = "black"
		else if (SStroke[s] == 0)
			SStroke[s] = ""
		else if ((SStroke[s]+0) == SStroke[s])
			SStroke[s] = CLA[SStroke[s]]

		if (SStroke[s] in Cs)
			SStroke[s] = Cs[SStroke[s]]

		if (SStroke[s] != "")
			SStroke[s] = "rgb" SStroke[s]

		SMarker[s] = A[3]
		if (SMarker[s] == "0")
			SMarker[s] = ""
		else if (SMarker[s] == "" && Mode != "scat")
			SMarker[s] = "b"
		else if (SMarker[s] == "")
		{
			SMarker[s] = MLA[cM]
			if (++cM > NML)
				cM = 1
		}
		else if ((SMarker[s]+0) == SMarker[s])
			SMarker[s] = MLA[SMarker[s]]


		if (substr(A[4],1,1) == "s" || substr(A[4],1,1) == "S")
		{
			SStep[s] = substr(A[4],1,1) == "s"  ? 1 : 2
			SLFmt[s] = substr(A[4],2)
		}
		else
			SLFmt[s] = A[4]

		if (SLFmt[s] == "0" || (SLFmt[s] == "" && Mode != "scat"))
			SLFmt[s] = ""
		else if (SLFmt[s] == "")
			SLFmt[s] = LLA[1]
		else if ((SLFmt[s]+0) == SLFmt[s])
			SLFmt[s] = LLA[SLFmt[s]]
		else
			SLFmt[s] = SLFmt[s]
		
		SLegend[s] = A[5]
		if (SLegend[s] == "")
			SLegend[s] = 1

		if (SLegend[s])
		{
			MaxSWidth = max(strwidth(SName[s], SFontSize), MaxSWidth)
			NumLegendS++
		}
	}

	################# CONVERT VALUES TO DIMENSIONS
	VMin[1] = XMin; VMax[1] = XMax; VLog[1] = XLog
	VNumTick[1] = NumXTick
	VMinTickMargin[1] = XMinTickMargin; VMaxTickMargin[1] = XMaxTickMargin
	BoxSize[1] = BoxWidth 
	AFmt[1] = XAxisFmt

	VMin[2] = YMin; VMax[2] = YMax; VLog[2] = YLog
	VNumTick[2] = NumYTick
	VMinTickMargin[2] = YMinTickMargin; VMaxTickMargin[2] = YMaxTickMargin 
	BoxSize[2] = BoxHeight
	AFmt[2] = YAxisFmt

	if (Mode == "scat")
	{
		# automatically determine BoxSize when - by maintaining aspect ratio
		for (i=1; i<=2; i++)
			if (BoxSize[i] == "-")
				BoxSize[i] = BoxSize[3-i] * abs((VMax[i]-VMin[i]) / (VMax[3-i]-VMin[3-i]))
	}
	else
		BoxSize[(Mode !~ /barh/) ? 1 : 2] = (BarWidth * ((Mode ~ /sbar/) ? 1 : NS) + BarSep) * NX + BarMargin * 2

	if (Mode ~ /barh/)
	{
		temp = XTitle
		XTitle = YTitle
		YTitle = temp
		
		VMin[1] = VMin[2]
		VMax[1] = VMax[2]
		XGrid = YGrid
	}

	for (d=1; d<=2; d++)
	{
		NumVGrid[d] = split(d==1 ? XGrid : YGrid, A, / /)
		for (i=1; i<=NumVGrid[d]; i++)
			VGrid[d,i] = A[i]
	}


	################# REPLACE %%[XY](Min|Max)%%
	ValRepl["%%XMin%%"] = VMin[1]; ValRepl["%%YMin%%"] = VMin[2]; ValRepl["%%Min%%"] = min(VMin[1],VMin[2]);
	ValRepl["%%XMax%%"] = VMax[1]; ValRepl["%%YMax%%"] = VMax[2]; ValRepl["%%Max%%"] = max(VMax[1], VMax[2]);
	for (i in ValV)
		if (ValV[i] in ValRepl)
			ValV[i] = ValRepl[ValV[i]]
	for (i in ValE)
		if (ValE[i] in ValRepl)
			ValE[i] = ValRepl[ValE[i]]
	delete ValRepl

	################# TRANSFORM INTO LOG SPACE, IF NECESSARY
	# will also need to adjust the tick values
	for (d=1; d<=3; d++)
		if ((d == 2 || (Mode == "scat" && d == 1) || (Mode != "scat" && d == 3)) && VLog[min(d,2)])
		{
			for (v=1; v<=NumVal; v++)
			{
				ValV[v, d] = log(ValV[v, d])
				
				for (i=1; i<=NumValE[v,d]; i++)
					ValE[v,d,i] = log(ValE[v,d,i])
			}

			if (d != 3)
			{
				VMin[d] = log(VMin[d]) 
				VMax[d] = log(VMax[d])

				for (i=1; i<=NumVGrid[d]; i++)
					VGrid[d,i] = log(VGrid[d,i])
			}
		}

	################# SETUP SOME POSITIONS WITHIN THE SVG
	# this is the margin on the slide of the plot with the max/min plotted value
	# note: for Y values the max is plotted at the bottom
	if (Mode ~ /barh/)
	{
		for (x=1; x<=NX; x++)
			if (MaxYTitleLen == "" || strwidth(XStr[x], SFontSize) > MaxYTitleLen)
				MaxYTitleLen = strwidth(XStr[x], SFontSize)
	}
	else
		for (i=0; i<VNumTick[2]; i++)
		{
			iv = tickv(VMin[2], VMax[2], VMinTickMargin[2], VMaxTickMargin[2], VNumTick[2], i)
			io = prettynum(iv, VLog[2], AFmt[2])
			if (MaxYTitleLen == "" || strwidth(io, SFontSize) > MaxYTitleLen)
				MaxYTitleLen = strwidth(io, SFontSize)
		}

	if (!XValTextFlip)
		MaxXTitleLen = SFontSize
	else if (Mode ~ /bar$/)
	{
		for (x=1; x<=NX; x++)
			if (MaxXTitleLen == "" || strwidth(XStr[x], SFontSize) > MaxXTitleLen)
				MaxXTitleLen = strwidth(XStr[x], SFontSize)
	}
	else
		for (i=0; i<VNumTick[1]; i++)
		{
			iv = tickv(VMin[1], VMax[1], VMinTickMargin[1], VMaxTickMargin[1], VNumTick[1], i)
			io = prettynum(iv, VLog[1], AFmt[1])
			if (MaxXTitleLen == "" || strwidth(io, SFontSize) > MaxXTitleLen)
				MaxXTitleLen = strwidth(io, SFontSize)
		}

	MinMargin[1] = strheight(YTitle) * LFontSize + MaxYTitleLen + 20   # at the left
	MaxMargin[1] = 10                                                  # at the right
	MaxMargin[2] = strheight(Title) * HFontSize + 20                   # at the top
	MinMargin[2] = strheight(XTitle) * LFontSize + MaxXTitleLen + 20   # at the bottom

	# these are the positions corresponding to the smallest/biggest values on the plot
	BMin[1] = MinMargin[1]
	BMax[1] = MinMargin[1] + BoxSize[1]

	BMin[2] = BoxSize[2] + MaxMargin[2]
	BMax[2] = MaxMargin[2]
	
	SvgWidth = BMax[1] + MaxMargin[1]
	SvgHeight = BMin[2] + MinMargin[2]

	if (NumLegendS == 0)
		DoLegend = 0

	if (DoLegend != 0)
	{
		LegendWidth = MaxSWidth + 55 + LegendWidthDelta
		LegendHeight = 10 + NumLegendS * (5 + SFontSize)
	}

	if (DoLegend == 1)
	{
		SvgWidth += LegendWidth + 10
		SvgHeight = max(SvgHeight, BMax[2] + 10 + LegendHeight) 
	}
	else
		SvgWidth += 20 # add some space for axis label

	print "<?xml version=\"1.0\" standalone=\"no\"?><!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
	print "<svg width=" Q(SvgWidth) " height=" Q(SvgHeight) " version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">"

	################# PRODUCE GRID LINES
	for (d=1; d<=2; d++)
		if (!(d==1 && Mode~/bar$/) && !(d==2 && Mode~/barh/))
			for (i=1; i<=NumVGrid[d]; i++)
			{
				val = map(d, VGrid[d,i])
				print poly(d, val " " val, BMin[3-d] " " BMax[3-d], "", "black")
			}
		else if (Mode ~ /bar/ && BarSepLines)
			for (i=1; i<NX; i++)
			{
				val = ((d == 1) ? BMin[1] : BMax[2]) + BarMargin + i * (((Mode ~ /sbar/) ? 1 : NS) * BarWidth + BarSep)
				print poly(d, val " " val, BMin[3-d] " " BMax[3-d], "", "black")
			}

	################# PLOT VALUES
	if (Mode == "scat")
	{
		################# PLOT THE INDIVIDUAL POINTS (+error bars) AND LINES BETWEEN THEM 
		# start error bars first so they go under the plots
		for (v=1; v<=NumVal; v++)
		{
			sr = ValS[v]

			for (d=1; d<=3; d++)
				b[d] = map(d, ValV[v,d])

			# error bars
			for (d=1; d<=2; d++)
				for (i=1; i<=NumValE[v,d]; i++)
				{
					eb = map(d,ValE[v,d,i])
					x = line(d, b[d], eb, b[3-d], b[3-d], "black", 0.5)

					if (x != "")
					{	
						print x
						print line(d, eb, eb, (b[3-d]-4), (b[3-d]+4), "black", 0.5)
					}
				}
		}

		# lines between points
		for (v=1; v<=NumVal; v++)
			if (SValNum[v] > 1)
			{
				sr = ValS[v]

				for (d=1; d<=2; d++)
					b[d] = map(d, ValV[v,d])

				lv = SVal[sr,SValNum[v]-1]

				for (d=1; d<=2; d++)
					lb[d] = map(d, ValV[lv,d])

				if (SStep[sr] == 1)
				{
					# step line (over then up)
					print line(1, lb[1], b[1], lb[2], lb[2], SColor[sr], SLFmt[sr], PLen[sr])
					PLen[sr] += abs(b[1]-lb[1])
					print line(1, b[1], b[1], lb[2], b[2], SColor[sr], SLFmt[sr], PLen[sr])
					PLen[sr] += abs(b[2]-lb[2])
				}
				else if (SStep[sr] == 2)
				{
					# step line (up then over)
					print line(1, lb[1], lb[1], lb[2], b[2], SColor[sr], SLFmt[sr], PLen[sr])
					PLen[sr] += abs(b[2]-lb[2])
					print line(1, lb[1], b[1], b[2], b[2], SColor[sr], SLFmt[sr], PLen[sr])
					PLen[sr] += abs(b[1]-lb[1])
				}
				else
				{
					print line(1, lb[1], b[1], lb[2], b[2], SColor[sr], SLFmt[sr], PLen[sr])

					# keep track of the current length of this line for the dash array
					PLen[sr] += sqrt((b[1]-lb[1])^2 + (b[2]-lb[2])^2)
				}
			}

		# the points themselves 
		for (v=1; v<=NumVal; v++)
		{
			sr = ValS[v]

			for (d=1; d<=2; d++)
				b[d] = map(d, ValV[v,d])

			print point(b[1], b[2], SMarker[sr], SColor[sr], SStroke[sr], ValZ[v])

			# text
			if (ValT[v] != "")
				print text(d, ValT[v], b[1] + 5, b[2], "", "", SFontSize)
		}

	}
	else
	{
		d = (Mode ~ /barh/ ) ? 1 : 2

		# lines between bars/points
		for (v=1; v<=NumVal; v++)
			if (SValNum[v] > 1)
			{
				sr = ValS[v]
				x = BarX(v)
				y = map(d, ValV[v,2])
				lv = SVal[sr,SValNum[v]-1]
				lx = BarX(lv)
				ly = map(d, ValV[lv,2])

				print line(3-d, lx, x, ly, y, SColor[sr], SLFmt[sr], PLen[sr])
				PLen[sr] += sqrt((lx-x)^2 + (ly-y)^2)
			}

		# plot the bars/points/errorbars themselves
		for (v=1; v<=NumVal; v++)
		{
			sr = ValS[v]
			xm = BarX(v)
			x1 = xm - 0.5 * BarWidth
			x2 = xm + 0.5 * BarWidth

			if (SMarker[sr] != "b")
			{
				# print a point
				y = map(d, ValV[v,2])
				print point(d==2 ? xm : y, d==2 ? y : xm, SMarker[sr], SColor[sr], SStroke[sr], ValZ[v])
			}
			else if (ValV[v,2] != "")
			{
				y1 = trunc(map(d, ValV[v,2]), BMin[d], BMax[d])
				y2 = trunc(map(d, ValV[v,3] == "" ? BarZero : ValV[v,3]), BMin[d], BMax[d])
				
				if (!Quiet && (y1 != map(d, ValV[v, 2]) || (ValV[v,3] != "" && y2 != map(d, ValV[v, 3]))))
					print "Bar truncated!" > "/dev/stderr"

				print poly(3-d, x1 " " x1 " " x2 " " x2 , y1 " " y2 " " y2 " " y1, SColor[sr], SStroke[sr])
			}

			# error bars
			for (ed=2; ed<=3; ed++)
				for (i=1; i<=NumValE[v,ed]; i++)
				{
					eb = map(d,ValE[v,ed,i])
					xt = line(d, ed==2 ? y1 : y2, eb, xm, xm, "black", 0.5)

					if (xt != "")
					{
						print xt
						print line(d, eb, eb, (xm-4), (xm+4), "black", 0.5)
					}
				}

			# text
			if (ValT[v] != "")
			{
				# to/ta assume barh... flipped for bar
				t = ValT[v]

				if (ValTO[v] ~ /P/)
					t = prettynum(t, 0, "%s")

				# text orientation (1 -- vertical, "" -- horizontal)
				to = ValTO[v] ~ /S/ ? 1 : ""

				# L/B/Z/C/M/X/A/R (plot left, before bar, zero, plot center, bar center, bar end, after bar, right)
				tx = (ValTO[v]~/L/) ? BMin[d] : (ValTO[v]~/[BZ]/) ? y2 : (ValTO[v]~/C/) ? ((BMin[d] + BMax[d])/2) : (ValTO[v]~/M/) ? ((y1 + y2)/2) : (ValTO[v]~/[AX]/) ? y1 : BMax[d]

				# text alignment within bar ("" -- left, 1 -- right; 2 -- center)
				ta = ValTO[v] ~ /[CM]/ ? 2 : (ValTO[v] ~ /[LZA]/ ? "" : 1)
				tc = (ValTO[v] ~ /^\([0-9]*,[0-9]*,[0-9]*\)/) ? ("rgb" gensub(/\).*/, ")", 1, ValTO[v])) : "black"
				ts = (gensub(/^\([0-9]*,[0-9]*,[0-9]*\)/, "", 1, ValTO[v]) ~ /^[0-9.]/) ? (gensub(/^\([0-9]*,[0-9]*,[0-9]*\)([0-9.]*)[A-Za-z]*/, "\\1", 1, ValTO[v]) + 0) : SFontSize
				tx += (to ? (ts/2) : 5) * (d==1 ? 1 : -1) * (!ta ? 1 : ta == 1 ? -1 : 0)
				print text(d, t, tx, xm, to, to ? 2 : d==2 ? (ta ? "" : 1) : ta, ts, tc)
			}
		}
	}

	################# DRAW BOX AROUND PLOTTING AREA
	print poly(1, BMin[1] " " BMin[1] " " BMax[1] " " BMax[1], BMax[2] " " BMin[2] " " BMin[2] " " BMax[2], "", "black")

	# titles around boxes
	print text(1, Title, (BMin[1] + BMax[1])/2, BMax[2]/2, "", 2, HFontSize, "black")
	print text(1, XTitle, (BMin[1] + BMax[1])/2, SvgHeight - strheight(XTitle) * LFontSize/2 - 10, "", 2, LFontSize, "black")
	print text(2, YTitle, (BMin[2] + BMax[2])/2, strheight(YTitle) * LFontSize / 2 + 10, "", 2, LFontSize, "black")

	# and the tick marks
	for (d=1; d<=2; d++)
		if (!(d==1 && Mode~/bar$/) && !(d==2 && Mode~/barh/))
		{
			vd = (Mode~/barh/) ? 1 : d
			for (i=0; i<VNumTick[vd]; i++)
			{
				iv = tickv(VMin[vd], VMax[vd], VMinTickMargin[vd], VMaxTickMargin[vd], VNumTick[vd], i)

				tick = map(d,iv)
				print poly(d, tick " " tick, (BMin[3-d] - (d==2 ? -1 : 1) * 4) " " BMin[3-d], "", "black")
				print poly(d, tick " " tick, (BMax[3-d] + (d==2 ? -1 : 1) * 4) " " BMax[3-d], "", "black")

				if (d == 2)
					print text(d, prettynum(iv, VLog[vd], AFmt[vd]), tick, BMin[3-d] - 5, 1, "", SFontSize, "black")
				else if (!XValTextFlip)
					print text(d, prettynum(iv, VLog[vd], AFmt[vd]), tick, BMin[3-d] + (SFontSize/2) + 4, "", 2, SFontSize, "black")
				else
					print text(d, prettynum(iv, VLog[vd], AFmt[vd]), tick, BMin[3-d] + (SFontSize/2) + 4, 1, 1, SFontSize, "black")
			}
		}

	# X/Y titles for bar plots
	if (Mode~/bar$/)
		for (x=1; x<=NX; x++)
		{
			if (!XValTextFlip)
				print text(1, XStr[x], BMin[1] + BarMargin + (x-0.5) * (((Mode ~ /sbar/) ? 1 : NS) * BarWidth + BarSep), BMin[2] + (SFontSize/2) + 4, "", 2, SFontSize, "black")
			else
				print text(1, XStr[x], BMin[1] + BarMargin + (x-0.5) * (((Mode ~ /sbar/) ? 1 : NS) * BarWidth + BarSep), BMin[2] + (SFontSize/2) + 4, 1, 1, SFontSize, "black")
		}

	if (Mode~/barh/)
		for (x=1; x<=NX; x++)
			print text(2, XStr[x], BMax[2] + BarMargin + (x-0.5) * (((Mode ~ /sbar/) ? 1 : NS) * BarWidth + BarSep), BMin[1] - 5, 1, "", SFontSize, "black")

	################# PRODUCE THE LEGEND
	if (DoLegend != 0)
	{
		if (DoLegend == 1)
		{
			LegendTop = BMax[2]
			LegendLeft = BMax[1] + 10
		}

		if (DoLegend == 2 || DoLegend == 3)
			LegendTop = BMax[2] + 10

		if (DoLegend == 4 || DoLegend == 5)
			LegendTop = BMin[2] - LegendHeight - 10

		if (DoLegend == 2 || DoLegend == 4)
			LegendLeft = BMin[1] + 10

		if (DoLegend == 3 || DoLegend == 5)
			LegendLeft = BMax[1] - LegendWidth - 10

		print poly(1, LegendLeft " " LegendLeft " " (LegendLeft + LegendWidth) " " (LegendLeft + LegendWidth), LegendTop " " (LegendTop + LegendHeight) " " (LegendTop + LegendHeight) " " LegendTop, "", "black")

		nl = 0
		for (sr=1; sr<=NS; sr++)
			if (SLegend[sr])
			{
				nl++

				y = LegendTop + 5 + (nl-0.5) * (5 + SFontSize)

				if (SMarker[sr] == "b")
					print poly(1, (LegendLeft+10) " " (LegendLeft+10) " " (LegendLeft+35) " " (LegendLeft+35) , (y+5) " " (y-5) " " (y-5) " " (y+5), SColor[sr], SStroke[sr])
				print line(1, (LegendLeft+10), (LegendLeft+35), y, y, SColor[sr], SLFmt[sr], 0, 1)
				if (SMarker[sr] != "b")
					print point(LegendLeft+22.5, y, SMarker[sr], SColor[sr], SStroke[sr], "", 1)

				print text(1, SName[sr], LegendLeft+45, y-1, "", "", SFontSize, "black")
			}
	}

	
	################# END OF SVG
	print "</svg>"
}

# get the middle x-value for a bar plot
function BarX (v)
{
	if (v < 1 || v > NumVal)
		return

	if (Dither && !(v in DitherX))
		DitherX[v] = rand() * 0.8 - 0.4

	if (Mode ~ /sbar/)
		return ((d == 2) ? BMin[1] : BMax[2]) + 0.5 * BarSep + (ValX[v]-1) * (BarWidth + BarSep) + BarMargin + (0.5 + DitherX[v]) * BarWidth
	else
		return ((d == 2) ? BMin[1] : BMax[2]) + 0.5 * BarSep + (ValX[v]-1) * (NS * BarWidth + BarSep) + BarMargin + (ValS[v] - 0.5 + DitherX[v]) * BarWidth
}

# add commas to values
# also in TT.awk
function AC(x,   A, LN, L)
{
	split(x~/^-/ ? substr(x,2) : x,A,/[.]/)
	LN=length(A[1])
	while (LN>3)
	{ 
		L = "," substr(A[1], LN-2) L
		A[1] = substr(A[1],1,LN-3)
		LN -= 3
	}
	
	L = (x~/^-/ ? "-" : "") A[1] L ((2 in A) ? "." A[2] : "")
	
	return L
}

function min(a,b) { return a < b ? a : b }

# if x is outside of min->max, returns min/max (whichever is closer)
# note: min:max is a range and min>max
function trunc (x, min, max)
{
	if (min < max)
	{
		if (x < min)
			return min
		if (x > max)
			return max
	}
	else
	{
		if (x > min)
			return min
		if (x < max)
			return max
	}

	return x
}

# produces a line clipped to the plotting area
function line (dim, x1, x2, y1, y2, str, lfmt, prelen, nocheck,   A)
{
	if (lfmt == "") return

	# force dim to 1
	if (dim == 2)
		return line(1, y1, y2, x1, x2, str, lfmt, prelen)

	if (!nocheck)
	{
		# if any point is blank, or if both of the same dimension are outside, we are done
		if (x1 == "" || x2 == "" || y1 == "" || y2 == "" || (x1 < BMin[1] && x2 < BMin[1]) || (x1 > BMax[1] && x2 > BMax[1]) || (y1 > BMin[2] && y2 > BMin[2]) || (y1 < BMax[2] && y2 < BMax[2]))
			return

		# if we intersect any of the four edges, clip the one outside to the inside
		if ((x1 >= BMin[1]) != (x2 >= BMin[1]))
		{
			if (x1 > x2)
				return line(dim, x2, x1, y2, y1, str, lfmt, prelen)

			return line(dim, BMin[1], x2, y1 + (y2-y1) / (x2 - x1) * (BMin[1] - x1), y2, str, lfmt, prelen)
		}

		if ((x1 <= BMax[1]) != (x2 <= BMax[1]))
		{
			if (x1 < x2)
				return line(dim, x2, x1, y2, y1, str, lfmt, prelen)

			return line(dim, BMax[1], x2, y1 + (y2-y1) / (x2 - x1) * (BMax[1] - x1), y2, str, lfmt, prelen)
		}

		if ((y1 >= BMax[2]) != (y2 >= BMax[2]))
		{
			if (y1 > y2)
				return  line(dim, x2, x1, y2, y1, str, lfmt, prelen)

			return line(dim, x1 + (x2-x1)/(y2-y1) * (BMax[2] - y1), x2, BMax[2], y2, str, lfmt, prelen)
		}

		if ((y1 <= BMin[2]) != (y2 <= BMin[2]))
		{
			if (y1 < y2)
				return  line(dim, x2, x1, y2, y1, str, lfmt, prelen)

			return line(dim, x1 + (x2-x1)/(y2-y1) * (BMin[2] - y1), x2, BMin[2], y2, str, lfmt, prelen)
		}
	}

	split(lfmt, A, /:/)
	return "<line x1=" Q(x1) " x2=" Q(x2)" y1=" Q(y1)" y2=" Q(y2) " stroke=" Q(str) " stroke-width=" Q(X(A[1],A[1],1)) X(A[2], " style=\"stroke-dasharray: " A[2] "; stroke-dashoffset: " (prelen+0) ";\"", "") "/>"
}



function Q(x) { return "\"" x "\"" }
function X(x,y,def) { return x != "" ? y : def }

# replaces %%tag%%text%%/tag%% with replace (which can include \\1 for the text)
function tagrep (t, tag, replace) {
	# first replace the end with SUBSEP (i.e. something that should never be in the text)
	gsub("%%/" tag "%%", SUBSEP, t)

	# now do the replacement
	return gensub("%%" tag "%%([^" SUBSEP "]*)" SUBSEP, replace, "g", t)
}

# like tagrep except tags of form %%tag=val%%text%%/tag%% with replace having \\1 and \\2
function tbgrep  (t, tag, replace) {
	gsub("%%/" tag "%%", SUBSEP, t)
	return gensub("%%" tag "=([^%]*)%%([^" SUBSEP "]*)" SUBSEP, replace, "g", t)
}

function strheight(x, A)
{
	return split(x, A, /%%Newline%%/)
}

function strwidth(x, size,    i, l, j, mx, A)
{
	gsub(/%%Semicolon%%/, ";", x)
	gsub(/%%VBar%%/, "|", x)
	split(x, A, /%%Newline%%/)

	# compute max over all lines
	for (j in A)
	{
		l = 0
		gsub(/%%[^%]*%%/, "", A[j])
		gsub(/&#[^;]*;/, "W", A[j])
		for (i=1; i<=length(A[j]); i++)
			l += CharWidth[substr(A[j],i,1)]
		if (l >= mx)
			mx = l
	}

	return mx * (size == "" ? FontSize : size)/50
}

# produce axis labels in nice format
function prettynum (v, lg, fmt,  A)
{
	# TODO: this needs to be designed so that a variety of formats
	#       can be used
	if (fmt ~ /^#/)
	{
		v = substr(fmt,2) "%%Super%%" AC(log(lg ? exp(v) : v)/log(substr(fmt,2)))
		v = v "%%/Super%%"
	}
	else
	{
		v = sprintf(fmt, lg ? exp(v) : v)

		if (v ~ /^[0-9.-]*e[-+0-9]*$/)
		{
			split(v, A, /e/)

			v = ""

			if (A[1] ~ /^-/)
			{
				v = "-"
				A[1] = substr(A[1],2)
			}

			sub (/[.]?0*$/, "", A[1])

			if (A[1] != "1")
				v = v A[1]
				
			if (fmt ~ /e/ || (A[2]+0)!=0)
			{
				if (v != "" && v != "-")
					v = v " x "
				v = v "10%%Super%%" AC(A[2]+0) "%%/Super%%"
			}

			if (v == "")
				v = 1
			else if (v == "-")
				v = "-1"
		}
		else 
			v = AC(v)
	}

	return v
}

########################## FUNCTIONS THAT ACTUALLY WRITE SVG
# align "" --> left
# align 1  --> right
# align 2  --> center
function text (dim, t, x, y, rot, align, size, txtcolor,  temp) 
{
	gsub(/%%Semicolon%%/, ";", t)
	gsub(/%%VBar%%/, "|", t)
	if (dim == 2)
	{
		temp = x
		x = y
		y = temp
		rot = rot ? "" : 1
		align = align == 2 ? align : (align == "" ? 1 : "")
	}

	if (size == "") size = FontSize
	if (rot)
		x += (size+1)/3
	else
		y += (size+1)/3

	gsub(/</, "\\&lt;", t)
	gsub(/>/, "\\&gt;", t)

	t = tagrep(t, "Super", "<tspan font-size=" Q("70%") " baseline-shift="Q("super")">\\1</tspan>")
	t = tagrep(t, "Sub", "<tspan font-size=" Q("70%") " baseline-shift="Q("sub")">\\1</tspan>")
	t = tagrep(t, "Monospace", "<text font-family=" Q("Monospace") ">\\1</text>")

	t = tagrep(t, "Bold", "<text font-weight=" Q("bold") ">\\1</text>")
	t = tagrep(t, "Underline", "<text text-decoration=" Q("underline") ">\\1</text>")
	t = tagrep(t, "Italic", "<text font-style=" Q("italic") ">\\1</text>")

	t = tbgrep(t, "Size", "<tspan font-size=" Q("\\1") ">\\2</tspan>")
	t = tbgrep(t, "Color", "<tspan fill=" Q("\\1") ">\\2</tspan>")

	# deal with multiple lines
	temp = strheight(t)
	if (temp > 1)
	{
		# add a non-breaking space between every pair of newlines
		gsub(/%%Newline%%%%Newline%%/, "%%Newline%%\\&#xa0;%%Newline%%", t)

		# also add before and after newline if at start/end
		if (t ~ /^ *%%Newline%%/)
			t = "&#xa0;" t

		if (t ~ /%%Newline%% *$/)
			t = t "&#xa0;" 

		# shifting the whole text up half the total height
		t = "<tspan dy=" Q(-0.5*(size * (temp-1))) ">" t "</tspan>"

		# replace newlines with an appropriate code
		gsub(/%%Newline%%/, "</tspan><tspan x=" Q(x) " dy=" Q(size) ">", t)
	}

	return X(t,"<text xml:space=\"preserve\" x=" Q(x) " y=" Q(y) " font-family=" Q(Font) " font-size=" Q(size) X(txtcolor, " fill=" Q(txtcolor)) X(rot, " transform=" Q("rotate(270 " x " " y ")")) X(align, " text-anchor=" Q(align==1 ? "end" : (align==2 ? "middle" : "start"))) ">" t "</text>")
}

function map (d, val)
{
	if (val == "" || VMax[d] == VMin[d])
		return ""
	return (val - VMin[d]) / (VMax[d] - VMin[d]) * (BMax[d] - BMin[d]) + BMin[d]
}

# x, y are in svg coordinates
function point (x, y, marker, cl, str, r, nocheck)
{
	if (r == "")
		r = PSize

	if (nocheck || ((x >= BMin[1]) && (x <= BMax[1]) && (y <= BMin[2]) && (y >= BMax[2])))
	{
		if (marker == "c")
			return X(cl str, "<circle cx=" Q(x) " cy=" Q(y) " r=" Q(r) " fill=" Q(X(cl,cl,"none")) X(str, " stroke=" Q(str) " stroke-width=" Q(PointStroke)) "/>")
		else if (marker == "d")
			return poly(1, (x-r) " " x " " (x+r) " " x, y " " (y+r) " " y " " (y-r), cl, str, PointStroke)
		else if (marker == "t")
			return poly(1, x " " (x-sqrt(3)/2*r) " " (x+sqrt(3)/2*r), (y+r) " " (y-r/2) " " (y-r/2), cl, str, PointStroke)
		else if (marker == "u")
			return poly(1, x " " (x-sqrt(3)/2*r) " " (x+sqrt(3)/2*r), (y-r) " " (y+r/2) " " (y+r/2), cl, str, PointStroke)
		else if (marker == "s")
			return poly(1, (x - r/sqrt(2)) " " (x + r/sqrt(2)) " " (x + r/sqrt(2)) " " (x - r/sqrt(2)), (y + r/sqrt(2)) " " (y + r/sqrt(2)) " " (y - r/sqrt(2)) " " (y - r/sqrt(2)), cl, str, PointStroke)
		else if (marker == "")
			return
	}
	
	if (!Quiet)
		print "Point outside plotting area!" > "/dev/stderr"
}

function poly (dim,  xs, ys, cl, str,  strw,      Xa, Ya, n, i, L)
{
	n = split(dim==2 ? ys : xs, Xa, / /)
	    split(dim==2 ? xs : ys, Ya, / /)

	for (i=1; i<=n; i++)
		L = L " " Xa[i] "," Ya[i]
	
	return X(cl str, "<polygon points=" Q(substr(L,2)) " fill=" Q(X(cl,cl,"none")) X(str, " stroke=" Q(str) " stroke-width=" Q(X(strw,strw,StrokeW))) "/>")
}

function tickv (mn, mx, mnmarg, mxmarg, ntick, i,   rev)
{
	rev = (mn > mx) ? -1 : 1

	mn += mnmarg * rev
	mx -= mxmarg * rev

	return (mn * (ntick-1-i) + mx * i) / (ntick-1)
}

function max(a,b) { return (a > b) ? a : b }
function abs(a,b) { return (a >= 0) ? a : (-a) }

########################## READ IN STDIN
{
	if (!($2 in SNum))
	{
		SNum[$2] = ++NS
		SStr[NS] = $2
	}

	if (Mode == "scat" && $1 == "" && $3 == "")
		next

	if (Mode != "scat")
	{
		# keep track of all the unique x values
		if (!($1 in XNum))
		{
			XNum[$1] = ++NX
			XStr[NX] = $1
		}

		if ($3 == "")
			next
	}

	++NumVal

	split($1, A, /[|]/)
	split($3, B, /[|]/)

	ValS[NumVal] = SNum[$2]
	ValV[NumVal,1] = A[1]
	ValV[NumVal,2] = B[1]

	# custom size of point
	# note: in /bar/ mode this is usually not ValZ
	#       except for when marker is defined
	#       and we are making a point
	ValZ[NumVal] = B[3]

	# keep track of the values from a series (for drawing lines between them)
	SVal[SNum[$2],++NumSVal[SNum[$2]]] = NumVal
	SValNum[NumVal] = NumSVal[SNum[$2]]

	if (Mode == "scat")
	{
		ValT[NumVal] = B[4]
		ValTO[NumVal] = B[5]
	}
	else
	{
		# keep track of all the values for the XVal (which is a name)
		XVal[XNum[$1], ++NumXVal[XNum[$1]]] = NumVal
		ValX[NumVal] = XNum[$1]

		# lower range of bars
		ValV[NumVal,3] = B[3]

		split(B[4], BB, /[:]/)
		for (i=1; i in BB; i++)
			ValE[NumVal,3,++NumValE[NumVal,3]] = BB[i]

		# text and textopts for bar
		ValT[NumVal] = B[5]
		ValTO[NumVal] = B[6]
	}

	# error bars (both may not be used for bar plots)
	split(A[2], AA, /[:]/)
	for (i=1; i in AA; i++)
		ValE[NumVal,1,++NumValE[NumVal,1]] = AA[i]
	split(B[2], BB, /[:]/)
	for (i=1; i in BB; i++)
		ValE[NumVal,2,++NumValE[NumVal,2]] = BB[i]
}

