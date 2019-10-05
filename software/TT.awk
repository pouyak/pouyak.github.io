#!/usr/bin/awk -f
# Written by Pouya Kheradpour (pouyak@gmail.com)
# Copyright (c) 2005-2015 Massachusetts Institute of Technology
# Released under the MIT license
# http://compbio.mit.edu/pouyak/software/LICENSE.mit

# version 20140120

# performs various operations on tables

# PARAMETERS
#
#   IF/OF          Input/output format  (default: 3)
#                  3      3 column format (X Y Val)
#                  m      matrix format
#                 
#   OF             Special output only formats
#                  mhtml  matrix format (html) 
#                  mtex   matrix format (latex tabular)
#               
#   D              Default value (default: "")
#   OFmt           Output format to apply (default: %s)
#                  Custom type %C: add's commas to float/int data types 
#   OY             FS separated column names to print in order (default: all)
#   OX             FS separated row names to print in order (default: all)
#               
#   Ops            Space separated list of operations to perform on data (in order)
#                  X - rows, Y - columns, M - matrix
#                  Scale01[XYM]  Scale from 0->1
#                  NormMStd[XYM] Subtract the mean, divide by stdev
#                  SubMean[XYM]  Subtract the mean
#                  DivSum[XYM]   Divide by sum
#                  DivSumB       Divide each value by the product of the row/column sums
#                  DivAMax[XYM]  Divide by the maximum absolute value
#                  Transpose     Switch X/Y
#                  Flip[XY]      Flip the order of the row/columns
#                  Sum[XY]       Sum each row/column
#                  Mean[XY]      Mean each row/column
#                  Stats[XY]     N, Mean and stdev for each row/column
#                  Max[XY]       Max for each row/column
#                  Min[XY]       Max for each row/column
#                  Rank[XY]      Replace values with their ranks rows/column
#                  Quant[XY]     Perform quantile normalization (to median) on rows/column
#
#   OpsPre         Prefix to add to each op name (useful for making unique)
#               
#   for IF==3   
#      XC/YF/VC    columns for X, Y, Value (default: 1/2/3)
#               
#   for IF==m   
#      HasXH       whether input has a names for rows (default: 1)
#      HasYH       whether input has a names for columns (default: 1)
#      HasN1       whether there is an N1 (only meaningful when HasYH=1) (default: 1)
#      IX          FS separated header row names
#      IY          FS separated header column names
#               
#   for OF==3   
#      OutM        Output missing values using D (default: 1)
#                  Note: most operations will effectively override this to 1
#               
#   for OF==m   
#      N1          value in upper left (default: "")
#      OutXH       print out header column (default: 1 if available, otherwise 0)
#      OutYH       print out header row (default: 1 if available, otherwise 0)
#      OutN1       print out N1 with OutXH and OutYH (default: 1)
#
#   for IF==m and OF==m
#      SaveN1      Save the N1 value from the input (default: 0)
#
#   for OF==mhtml
#      HtmlBorder  Border width (default: 2)
#

BEGIN{
	if (IF == "")
		IF = 3

	if (OF == "")
		OF = 3

	if (XC == "")
		XC = 1

	if (YC == "")
		YC = 2

	if (VC == "")
		VC = 3

	if (N1 == "")
		N1 = ""

	if (D == "")
		D = ""

	if (OFmt == "")
		OFmt = "%s"

	if (HasYH == "")
		HasYH = 1

	if (HasXH == "")
		HasXH = 1

	if (IX != "")
		nX = split(IX, X)

	if (IY != "")
		nY = split(IY, Y)

	if (SaveN1 == "")
		SaveN1 = 0

	if (HasN1 == "")
		HasN1 = 1

	if (OutYH == "")
		OutYH = (HasYH || IF==3 || IY != "") ? 1 : 0

	if (OutXH == "")
		OutXH = (HasXH || IF==3 || IX != "") ? 1 : 0

	if (OutN1 == "")
		OutN1 = 1

	if (OutM == "")
		OutM = 1

	if (HtmlBorder == "")
		HtmlBorder = 2
}

IF == 3{
	if (!($XC in Xi))
	{
		Xi[$XC] = ++nX
		X[nX] = $XC
	}

	if (!($YC in Yi))
	{
		Yi[$YC] = ++nY
		Y[nY] = $YC
	}

	V[Xi[$XC], Yi[$YC]] = $VC	
}

IF == "m"{
	if (HasYH && NR == 1)
	{
		if (SaveN1 && HasXH && HasN1)
			N1 = $1

		if (IY == "")
			for (i=(HasXH && HasN1) ? 2 : 1; i<=NF; i++)
				Y[++nY] = $i
	}
	else
	{
		x = NR - (HasYH ? 1 : 0)

		if (!HasYH)
			while((NF - (HasXH ? 1 : 0)) > nY)
			{	
				++nY
				Y[nY] = nY
			}

		if (IX == "")
		{
			X[x] = HasXH ? $1 : x
			nX = x
		}

		for (i=HasXH ? 2 : 1; i<=NF; i++)
			V[x, i - (HasXH ? 1 : 0)] = $i
	}
}

# transpose data
function T(  i, temp1, temp2)
{
	# swap V
	for (i in V) temp1[i] = V[i]
	delete V
	for (i=1; i<=nX; i++)
		for (j=1; j<=nY; j++)
			if ((i,j) in temp1)
				V[j,i] = temp1[i,j]

	# swap Y/X
	for (i in X) temp2[i] = X[i]
	delete X
	for (i in Y) X[i] = Y[i]
	delete Y
	for (i in temp2) Y[i] = temp2[i]

	# swap nX/nY
	i = nX
	nX = nY
	nY = i
}

function SelRows (OX,   Xi, Vn, i, j) {
	# generate Xi (we ignore the one in the higher scope)
	for (i=1; i<=nX; i++)
		Xi[X[i]] = i

	# overwrite the old X
	nX = split(OX, X)

	for (i=1; i<=nX; i++)
		for (j=1; j<=nY; j++)
			if ((Xi[X[i]],j) in V)
				Vn[i,j] = V[Xi[X[i]],j]
	delete V

	for (i in Vn) 
		V[i] = Vn[i]
}


# reshapes matrix (does NOT update X, Y)
# nXnew must evenly divide nX * nY
function ReShape(nXnew,   nYnew, Vn, i, j, ni, nj) {
	nYnew = nX * nY / nXnew
	ni = 1
	for (i=1; i<=nX; i++)
		for (j=1; j<=nY; j++)
		{
			nj++

			if (nj > nYnew)
			{
				nj = 1
				ni++
			}
			Vn[ni, nj] = V[i, j]
		}
	delete V
	for (i in Vn) 
		V[i] = Vn[i]
	nX = nXnew
	nY = nYnew
}

END{
	######################### REORDER/RESELECT VALUES
	if (OX != "")
		SelRows(OX)

	if (OY != "")
	{
		T()
		SelRows(OY)
		T()
	}

	######################### PERFORM OPERATIONS
	nOp = split(Ops, OpA, / /)

	# incorporate default value before Ops
	if (nOp > 0)
		for (i=1; i<=nX; i++)
			for (j=1; j<=nY; j++)
				if (!((i,j) in V))
					V[i,j] = D

	for (op=1; op<=nOp; op++)
		if (OpA[op] == "Transpose")
			T()
		else if (OpA[op] ~ /[XYBM]$/)
		{
			# if operation should have been done in columns
			# transpose before and after
			if (OpA[op] ~ /Y$/) T()
			if (OpA[op] ~ /M$/) {nXsaved = nX; ReShape(1)}

			if (OpA[op] ~ /^Scale01[XYM]$/ || OpA[op] ~ /^DivAMax[XYM]$/)
			{
				# scale rows to go from 0->1
				for (i=1; i<=nX; i++)
				{
					min = max = V[i,1]

					for (j=1; j<=nY; j++)
					{
						if (V[i,j] > max)
							max = V[i,j]
						if (V[i,j] < min)
							min = V[i,j]
					}
					
					if (OpA[op] ~ /^DivAMax[XYM]$/)
					{
						div = (max > -min) ? max : -min
						min = 0
					}
					else
						div = max - min

					for (j=1; j<=nY; j++)
						if (div != 0)
							V[i,j] = (V[i,j] - min) / div
						else
							delete V[i,j]
				}
			}
			else if (OpA[op] ~ /^DivSum[XYMB]$/)
			{
				# get the column sums if necessary
				if (OpA[op] ~ /^DivSumB$/)
					for (i=1; i<=nX; i++)
						for (j=1; j<=nY; j++)
							ColSum[j] += V[i,j]

				# scale rows dividing by the sum
				for (i=1; i<=nX; i++)
				{
					sum = 0

					for (j=1; j<=nY; j++)
						sum += V[i,j]
					
					for (j=1; j<=nY; j++)
						if (sum != 0)
							V[i,j] = V[i,j]/sum
						else
							delete V[i,j]
				}

				# apply the column sums if necessary
				if (OpA[op] ~ /^DivSumB$/)
				{
					for (i=1; i<=nX; i++)
						for (j=1; j<=nY; j++)
							if (ColSum[j] != 0)
								V[i,j] = V[i,j] / ColSum[j]
							else
								delete V[i,j]

					delete ColSum
				}
			}
			else if (OpA[op] ~ /^NormMStd[XYM]$/)
			{
				for (i=1; i<=nX; i++)
				{
					mean = std = sum = sum2 = num = 0

					for (j=1; j<=nY; j++)
						if (V[i,j] != "")
						{
							sum += V[i,j]
							sum2 += V[i,j] * V[i,j]
							num++
						}
					
					if (num > 1)
					{
						mean = sum / num

						if (sum2 > (sum*mean))
							std = sqrt((sum2 - sum * mean) / (num-1))
					}

					for (j=1; j<=nY; j++)
						if (std != 0)
							V[i,j] = (V[i,j] - mean)/std
						else
							delete V[i,j]
				}
			}
			else if (OpA[op] ~ /^SubMean[XYM]$/)
			{
				for (i=1; i<=nX; i++)
				{
					mean = sum = num = 0

					for (j=1; j<=nY; j++)
						if (V[i,j] != "")
						{
							sum += V[i,j]
							num++
						}
					
					if (num > 0)
						mean = sum / num

					for (j=1; j<=nY; j++)
						V[i,j] = (V[i,j] - mean)
				}
			}
			else if (OpA[op] ~ /^Quant[XYM]$/ || OpA[op] ~ /^Rank[XYM]$/)
			{
				SC = "sort -k1,1g" 
				for (i=1; i<=nX; i++)
					for (j=1; j<=nY; j++)
						print (0  + (V[i,j] != "" ? V[i, j] : D)) FS i FS j |& SC

				close(SC, "to")

				while (SC |& getline)
				{	
					V[$2, $3] = ++Rank[$2]
					
					NumSeen[Rank[$2]]++

					if ((nX%2 == 1) && ((nX+1) == (2*NumSeen[Rank[$2]])))
						MedVal[Rank[$2]] = $1
					else if ((nX%2 == 0) && ((nX == (2*NumSeen[Rank[$2]])) || (2+nX) == (2*NumSeen[Rank[$2]])))
						MedVal[Rank[$2]] += $1 / 2
				}

				close(SC)

				if (OpA[op] ~ /^Quant[XYM]$/)
					for (i=1; i<=nX; i++)
						for (j=1; j<=nY; j++)
							V[i,j] = MedVal[V[i,j]]

				delete NumSeen
				delete MedVal
				delete Rank
			}
			else if (OpA[op] ~ /^Min[XY]$/)
			{
				Y[++nY] = OpsPre "Min"
				for (i=1; i<=nX; i++)
					for (j=1; j<nY; j++)
						if (!((i,nY) in V) || (V[i,j]!="" &&  V[i,j] < V[i,nY]))
							V[i,nY] = V[i,j]
			} 
			else if (OpA[op] ~ /^Max[XY]$/)
			{
				Y[++nY] = OpsPre "Max"
				for (i=1; i<=nX; i++)
					for (j=1; j<nY; j++)
						if (!((i,nY) in V) || (V[i,j]!="" && V[i,j] > V[i,nY]))
							V[i,nY] = V[i,j]
			} 
			else if (OpA[op] ~ /^Sum[XY]$/)
			{
				Y[++nY] = OpsPre "Sum"
				for (i=1; i<=nX; i++)
					for (j=1; j<nY; j++)
						V[i,nY] += V[i,j]
			} 
			else if (OpA[op] ~ /^Stats[XY]$/)
			{
				for (i=1; i<=nX; i++)
				{
					sum = sum2 = 0

					for (j=1; j<=nY; j++)
						if (V[i,j] != "")
						{
							sum += V[i,j]
							sum2 += V[i,j] * V[i,j]
							V[i,nY+1]++
						}

					if (V[i,nY+1] > 0)
					{
						V[i, nY+2] = sum / V[i,nY+1]
						V[i, nY+3] = sqrt((sum2 - sum * V[i, nY+2]) / (V[i,nY+1]-1))
					}
				}

				Y[++nY] = OpsPre "N"
				Y[++nY] = OpsPre "Mean"
				Y[++nY] = OpsPre "STDev"
			} 
			else if (OpA[op] ~ /^Mean[XY]$/)
			{
				for (i=1; i<=nX; i++)
				{
					num = sum = 0

					for (j=1; j<=nY; j++)
						if (V[i,j] != "")
						{
							sum += V[i,j]
							num++
						}

					if (num > 0)
						V[i, nY+1] = sum / num
				}
				Y[++nY] = OpsPre "Mean"
			} 
			else if (OpA[op] ~ /^Flip[XY]$/)
			{
				new_x = ""
				for (i=1; i<=nX; i++)
					new_x = X[i] (i==1 ? "" : FS) new_x
				SelRows(new_x)
			}

			if (OpA[op] ~ /Y$/) T()
			if (OpA[op] ~ /M$/) ReShape(nXsaved)
		}



	######################### OUTPUT VALUES
	if (OF == 3)
	{
		for (i=1; i<=nX; i++)
			for (j=1; j<=nY; j++)
				if (OutM || ((i,j) in V))
					print X[i], Y[j], OutV(i,j)
	}
	else if (OF == "m")
	{
		if (OutYH)
		{
			L = (OutXH && OutN1) ? N1 : ""
			for (i=1; i<=nY; i++)
				L = L OFS Y[i]
			print substr(L, (OutXH && OutN1) ? 1 : (length(OFS)+1))
		}

		for (i=1; i<=nX; i++)
		{
			L = OutXH ? X[i] : ""
			for (j=1; j<=nY; j++)
				L = L OFS OutV(i,j)
			print substr(L, OutXH ? 1 : (length(OFS)+1))
		}
	}
	else if (OF == "mhtml")
	{
		print "<table border=\"" HtmlBorder "\">"

		if (OutYH)
		{
			print "<tr>"
			if (OutXH)
				print htmltd("<b>"N1"</b>")
			for (i=1; i<=nY; i++)
				print htmltd("<b>"Y[i]"</b>")
			print "</tr>"
		}

		for (i=1; i<=nX; i++)
		{
			print "<tr>"
			if (OutXH)
				print htmltd(X[i])
			for (j=1; j<=nY; j++)
				print htmltd(OutV(i,j))
			print "</tr>"
		}

		print "</table>"
	}
	else if (OF == "mtex")
	{
		print "\\begin{tabular}{" repeat("l", nY + (OutXH?1:0)) "}"
		if (OutYH)
		{
			if (OutXH)
				print N1
			for (i=1; i<=nY; i++)
			{
				if (OutXH || i>1)
					print "&"
				print Y[i]
			}
			print "\\\\"
		}

		for (i=1; i<=nX; i++)
		{
			if (OutXH)
				print X[i]
			for (j=1; j<=nY; j++)
			{
				if (OutXH || i>1)
					print "&"
				print OutV(i,j)
			}
			print "\\\\"
		}

		print "\\end{tabular}"
	}
}

function repeat(t, n,   l) 
{
	for (; n>0; n--)
		l = l t
	return l
}

function htmltd(x) 
{ 
	return "<td>" (x != "" ? x : "&nbsp;") "</td>" 
}

function OutV(i, j,   x)
{
	x = ((i,j) in V) ? V[i,j] : D

	if (OFmt ~ "%C")
		return gensub(/%C/, AC(x), "g", OFmt)

	return sprintf(OFmt, x)
}

# add commas to values
# also in Plot.awk
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

