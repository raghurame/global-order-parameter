import lzma as l
import time as t
import argparse as a
import decimal as d
import numpy as n
import math as m
import os
from tqdm import tqdm

def extract_numbers (line):
	lineString = line.replace ("\t", ",").replace (" ", ",").replace ("\n", "")
	for item in lineString.split (","):
		try:
			yield d.Decimal (item)
		except:
			pass

def unwrapCoordinates (rawCoords, image, loBounds, hiBounds):
	if (image > 0):
		return hiBounds + ((abs (image) - 1) * (hiBounds - loBounds)) + (rawCoords - loBounds)
	elif (image < 0):
		return loBounds + ((abs (image) - 1) * (hiBounds - loBounds)) + (hiBounds - rawCoords)
	else:
		return rawCoords

def zeroListMaker(n):
	listOfZeroes = [0] * n
	return listOfZeroes

def computeOrderParameter (mainEntry, xlow, xhigh, ylow, yhigh, zlow, zhigh, cxlow, cxhigh, cylow, cyhigh, czlow, czhigh):

	orderParameter_array = []

	for x in range (1, len (mainEntry) - 9, 1):

		try:
			xCOM = (mainEntry[x][3] + mainEntry[x + 9][3]) / 2
			yCOM = (mainEntry[x][4] + mainEntry[x + 9][4]) / 2
			zCOM = (mainEntry[x][5] + mainEntry[x + 9][5]) / 2

			if (xCOM < cxhigh and xCOM > cxlow and yCOM > cylow and yCOM < cyhigh and zCOM < czhigh and zCOM > czlow):
				
				x1 = unwrapCoordinates (float (mainEntry[x][3]), float (mainEntry[x][9]), xlow, xhigh)
				y1 = unwrapCoordinates (float (mainEntry[x][4]), float (mainEntry[x][9]), ylow, yhigh)
				z1 = unwrapCoordinates (float (mainEntry[x][5]), float (mainEntry[x][9]), zlow, zhigh)
				x2 = unwrapCoordinates (float (mainEntry[x + 9][3]), float (mainEntry[x + 9][9]), xlow, xhigh)
				y2 = unwrapCoordinates (float (mainEntry[x + 9][4]), float (mainEntry[x + 9][9]), ylow, yhigh)
				z2 = unwrapCoordinates (float (mainEntry[x + 9][5]), float (mainEntry[x + 9][9]), zlow, zhigh)

				dotProduct = float ((x2 - x1) * 0 + (y2 - y1) * 0 + (z2 - z1) * 1)
				mag1 = float (m.sqrt (pow ((x2 - x1), 2) + pow ((y2 - y1), 2) + pow ((z2 - z1), 2)))
				mag2 = float (1)
				cosTheta = dotProduct / (mag1 * mag2)
				theta = m.acos (dotProduct / (mag1 * mag2))
				orderParameter = ((3 * cosTheta * cosTheta) - 1) / 2
				orderParameter_array.append (orderParameter)
		except:
			pass

	return orderParameter_array

def findPopDist (orderParameter, orderParameter_popDist, binwidth):

	for value in orderParameter:
		for x in range (0, int (1.5/binwidth)):
			binLow = (x * binwidth) - 0.5
			binHigh = binLow + binwidth
			if (value < binHigh and value > binLow):
				orderParameter_popDist[x] += 1

	return orderParameter_popDist

def readXZ (inputFile, dirs, cxlow, cxhigh, cylow, cyhigh, czlow, czhigh, binwidth, finalFrame, outputDist, outputAvg):
	lineNumber = 0
	entryDict = {}
	mainEntry = []
	nAtoms = []
	orderParameter = []
	orderParameter_timeData = []
	orderParameter_popDist = zeroListMaker(int (1.5/binwidth))
	nTimeframe = 0
	isDimX = 1
	isDimY = 1
	isDimZ = 1
	
	with l.open (inputFile, mode = "rt") as file:
		for line in tqdm (file):
			lineNumber += 1

			if (lineNumber == 4):
				nAtoms = list (extract_numbers (line))
				nTimeframe += 1

			if (nTimeframe == finalFrame):
				break

			if (isDimX and lineNumber == 6):
				print ("Reading boundary...")
				xlow, xhigh = list (extract_numbers (line))
				isDimX = 0
				print ("xlow: {}\txhigh: {}".format (xlow, xhigh))
			if (isDimY and lineNumber == 7):
				ylow, yhigh = list (extract_numbers (line))
				isDimY = 0
				print ("ylow: {}\tyhigh: {}".format (ylow, yhigh))
			if (isDimZ and lineNumber == 8):
				zlow, zhigh = list (extract_numbers (line))
				isDimZ = 0
				print ("zlow: {}\tzhigh: {}".format (zlow, zhigh))

			if (lineNumber > 9 and lineNumber <= nAtoms[0] + 9):
				mainEntry.append (list (extract_numbers (line)))

				if (lineNumber == nAtoms[0] + 9):
					orderParameter = computeOrderParameter (mainEntry, xlow, xhigh, ylow, yhigh, zlow, zhigh, cxlow, cxhigh, cylow, cyhigh, czlow, czhigh)
					orderParameter_timeData.append (n.average (orderParameter))
					orderParameter_popDist = findPopDist (orderParameter, orderParameter_popDist, binwidth)

					lineNumber = 0
					mainEntry = []


	with open (os.path.join (dirs, outputDist), "w") as file:
		file.write (str(orderParameter_popDist))
	with open (os.path.join (dirs, outputAvg), "w") as file:
		file.write ("Average: {}\nStandard Deviation: {}".format (n.average (orderParameter_timeData), n.std (orderParameter_timeData)))

	print ("\n\nOrder parameter frequency")
	print (orderParameter_popDist)
	print ("\nAverage: {}\nStandard Deviation: {}\n".format (n.average (orderParameter_timeData), n.std (orderParameter_timeData)))

if __name__ == '__main__':
	parser = a.ArgumentParser (description = "Compute global order parameter from input LAMMPS dump file in XZ format.")
	parser.add_argument ("--input", "-i", required = True, type = str, help = "Input LAMMPS dump file name in compressed XZ format [required arg]")
	parser.add_argument ("--cxlow", "-cxl", default = -99999999, type = float, help = "Lower bounds in X axis for analysis [optional arg]")
	parser.add_argument ("--cxhigh", "-cxh", default = 99999999, type = float, help = "Upper bounds in X axis for analysis [optional arg]")
	parser.add_argument ("--cylow", "-cyl", default = -99999999, type = float, help = "Lower bounds in Y axis for analysis [optional arg]")
	parser.add_argument ("--cyhigh", "-cyh", default = 99999999, type = float, help = "Upper bounds in Y axis for analysis [optional arg]")
	parser.add_argument ("--czlow", "-czl", default = -99999999, type = float, help = "Lower bounds in Z axis for analysis [optional arg]")
	parser.add_argument ("--czhigh", "-czh", default = 99999999, type = float, help = "Upper bounds in Z axis for analysis [optional arg]")
	parser.add_argument ("--binwidth", "-bw", default = 0.1, type = float, help = "Bin width for calculating order parameter frequency [optional arg; default value = 0.1]")
	parser.add_argument ("--finalframe", "-ff", default = -1, type = int, help = "Timeframe to trigger a shutdown [optional arg; default value = last timeframe]")
	parser.add_argument ("--outputDist", "-od", default = "globalOrderParameter.dist", type = str, help = "Timeframe to trigger a shutdown [optional arg; default value = last timeframe]")
	parser.add_argument ("--outputAvg", "-oa", default = "globalOrderParameter.average", type = str, help = "Timeframe to trigger a shutdown [optional arg; default value = last timeframe]")
	args = parser.parse_args()

	for dirs, dirname, files in os.walk("."):
		for file in files:
			if args.input in file:
				filePath = dirs + "/" + file
				print (filePath)
				readXZ (filePath, dirs, args.cxlow, args.cxhigh, args.cylow, args.cyhigh, args.czlow, args.czhigh, args.binwidth, args.finalframe, args.outputDist, args.outputAvg)
