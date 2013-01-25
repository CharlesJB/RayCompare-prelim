#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-01-23

"""
This script compares multiple assemblies by k-merizing and coloring them, simulating
what will eventually be done with greater efficiency with RayCompare plugin.

Usage:
./AssemblyCompare.py <inputFile> <kmer-length>
	inputFile: A file containing the path to every file to compare. One path per line.
		   Files must be in fasta format.
	kmer-length: Size of k-mers to produce.
"""

class KmerParser:
	def __init__(self, inputFile, kmerLength):
		self.m_fileList = []
		self._parseInputFile(inputFile)
		self.m_kmerLength = kmerLength
		self.m_kmerList = {}

	def parseFiles(self):
		for filename in self.m_fileList:
			self._parseFile(filename)

	def printSummary(self):
		summary = {}
		# Compute count for every virtual color
		for key in self.m_kmerList:
			color = self.m_kmerList[key]
			if color not in summary:
				summary[color] = 1
			else:
				summary[color] = summary[color] + 1
		# Calculate total kmer count
		totalKmerCount = sum(summary.values())
		# Print summary file
		print "Filenames\tValues"
		for filename in self.m_fileList:
			print filename + "\t" + self._printColor(self._convertFilenameToColorValue(filename))
		print ""
		print "VirtualColor\tCount\tPercentage"
		for key in summary:
			percentage = 0
			if totalKmerCount != 0:
				percentage = (float(summary[key]) / totalKmerCount) * 100.0
			print self._printColor(key) + "\t" + str(summary[key]) + "\t" + str(percentage)
		print ""
		print "Unique k-mers total count\t" + str(totalKmerCount)

	def _parseFile(self, filename):
		currentFilenameColorValue = self._convertFilenameToColorValue(filename)
		entrySequence = ""
		nucleotideCount = 0
		for line in open(filename):
			if line[0] == '>':
				self._kmerizeSequence(entrySequence.lower(), currentFilenameColorValue)
				entrySequence = ""
			else:
				entrySequence = entrySequence + line.strip()
				nucleotideCount += len(line.strip())
				if nucleotideCount > 1000: # To avoid using too much memory when parsing a sequence, we kmerize it once in a while
					self._kmerizeSequence(entrySequence.lower(), currentFilenameColorValue)
					# We keep what we know for now about the next k-mer
					entrySequence = entrySequence[len(entrySequence)-self.m_kmerLength+1]
					nucleotideCount = len(entrySequence)
		self._kmerizeSequence(entrySequence.lower(), currentFilenameColorValue)

	def _kmerizeSequence(self, sequence, currentFilenameColorValue):
		if len(sequence) >= self.m_kmerLength:
			start = 0
			end = start + self.m_kmerLength - 1
			while end < len(sequence):
				self._addKmer(sequence[start:end], currentFilenameColorValue)
				start += 1
				end += 1

	def _addKmer(self, kmer, colorValue):
		if kmer in self.m_kmerList:
			lastColorValue = self.m_kmerList[kmer]
			self.m_kmerList[kmer] = self._updateVirtualColor(colorValue, lastColorValue)
		else:
			self.m_kmerList[kmer] = colorValue

	def _parseInputFile(self, inputFile):
		for line in open(inputFile):
			self.m_fileList.append(line.strip())
	
	def _convertFilenameToColorValue(self, filename):
		filename_index = self.m_fileList.index(filename)
		return 10**filename_index

	def _updateVirtualColor(self, filenameColorValue, lastColorValue = 0):
		newValue = lastColorValue
		if ((lastColorValue / filenameColorValue) % 2) == 0:
			newValue = lastColorValue + filenameColorValue
		return newValue

	def _printColor(self, colorValue):
		toPrint = ""
		fileCount = len(self.m_fileList)
		for i in range(1, fileCount + 1):
			if (colorValue / (10**(fileCount - i)) % 2) == 0:
				toPrint = toPrint + "0"
			else:
				toPrint = toPrint + "1"
			if i != (fileCount):
				toPrint = toPrint + "-"
		return toPrint

import sys

if __name__=="__main__":
        if len(sys.argv)!=3:
                print __doc__
                sys.exit(1)

	inputFile = sys.argv[1]
	kmerLength = int(sys.argv[2])
	kmerParser = KmerParser(inputFile, kmerLength)
	kmerParser.parseFiles()
	kmerParser.printSummary()
