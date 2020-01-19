#!/usr/bin/python
#
import math
import os
import psutil
#
###############
# PERSISTANCE #
###############
#
# flushBuffers - Transfer all data stored in the script buffers to disk, then free the memory in the buffers
#
# buffers - dictionary object with filenames as the keys, and a list of strings as the values
#
def flushBuffers(directory, buffers):

	for bufferKey in buffers:
		file = open(os.path.join(directory, bufferKey), "a+")
		file.write(''.join(buffers[bufferKey]))
		file.close()
	buffers.clear()
	
#
# createDirectory - Create a directory in the current working directory
#
# directory - Name of the directory to create
#
def createDirectory(directory):
	try:  
    		os.mkdir(directory)
	except OSError:  
    		print ("Directory already present [{}]".format(directory))
	else:  
  			print ("Successfully created the directory [{}]".format(directory))

#
# clearFilesMatchingFilter - Clean up any files that didn't meet filter criteria, having less reads than the minimum reads
#
# directory - Directory to clean up files from
# filter - Lambda predicate function to test whether or not to 
#
def clearFilesMatchingFilter(directory, filter):
	
	for (dirpath, dirnames, filenames) in os.walk(directory):
		for filename in filenames:
			if filter(filename):
				os.remove(os.path.join(directory, filename))
			
#
# countFiles - Count the number of files in a given directory matching a specified filter
#
# directory - Directory to count files in
# filter - Lambda predicate function to test whether or not to include the file in the count
#		
def countFilesMatchingFilter(directory, filter):

	count = 0
	for (dirpath, dirnames, filenames) in os.walk(directory):
		for filename in filenames:
			if filter(filename):
				count = count + 1
	return count
	
#
##########
# MEMORY #
##########
#
# analyzeMemoryUsage - Analyze the current memory usage of the script.  Return RSS memory used in bytes
#
def analyzeMemoryUsage():

	pid = os.getpid()
	py = psutil.Process(pid)

	rss = py.memory_info().rss				
	return rss	
#
# bytesToDisplay - Convert number of bytes to a human-readable string
#
# bytes - number of bytes
#
def bytesToDisplay(bytes):
	if bytes == 0:
		return "0B"
	units = ("B", "KB", "MB", "GB")
	i = int(math.floor(math.log(bytes, 1024)))
	p = math.pow(1024, i)
	s = round(bytes / p, 2)
	return "%s %s" % (s, units[i])	
#
# mbToBytes - Convert input in MB to bytes
#
# value - Value to convert, in megabytes
#				
def mbToBytes(value):
	mb = int(value)
	return mb * 1024 * 1024		

#
#####################
# DATA MANIPULATION #
#####################
#
# addToDictionaryList - Either append new value to an existing dictionary value list or initialize a new value list
#
# dictionary - Dictionary data structure to populate
# key - Key in the dictionary, will be created if it does not previously exist
# value - Value to add to the list at the given key, create a new list if key did not exist previously
#
def addToDictionaryList(dictionary, key, value):

	if key in dictionary:
		dictionary[key].append(value)
	else:
		dictionary[key] = []
		dictionary[key].append(value)
		
#
# addToDictionarySet - Either append new value to an existing dictionary value set or initialize a new value set
#
# dictionary - Dictionary data structure to populate
# key - Key in the dictionary, will be created if it does not previously exist
# value - Value to add to the set for the given key, create a new set if key did not exist previously
#
def addToDictionarySet(dictionary, key, value):

	if key in dictionary:
		dictionary[key].add(value)
	else:
		dictionary[key] = set()
		dictionary[key].add(value)

		