#!/usr/bin/python3

import  sys, os
from glob import glob

# mlsqlog-scan-sips.py
# Prints informative table of how chunks of mlsl2 job are progressing
# run at sips
# Based on Alyn Lambert's clever work

# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# "$Id$"

dirin  = sys.argv[1] + '/'
src = sys.argv[2]
l2cf =  sys.argv[3]
masterlog = sys.argv[4]

# version says whether we can expect to find 'Chunk number' in the chunk's own
# log file (version 2 or greater)
# version = 1.52
version = sys.argv[5]

ChunkDict = {}
TidDict = {}
PhaseDict = {'Completed':('Completed',0,0)}
StatusDict = {'Completed':'XXXX', 'Total':'XXXX', 'Underway':'XXXX', 'Abandoned':'XXXX', 'Remaining':'XXXX'}


lines = []
with open(dirin + '*.log', 'r') as f:
    for line in f:
        if 'died, try again' in line:
            lines.append(line)
for line in lines:
    line = line[:-1]
    chunkno = line.split("run of chunk ")[1].split(" ")[0]
    ChunkDict[chunkno] = 'Died'

lines = []
#print "Phases"
with open(l2cf, 'r') as f:
    for lineNumber, line in enumerate(f, start=1):
        if '= Phase:' in line:
            lines.append(lineNumber + ': ' + line)
for line in lines:
    line = line[:-1]
#    print (line.split(':'))[0], (line.split(':'))[2]
    PhaseDict[(line.split(':'))[0]] = ((line.split(':'))[2], 0, 0)

#print PhaseDict

# search master log
f = open(masterlog)
lines = f.readlines()

w = ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']
for line in lines:
    line = line[:-1]
    if line.find('chunks completed') != -1:
        if line.find(':') != -1:
            v = ((line.split(':'))[1]).split()
        else :
            v = line.split()
        # print v, len(v)
        w[0:len(v)-1] = v
        StatusDict = {'Completed':w[0], 'Total':w[2], 'Underway':w[5], 'Abandoned':w[7], 'Remaining':w[9]}

# print masterlog
# print len(lines)
c = ['0', '0', '0', '0', '0']
for line in lines:
    line = line[:-1]
    if line.find('Launched chunk') != -1 :
        b = line.split('Launched chunk')[1].split()
        c[0:len(b)] = b
        chunkno = c[0]
        tidno = c[4]
        TidDict[tidno] = chunkno
        # print c, tidno, chunkno

f.close()
        

#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#   #start search in log file directory
fileList = glob(dirin + '/*/*.log')
#   
#   nlogfiles = len(os.popen(cmd).readlines())
#   
#   if nlogfiles == 1:
#       file = (os.popen(cmd).readlines())[0]
#       filein = file[:-1] 
#   
#       print filein
#   
#       f = open(filein)
#       lines = f.readlines()
#       
#        for line in lines:
#           line = line[:-1]      
#       
#           if line.find('Slave task') != -1 :
#    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx           


# print 'xxxxxxxxxxxxx'

for file in fileList:
    filein = file[:-1] 
#    print filein
    f = open(filein)
    chunkno="-999"
    time=-999
    lno="-999"
    lines = f.readlines()
    lastline = len(lines)
    completed = False
    
    for line in lines:
        line = line[:-1]      
                 
        if version > 1.9 :
            if line.find('Chunk number') != -1 :
                chunkno = (((line.split('Chunk number'))[1]).split(' '))[0]
        else :
            if line.find('Task ID:') != -1 :
                tidno = (((line.split('Task ID:'))[1]).split())[0]
                if tidno in TidDict :
                    chunkno = TidDict[tidno]
                # print line, (((line.split('Task ID:'))[1]).split())[0], chunkno

        if line.find('.Enter') != -1 and line.find('line') and line.find('column') != -1 :
            lno = (((line.split('line'))[1].split('column'))[0]).strip()
            # In case the lno has a trailing comma
            lno = lno.split(',')[0]

        if line.find('ending mlsl2') != -1 :
            completed = True

        died = ChunkDict.has_key(chunkno)
            
    f.close()
    # print chunkno
    for k,v in PhaseDict.iteritems() :
        # print k, lno, died, completed, lno.isdigit()
        if k == 'Completed' :
            if completed :  PhaseDict[k] = (PhaseDict[k][0],PhaseDict[k][1]+1,PhaseDict[k][2])
            if died      :  PhaseDict[k] = (PhaseDict[k][0],PhaseDict[k][1],PhaseDict[k][2]+1)
        elif lno.isdigit() :
            if int(lno) >= int(k) :
                if died :
                    PhaseDict[k] = (PhaseDict[k][0],PhaseDict[k][1]+1,PhaseDict[k][2]+1)
                else :
                    PhaseDict[k] = (PhaseDict[k][0],PhaseDict[k][1]+1,PhaseDict[k][2])
                    
print
print( 'Total'.rjust(10),    StatusDict['Total'].rjust(10))
print('Underway'.rjust(10),  StatusDict['Underway'].rjust(10))
print('Remaining'.rjust(10), StatusDict['Remaining'].rjust(10))
print('Abandoned'.rjust(10), StatusDict['Abandoned'].rjust(10))
print('Completed'.rjust(10), StatusDict['Completed'].rjust(10))
print(" ")

n=-1
PhaseStartLineNos=PhaseDict.keys()
PhaseStartLineNos.sort()
# print PhaseStartLineNos
print('LineNo'.rjust(10), 'Phase Name'.rjust(20), '#Chunks'.rjust(10))
for i in PhaseStartLineNos :
    n = n + 1
    print('%s %s %10d %10d' % (i.rjust(10), PhaseDict[i][0].rjust(20), PhaseDict[i][1], PhaseDict[i][2]), end=" ")
    if n != len(PhaseStartLineNos)-1 :
        print('%10d %10d' % (PhaseDict[PhaseStartLineNos[n]][1] - PhaseDict[PhaseStartLineNos[n+1]][1], PhaseDict[PhaseStartLineNos[n]][2] - PhaseDict[PhaseStartLineNos[n+1]][2]))
    else :
        print(" ")
                                                          
print

PhaseStartLineNos=PhaseDict.keys()
PhaseStartLineNos.sort()
print('LineNo'.rjust(10), 'Phase Name'.rjust(20), '#Chunks'.rjust(10))
for i in PhaseStartLineNos :
    print('%s %s %10d %10d' % (i.rjust(10), PhaseDict[i][0].rjust(20), PhaseDict[i][1], PhaseDict[i][2]))

    
print("LineNumbers")
for file in fileList:
    filein = file[:-1] 
    
    f = open(filein)
    chunkno=-999
    time=-999
    lno="-999"
    PhaseName=-999
    lines = f.readlines()
    lastline = len(lines)
    for line in lines:
        line = line[:-1]
        if version > 1.9 :
            if line.find('Chunk number') != -1 :
                chunkno = (((line.split('Chunk number '))[1]).split(' '))[0]
        else :
            if line.find('Task ID:') != -1 :
                tidno = (((line.split('Task ID:'))[1]).split())[0]
                if tidno in TidDict :
                    chunkno = TidDict[tidno]
        
        if line.find('.Enter') != -1 and line.find('line') != -1 :
            time = (((line.split('at '))[1]).split(' '))[0]
            lno =(((line.split('line'))[1].split('column'))[0]).strip()
            #print chunkno, time, lno
            # In case the lno has a trailing comma
            lno = lno.split(',')[0]

        if lno.isdigit() : 
            if (line.lower()).find(src.lower()) !=-1 :
                for i in PhaseStartLineNos[0:-1] :
                    if int(lno) >= int(i) : PhaseName = PhaseDict[i][0]
                    print('>>>>>>>>>>>', PhaseName, chunkno, time, lno)
                    print(line)
    f.close()
    

# $Log$
# Revision 1.1  2006/10/19 18:29:10  pwagner
# First commit
#
