from lofar.pythondppp import DPStep
from lofar.parameterset import parameterset
import numpy as np
import math
from stefcal import stefcal
from time import time as systime

class DPStefcal(DPStep):
    def __init__(self, parsetDict):
        # The constructor gets the subset of the NDPPP parset containing
        # all keys-value pairs for this step.
        # Note: the superclass constructor MUST be called.
        DPStep.__init__(self, parsetDict)
        parset = parameterset(parsetDict)

        self.itsTimeFill = 0.
        self.itsTimeFlag = 0.
        self.itsTimeReorder = 0.
        self.itsTimeSolve = 0.

        self.itsSols=[]
        self.itsTimeSlot=0
        self.itsMinBlPerAnt=parset.getInt('minBlPerAnt',4)
        self.itsSolInt=parset.getInt('solint',1)
        if self.itsSolInt>1:
          raise("SolInt>1 is not yet supported")
        self.itsNChan=parset.getInt('nchan',0)


    def updateInfo(self, dpinfo):
        # This function must be implemented.
        self.itsInfo = dpinfo
        # Make the arrays that will get the input buffer data from
        # the getData, etc. calls in the process function.
        self.itsData = self.makeArrayDataIn()
        self.itsModelData = self.makeArrayDataIn()
        self.itsFlags = self.makeArrayFlagsIn()
        self.itsWeights = self.makeArrayWeightsIn()
        self.itsUVW = self.makeArrayUVWIn()

        self.itsAntUsedNames=[self.itsInfo['AntNames'][i] for i in self.itsInfo['AntUsed']]
#        ['AntNames', 'RefFreq', 'AutoCorrIndex', 'MSName', 'NTime', 'OrigNChan', 'WriteData', 'NeedVisData', 'ChanFreqs', 'Ant1', 'Ant2', 'ChanAvg', 'TimeAvg', 'AntDiam', 'ChanWidths', 'TimeInterval', 'AntennaSet', 'MetaChanged', 'AntUsed', 'BLength', 'WriteFlags', 'EffectiveBW', 'NCorr', 'AntMap', 'NChan', 'WriteWeights', 'StartTime', 'TotalBW', 'Resolutions', 'StartChan']

        nStUsed=len(self.itsAntUsedNames)

        if self.itsNChan==0:
          self.itsNChan=self.itsInfo['NChan']

        if self.itsNChan%self.itsInfo['NChan']!=0:
          raise("nchan needs to divide total number of channels")

        nCh=self.itsInfo['NChan']/self.itsNChan # number of 'channel' solutions
        nDb=self.itsSolInt*self.itsNChan # number of data instances
        
        self.itsDataCube=np.zeros((nStUsed,2,nDb,nCh,2,nStUsed),
                                  dtype=complex,order='F')
        self.itsModelDataCube=np.zeros((nStUsed,2,nDb,nCh,2,nStUsed),
                                       dtype=complex,order='F')

        # Return the dict with info fields that change in this step.
        return {}

    def findFlaggedStations(self, flags):
        stCount={}
        flaggedStations=[]

        for st in self.itsInfo['AntUsed']:
          stCount[st]=0

        for bl,(st1,st2) in enumerate(zip(self.itsInfo['Ant1'], self.itsInfo['Ant2'])):
          stCount[st1]+=flags[bl,:,:].size-np.sum(flags[bl,:,:])  # adds number of Falses (unflagged datapoints)
          stCount[st2]+=flags[bl,:,:].size-np.sum(flags[bl,:,:])

        for st in stCount:
          if stCount[st]<self.itsMinBlPerAnt:
            flaggedStations.append(st)

        return flaggedStations

    def fillDataCubes(self, data, modeldata, flags, weights):
        self.itsDataCube*=0
        self.itsModelDataCube*=0
        nSt=len(self.itsAntUsedNames)

        for bl,(st1,st2) in enumerate(zip(self.itsInfo['Ant1'], self.itsInfo['Ant2'])):
          for ch in range(data.shape[1]):
            chset=ch/self.itsNChan
            chrest=ch%self.itsNChan
            for pol in range(data.shape[2]):
              if not flags[bl,ch,pol]:
                s1=self.itsInfo['AntMap'][st1]
                s2=self.itsInfo['AntMap'][st2]
                pol1=pol%2
                pol2=pol/2
                self.itsDataCube[s1,pol1,]=math.sqrt(weights[bl,ch,pol])*data[bl,ch,pol]
                self.itsModelDataCube[s1,pol1,]=math.sqrt(weights[bl,ch,pol])*modeldata[bl,ch,pol]
 
    def process(self, time, exposure):
        # This function must be implemented.
        # First get the data arrays needed by this step.
        time0=systime()
        self.getData (self.itsData);
        self.getModelData (self.itsModelData);
        self.getFlags (self.itsFlags);
        self.getWeights (self.itsWeights);
        self.getUVW (self.itsUVW);
        self.itsTimeFill += systime() - time0
        # Process the data.
        #print "process tPythonStep", time-4.47203e9, exposure, self.itsData.sum(), self.itsFlags.sum(), self.itsWeights.sum(), self.itsUVW.sum()
        #print self.itsData.shape
        # Execute the next step in the DPPP pipeline. TIME,UVW are changed.

        time0=systime()
        flaggedStations=self.findFlaggedStations(self.itsFlags)
        self.itsTimeFlag += systime() - time0

        if len(flaggedStations)>0:
          raise Exception('Yikes, cannot handle totally flagged station, which is the case for station '+str(flaggedStations[0])+' ('+self.itsInfo['AntNames'][flaggedStations[0]]+')')

        time0=systime()
        self.fillDataCubes(self.itsData, self.itsModelData, self.itsFlags, self.itsWeights)
        self.itsTimeReorder += systime() - time0
        self.itsTimeSlot+=1

        time0=systime()
        self.itsSols.append(stefcal(self.itsDataCube, self.itsModelDataCube))
        self.itsTimeSolve += systime() - time0

        return self.processNext ({'DATA': self.itsData})

    def finish(self):
        # Finish the step as needed.
        # This function does not need to be implemented.
        # Note: finish of the next step is called by the C++ layer.
        return
        #print self.itsSols

    def showCounts(self):
        # Show the counts of this test.
        # This function does not need to be implemented.
        timestatus = "    Time spent in getting data:    "+str(self.itsTimeFill)+'\n'
        timestatus+= "    Time spent in checking flags:  "+str(self.itsTimeFlag)+'\n'
        timestatus+= "    Time spent in reordering data: "+str(self.itsTimeReorder)+'\n'
        timestatus+= "    Time spent in actual stefcal:  "+str(self.itsTimeSolve)
        return timestatus
