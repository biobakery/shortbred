#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 19:30:03 2018

@author: jim
"""


class SBMarker:
    def __init__(self,strFamily,iCount,strType,strSeq):
        self.strFamily = strFamily
        self.iCount = iCount
        self.strType = strType
        self.strSeq = strSeq
        
    def __repr__(self):
        return "Family:" + self.strFamily + "\nNumber:" + str(self.iCount) + \
    "\nType:" + str(self.strType) + "\nSequence:" + str(self.strSeq) 
    

class SBFamily: 
    def __init__(self,strFamily,strConsensusSeq,aSBMarkers=None,
                 dictMemberSeqs=None,aiOverlapGOI=None,atupOverlapGOI=None,
                 aiOverlapRef=None,atupOverlapRef=None,aiOverlapTotal=None):
        self.strFamily = strFamily
        self.strConsensusSeq = strConsensusSeq
        
        if aSBMarkers is None:
            self.aSBMarkers = []
        else:
            self.aSBMarkers = list(aSBMarkers)
            
        if dictMemberSeqs is None:
            self.dictMemberSeqs = {}
        else:
            self.dictMemberSeqs = dictMemberSeqs
 
        if aiOverlapGOI is None:
            self.aiOverlapGOI = [0]*len(self.strConsensusSeq)
        else:
            self.aiOverlapGOI = aiOverlapGOI
            
        if atupOverlapGOI is None:
            self.atupOverlapGOI = {}
        else:
            self.atupOverlapGOI = atupOverlapGOI

 
        if aiOverlapRef is None:
            self.aiOverlapRef = [0]*len(self.strConsensusSeq)
        else:
            self.aiOverlapRef = aiOverlapRef
            
        if atupOverlapRef is None:
            self.atupOverlapRef = {}
        else:
            self.atupOverlapRef = atupOverlapRef
            
        if aiOverlapTotal is None:
            self.aiOverlapTotal = [iGOI + iRef for iGOI,iRef in zip(self.aiOverlapGOI,self.aiOverlapRef)]
        else:
            self.aiOverlapTotal = aiOverlapTotal
           
        self.iTM= 0
        self.iJM= 0
        self.iQM = 0
 
    def __repr__(self):
        return "Name:" + self.strFamily + "\nConsensus Sequence:" + str(self.strConsensusSeq) + \
    "\nGOI_Overlap:" + str(self.aiOverlapGOI) + "\nRef_Overlap:" + str(self.aiOverlapRef) + \
    "\nTotal_Overlap:" + str(self.aiOverlapTotal) + "\nTrue Markers:" + str(self.iTM) + \
    "\nJunction Markers:" + str(self.iJM) + "\nQuasi Markers:" + str(self.iQM)
            
