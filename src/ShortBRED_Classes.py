#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 19:30:03 2018

@author: jim
"""


class SBMarker:
    def __init__(self,strFamily,iCount,strType,strSeq,aiOL=None,dictOL=None):
        self.strFamily = strFamily
        self.iCount = iCount
        self.strType = strType
        self.strSeq = strSeq
        
        if aiOL is None:
            self.aiOL = []
        else:
            self.aiOL = aiOL
            
        if dictOL is None:
            self.dictOL = {}
        else:
            self.dictOL = dictOL
    

class SBFamily: 
    def __init__(self,strFamily,strConsensusSeq,aSBMarkers=None,
                 dictMemberSeqs=None):
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
            
    def __repr__(self):
        return "Name:" + self.strFamily + "\n" + str(self.strConsensusSeq)
            
