#!/bin/bash
import rat
import numpy as np 
import time
from math import exp, pi 
import os.path
from os import path
#import matplotlib.pyplot as plt 
#from mpl_toolkits import mplot3d
#from matplotlib.pyplot import colormaps as cmaps
import ROOT
"""
Script to implement the MEGACUBE (TM) multi-site analysis method. Fitted vertex forms the centre of a 
cube, of side length 1m. This MEGACUBE is filled with minicube cubes, side length = 1/2 pos. resolution. 
For each mini-cube, a line is drawn to every hit PMT. ET1D PDF is used to obtain probability light originated 
from that mini-cube. The probability is stored within the mini-cube, and an image of the 3D probability 
field is created. Funky visuals. 0vbb and Co-60 (single vs multi-site) events will look different, and an ML 
should be able to discriminate between these events with ease.  
"""

class Minicube(object):
    def __init__(self, miniSidel, calPMTS, eventTime, times, probs):
        self.probability = 0
        self.position = np.zeros(3)    # center of the Minicube
        self.sidel = miniSidel         # side length of Minicube   
        self.calPMTS = calPMTS         # all the triggered PMTs 
        self.eventTime = eventTime     # global fitted time 
        self.times = times 
        self.probs = probs 

    def obtain_probability(self):
        """
        Function draws line to each hit PMT from center of minicube. Time of 
        flight is obtained. Emission time is obtained via: 

        T_em = T_trig - T_flight - T_event  

        ET1D "effective time 1D" PDF is used to return probability photon 
        originated in that minicube. Photon attenuation length and solid angle to 
        each PMT is also calculated. Probability is updated. 
        """           

        light_path = rat.utility().GetLightPathCalculator()
        group_velocity = rat.utility().GetGroupVelocity()
        pmt_info = rat.utility().GetPMTInfo()

        # have to convert np to TVector for ROOT/RAT functions 
        position = ROOT.std.TVector3(self.position[0], self.position[1], self.position[2])
        #print("There are {} triggered PMTs.".format(self.calPMTS.GetCount()))
        #start = time.time()
        for ipmt in range(self.calPMTS.GetCount()):
            # obtain PMT position 
            pmt_cal = self.calPMTS.GetPMT(ipmt)
            pmtPos = pmt_info.GetPosition(pmt_cal.GetID())

            # find lightpath from minicube to triggered pmt 
            light_path.CalcByPosition(position, pmtPos)
            inner_av_distance = light_path.GetDistInInnerAV()
            av_distance = light_path.GetDistInAV()
            water_distance = light_path.GetDistInWater()
            totalDistance = inner_av_distance + av_distance + water_distance
            #print("inner distance is: ", inner_av_distance)
            # obtain time of flight from minicube to triggered pmt 
            transit_time = group_velocity.CalcByDistance(inner_av_distance, av_distance, water_distance)

            # obtain emission time 
            tEm = np.array([pmt_cal.GetTime() - transit_time - eventTime])

            # obtain probability from ET1D PDF 
            binIdx = np.digitize(tEm, bins = self.times)

            # probability incorporates photon attenuation, ET1D PDF and solid angle 
            self.probability += self.probs[binIdx] * exp(-totalDistance/40000) * (1/totalDistance**2) #((72361*pi/40) / totalDistance**2) 
        
        #stop = time.time()
        #print("pmt loop time: ", stop - start)
        
class Megacube(object):
    def __init__(self, eventPosition, megaSidel, numMiniCubes, calPMTs, eventTime, path, count, type):
        self.numMiniCubes   = numMiniCubes              # num minicubes on each axis
        self.eventPosition  = eventPosition             # center of the Megacube (fit vertex) 
        self.eventTime      = eventTime                 # vertex fitted time of event
        self.megaSidel      = megaSidel                 # side length of Megacube defining volume around fit vertex  
        self.miniSidel      = megaSidel / numMiniCubes  # side length of minicube 
        self.calPMTS        = calPMTs                   # all PMTs triggered by event
        self.eventTime      = eventTime
        self.probs          = None                      # ET1D PDF probabilities (y axis)
        self.times          = None                      # ET1D PDF times (x axis)
        self.savePath       = path                      # path to save output arrays 
        self.count          = count                     # ID of specific Megacube used in save function 
        self.type           = type                      # type of event creating cube for saving

        # load the ET1D PDF 
        s1 = time.time()
        self.load_ET1D()
        e1 = time.time()
        print("Loading PDF took: ", e1-s1)
        
        # create all Minicube objects contained by Megacube 
        s2 = time.time()
        self.miniCubes = np.array([Minicube(self.miniSidel, self.calPMTS, self.eventTime,
                         self.times, self.probs) for i in range(numMiniCubes**3)]) 
        e2 = time.time()
        print("Creating mini-cubes took: ", e2-s2)
        
        # Megacube side length / 2 either side of center, split into minicubes
        min_x = self.eventPosition[0] - self.megaSidel/2
        max_x = self.eventPosition[0] + self.megaSidel/2
        min_y = self.eventPosition[1] - self.megaSidel/2
        max_y = self.eventPosition[1] + self.megaSidel/2
        min_z = self.eventPosition[2] - self.megaSidel/2
        max_z = self.eventPosition[2] + self.megaSidel/2

        xAxis = np.linspace(min_x, max_x, self.numMiniCubes)
        yAxis = np.linspace(min_y, max_y, self.numMiniCubes)
        zAxis = np.linspace(min_z, max_z, self.numMiniCubes)

        # meshgrid this 
        xx, yy, zz = np.meshgrid(xAxis, yAxis, zAxis)

        # flatten meshed arrays 
        xx = np.ravel(xx)
        yy = np.ravel(yy)
        zz = np.ravel(zz)
        done = 0 

        s3 = time.time()
        for cube in self.miniCubes:
            # populate position attribute for every Minicube object 
            cube.position = np.array([xx[done], yy[done], zz[done]])
        
            # populate probability attribute of every Minicube  
            #self.fill_probabilities(cube)
            cube.obtain_probability()
            
            done += 1 
        e3 = time.time()
        print("Filling probability took: ", e3-s3)

        # normalise the probability 
        s4 = time.time()
        total_prob = sum([cube.probability for cube in self.miniCubes])
        for cube in self.miniCubes:
            cube.probability /= total_prob
        e4 = time.time()
        print("Normalising took: ", e4-s4)

        s5 = time.time()
        # save the Megacube 
        self.save_cube()
        e5 = time.time()
        print("Saving cube took: ", e5-s5)
    
    # def fill_probabilities(self, cube):
    #     """
    #     For each Minicube, fill with probability photons originated from its position.
    #     """
        
    #     #done = 0 
    #     #for cube in self.miniCubes:
    #     #    cube.obtain_probability()
    #     #    done += 1
    #     #    print("done: {}/{}".format(done, self.numMiniCubes**3)) 

    #     cube.obtain_probability()

    def load_ET1D(self):
        """
        Loads the probability and timing information and saves 
        them as global arrays for use by the minicubes. 
        """

        self.times = np.load("/home/hunt-stokes/times_ET1D.npy")
        self.probs = np.load("/home/hunt-stokes/probs_ET1D.npy")

    def save_cube(self):
        """
        Plots the mini cube positions.
        """

        xx = np.zeros(self.numMiniCubes**3)
        yy = np.zeros(self.numMiniCubes**3)
        zz = np.zeros(self.numMiniCubes**3)
        probs = np.zeros(self.numMiniCubes**3)
        for cubeIdx in range(np.size(self.miniCubes)):
            xx[cubeIdx] = self.miniCubes[cubeIdx].position[0]
            yy[cubeIdx] = self.miniCubes[cubeIdx].position[1]
            zz[cubeIdx] = self.miniCubes[cubeIdx].position[2]
            probs[cubeIdx] = self.miniCubes[cubeIdx].probability
        
        # normalise the probability 
        # self.probs = probs / np.sum(probs)
        
        # save 
        if path.isdir(self.savePath + "/" + self.type) == False:
            os.mkdir(self.savePath + "/" + self.type)
        np.save(self.savePath + "/" + self.type + "/x_{}ttest".format(self.count), xx)
        np.save(self.savePath + "/" + self.type + "/y_{}ttest".format(self.count), yy)
        np.save(self.savePath + "/" + self.type + "/z_{}ttest".format(self.count), zz)
        np.save(self.savePath + "/" + self.type + "/probs_{}ttest".format(self.count), probs)

if __name__ == "__main__":

    # load the data 
    dsEntry = rat.dsreader("/data/snoplus/hunt-stokes/multisite/cobalt_test_extras/0/cobalt_test_53_0.root")
    count = 0 
    for entry, _ in dsEntry:

        # ignore re-triggers 
        if entry.GetEVCount() > 0:
            event = entry.GetEV(0)

            # check fit exists & is valid 
            if not event.FitResultExists("scintFitter"):
                # fit result does not exist 
                continue 
            vertex = event.GetFitResult("scintFitter").GetVertex(0)
            if (not vertex.ContainsPosition() or 
                not vertex.ContainsTime() or 
                not vertex.ValidPosition() or 
                not vertex.ValidTime()):
                # didn't fit correctly 
                continue 
                
            # apply the radial cut at 3.5m 
            eventPosition = vertex.GetPosition()
            radius = eventPosition.Mag()
            if radius < 3500:
                eventTime = vertex.GetTime()
                calibratedPMTs = event.GetCalPMTs()
                eventTime = vertex.GetTime()
                # create the MEGACUBE 
                start = time.time()
                x = Megacube(np.array([eventPosition[0],eventPosition[1],eventPosition[2]]), 
                                       1000, 10, calibratedPMTs, eventTime, "/home/hunt-stokes/megacube/data", 
                                       count, "tester")
                end = time.time()
                print("Total Runtime: ", end - start) 
                print("##### NEW CUBE #####")
                count += 1                
            
                if count == 1:
                    break 