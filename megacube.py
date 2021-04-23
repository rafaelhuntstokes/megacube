#!/bin/bash
import rat
import numpy as np 
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
class Metacube(object):
    """
    Parent class of Megacube and Minicube objects. Contains general, event specific, attributes 
    accessible by both daughter classes.  
    """

    def __ini__(self, eventPosition, sideLength, numMiniCubes, calPMTs, eventTime, count):
        self.eventPosition  = eventPosition   # fitted vertex position, 3D np array (x,y,z) 
        self.eventTime      = eventTime       # fitted vertex time 
        self.sideLength     = sideLength      # side length of Megacube drawn around vertex 
        self.numMiniCubes   = numMiniCubes    # number of Minicubes along 1 axis of Megacube
        self.calPMTs        = calibratedPMTs  # all the hit PMT objects in event  
        self.count          = count           # event counter to track saved output 

        # obtain side length of each Minicube
        self.miniSideLength = sideLength / numMiniCubes

        # obtain effective hit time PDFs for use in analysis 
        self.times = np.load("/home/hunt-stokes/times_ET1D.npy")
        self.probs = np.load("/home/hunt-stokes/probs_ET1D.npy")

class Megacube(Metacube):
    """
    Megacube object contains each Minicube as an attribute and methods to run the analysis and 
    extract the outputted probability density field.  
    """

    def populate(self):
        """
        Fills an array with Minicubes. 
        """
        
        self.miniCubes    = np.array([Minicube() 
                            for i in range(self.numMiniCubes**3)]) # 1m^3 megacube filled with side length 2cm mincubes  
        
        self.fill_positions()                          # populate minicube positions 
        done = 0 
        for cube in self.miniCubes:
            cube.obtain_probability()
            done += 1
            print("done: {}/{}".format(done, self.numMiniCubes**3)) 
        self.plot_cubes()
        self.flatten_cube(count)
        

    def plot_cubes(self):
        """
        Plots the mini cube positions.
        """

        # fig = plt.figure()
        # ax = plt.axes(projection='3d')
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
        self.probs = probs / np.sum(probs)
        
        #scat = ax.scatter(xx,yy,zz, c = probs, marker = "s")
        #print(max(probKeep), min(probKeep))
        #ax.scatter(xKeep, yKeep, zKeep, c = probKeep)
        np.save("x_co_2ns.npy", xx)
        np.save("y_co_2ns.npy", yy)
        np.save("z_co_2ns.npy", zz)
        np.save("probs_co_2ns.npy", probs)
        
        #fig.colorbar(scat)
        #plt.show()

    def flatten_cube(self, count):
        """
        Takes a 3D cube and produces a flattened image. 
        """

        numCubes = self.numMiniCubes**3
        print(numCubes)
        print(int(np.sqrt(numCubes)) * int(np.sqrt(numCubes)))
        #flatmap = np.zeros(numCubes)
        # for cubeIdx in range(self.numMiniCubes**3):
        #     flatmap[cubeIdx] = self.miniCubes[cubeIdx].probability
        
        # reshape the map to a 2D array 
        flatmap = np.reshape(self.probs, (50,50,50))
        # imshow the flatmap 
        #plt.imshow(flatmap)
        #plt.savefig("./multi_site_big{}".format(0))
        np.save("prob_co_cube_2ns.npy", flatmap)
        
    def fill_positions(self):
        """
        Given Megacube central position, populate the mini-cube central positions. 
        """

        # Megacube side length / 2 either side of center, split into minicubes
        min_x = self.eventPosition[0] - self.sideLength/2
        max_x = self.eventPosition[0] + self.sideLength/2
        min_y = self.eventPosition[1] - self.sideLength/2
        max_y = self.eventPosition[1] + self.sideLength/2
        min_z = self.eventPosition[2] - self.sideLength/2
        max_z = self.eventPosition[2] + self.sideLength/2

        xAxis = np.linspace(min_x, max_x, self.numMiniCubes)
        yAxis = np.linspace(min_y, max_y, self.numMiniCubes)
        zAxis = np.linspace(min_z, max_z, self.numMiniCubes)

        idx = 0 
        for x in xAxis: 
            for y in yAxis:
                for z in zAxis:
                    self.miniCubes[idx].position[0] = x
                    self.miniCubes[idx].position[1] = y
                    self.miniCubes[idx].position[2] = z
                    idx += 1

class Minicube(Megacube):
    def __init__(self):
        self.probability = 0
        self.position = np.zeros(3)  


    def obtain_probability(self):
        """
        Function draws line to each hit PMT from center of minicube. Time of 
        flight is obtained. Emission time is obtained via: 

        T_em = T_trig - T_flight - T_event  

        ET1D "effective time 1D" PDF is used to return probability photon 
        originated in that minicube. Probability is updated. 
        """           

        light_path = rat.utility().GetLightPathCalculator()
        group_velocity = rat.utility().GetGroupVelocity()
        pmt_info = rat.utility().GetPMTInfo()

        # have to convert np to TVector for ROOT/RAT functions 
        position = ROOT.std.TVector3(self.position[0], self.position[1], self.position[2])
        #print("There are {} triggered PMTs.".format(self.calPMTS.GetCount()))
        for ipmt in range(self.calPMTs.GetCount()):
            # obtain PMT position 
            pmt_cal = self.calPMTs.GetPMT(ipmt)
            pmtPos = pmt_info.GetPosition(pmt_cal.GetID())

            # find lightpath from minicube to triggered pmt 
            light_path.CalcByPosition(position, pmtPos)
            inner_av_distance = light_path.GetDistInInnerAV()
            av_distance = light_path.GetDistInAV()
            water_distance = light_path.GetDistInWater()

            # obtain time of flight from minicube to triggered pmt 
            transit_time = group_velocity.CalcByDistance(inner_av_distance, av_distance, water_distance)

            # obtain emission time 
            tEm = np.array([pmt_cal.GetTime() - transit_time - eventTime])

            # obtain probability from ET1D PDF 
            binIdx = np.digitize(tEm, bins = self.times)
            self.probability += self.probs[binIdx]            

if __name__ == "__main__":

    # load the data 
    dsEntry = rat.dsreader("/data/snoplus/hunt-stokes/multisite/cobalt_test_extras/0/cobalt_test_53_0.root")
    count = 0 
    for entry, _ in dsEntry:
        # light_path = rat.utility().GetLightPathCalculator()
        # group_velocity = rat.utility().GetGroupVelocity()
        # pmt_info = rat.utility().GetPMTInfo()

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
                x = Megacube(np.array([eventPosition[0],eventPosition[1],eventPosition[2]]), 
                                       1000, 5, calibratedPMTs, eventTime, count)
                count += 1                
            
                break