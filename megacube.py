import rat
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d
from matplotlib.pyplot import colormaps as cmaps
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
        originated in that minicube. Probability is updated. 
        """           

        light_path = rat.utility().GetLightPathCalculator()
        group_velocity = rat.utility().GetGroupVelocity()
        pmt_info = rat.utility().GetPMTInfo()

        # have to convert np to TVector for ROOT/RAT functions 
        position = ROOT.std.TVector3(self.position[0], self.position[1], self.position[2])
        print("There are {} triggered PMTs.".format(self.calPMTS.GetCount()))
        for ipmt in range(self.calPMTS.GetCount()):
            # obtain PMT position 
            pmt_cal = self.calPMTS.GetPMT(ipmt)
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

class Megacube(object):

    def __init__(self, position, megaSidel, numMiniCubes, calPMTs, eventTime, count):
        self.numMiniCubes = numMiniCubes               # num minicubes on each axis
        self.position     = position                   # center of the Megacube (fit vertex) 
        self.megaSidel    = megaSidel                  # side length of Megacube defining volume around fit vertex  
        self.miniSidel    = megaSidel / numMiniCubes   # side length of minicube 
        print("Side Length: {}".format(self.miniSidel))
        self.calPMTS      = calPMTs                    # all PMTs triggered by event
        self.eventTime    = eventTime
        self.load_ET1D()
        self.miniCubes    = np.array([Minicube(self.miniSidel, self.calPMTS, self.eventTime,
                            self.times, self.probs) 
                            for i in range(numMiniCubes**3)]) # 1m^3 megacube filled with side length 2cm mincubes  
        self.eventTime    = eventTime
        self.fill_positions()                          # populate minicube positions 
        done = 0 
        for cube in self.miniCubes:
            cube.obtain_probability()
            done += 1
            print("done: {}/{}".format(done, numMiniCubes**3)) 
        #self.plot_cubes()
        self.flatten_cube(count)
        
    def load_ET1D(self):
        """
        Loads the probability and timing information and saves 
        them as global arrays for use by the minicubes. 
        """

        self.times = np.load("/home/hunt-stokes/times_ET1D.npy")
        self.probs = np.load("/home/hunt-stokes/probs_ET1D.npy")


    def plot_cubes(self):
        """
        Plots the mini cube positions.
        """

        fig = plt.figure()
        ax = plt.axes(projection='3d')
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
        probs = probs / np.sum(probs)
        # xKeep = []
        # yKeep = []
        # zKeep = []
        # probKeep = []
        # for cubeIdx in range(np.size(self.miniCubes)):
        #     if probs[cubeIdx]: 
        #         xKeep.append(xx[cubeIdx])
        #         yKeep.append(yy[cubeIdx])
        #         zKeep.append(zz[cubeIdx])
        #         probKeep.append(probs[cubeIdx])
        
        scat = ax.scatter(xx,yy,zz, c = probs, marker = "s")
        #print(max(probKeep), min(probKeep))
        #ax.scatter(xKeep, yKeep, zKeep, c = probKeep)
        np.save("./x.npy", xx)
        np.save("./y.npy", yy)
        np.save("./z.npy", zz)
        np.save("./probs.npy", probs)
        
        fig.colorbar(scat)
        plt.show()

    def flatten_cube(self, count):
        """
        Takes a 3D cube and produces a flattened image. 
        """

        numCubes = self.numMiniCubes**3
        print(numCubes)
        print(int(np.sqrt(numCubes)) * int(np.sqrt(numCubes)))
        flatmap = np.zeros(numCubes)
        for cubeIdx in range(self.numMiniCubes**3):
            flatmap[cubeIdx] = self.miniCubes[cubeIdx].probability

        # reshape the map to a 2D array 
        flatmap = np.reshape(flatmap, (250,500))
        print(np.shape(flatmap))
        # imshow the flatmap 
        plt.imshow(flatmap)
        plt.savefig("./multi_site_big{}".format(0))
        



    def fill_positions(self):
        """
        Given Megacube central position, populate the mini-cube central positions. 
        """

        # Megacube side length / 2 either side of center, split into minicubes
        min_x = self.position[0] - self.megaSidel/2
        max_x = self.position[0] + self.megaSidel/2
        min_y = self.position[1] - self.megaSidel/2
        max_y = self.position[1] + self.megaSidel/2
        min_z = self.position[2] - self.megaSidel/2
        max_z = self.position[2] + self.megaSidel/2

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

            

if __name__ == "__main__":

    # load the data 
    dsEntry = rat.dsreader("/data/snoplus/hunt-stokes/multisite/cobalt_test/0/cobalt_test_4_0.root")
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
                                       1000, 50, calibratedPMTs, eventTime, count)
                count += 1                
            
                break