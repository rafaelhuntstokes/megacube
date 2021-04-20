#import rat
import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d

"""
Script to implement the MEGACUBE (TM) multi-site analysis method. Fitted vertex forms the centre of a 
cube, of side length 1m. This MEGACUBE is filled with minicube cubes, side length = 1/2 pos. resolution. 
For each mini-cube, a line is drawn to every hit PMT. ET1D PDF is used to obtain probability light originated 
from that mini-cube. The probability is stored within the mini-cube, and an image of the 3D probability 
field is created. Funky visuals. 0vbb and Co-60 (single vs multi-site) events will look different, and an ML 
should be able to discriminate between these events with ease.  
"""
class Minicube(object):
    def __init__(self, miniSidel):
        self.probability = 0
        self.position = np.zeros(3)    # center of the Minicube
        self.sidel = miniSidel  # side length of Minicube   

    def obtain_probability(self, hitPMTS, ET1D):
        """
        Function draws line to each hit PMT from center of minicube. Time of 
        flight is obtained. Emission time is obtained via: 

        T_em = T_trig - T_flight 

        ET1D PDF is used to return probability photon originated in that minicube. 
        Probability is updated. 
        """           




class Megacube(object):

    def __init__(self, position, megaSidel, numMiniCubes):
        self.numMiniCubes = numMiniCubes               # num minicubes on each axis
        self.position     = position                   # center of the Megacube (fit vertex) 
        self.megaSidel    = megaSidel                  # side length of Megacube defining volume around fit vertex  
        self.miniSidel    = megaSidel / numMiniCubes   # side length of minicube 
        self.miniCubes    = np.array([Minicube(self.miniSidel) 
                            for i in range(numMiniCubes**3)]) # 1m^3 megacube filled with side length 2cm mincubes  
        
        self.fill_positions()                          # populate minicube positions 

        self.plot_cubes()

        

    def plot_cubes(self):
        """
        Plots the mini cube positions.
        """
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        xx = np.zeros(self.numMiniCubes**3)
        yy = np.zeros(self.numMiniCubes**3)
        zz = np.zeros(self.numMiniCubes**3)

        for cubeIdx in range(np.size(self.miniCubes)):
            xx[cubeIdx] = self.miniCubes[cubeIdx].position[0]
            yy[cubeIdx] = self.miniCubes[cubeIdx].position[1]
            zz[cubeIdx] = self.miniCubes[cubeIdx].position[2]
        ax.scatter(xx,yy,zz)
        plt.show()

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

            

#if __name__ == "__main__":
    # load the data 
    # dsEntry = rat.dsreader("/data/snoplus/hunt-stokes/multisite/tts_3.7/0/*.root")

    # for entry, _ in dsEntry:
    #     light_path = rat.utility().GetLightPathCalculator()
    #     group_velocity = rat.utility().GetGroupVelocity()
    #     pmt_info = rat.utility().GetPMTInfo()

    #     # ignore re-triggers 
    #     if entry.GetEVCount() > 0:
    #         event = entry.GetEV(0)

    #         # check fit exists & is valid 
    #         if not event.FitResultExists("scintFitter"):
    #             # fit result does not exist 
    #             continue 
    #         vertex = event.GetFitResult("scintFitter").GetVertex(0)
    #         if (not vertex.ContainsPosition() or 
    #             not vertex.ContainsTime() or 
    #             not vertex.ValidPosition() or 
    #             not vertex.ValidTime()):
    #             # didn't fit correctly 
    #             continue 
                
    #         # apply the radial cut at 3.5m 
    #         eventPosition = vertex.GetPosition()
    #         radius = eventPosition.Mag()
    #         if radius < 3500:
    #             eventTime = vertex.GetTime()
    #             calibratedPMTs = event.GetCalPMTs()
                
                # create the MEGACUBE 

x = Megacube(np.array([0,0,0]), 1000, 10)