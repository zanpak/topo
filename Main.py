
from utilities import Modules
from utilities import DataPoints
from utilities import Methods



'''Programm'''
M = Modules.pd.read_pickle(r'C:\Users\arvid\Desktop\TDA\DATAFILES\Preem\data\preemdataFIXED0629')
M = M.values
M = M[0:20157,:]
#A = M[9500:13000]
A = M[3000:13000]

#M = M[::-1]
#M = Modules.np.load(r'C:\Users\arvid\Desktop\TDA\DATAFILES\CrimeData.npy')#Modules.pdb.set_trace()

variableValues = Modules.np.zeros(len(M[:,0]))
for i in range(0, len(M), 1):
    #variableValues[i] = M[i, 0]
    variableValues[i] = i
    #variableValues[i] = Modules.rng.random()
#for i in range(0, 400, 1):
#    M[i + 1000, 3] = 10000
M2 = M[0:len(M), 1:len(M[0])]


#frameClusters, intervals, outliers = Methods.__GetClustersFromFrames__(M=M2, n=1900, timeOverlap=.9, nFrames=20,
#                                                eps=2, min_samples=20, variableValues = variableValues) #For crimedata




#1500-2000
frameClusters, intervals, outliers = Methods.__GetClustersFromFrames__(M=A, n=10000, timeOverlap=0, nFrames=100,
                                                eps=2, min_samples=2, variableValues = variableValues)
evo = Methods.__CompareKernels__(frameClusters)
con = Methods.__Connections__(evo, 2.25) #delta
indx = Methods.__ConnectedFrames__(intervals)

for i in range(0, len(indx), 1):
    indx[i][0] = 0

nds, wnd, edg, wdth, szstmp, val, fnds = Methods.__MakeFrameGraphs__(evo, con, indx, frameClusters)
clrstmp = Methods.__graphColorAssignment__(nds, val)
szs, clrs = Methods.__GraphProperties__(wnd, szstmp, clrstmp)
Methods.__MakeGraphs__(edg, wnd, wdth, szs, clrs, False)
Methods.__PlotCovers__(intervals, A)
Modules.pdb.set_trace()
Modules.mp.figure()
Methods.__windowIntervals__(indx, intervals)

Modules.mp.show()
Modules.pdb.set_trace()

''' #Nice permutations...

frameClusters, intervals, outliers = Methods.__GetClustersFromFrames__(M=A, n=10000, timeOverlap=0, nFrames=100,
                                                eps=2, min_samples=2, variableValues = variableValues)
evo = Methods.__CompareKernels__(frameClusters)
con = Methods.__Connections__(evo, 2.25) #delta
'''

#100 is last frame