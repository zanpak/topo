
from utilities import Modules
from utilities import DataPoints
from utilities import Methods



'''Programm'''
M = Modules.np.load(r'C:\Users\arvid\Desktop\TDA\DATAFILES\Preem\MATRIX.npy')
#Modules.pdb.set_trace()
M2 = M[0:len(M)]
frameClusters, intervals = Methods.__GetClustersFromFrames__(M=M2, n=20000, timeOverlap=.99, nFrames=20,
                                                eps=1, min_samples=40, variable=-1)

evo = Methods.__CompareKernels__(frameClusters)
con = Methods.__Connections__(evo, 2)
indx = Methods.__ConnectedFrames__(intervals)
nds, wnd, edg, wdth, szstmp, val, fnds = Methods.__MakeFrameGraphs__(evo, con, indx, frameClusters)
clrstmp = Methods.__graphColorAssignment__(nds, val)
szs, clrs = Methods.__GraphProperties__(wnd, szstmp, clrstmp)
Methods.__MakeGraphs__(edg, wnd, wdth, szs, clrs, False)


Modules.mp.show()

Modules.pdb.set_trace()





