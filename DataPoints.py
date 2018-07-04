from utilities import Modules





'''
    Point object has position and value.

    vars:
    pos: vector (array) with chartesian pos of point
    value: value of the point could e.g. be the value of a function in that point

    methods:
    __init__: constructs point
'''
class Point:

    def __init__(self, pos):
        self.pos = pos
        self.value = []

    def __dist__(self, other): #THIS USES EUCLIDIAN
        return Modules.np.linalg.norm(Modules.np.subtract(self.pos, other.pos))



'''
    Points object contains list of points and a number that tells us which cover it belongs to.

    vars:
    pList: Llist containing point objects

    methods:
    __init__: constructor
    __add__: Appends point to pList
    __getPosMat__: Returns matrix containing positions of points in pList
    __getValue__: Returns mean of "value" among point objects in pList

'''
class Points:

    def __init__(self):
        self.pList = []
        self.kernel = []
        self.val = []

    def __dataToPoints__(self, M, frm, to):
        M = M[len(M) + frm:len(M) + to][:]  # here we assume that we want to use the points up until the last timestamp, fix later

        for i in range(0, len(M[:, 0]), 1):
            tmpPoint = Point(M[i, :])
            tmpPoint.value = i + frm  # value is set here, fix so this is outside the method or smth
            self.__add__(tmpPoint)

    def __add__(self, p):
        self.pList.append(p)

    def __getPosMat__(self):

        if len(self.pList) > 0:
            M = Modules.np.zeros((len(self.pList), len(self.pList[0].pos)))
            for i in range(0, len(self.pList), 1):
                for j in range(0, len(self.pList[0].pos), 1):
                    M[i][j] = self.pList[i].pos[j]
            return M
        else:
            return 0

    def __getValue__(self):
        sum = 0
        for i in range(0, len(self.pList), 1):
            sum = sum + self.pList[i].value
        return sum / (len(self.pList))

    def __getKernel__(self):
        return Modules.np.mean(self.__getPosMat__(), 0)

    def __dist__(self, other): #THIS USES EUCLIDIAN
        return Modules.np.linalg.norm(Modules.np.subtract(self.__getKernel__(), other.__getKernel__()))



'''
    Groups object has array of Points objects

    vars:
    groupArray: List containing points objects

    methods:
    __init__: Constructor
    __add__: Adds Points object to groupArray
    __getValues__: returns array with values for groups
    __getKernelsList__: returns array with kernels for groups
'''
class Groups:

    def __init__(self):
        self.groupArray = []
        self.kernels = []
        self.values = []

    def __add__(self, newGroup):
        self.groupArray.append(newGroup)

    def __getValues__(self):
        tmp = Modules.np.zeros(len(self.groupArray))
        for i in range(0, len(self.groupArray), 1):
            tmp[i] = self.groupArray[i].val
        self.values = tmp
        return tmp

    def __getKernelsList__(self):
        tmp = Modules.np.zeros([len(self.groupArray), len(self.groupArray[0].pList[0].pos)])
        for i in range(0, len(self.groupArray), 1):
            tmp2 = self.groupArray[i].__getKernel__()
            for j in range(0, len(self.groupArray[0].pList[0].pos), 1):
                tmp[i][j] = tmp2[j]
        return tmp
