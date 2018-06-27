from utilities import Modules
from utilities import DataPoints
'''
    __Initialize__ Returns "dataPoints", this is a Points object that contains all the points from the dataset.

    vars:
    M: Data in matrix form
    nPoints: Number of points we want to look at

    returns: 
    dataPoints: input data represented as Points object
'''




'''
   params:
   M: data matrix
   n: nbr of points in total from data
   timeOverlap: % overlap between frames in time
   nFrames: number of frames
   eps, min_samples: params for clustering algorithm
   variable: what variable is used for coloring
   
   returns:
    frameClusters: array with nFrames elements, each element is a Groups object, containing the clusters
                  found in the corresponding frame.
    intervalse: array with start and stop for frame intervals
'''
def __GetClustersFromFrames__(M, n, timeOverlap, nFrames, eps, min_samples, variable):

    intervals = Modules.np.zeros([nFrames, 2])
    covlen = int(n/(nFrames-(nFrames-1)*timeOverlap))
    frameClusters = []
    for i in range(0, len(M[0, :]), 1):
        M[:, i] = M[:, i] / Modules.np.std(M[:, i])
    for i in range(0, nFrames, 1):
        frm2 = int(-covlen * i + timeOverlap*covlen * i)
        to2 = frm2 - covlen
        intervals[i,0] = to2
        intervals[i,1] = frm2
        dataPoints = DataPoints.Points()
        dataPoints.__DataToPoints__(M, to2, frm2) #Frm2 is closest in time
        clusters, labelArray, ntot, outliers = __ClusterDBSCAN__(points=dataPoints, eps=eps, min_samples=min_samples)
        temp = clusters.__getKernelsList__()
        clusters.kernels = temp

        for j in range(0, len(clusters.groupArray), 1): #values are set here, move this pls
            clusters.groupArray[j].val = i+1

        frameClusters.append(clusters)
        print("FRAME " + str(nFrames-(i + 1)) + ": INTERVAL: " + str(to2) + " TO " + str(frm2) +
              ", HAS " + str(len(labelArray[0])) + " PTS, "
              + str(len(clusters.groupArray)) + " CLUSTERS")
    frameClusters = frameClusters[::-1]

    return frameClusters, intervals



'''
    params:
    points: pts to be clustered
    eps: radius for "checking adjacent points"
    min_samples: min nr of pts in a cluster
    
    returns: 
    clusters: groups object with point groups that are seen as clusters 
    labelsArray: array with information about which points belong to which cluster
    ntot: nbr of clusters
    outliers: outlier points in points object
'''
def __ClusterDBSCAN__(points, eps, min_samples):
    clusters = DataPoints.Groups()
    outliers = DataPoints.Points()
    labelsArray = []
    ntot = 0
    if isinstance(points.__getPosMat__(), Modules.numbers.Number) == False:
        X = []
        posMat = points.__getPosMat__()
        X = Modules.StandardScaler().fit_transform(posMat)
        db = Modules.DBSCAN(eps=eps, min_samples=min_samples).fit(X)
        core_samples_mask = Modules.np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        labelsArray.append(labels)
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        ntot = ntot + n_clusters_
        for i in range(-1, n_clusters_, 1):
            pointsTmp = DataPoints.Points()
            for j in range(0, len(labels) - 1, 1):
                if labels[j] == i:
                    pointsTmp.__add__((points.pList[j]))
            if i != -1:
                clusters.__add__(pointsTmp)
            else:
                for k in range(0, len(pointsTmp.pList), 1):
                    outliers.__add__(pointsTmp.pList[k])
    return clusters, labelsArray, ntot, outliers



'''
    params:
    frameClusters: Frames and their corresponding clusters
    
    returns:
    EVOLUTION: Array with len(frameclusters)-1 elements, each element is an mxn matrix, where each element in the matrix
               tells us the distance between kernel i in the frame with index in the array EVO + 1, and kernel j in the 
               previous frame.
'''
def __CompareKernels__(frameClusters):
    EVOLUTION = []
    for i in range(1, len(frameClusters), 1):
        EVTMP = Modules.np.zeros([len(frameClusters[i].kernels),len(frameClusters[i-1].kernels)])
        for j in range(0, len(frameClusters[i].kernels), 1):
            for k in range(0, len(frameClusters[i-1].kernels), 1):
                tmp = Modules.np.linalg.norm(Modules.np.subtract(frameClusters[i].kernels[j]
                                                                 , frameClusters[i-1].kernels[k]))
                EVTMP[j, k] = tmp
        EVOLUTION.append(EVTMP)
    return EVOLUTION



'''
    params:
    intervals: frame intervals
    
    returns:
    fromTo: array containing information about what frames are included in each window
'''
def __ConnectedFrames__(intervals):
    fromTo = []
    for i in range(0, len(intervals), 1):
        LOWER = Modules.np.abs(intervals[i][0])
        set = False
        for j in range(0, len(intervals), 1):
            CURRENT = Modules.np.abs(intervals[j][1])

            if CURRENT > LOWER and set == False:
                fromTo.append([j-1, i])
                set = True
        if set == False:
            fromTo.append([len(intervals)-1, i])
    fromTo = fromTo[::-1]

    for i in range(0, len(fromTo), 1):
        fromTo[i][0] = Modules.np.abs(fromTo[i][0]-len(fromTo)+1)
        fromTo[i][1] = Modules.np.abs(fromTo[i][1]-len(fromTo)+1)

    return fromTo



'''
    params:
    EVOLUTION: Contains information about distances between kernels from different frames
    tolerance: radius within which kernels are connected
    
    returns:
    CONNECTIONS: Array with len(frameclusters)-1 elements, each element is an mxn matrix, where each element in the matrix
                 tells wether kernel i in the frame with index in the array EVO + 1, and kernel j in the 
                 previous frame are to be connected with an edge with respect to the given tolerance.
                
'''
def __Connections__(EVOLUTION, tolerance):
    CONNECTIONS = []
    for i in range(0, len(EVOLUTION), 1):
        tmp = Modules.np.zeros([len(EVOLUTION[i][:,0]), len(EVOLUTION[i][0,:])])
        for j in range(0, len(EVOLUTION[i][:,0]), 1):
            for k in range(0, len(EVOLUTION[i][0,:]), 1):
                if EVOLUTION[i][j, k] <= tolerance:
                    tmp[j, k] = 1
                else:
                    tmp[j, k] = 0
        CONNECTIONS.append(tmp)
    return CONNECTIONS


'''
    params:
    evo: evolution array, check __CompareKernels__ method for more detailed explanation
    con: conn array, check __Connections__ method foro more detailed explanation
    indx: ... check Connectedframes ...
    frameClusters: ... check GetClustersFromFrame ...
    
    returns:
    NODES: All nodes, eg. 1, 2, ... , n
    WINDOWNODES: All nodes that are contained in each window.
    WINDOWEDGES: All edges that are contained within each window
    SIZES: Nbr of pts in the clusters corresponding to each node
    VALUES: values for clusters corresponding to nodes
    FRAMENODES: Nodes within each frame
'''
def __MakeFrameGraphs__(evo, con, indx, frameClusters):
    NODES = []
    FRAMENODES = []
    FRAMEEDGES = []
    WINDOWEDGES = []
    SIZES = []
    VALUES = []
    FRAMEWIDTH = []
    WINDOWWIDTH = []

    for i in range(0, len(frameClusters), 1):
        for j in range(0, len(frameClusters[i].groupArray)):
            SIZES.append(len(frameClusters[i].groupArray[j].pList))
            VALUES.append(frameClusters[i].groupArray[j].__getValue__())

    #THIS PART CREATES NODES
    for i in range(0, len(evo), 1):
        if i == 0:
            tmp = Modules.np.zeros(len(evo[0][0]))
            for j in range(1, len(tmp)+1, 1):
                NODES.append(j)
                tmp[j-1] = j
            FRAMENODES.append(tmp)
        tmp2 = Modules.np.zeros(len(evo[i]))
        tmp3 = NODES[len(NODES)-1]
        for j in range(1, len(tmp2) + 1, 1):
            tmp2[j-1] = tmp3 + j
            NODES.append(NODES[len(NODES)-1]+1)
        FRAMENODES.append(tmp2)

    #EDGES IN FRAMES
    for i in range(0, len(con), 1):
        tmp = []
        tmp2 = []
        for j in range(0, len(con[i][0]), 1):
            for k in range(0, len(con[i]), 1):
                if con[i][k][j] == 1:
                    tmp2.append(evo[i][k][j])
                    tmp.append([int(FRAMENODES[i][j]), int(FRAMENODES[i+1][k])])
        FRAMEWIDTH.append(tmp2)
        FRAMEEDGES.append(tmp)

    WINDOWNODES = []
    out = indx[len(indx)-1][1]
    #EDGES IN WINDOWS
    for i in range(0, len(indx), 1):
        tmp = []
        tmp2 = []
        for j in range(indx[i][0], indx[i][1], 1):
            for k in range(0, len(FRAMEEDGES[j]), 1):
                tmp2.append(FRAMEWIDTH[j][k])
                tmp.append(FRAMEEDGES[j][k])
        WINDOWEDGES.append(tmp)
        WINDOWWIDTH.append(tmp2)
        WINDOWNODES = []


    for i in range(0, len(indx), 1):
        tmp = []
        for j in range(indx[i][0], indx[i][1]+1, 1):
            for k in range(0, len(FRAMENODES[j]), 1):
                tmp.append(int(FRAMENODES[j][k]))
        WINDOWNODES.append(tmp)

    return NODES, WINDOWNODES, WINDOWEDGES, WINDOWWIDTH, SIZES, VALUES, FRAMENODES


'''
    params:
    nodes: Nodes..
    values: values for the nodes
    
    returns:
    colors_sorted: colors for the nodes based on their values
'''
def __graphColorAssignment__(nodes, values):
    '''Because the "Colours" module works so poorly with networkx getting colors to allign with values for the
        clusters is abit of a hassle. It includes adding together multiple lists of colors, removing colors from the list
        that are not compatible with the networkx draw fcn and lastly actually assigning colors to the clusters//nodes based
        on which percentile of the value we look at it is in.'''
    intervals = 99
    intThird = 33
    breaks = Modules.np.zeros(intervals)
    for i in range(0, intervals, 1):
        breaks[i] = Modules.np.percentile(values, int(100 / intervals * i))

    from1 = Modules.Color("#a100ff")
    to1 = Modules.Color("#00a9ff")
    to2 = Modules.Color("#50ff00")
    to3 = Modules.Color("#ff0000")
    colors = list(from1.range_to(to1, intThird))
    colors.extend(list(to1.range_to(to2, intThird)))
    colors.extend(list(to2.range_to(to3, intThird)))
    colors_tmp = []
    for i in range(0, len(colors), 1):
        colors_tmp.append(str(colors[i]))
    for i in range(0, len(colors_tmp), 1):
        if len(colors_tmp[i]) < 5:
            colors_tmp[i] = colors_tmp[Modules.np.abs(i - 1)]

    colors_sorted = []
    for i in range(0, len(nodes), 1):
        set = 0
        for j in range(0, intervals, 1):
            lower = Modules.np.percentile(values, round(100 / intervals * j))
            upper = Modules.np.percentile(values, round(100 / intervals * (j + 1)))
            valTmp = values[i]
            if lower <= valTmp <= upper and set == 0:
                colors_sorted.append(colors_tmp[j])
                set = 1
        if set == 0:
            colors_sorted.append(str(Modules.Color("black")))

    return colors_sorted



'''
    params:
    wndtmp: nodes in windows
    szstmp: sizes for all nodes
    clrstmp: colors for allnodes
    
    returns:
    szs: sizes for nodes in certian window
    clrs: colors for nodes in certian window
'''
def __GraphProperties__(wndtmp, szstmp, clrstmp):
    szs = []
    clrs = []
    for i in range(0, len(wndtmp), 1):
        tmp1 = []
        tmp2 = []
        for j in range(int(wndtmp[i][0]-1), int(wndtmp[i][len(wndtmp[i])-1]), 1):
            tmp1.append(szstmp[j])
            tmp2.append(clrstmp[j])
        szs.append(tmp1)
        clrs.append(tmp2)
    return szs, clrs



'''
    params:
    edg: edges
    wnd: nodes in windows
    wdth: edge widths
    szs: node sizes
    clrs: node colors
    lbls: labels? true/false
    
    returns:
    none, but plots graphs
'''
def __MakeGraphs__(edg, wnd, wdth, szs, clrs, lbls):
    for i in range(len(wnd)-1, -1, -1):
        Modules.mp.figure()
        G = Modules.nx.Graph()
        G.add_nodes_from(wnd[i])
        G.add_edges_from(edg[i])
        pos = Modules.nx.spring_layout(G, dim=2, iterations=90)
        sizes = Modules.np.add(Modules.np.multiply(szs[i], .4), Modules.np.ones(len(szs[i])) * 50)  # resizers
        width = 1/(Modules.np.multiply(0.3, wdth[i])+Modules.np.multiply(.1, Modules.np.ones(len(wdth[i]))))
        Modules.nx.draw(G, with_labels=lbls, node_size=sizes, pos=pos, alpha=0.6, node_color=clrs[i], width=width)



'''
    params:
    node: the node you want the kernel for
    frameClusters: clusters in frames
    frameNodes: nodes in frames
    
    returns:
    kernel: the kernel corresponding to the node which was specified in the input
'''
def __getClusterKernel__(node, frameClusters, frameNodes):
    frame = 0
    frameIndx = 0
    for i in range(0, node, 1):
        if frameIndx == len(frameClusters[frame].groupArray):
            frameIndx = 0
            frame = frame + 1
        if frameNodes[frame][frameIndx] == node:
            Modules.pdb.set_trace()
            return frameClusters[frame].groupArray[frameIndx].__getKernel__()
        frameIndx = frameIndx + 1
    return "Node not found"








def __ClusterDBSCANOLD__(groups, eps, min_samples):
    # Clusters the points in the covers with DPSCAN, params are eps, min_samples
    # Also plots the clusters in different colors and saves clusters in array as "points" objects
    Modules.pdb.set_trace()
    clusters = DataPoints.Groups()
    outliers = DataPoints.Points()
    labelsArray = []
    ntot = 0
    #print("Cover nr:, points: and clusters:")
    for o in range(0, len(groups.groupArray), 1):
        pointsIndex = o
        if isinstance(groups.groupArray[o].__getPosMat__(), Modules.numbers.Number) == False:
            X = []
            posMat = groups.groupArray[o].__getPosMat__()
            X = Modules.StandardScaler().fit_transform(posMat)
            # mp.title('Clusters from DBSCAN, fr. cover ' + str(o + 1))

            db = Modules.DBSCAN(eps=eps, min_samples=min_samples).fit(X)
            core_samples_mask = Modules.np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_
            labelsArray.append(labels)
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
            #print(str(o + 1) + "    " + str(len(groups.groupArray[o].pList)) + "    " + str(n_clusters_))
            ntot = ntot + n_clusters_
            for i in range(-1, n_clusters_, 1):
                pointsTmp = DataPoints.Points(pointsIndex=pointsIndex)
                for j in range(0, len(labels) - 1, 1):
                    if labels[j] == i:
                        pointsTmp.__add__((groups.groupArray[o].pList[j]))
                if i != -1:
                    clusters.__add__(pointsTmp)
                else:
                    for k in range(0, len(pointsTmp.pList), 1):
                        outliers.__add__(pointsTmp.pList[k])
    return clusters, labelsArray, ntot, outliers




