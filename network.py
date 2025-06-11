import sys
sys.path.append("/booleanfs2/sahoo/Hegemon/")
sys.path = ["/booleanfs2/sahoo/BoNE/"] + sys.path

import bone

def T1():
    import bone
    hu = bone.hu
    nodes, edges = bone.BIGraph.readEqGraph("results/ad-net-1-res.txt")
    df = bone.BIGraph.pruneEqGraph(edges)
    df.to_csv("results/ad-net-1-eq.txt", sep="\t", header=False, index=False)
    thr = hu.getThrData(df[1])
    print(thr)
    cls = bone.BIGraph.getClusters(df)
    bone.BIGraph.saveClusters("results/ad-net-1-cls.txt", cls)
    cfile = "results/ad-net-1.rl"
    net = bone.BINetwork(cfile)
    net.init()
    cg = bone.BIGraph.getClustersGraph(net, cls)
    cg.to_csv("results/ad-net-1-eq-g.txt", sep="\t", header=False, index=False)

T1()
