import bone
pd = bone.pd
re = bone.re
hu = bone.hu
np = bone.np
os = bone.os

basedir = os.path.join(os.getcwd(), os.path.dirname(__file__)) + "/"

def getMSigDB(gs):
    url = "https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=" + gs + "&fileType=txt"
    df = pd.read_csv(url, sep="\t")
    df.columns.values[0] = 'ID'
    l1 = [list(df.ID[1:])]
    wt1 = [1]
    return wt1, l1
bone.getMSigDB = getMSigDB

def plotTitleBar(cval, atypes, params):
    dpi = 100
    if 'dpi' in params:
        dpi = params['dpi']
    w,h = (5, 0.8)
    if 'w' in params:
        w = params['w']
    if 'h' in params:
        h = params['h']
    color_sch1 = ["#3B449C", "#B2509E","#EA4824"]
    color_sch1 = ["#00CC00", "#EFF51A","#EC008C", "#F7941D", "#808285",
            'cyan', 'blue', 'black', 'green', 'red']
    if 'acolor' in params:
        color_sch1 = params['acolor']
    if 'cval' in params:
        cval = params['cval']

    ax = None
    if 'ax' in params:
        ax = params['ax']
    if ax is None:
        fig = bone.plt.figure(figsize=(w,h), dpi=dpi)
        ax = fig.add_subplot(1, 1, 1)
    nAt = len(cval[0])
    extent = [0, nAt, 0, 5]
    ax.axis(extent)
    cmap = bone.colors.ListedColormap(color_sch1)
    boundaries = range(len(color_sch1) + 1)
    norm = bone.colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    #ax.imshow(cval, interpolation='none', cmap=cmap, \
    #                  norm=norm, extent=extent, aspect="auto")
    y = [0, 5]
    x = bone.np.arange(nAt + 1)
    ax.pcolormesh(x, y, cval, cmap=cmap, norm=norm, zorder=-1.0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(top=False, left=False, bottom=False, right=False)
    ax.set_xticks(bone.np.arange(0, nAt, 1))
    ax.grid(which='major', alpha=0.2, linestyle='-', linewidth=0.5,
            color='black', zorder=2.0)
    for edge, spine in ax.spines.items():
                spine.set_visible(False)
    divider = bone.make_axes_locatable(ax)
    width = bone.axes_size.AxesX(ax, aspect=1./20)
    spaceAnn = 70
    widthAnn = 3
    tAnn = 1
    if 'spaceAnn' in params:
        spaceAnn = params['spaceAnn']
    if 'widthAnn' in params:
        widthAnn = params['widthAnn']
    if 'tAnn' in params:
        tAnn = params['tAnn']
    pad = bone.axes_size.Fraction(0.1, width)
    lax = divider.append_axes("top", size="100%", pad="20%", frame_on=False)
    lax.axison = False
    lax.axis(extent)
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.grid(False)
    lax.tick_params(top=False, left=False, bottom=False, right=False)
    if 'atypes' in params:
        atypes = params['atypes']
    bone.barTop(lax, atypes, color_sch1, params)
    return ax

bone.plotTitleBar = plotTitleBar

def plotTitleBarSingle(ana, ax=None, desc=None):
    expr = ana.f_ranks
    v1 = [expr[i-2] for i in ana.state[0]]
    v2 = [expr[i-2] for i in ana.state[1]]
    t, p = bone.ttest_ind(v1,v2, equal_var=False)
    expr = [expr[i-2] for i in ana.order]
    ana.cval = [[ana.aval[ana.order[i]] for i in bone.np.argsort(expr)]]
    ax = ana.printTitleBar({'w':4, 'h':0.5, 'spaceAnn': len(ana.order)/len(ana.atypes),
                            'tAnn': 1, 'widthAnn':1, 'ax':ax})
    if desc == None:
        desc = ana.h.getSource() + f"({len(ana.state[0])},{len(ana.state[1])})"
    ax.text(-1, 2, desc, horizontalalignment='right', verticalalignment='center')
    actual = ana.cval[0]
    score = [expr[i] for i in bone.np.argsort(expr)]
    fpr, tpr, thrs = bone.roc_curve(actual, score, pos_label=1)
    auc = bone.auc(fpr, tpr)
    roc_auc = "%.2f %.3g" % (auc, p)
    print(roc_auc)
    ax.text(len(ana.cval[0]), 2, roc_auc)
    return [auc, p, t]
bone.plotTitleBarSingle = plotTitleBarSingle

def adjustFigure(ana, fig, xlabel="Score", wfactor=1):
    rocauc = ana.getROCAUC()
    ax = fig.axes[1]
    #ax.set_xlim([-5, 10])
    ax.set_title(f"ROC-AUC = {rocauc}", fontsize=16)
    ax.set_xlabel(xlabel, fontsize=20)
    ax.tick_params(axis='x', which='major', labelsize=12)
    ax.tick_params(axis='y', which='major', labelsize=20)
    #fig.patch.set_facecolor('#FFF7CD')
    #ax.set_facecolor('#FFF7CD')
    for tobj in ax.texts:
        tobj.set_fontsize(16)
        #tobj.set_color('black')
    #for tobj in fig.findobj(plt.Text):
    #    tobj.set_color('black')        
    fw, fh = fig.get_size_inches()
    bbox = fig.axes[2].get_position()
    aw, ah = bbox.width * fw, bbox.height * fh
    xlim = fig.axes[2].get_xlim()
    ylim = fig.axes[2].get_ylim()
    for p in fig.axes[2].patches:
        bbox1 = p.get_bbox()
        pw = bbox1.width/(xlim[1]-xlim[0]) * aw * 72
        ph = bbox1.height/(ylim[1]-ylim[0]) * ah * 72
        sw = ph/72/aw*(xlim[1]-xlim[0])/wfactor
        p.set_width(sw)
    return
bone.adjustFigure = adjustFigure

def getADGeneSet(index=0):
    wt1, l1, name = [], [], "Unknown"
    if index == 0 or index == "BoNE":
        df = pd.read_csv("results/model-1.txt", sep="\t")
        wt1 = list(pd.to_numeric(df.columns))
        l1 = [list(df[k].dropna()) for k in df.columns]
        name = "BoNE"
    if index == 1 or index == 'PMID29656362-Kim':
        up = bone.getEntries(basedir + "database/PMID29656362-Kim-up.txt", 0)
        down = bone.getEntries(basedir + "database/PMID29656362-Kim-down.txt", 0)
        wt1, l1 = [1, -1], [up, down]
        name = 'PMID29656362-Kim'
    if index == 2 or index == 'PMID25710473-IL15/RA':
        wt1, l1 = [1], [['IL15', 'IL15RA']]
        name = 'PMID25710473-IL15/RA'
    if index == 3 or index == 'PMID33339153-IL6/R':
        wt1, l1 = [1], [['IL6', 'IL6R']]
        name = 'PMID33339153-IL6/R'
    if index == 4 or index == 'PMID34410590-GWAS':
        wt1, l1 = [-1], [["PCDH11Y" , "PCDH11X", "TPTE2", "LOC107985902", "MUC16" , "LINC01621"]]
        name = 'PMID34410590-GWAS'
    if index == 5 or index == 'PMID33057964-Sanfilippo':
        wt1, l1 = [-1], [ ["GNG13", "PCP2", "GPSM4", "PCP4", "PEP19", "CBLN", "CBLN3"]]
        name = 'PMID33057964-Sanfilippo'
    if index == 6 or index == 'PMID34234395-Madar':
        wt1, l1 = [1], [["CNPY3", "GPR84", "HIST1H2AB", "HIST1H2AE", "IFNAR1", "LMO3",
                     "MYO18A", "N4BP2L1", "PML", "SLC4A4", "ST8SIA4", "TLE1" , "N4BP2L1"] ]
        name = 'PMID34234395-Madar'
    if index == 7 or index == 'PMID34384464-Bai':
        up = bone.getEntries(basedir + "database/PMID34384464-Bai-up.txt", 0)
        down = bone.getEntries(basedir + "database/PMID34384464-Bai-down.txt", 0)
        wt1, l1 = [1, -1], [up, down]
        name = 'PMID34384464-Bai'
    if index == 8 or index == 'PMID34194312-Lasso-Yu':
        lassof=basedir + "database/PMID34194312-Lasso-Yu.txt"
        lasso = bone.getEntries(lassof, 0)
        ddeg_yuf=basedir + "database/PMID34194312-DEG-Yu-down.txt"
        udeg_yuf=basedir + "database/PMID34194312-DEG-Yu-up.txt"
        ddeg_yu = bone.getEntries(ddeg_yuf, 0)
        udeg_yu = bone.getEntries(udeg_yuf, 0)
        wt1, l1 = [-1, 1], [set(lasso).intersection(ddeg_yu), set(lasso).intersection(udeg_yu)]
        name = 'PMID34194312-Lasso-Yu'
    if index == 9 or index == 'PMID34194312-DEG-Yu':
        ddeg_yuf=basedir + "database/PMID34194312-DEG-Yu-down.txt"
        udeg_yuf=basedir + "database/PMID34194312-DEG-Yu-up.txt"
        ddeg_yu = bone.getEntries(ddeg_yuf, 0)
        udeg_yu = bone.getEntries(udeg_yuf, 0)
        wt1, l1 = [1, -1], [udeg_yu, ddeg_yu]
        name = 'PMID34194312-DEG-Yu'
    if index == 10 or index == 'PMID34194312-HUB-Yu':
        dhup_yuf=basedir + "database/PMID34194312-HUB-Yu-down.txt"
        uhub_yuf=basedir + "database/PMID34194312-HUB-Yu-down.txt"
        dhup_yu = bone.getEntries(dhup_yuf, 0)
        uhub_yu = bone.getEntries(uhub_yuf, 0)
        wt1, l1 = [1, -1], [uhub_yu, dhup_yu]
        name = 'PMID34194312-HUB-Yu'
    if index == 11 or index == 'PMID32699331-Zakeri': # Blood based
        final_list = ["MBOAT1", "ARMC7", "RABL2B", "HNRNPUL1", "LAMTOR1",
                      "PLAGL2", "CREBRF", "LCOR",  "MRI1"]
        wt1, l1 = [1], [ final_list]
        name = 'PMID32699331-Zakeri'
    if index == 12 or index == 'PMID32699331-Zakeri-deg': # Blood based
        deg_z = bone.getEntries(basedir + "database/PMID32699331-Zakeri-deg.txt", 0)
        wt1, l1 = [1], [ deg_z]
        name = 'PMID32699331-Zakeri-deg'
    if index == 13 or index == 's41048-019-0086-2-Wu':
        up = bone.readList(basedir + "database/10.1007_s41048-019-0086-2-Wu-up.txt")
        down = bone.readList(basedir + "database/10.1007_s41048-019-0086-2-Wu-down.txt")
        wt1, l1 = [1, -1], [up, down]
        name = 's41048-019-0086-2-Wu'
    if index == 14 or index == 'PMID31984950-Lui':
        up_hub = [ "ENO2", "NFKBIA",  "ACACB"]
        down_hub =  [ "CCT2", "CALM2", "ATP5B", "MDH1", "PPP2CA", "PSMD14", "ATP5C1",
                     "LDHA", "CCT4", "COPS5", "TXN", "PGK1", "NDUFA4", "ACTR10", "IMMT",
                     "ATP5F1", "NDUFAB1", "CUL1", "PSMB7", "SKP1", "MTHFD1", "TUBA4A",
                     "PSMA5", "TUBA1B", "TUBA1C", "HINT1", "NME1", "PSMA1", "TUBB4B"]
        wt1, l1 = [1, -1], [up_hub, down_hub]
        name = 'PMID31984950-Lui'
    if index == 15 or index == 'PMID32255787-Klein':
        up = bone.readList(basedir + "database/PMID32255787-Klein-up.txt")
        down = bone.readList(basedir + "database/PMID32255787-Klein-down.txt")
        wt1, l1 = [1, -1], [up, down]
        name = 'PMID32255787-Klein'
    if index == 16 or index == 'PMID29802388-S3m109':
        ff= bone.getEntries(basedir + "database/PMID29802388-S3m109.txt", 0)
        wt1, l1 = [1], [ ff]
        name = 'PMID29802388-S3m109'
    if index == 17 or index == 'PMID27799057-S7':
        up = bone.readList(basedir + "database/PMID27799057-S7-up.txt")
        down = bone.readList(basedir + "database/PMID27799057-S7-down.txt")
        wt1, l1 = [1, -1], [up, down]
        name = 'PMID27799057-S7'
    if index == 18 or index == 'PMID27989508-S8M14':
        ff= bone.getEntries(basedir + "database/PMID27989508-S8M14.txt" , 0)
        wt1, l1 = [1], [ ff]
        name = 'PMID27989508-S8M14'
    if index == 19 or index == 'PMID29110684-S5D2':
        up = bone.readList(basedir + "database/PMID29110684-S5D2-up.txt")
        down = bone.readList(basedir + "database/PMID29110684-S5D2-down.txt")
        wt1, l1 = [1, -1], [up, down]
        name = 'PMID29110684-S5D2'
    if index == 20 or index == 'PMID29107053-S3':
        up = bone.readList(basedir + "database/PMID29107053-S3-up.txt")
        down = bone.readList(basedir + "database/PMID29107053-S3-down.txt")
        wt1, l1 = [1, -1], [up, down]
        name = 'PMID29107053-S3'
    if index == 21 or index == 'PMID29358399-S02':
        up = bone.readList(basedir + "database/PMID29358399-S02-up.txt")
        down = bone.readList(basedir + "database/PMID29358399-S02-down.txt")
        wt1, l1 = [1, -1], [up, down]
        name = 'PMID29358399-S02'
    if index == 22 or index == 'PMID30136084-S1':
        up = bone.readList(basedir + "database/PMID30136084-S1-up.txt")
        down = bone.readList(basedir + "database/PMID30136084-S1-down.txt")
        wt1, l1 = [1, -1], [up, down]
        name = 'PMID30136084-S1'
    if index == 23 or index == 'PMID29598827-S8a':
        ff= bone.getEntries(basedir + "database/PMID29598827-S8a.txt" , 0)
        wt1, l1 = [1], [ ff]
        name = 'PMID29598827-S8a'
    if index == 24 or index == 'PMID30297968-S8':
        up = bone.readList(basedir + "database/PMID30297968-S8-up.txt")
        down = bone.readList(basedir + "database/PMID30297968-S8-down.txt")
        wt1, l1 = [1, -1], [up, down]
        name = 'PMID30297968-S8'
    return wt1, l1, name
bone.getADGeneSet = getADGeneSet

def ADComparison(ana, desc=''):
    def getL(l1):
        return '(' + ",".join([str(len(k)) for k in l1]) +')'
    
    res = []

    for i in range(25):
        wt1, l1, name = getADGeneSet(i)
        ann = [re.sub(" .*", "", ana.h.getSource()),'Hs',name, getL(l1)]
        res += [ana.getStats(l1, wt1, ann)]
    
    cols = ['GSEID', 'ROC-AUC', 'pvalue', '#Cont', '#Expt',
            'Series', 'Species', 'Signature', '#Genes']
    df = pd.DataFrame(res, columns=cols)
    df['Condition'] = desc
    return df
bone.ADComparison = ADComparison

def ADResults(arg, desc=''):
    ana, wt1, l1, signame = arg
    def getL(l1):
        return '(' + ",".join([str(len(k)) for k in l1]) +')'

    res = []
    ann = [re.sub(" .*", "", ana.h.getSource()),'Hs',signame, getL(l1)]
    res += [ana.getStats(l1, wt1, ann) + [desc] ]

    return res
bone.ADResults = ADResults

def ADSummary(ana, desc=''):
    src = re.sub(" .*", "", ana.source)
    num = '(' + ",".join([str(len(k)) for k in ana.state]) +')'
    res = [src, num, re.sub(" .*", "", ana.h.getSource()), desc]
    return res
bone.ADSummary = ADSummary

def ADTraining(wt1, l1, signame='BoNE'):
    ana = bone.IBDAnalysis()
    arg = (ana, wt1, l1, signame)
    res = []
    ana.getFriedman2017()
    res += ADResults(arg, desc=ana.h.getSource() + "-fusiform gyrus")
    ana.getPatel2019(2, 0)
    res += ADResults(arg, desc=ana.h.getSource() + "-Entorhinal_Cortex")
    ana.getWebster2009()
    res += ADResults(arg, desc=ana.h.getSource() + "-cortical")
    cols = ['GSEID', 'ROC-AUC', 'pvalue', '#Cont', '#Expt',
            'Series', 'Species', 'Signature', '#Genes', 'Condition']
    df = pd.DataFrame(res, columns=cols)
    return df
bone.ADTraining = ADTraining

def ADTrainingAll():
    ana = bone.IBDAnalysis()
    res = []
    ana.getFriedman2017()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-fusiform gyrus")]
    ana.getPatel2019(2, 0)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Entorhinal_Cortex")]
    ana.getWebster2009()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-cortical")]
    df = pd.concat(res, sort=True)
    return df
bone.ADTrainingAll = ADTrainingAll

def ADValidations(wt1, l1, signame='BoNE'):
    ana = bone.IBDAnalysis()
    arg = (ana, wt1, l1, signame)
    res = []
    ana.getPatel2019(2, 1)
    res += ADResults(arg, desc=ana.h.getSource() + "-Temporal_Cortex")
    ana.getPatel2019(2, 2)
    res += ADResults(arg, desc=ana.h.getSource() + "-Frontal_Cortex")
    ana.getLiang2007()
    res += ADResults(arg, desc=ana.h.getSource() + "-all")
    ana.getLiang2007(2, 0)
    res += ADResults(arg, desc=ana.h.getSource() + "-MTG")
    ana.getLiang2007(2, 1)
    res += ADResults(arg, desc=ana.h.getSource() + "-HIP")
    ana.getLiang2007(2, 2)
    res += ADResults(arg, desc=ana.h.getSource() + "-SFG")
    ana.getLiang2007(2, 3)
    res += ADResults(arg, desc=ana.h.getSource() + "-PVC")
    ana.getLiang2007(2, 5)
    res += ADResults(arg, desc=ana.h.getSource() + "-EC")
    ana.getLiang2007(3)
    res += ADResults(arg, desc=ana.h.getSource() + "-PCC")
    ana.getWang2016(3)
    res += ADResults(arg, desc=ana.h.getSource() + "-all")
    ana.getWang2016(4, 0)
    res += ADResults(arg, desc=ana.h.getSource() + "-Amygdala")
    ana.getWang2016(4, 1)
    res += ADResults(arg, desc=ana.h.getSource() + "-Nucleus Accumbens")
    ana.getWang2016II(3)
    res += ADResults(arg, desc=ana.h.getSource() + "-II-all")
    ahash = {'Inferior Temporal Gyrus':0, 'Parahippocampal Gyrus':1,
            'Middle Temporal Gyrus':2, 'Occipital Visual Cortex':3,
            'Prefrontal Cortex':4, 'Hippocampus':5, 'Caudate Nucleus':6,
            'Frontal Pole':7, 'Precentral Gyrus':8,
            'Posterior Cingulate Cortex':9, 'Superior Temporal Gyrus':10,
            'Superior Parietal Lobule':11, 'Temporal Pole':12,
            'Anterior Cingulate':13, 'Inferior Frontal Gyrus':14,
            'Dorsolateral Prefrontal Cortex':15, 'Putamen':16}
    for k in ahash:
        ana.getWang2016II(4, ahash[k])
        res += ADResults(arg, desc=ana.h.getSource() + "-" + k)
    ana.getBerchtold2014()
    res += ADResults(arg, desc=ana.h.getSource() + "-all")
    ana.getMarttinen2019()
    res += ADResults(arg, desc=ana.h.getSource() + "-Braak 0 vs")
    ana.getMarttinen2019(2)
    res += ADResults(arg, desc=ana.h.getSource() + "-Braak 0-2 vs")
    ana.getLu2014()
    res += ADResults(arg, desc=ana.h.getSource() + "-Ageing")
    ana.getLu2004()
    res += ADResults(arg, desc=ana.h.getSource() + "-Ageing")
    ana.getBerchtold2019(2)
    res += ADResults(arg, desc=ana.h.getSource() + "-Activity")
    ana.getSrinivasan2020()
    res += ADResults(arg, desc=ana.h.getSource() + "-fusiform gyrus")
    ahash = {'myeloid':0, 'astrocyte':1, 'neuron':2, 'endothelial':3}
    for k in ahash:
        ana.getSrinivasan2020II(2, ahash[k])
        res += ADResults(arg, desc=ana.h.getSource() + "-" + k)
    ana.getFriedman2018()
    res += ADResults(arg, desc=ana.h.getSource() + "-fusiform gyrus")
    ana.getMiyashita2014(2)
    res += ADResults(arg, desc=ana.h.getSource() + "-ADC")
    ana.getMcKay2019(2)
    res += ADResults(arg, desc=ana.h.getSource() + "-ADC")
    ana.getLow2021(2)
    res += ADResults(arg, desc=ana.h.getSource() + "-BA9")
    ana.getDunckley2006()
    res += ADResults(arg, desc=ana.h.getSource() + "-EC")
    ana.getBlalock2011(2)
    res += ADResults(arg, desc=ana.h.getSource() + "-HIP")
    cols = ['GSEID', 'ROC-AUC', 'pvalue', '#Cont', '#Expt',
            'Series', 'Species', 'Signature', '#Genes', 'Condition']
    df = pd.DataFrame(res, columns=cols)
    return df

def ADValidationsII(wt1, l1, signame='BoNE'):
    ana = bone.IBDAnalysis()
    arg = (ana, wt1, l1, signame)
    res = []
    ana.getNativio2020(2)
    res += ADResults(arg, desc=ana.h.getSource() + "-HIP")
    ana.getPascal2011(2, 0, 0)
    res += ADResults(arg, desc=ana.h.getSource() + "-EC")
    ana.getNancy2020()
    res += ADResults(arg, desc=ana.h.getSource() + "-PFC")
    ahash = {'unID':0, 'endo':1, 'oligo':2, 'astro':3, 'OPC':4, 'doublet':5,
             'neuron':6, 'mg':7}
    for k in ahash:
        ana.getChew2019(2, ahash[k])
        res += ADResults(arg, desc=ana.h.getSource() + "-" + k)
    ana.getADPooledDyn(2, 0)
    res += ADResults(arg, desc="Pooled-no ADC")
    ana.getADPooledDyn(2, 1)
    res += ADResults(arg, desc="Pooled-ADC")
    ahash = {'TCx':0, 'PCx':1, 'HIP':2, 'FWM':3}
    for k in ahash:
        ana.getMiller2017(2, ahash[k])
        res += ADResults(arg, desc=ana.h.getSource() + "-" + k)
    ana.getZhang2013("MAC39.1")
    res += ADResults(arg, desc=ana.h.getSource() + "-CR")
    ana.getZhang2013("MAC39.2")
    res += ADResults(arg, desc=ana.h.getSource() + "-PFC")
    ana.getZhang2013("MAC39.3")
    res += ADResults(arg, desc=ana.h.getSource() + "-VC")
    cols = ['GSEID', 'ROC-AUC', 'pvalue', '#Cont', '#Expt',
            'Series', 'Species', 'Signature', '#Genes', 'Condition']
    df = pd.DataFrame(res, columns=cols)
    return df

def ADValidationAll():
    ana = bone.IBDAnalysis()
    res = []
    ana.getPatel2019(2, 1)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Temporal_Cortex")]
    ana.getPatel2019(2, 2)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Frontal_Cortex")]
    ana.getLiang2007()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-all")]
    ana.getLiang2007(2, 0)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-MTG")]
    ana.getLiang2007(2, 1)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-HIP")]
    ana.getLiang2007(2, 2)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-SFG")]
    ana.getLiang2007(2, 3)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-PVC")]
    ana.getLiang2007(2, 5)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-EC")]
    ana.getLiang2007(3)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-PCC")]
    ana.getWang2016(3)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-all")]
    ana.getWang2016(4, 0)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Amygdala")]
    ana.getWang2016(4, 1)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Nucleus Accumbens")]
    ana.getWang2016II(3)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-all")]
    ahash = {'Inferior Temporal Gyrus':0, 'Parahippocampal Gyrus':1,
            'Middle Temporal Gyrus':2, 'Occipital Visual Cortex':3,
            'Prefrontal Cortex':4, 'Hippocampus':5, 'Caudate Nucleus':6,
            'Frontal Pole':7, 'Precentral Gyrus':8,
            'Posterior Cingulate Cortex':9, 'Superior Temporal Gyrus':10,
            'Superior Parietal Lobule':11, 'Temporal Pole':12,
            'Anterior Cingulate':13, 'Inferior Frontal Gyrus':14,
            'Dorsolateral Prefrontal Cortex':15, 'Putamen':16}
    for k in ahash:
        ana.getWang2016II(4, ahash[k])
        res += [ADComparison(ana, desc=ana.h.getSource() + "-" + k)]
    ana.getBerchtold2014()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-all")]
    ana.getMarttinen2019()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Braak 0 vs")]
    ana.getMarttinen2019(2)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Braak 0-2 vs")]
    ana.getLu2014()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Ageing")]
    ana.getLu2004()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Ageing")]
    ana.getBerchtold2019(2)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-Activity")]
    ana.getSrinivasan2020()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-fusiform gyrus")]
    ahash = {'myeloid':0, 'astrocyte':1, 'neuron':2, 'endothelial':3}
    for k in ahash:
        ana.getSrinivasan2020II(2, ahash[k])
        res += [ADComparison(ana, desc=ana.h.getSource() + "-" + k)]
    ana.getFriedman2018()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-fusiform gyrus")]
    ana.getMiyashita2014(2)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-ADC")]
    ana.getMcKay2019(2)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-ADC")]
    ana.getLow2021(2)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-BA9")]
    ana.getDunckley2006()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-EC")]
    ana.getBlalock2011(2)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-HIP")]
    # Set II
    ana.getNativio2020(2)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-HIP")]
    ana.getPascal2011(2, 0, 0)
    res += [ADComparison(ana, desc=ana.h.getSource() + "-EC")]
    ana.getNancy2020()
    res += [ADComparison(ana, desc=ana.h.getSource() + "-PFC")]
    ahash = {'unID':0, 'endo':1, 'oligo':2, 'astro':3, 'OPC':4, 'doublet':5,
             'neuron':6, 'mg':7}
    for k in ahash:
        ana.getChew2019(2, ahash[k])
        res += [ADComparison(ana, desc=ana.h.getSource() + "-" + k)]
    ana.getADPooledDyn(2, 0)
    res += [ADComparison(ana, desc="Pooled-no ADC")]
    ana.getADPooledDyn(2, 1)
    res += [ADComparison(ana, desc="Pooled-ADC")]
    ahash = {'TCx':0, 'PCx':1, 'HIP':2, 'FWM':3}
    for k in ahash:
        ana.getMiller2017(2, ahash[k])
        res += [ADComparison(ana, desc=ana.h.getSource() + "-" + k)]
    ana.getZhang2013("MAC39.1")
    res += [ADComparison(ana, desc=ana.h.getSource() + "-CR")]
    ana.getZhang2013("MAC39.2")
    res += [ADComparison(ana, desc=ana.h.getSource() + "-PFC")]
    ana.getZhang2013("MAC39.3")
    res += [ADComparison(ana, desc=ana.h.getSource() + "-VC")]
    df = pd.concat(res, sort=True)
    return df
bone.ADValidationAll = ADValidationAll

def ADDatasets():
    ana = bone.IBDAnalysis()
    res = []
    ana.getFriedman2017()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-fusiform gyrus")]
    ana.getPatel2019(2, 0)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Entorhinal_Cortex")]
    ana.getWebster2009()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-cortical")]
    ana.getPatel2019(2, 1)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Temporal_Cortex")]
    ana.getPatel2019(2, 2)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Frontal_Cortex")]
    ana.getLiang2007()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-all")]
    ana.getLiang2007(2, 0)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-MTG")]
    ana.getLiang2007(2, 1)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-HIP")]
    ana.getLiang2007(2, 2)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-SFG")]
    ana.getLiang2007(2, 3)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-PVC")]
    ana.getLiang2007(2, 5)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-EC")]
    ana.getLiang2007(3)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-PCC")]
    ana.getWang2016(3)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-all")]
    ana.getWang2016(4, 0)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Amygdala")]
    ana.getWang2016(4, 1)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Nucleus Accumbens")]
    ana.getWang2016II(3)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-all")]
    ahash = {'Inferior Temporal Gyrus':0, 'Parahippocampal Gyrus':1,
            'Middle Temporal Gyrus':2, 'Occipital Visual Cortex':3,
            'Prefrontal Cortex':4, 'Hippocampus':5, 'Caudate Nucleus':6,
            'Frontal Pole':7, 'Precentral Gyrus':8,
            'Posterior Cingulate Cortex':9, 'Superior Temporal Gyrus':10,
            'Superior Parietal Lobule':11, 'Temporal Pole':12,
            'Anterior Cingulate':13, 'Inferior Frontal Gyrus':14,
            'Dorsolateral Prefrontal Cortex':15, 'Putamen':16}
    for k in ahash:
        ana.getWang2016II(4, ahash[k])
        res += [ADSummary(ana, desc=ana.h.getSource() + "-" + k)]
    ana.getBerchtold2014()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-all")]
    ana.getMarttinen2019()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Braak 0 vs")]
    ana.getMarttinen2019(2)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Braak 0-2 vs")]
    ana.getLu2014()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Ageing")]
    ana.getLu2004()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Ageing")]
    ana.getBerchtold2019(2)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-Activity")]
    #ana.getSrinivasan2020()
    #res += [ADSummary(ana, desc=ana.h.getSource() + "-fusiform gyrus")]
    ahash = {'myeloid':0, 'astrocyte':1, 'neuron':2, 'endothelial':3}
    for k in ahash:
        ana.getSrinivasan2020II(2, ahash[k])
        res += [ADSummary(ana, desc=ana.h.getSource() + "-" + k)]
    ana.getFriedman2018()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-fusiform gyrus")]
    ana.getMiyashita2014(2)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-ADC")]
    ana.getMcKay2019(2)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-ADC")]
    ana.getLow2021(2)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-BA9")]
    ana.getDunckley2006()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-EC")]
    ana.getBlalock2011(2)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-HIP")]
    # Set II
    ana.getNativio2020(2)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-HIP")]
    ana.getPascal2011(2, 0, 0)
    res += [ADSummary(ana, desc=ana.h.getSource() + "-EC")]
    ana.getNancy2020()
    res += [ADSummary(ana, desc=ana.h.getSource() + "-PFC")]
    ahash = {'unID':0, 'endo':1, 'oligo':2, 'astro':3, 'OPC':4, 'doublet':5,
             'neuron':6, 'mg':7}
    for k in ahash:
        ana.getChew2019(2, ahash[k])
        res += [ADSummary(ana, desc=ana.h.getSource() + "-" + k)]
    ana.getADPooledDyn(2, 0)
    res += [ADSummary(ana, desc="Pooled-no ADC")]
    ana.getADPooledDyn(2, 1)
    res += [ADSummary(ana, desc="Pooled-ADC")]
    ahash = {'TCx':0, 'PCx':1, 'HIP':2, 'FWM':3}
    for k in ahash:
        ana.getMiller2017(2, ahash[k])
        res += [ADSummary(ana, desc=ana.h.getSource() + "-" + k)]
    ana.getZhang2013("MAC39.1")
    res += [ADSummary(ana, desc=ana.h.getSource() + "-CR")]
    ana.getZhang2013("MAC39.2")
    res += [ADSummary(ana, desc=ana.h.getSource() + "-PFC")]
    ana.getZhang2013("MAC39.3")
    res += [ADSummary(ana, desc=ana.h.getSource() + "-VC")]
    cols = ['GSEID', 'n', 'Series', 'Condition']
    df = pd.DataFrame(res, columns=cols)
    return df
bone.ADDatasets = ADDatasets

def DP1(dfs): # VERTICAL
    sns = bone.sns
    plt = bone.plt

    if len(dfs) <= 0:
        return None
    df1 = dfs[0]
    df1['Name'] = list(df1['Signature'])
    df1['Xl'] = list(df1['#Genes'])
    labels = [k['GSEID'][0] + ' ' + k['Condition'][0] + \
            "(%d,%d)" % (k['#Cont'][0], k['#Expt'][0]) \
            for k in dfs]
    n1 = df1.shape[0]
    rocauc = list(df1['ROC-AUC'])
    p = list(df1['pvalue'])
    y = [1] * n1
    for i in range(1, len(dfs)):
        rocauc += list(dfs[i]['ROC-AUC'])
        p += list(dfs[i]['pvalue'])
        y += [i+1] * n1
    df = pd.DataFrame()
    df['ROC-AUC'] = rocauc
    df['pvalue'] = p
    df['ROC-AUC'] = df['ROC-AUC'].apply(
               lambda x: max([float(k) for k in str(x).split(",")]))
    df['pvalue'] = df['pvalue'].apply(
               lambda x: min([float(k) for k in str(x).split(",")]))
    df['Y'] = y
    df['R'] = [(i - 0.5) if i !=200  else 200  for i in df['ROC-AUC']]      #df['ROC-AUC'] - 0.5 
    df['Ra'] = [abs(i)+ 0.5  if i<1  else 0  for i in df['R'] ]     # abs(df['R']) + 0.5 
    df['Ra1'] = [abs(i)+ 0.5  if i<1  else 0.51 if i==0 else 0.48  for i in df['R'] ] 
    
    df['AUC'] = ['Up' if i > 0 else 'Down' for i in df['R']]
    #['Up' if i > 0 else 'Down' if i!=0 else 'NP' for i in df['R']]
    df['code'] = [bone.getCode(i) for i in df['pvalue']]
    df['ROC-AUC'] = df['Ra1']
    sns.set()
    sns.set_style("white")
    sns.set_style({'xtick.color':'.5', 'ytick.color':'.5', 'axes.labelcolor': '.5'})
    sns.set_context("notebook")
    sns.set_palette([bone.adj_light(c, 0.7, 1) for c in ['red', 'blue']])
    x = [i + 1 for i in range(n1)] * len(labels)
    y = df['Y']
    fig, ax = plt.subplots(figsize=(6, len(dfs)*0.5+1), dpi=100)
    ax = sns.scatterplot(x=x, y=y, size="ROC-AUC", hue='AUC',
                         sizes = (0, 100), size_norm = (0.49, 1),
                         hue_order = ['Up', 'Down'], ax=ax, data=df);
    roc = list(df['Ra'])
    code = list(df['code'])
    for line in range(n1):
        ax.text(line + 1, len(labels) + .5, df1['Xl'][line],
                horizontalalignment='center', size='small', color='0.8',
                verticalalignment='bottom', rotation=90)
        for i in range(len(labels)):
            ax.text(line + 1, i + 0.9, "%.2f" % roc[line + n1 * i],
                    horizontalalignment='right', size='small', color='0.8',
                    verticalalignment='top',  rotation=90)
#             ax.text(line + 1.5, i + 0.9, code[line + n1 * i],
#                     horizontalalignment='right', size='small', color='0.8',
#                     verticalalignment='top',  rotation=90)

    x1 = [i + 1 for i in range(n1)]
    ax.set_yticks(range(1, len(labels) + 1))
    ax.set_yticklabels(labels)
    ax.set_xlim([0, len(x1)+1])
    ax.set_ylim([0, len(labels) + 2])
    ax.set_xticks(x1)
    ax.set_xticklabels(df1['Name'], rotation=90)
    ax.set_ylabel("")
    ax.grid(False)
    handles, labels = ax.get_legend_handles_labels()
    legend_sizes = [0.5, 0.6, 0.8, 1.0]
    labels_m = labels[:4] + [f"{s}" for s in legend_sizes] 
    handles_m = handles[:4] + [plt.scatter([], [], s=(s-0.5)*200, color='black')
                       for s in legend_sizes]
    legend = ax.legend(handles_m, labels_m, bbox_to_anchor=(1.3, 1))
    legend.get_frame().set_edgecolor('black')

    return df,ax,fig
bone.DP1 = DP1

def getADPooledDyn(self, tn=1):
    self.prepareData("AD7")
    atype = self.h.getSurvName('c AD specific');
    ahash = {'0':0, '1':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Disease State');
    atypes = ['N', 'AD']
    ahash = {'Normal':0, "Alzheimer's Disease":1, 'normal':0,
            'definite AD':0, 'Control':0}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getADPooledDyn = getADPooledDyn


def getLiang2007(self):
    self.prepareData("AD2")
    atype = self.h.getSurvName('c Disease State');
    atypes = ['N', 'AD']
    ahash = {'normal\xa0':0, "Alzheimer's Disease\xa0":1}
    ahash = asciiNorm(ahash)
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLiang2007 = getLiang2007


def getFriedman2017(self):
    self.prepareData("AD3")
    atype = self.h.getSurvName('c diagnosis');
    atypes = ['N', 'AD']
    ahash = {'control':0, "Alzheimer's disease":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFriedman2017 = getFriedman2017


def getBerchtold2018(self, tn = 1):
    self.prepareData("AZ3", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c physical activity tier');
    atypes = ['H', 'M', 'L']
    ahash = {'high activity':0,
            'low activity':2,
            'moderate activity':1}
    if tn == 2:
        atypes = ['H', 'L']
        ahash = {'high activity':0,
                'low activity':1,
                'moderate activity':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerchtold2018 = getBerchtold2018


def getSimpson2015(self, tn = 1):
    self.prepareData("AZ5", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c neuronal ddr');
    atypes = ['H', 'L']
    ahash = {'High DNA damage':0, 'Low DNA damage':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSimpson2015 = getSimpson2015


def getPiras2019(self, tn = 1):
    self.prepareData("AZ1", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    age = self.h.getSurvName('c expired_age');
    atype = self.h.getSurvName('c diagnosis');
    atypes = ['N', 'AD']
    ahash = {'ND':0, 'AD':1}
    if tn == 2:
        atype = [atype[i] if int(age[i]) > 90
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPiras2019 = getPiras2019


def getNarayanan2014(self, tn = 1):
    self.prepareData("AZ8", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c disease status');
    atypes = ['N', 'AD', 'HD']
    ahash = {'non-demented':0,
            "Alzheimer's disease":1,
            "Huntington's disease":2}
    if tn == 2:
        atypes = ['N', 'AD']
        ahash = {'non-demented':0,
                "Alzheimer's disease":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNarayanan2014 = getNarayanan2014


def getZhang2013m(self, dbid = 'AZ12', tn = 1):
    self.db = hu.Database("/Users/mgosztyl/public_html/Hegemon/explore.conf")
    self.dbid = dbid
    atype = self.h.getSurvName('c disease');
    atypes = ['N', 'A']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2013m = getZhang2013m


def getBerchtold2014(self, tn = 1):
    self.prepareData("AZ11", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c src1")
    atype = [str(i) for i in atype]
    res = []
    for k in atype:
        l1 = k.split(",")
        if (len(l1) != 4):
            res.extend([k])
        else:
            res.extend([l1[1].strip() + "_" + l1[2].strip()])
    atype = res
    atypes = ['N', 'AD']
    ahash = {'entorhinal cortex_male':0,
            'entorhinal cortex_male_AD':1,
            'entorhinal cortex_female':0,
            'entorhinal cortex_female_AD':1,
            'superior frontal gyrus_male':0,
            'superior frontal gyrus_male_AD':1,
            'superior frontal gyrus_female':0,
            'superior frontal gyrus_female_AD':1,
            'postcentral gyrus_male':0,
            'post-central gyrus_male_AD':1,
            'postcentral gyrus_female':0,
            'post-central gyrus_female_AD':1,
            'hippocampus_male':0,
            'hippocampus_male_AD':1,
            'hippocampus_female':0,
            'hippocampus_female_AD':1}
    if (tn == 2):
        ahash = {'entorhinal cortex_male':0,
                'entorhinal cortex_male_AD':1,
                'superior frontal gyrus_male':0,
                'superior frontal gyrus_male_AD':1,
                'postcentral gyrus_male':0,
                'post-central gyrus_male_AD':1,
                'hippocampus_male':0,
                'hippocampus_male_AD':1}
    if (tn == 3):
        ahash = {'entorhinal cortex_female':0,
                'entorhinal cortex_female_AD':1,
                'superior frontal gyrus_female':0,
                'superior frontal gyrus_female_AD':1,
                'postcentral gyrus_female':0,
                'post-central gyrus_female_AD':1,
                'hippocampus_female':0,
                'hippocampus_female_AD':1}
    if (tn == 4):
        ahash = {'entorhinal cortex_male':0,
                'entorhinal cortex_male_AD':1,
                'entorhinal cortex_female':0,
                'entorhinal cortex_female_AD':1}
    if (tn == 5):
        ahash = {'superior frontal gyrus_male':0,
                'superior frontal gyrus_male_AD':1,
                'superior frontal gyrus_female':0,
                'superior frontal gyrus_female_AD':1}
    if (tn == 6):
        ahash = {'postcentral gyrus_male':0,
                'post-central gyrus_male_AD':1,
                'postcentral gyrus_female':0,
                'post-central gyrus_female_AD':1}
    if (tn == 7):
        ahash = {'hippocampus_male':0,
                'hippocampus_male_AD':1,
                'hippocampus_female':0,
                'hippocampus_female_AD':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerchtold2014 = getBerchtold2014


def getWang2016(self, dbid="AD9", tn = 1):
    self.dbid = dbid
    atype = self.h.getSurvName('c brain region')
    ahash = {'Inferior Temporal Gyrus':3,
            'Parahippocampal Gyrus':4,
            'Middle Temporal Gyrus':5,
            'Occipital Visual Cortex':6,
            'Prefrontal Cortex':7,
            'Hippocampus':8,
            'Caudate Nucleus':9,
            'Frontal Pole':10,
            'Precentral Gyrus':11,
            'Posterior Cingulate Cortex':12,
            'Superior Temporal Gyrus':13,
            'Superior Parietal Lobule':14,
            'Temporal Pole':15,
            'Anterior Cingulate':16,
            'Inferior Frontal Gyrus':17,
            'Dorsolateral Prefrontal Cortex':18,
            'Putamen':19}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c neuropathological category')
    atypes = ['N', 'AD']
    ahash = {'definite AD':1, 'Possible AD':1,
            'Normal':0, 'Probable AD':1}
    if (tn >= 2):
        atypes = ['N', 'AD']
        ahash = {'definite AD':1, 'Normal':0}
    if (tn >= 3):
        atype = [atype[i] if rval[i] == tn else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2016 = getWang2016


def getBerchtold2014RMA(self, tn=1):
    self.prepareData("AD11")
    atype = self.h.getSurvName('c AD specific');
    ahash = {'0':0, '1':1}
    rval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName("c src1")
    atype = [str(i) for i in atype]
    res = []
    for k in atype:
        l1 = k.split(",")
        if (len(l1) != 4):
            res.extend([k])
        else:
            res.extend([l1[1].strip() + "_" + l1[2].strip()])
    atype = res
    atypes = ['N', 'AD']
    ahash = {'entorhinal cortex_male':0,
            'entorhinal cortex_male_AD':1,
            'entorhinal cortex_female':0,
            'entorhinal cortex_female_AD':1,
            'superior frontal gyrus_male':0,
            'superior frontal gyrus_male_AD':1,
            'superior frontal gyrus_female':0,
            'superior frontal gyrus_female_AD':1,
            'postcentral gyrus_male':0,
            'post-central gyrus_male_AD':1,
            'postcentral gyrus_female':0,
            'post-central gyrus_female_AD':1,
            'hippocampus_male':0,
            'hippocampus_male_AD':1,
            'hippocampus_female':0,
            'hippocampus_female_AD':1}
    if (tn == 4):
        ahash = {'entorhinal cortex_male':0,
                'entorhinal cortex_male_AD':1,
                'superior frontal gyrus_male':0,
                'superior frontal gyrus_male_AD':1,
                'postcentral gyrus_male':0,
                'post-central gyrus_male_AD':1,
                'hippocampus_male':0,
                'hippocampus_male_AD':1}
    if (tn == 5):
        ahash = {'entorhinal cortex_female':0,
                'entorhinal cortex_female_AD':1,
                'superior frontal gyrus_female':0,
                'superior frontal gyrus_female_AD':1,
                'postcentral gyrus_female':0,
                'post-central gyrus_female_AD':1,
                'hippocampus_female':0,
                'hippocampus_female_AD':1}
    if (tn == 6):
        ahash = {'entorhinal cortex_male':0,
                'entorhinal cortex_male_AD':1,
                'entorhinal cortex_female':0,
                'entorhinal cortex_female_AD':1}
    if (tn == 7):
        ahash = {'superior frontal gyrus_male':0,
                'superior frontal gyrus_male_AD':1,
                'superior frontal gyrus_female':0,
                'superior frontal gyrus_female_AD':1}
    if (tn == 8):
        ahash = {'postcentral gyrus_male':0,
                'post-central gyrus_male_AD':1,
                'postcentral gyrus_female':0,
                'post-central gyrus_female_AD':1}
    if (tn == 9):
        ahash = {'hippocampus_male':0,
                'hippocampus_male_AD':1,
                'hippocampus_female':0,
                'hippocampus_female_AD':1}
    if (tn == 2):
        atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerchtold2014RMA = getBerchtold2014RMA


def getGelman2012(self, tn = 1):
    self.prepareData("DE39", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c phenotype');
    atypes = ['N', 'HAND']
    ahash = {'HIV Infected':0,
            'HIV Infected with neurocognitive impairment (HAD: HIV-associated dementia)':1,
            'HIV Infected with HAD and HIV encephalitis (HIVE)':1,
            'normal (control)':0,
            'HIV Infected with HAD':1,
            'HIV Infected with HAD and encephalitis':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGelman2012 = getGelman2012

    
def getChenPlotkin2008(self, tn = 1):
    self.prepareData("DE1", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName("c src1")
    atype = [str(k).split(" ")[1] if len(str(k).split(" ")) > 1 else None
            for k in atype]
    atype = [str(k).split("-")[0] for k in atype]
    atypes = ['N', 'P', 'S']
    ahash = {'Normal':0,
            'Progranulin':1,
            'Sporadic':2}
    if tn == 2:
        atypes = ['N', 'FTD']
        ahash = {'Normal':0,
                'Progranulin':1}
    if tn == 3:
        atype = self.h.getSurvName('c disease and tissue')
        atypes = ['N', 'FTD']
        ahash = {'Normal hippocampus':0,
            'Progranulin hippocampus':1,
            'Sporadic hippocampus':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChenPlotkin2008 = getChenPlotkin2008

    
def getOlmosSerrano2016(self, tn = 1):
    self.prepareData("DE49", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c disease status');
    atypes = ['N', 'DS']
    ahash = {'CTL':0, 'DS':1}
    if tn == 2:
        atype = self.h.getSurvName('c disease and tissue');
        atypes = ['N', 'DS']
        ahash = {'CTL ITC':0, 'CTL STC':0, 'CTL HIP':0,
            'DS ITC':1, 'DS STC':1, 'DS HIP':1}
    if tn == 3:
        atype = self.h.getSurvName('c disease and tissue');
        atypes = ['N', 'DS']
        ahash = {'CTL CBC':0, 'DS CBC':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOlmosSerrano2016 = getOlmosSerrano2016

    
def getWilliams2009(self, tn = 1):
    self.prepareData("DE51", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c Disease State');
    atypes = ['N', 'MCI']
    ahash = {'Control':0, 'MCI':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWilliams2009 = getWilliams2009

    
def getBartolettiStella2019 (self, tn = 1):
    self.prepareData("DE52", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c condition');
    atypes = ['N', 'CJD']
    ahash = {'Control':0, 'sCJD affected':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBartolettiStella2019  = getBartolettiStella2019 

    
def getWes2014a (self, tn = 1):
    self.prepareData("AZ78", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genotype');
    atypes = ['N', 'rTg4510']
    ahash = {'Dbl Neg':0, 'Tg4510':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWes2014a  = getWes2014a 


def getWes2014b (self, tn = 1):
    self.prepareData("AZ77", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genotype');
    atypes = ['N', 'rTg4510']
    ahash = {'Dbl Neg':0, 'Tg4510':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWes2014b  = getWes2014b 

    
def getWes2014c (self, tn = 1):
    self.prepareData("AZ78", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genotype');
    atypes = ['N', 'rTg4510']
    ahash = {'Dbl Neg':0, 'Tg4510':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWes2014c  = getWes2014c 

    
def getHokama2013 (self, tn = 1):
    self.prepareData("AZ49", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genotype');
    atypes = ['N', '3xTg']
    ahash = {'non-Tg':0, '3xTg-Homo':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHokama2013  = getHokama2013 

    
def getWakutani2014 (self, tn = 1):
    self.prepareData("AZ60", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genetic background');
    atypes = ['N', 'TgCRND8']
    ahash = {'non-transgenic littermate mouse':0,
            'TgCRND8 transgenic mouse':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWakutani2014  = getWakutani2014 

    
def getMeyer2019 (self, dbid = "AZ29", tn = 1):
    self.db = hu.Database("/Users/mgosztyl/public_html/Hegemon/explore.conf")
    self.dbid = dbid
    atype = self.h.getSurvName('c diagnosis');
    atypes = ['N', 'AD']
    ahash = {'normal':0,
            'sporadic Alzheimer\'s disease':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMeyer2019  = getMeyer2019 

    
def getScheckel2019 (self, tn = 1):
    self.prepareData("AZ95", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c source_name (ch1)');
    atypes = ['N', 'AD']
    ahash = {'wiltype iPSC-derived neurons':0,
            'APP/PSEN1 double mutant iPSC-derived neurons':1}
    if tn == 2:
        ahash = {'wiltype iPSC-derived neurons':0,
            'APP mutant iPSC-derived neurons':1}
    if tn == 3:
        ahash = {'wiltype iPSC-derived neurons':0,
            'PSEN1 mutant iPSC-derived neurons':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getScheckel2019  = getScheckel2019 


def getPatel2019(self, tn=1, ta=0, tb=0):
    self.prepareData("AD8")
    atype = self.h.getSurvName('c tissue')
    ahash = {'Temporal_Cortex':1, 'Cerebellum':3,
             'Frontal_Cortex':2, 'Entorhinal_Cortex':0}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease state')
    atypes = ['N', 'A', "AD"]
    ahash = {'AsymAD':1, 'AD':2, 'control':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
        atypes = ['N', "AD"]
        ahash = {'AD':1, 'control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPatel2019 = getPatel2019


def getFriedman2017(self, tn=1, ta=0, tb=0):
    self.prepareData("AD3")
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['N', "AD"]
    ahash = {'control':0, "Alzheimer's disease":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFriedman2017 = getFriedman2017


def getLiang2007(self, tn=1, ta=0, tb=0):
    self.prepareData("AD2")
    atype = self.h.getSurvName('c organ region')
    atype1 = self.h.getSurvName('c Organ Region')
    atype = [atype[k] if atype[k] != '' else atype1[k]
                    for k in range(len(atype))]
    ahash = {'Medial Temporal Gyrus\xa0':0, 'hippocampus\xa0':1,
            'Superior Frontal Gyrus\xa0':2, 'Primary Visual Cortex\xa0':3,
            'Posterior Singulate\xa0':4, 'Entorhinal Cortex\xa0':5,
            'Posterior Cingulate\xa0':6, 'Posterior cingulate\xa0':7}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Disease State')
    atype1 = self.h.getSurvName('c disease state')
    atype = [atype[k] if atype[k] != '' else atype1[k]
            for k in range(len(atype))]
    atypes = ['N', "AD"]
    ahash = {'normal\xa0':0, "Alzheimer's Disease\xa0":1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        ah = {4:1, 6:1, 7:1}
        atype = [atype[i] if tval[i] in ah
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLiang2007 = getLiang2007


def getWang2016(self, tn=1, ta=0, tb=0):
    self.prepareData("AD9")
    atype = self.h.getSurvName('c brain region')
    ahash = {'Amygdala':0, 'Nucleus Accumbens':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c neuropathological category')
    atypes = ['N', "AD", 'pS', 'pB']
    ahash = {'definite AD':1, 'Possible AD':2, 'Normal':0, 'Probable AD':3}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['N', "AD"]
        ahash = {'definite AD':1, 'Normal':0}
    if (tn == 4):
        atypes = ['N', "AD"]
        ahash = {'definite AD':1, 'Normal':0}
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2016 = getWang2016


def getWang2016II(self, tn=1, ta=0, tb=0):
    self.prepareData("AD10")
    atype = self.h.getSurvName('c brain region')
    ahash = {'Inferior Temporal Gyrus':0, 'Parahippocampal Gyrus':1,
            'Middle Temporal Gyrus':2, 'Occipital Visual Cortex':3,
            'Prefrontal Cortex':4, 'Hippocampus':5, 'Caudate Nucleus':6,
            'Frontal Pole':7, 'Precentral Gyrus':8,
            'Posterior Cingulate Cortex':9, 'Superior Temporal Gyrus':10,
            'Superior Parietal Lobule':11, 'Temporal Pole':12,
            'Anterior Cingulate':13, 'Inferior Frontal Gyrus':14,
            'Dorsolateral Prefrontal Cortex':15, 'Putamen':16}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c neuropathological category')
    atypes = ['N', "AD", 'pS', 'pB']
    ahash = {'Possible AD':2, 'definite AD':1, 'Normal':0, 'Probable AD':3}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        atypes = ['N', "AD"]
        ahash = {'definite AD':1, 'Normal':0}
    if (tn == 4):
        atypes = ['N', "AD"]
        ahash = {'definite AD':1, 'Normal':0}
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2016II = getWang2016II


def getBerchtold2014(self, tn=1, ta=0, tb=0):
    self.prepareData("AD11")
    atype = self.h.getSurvName('c brain region')
    ahash = {'entorhinal cortex':0, 'postcentral gyrus':1, 'hippocampus':2,
            'superior frontal gyrus':3, 'post-central gyrus':4}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c AD specific')
    atypes = ['N', "AD"]
    ahash = {'0':0, '1':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerchtold2014 = getBerchtold2014


def getMarttinen2019(self, tn=1, ta=0, tb=0):
    self.prepareData("AD12")
    atype = self.h.getSurvName('c braak stage')
    atypes = ['N', "AD"]
    ahash = {'3':1, '5':1, '2':1, '6':1, '1':1, '0':0, '4':1}
    if (tn == 2):
        ahash = {'3':1, '5':1, '2':0, '6':1, '1':0, '0':0, '4':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMarttinen2019 = getMarttinen2019


def getWebster2009(self, tn=1, ta=0, tb=0):
    self.prepareData("AD14")
    atype = self.h.getSurvName('c src1')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['N', "AD"]
    ahash = { "Alzheimer's":1, 'neuropathologically':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWebster2009 = getWebster2009


def getLu2014(self, tn=1, ta=0, tb=0):
    self.prepareData("AD16")
    atype = self.h.getSurvName('c age')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['N', "AD"]
    ahash = {}
    for k in atype[2:]:
        if float(k) < 80:
            ahash[k] = 0
        else:
            ahash[k] = 1
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLu2014 = getLu2014


def getLu2004(self, tn=1, ta=0, tb=0):
    self.prepareData("AD18")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(" .*", "", str(k)) for k in atype]
    atypes = ['N', "AD"]
    ahash = {}
    for k in atype[2:]:
        if k == 'another':
            continue
        if float(k) < 80:
            ahash[k] = 0
        else:
            ahash[k] = 1
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLu2004 = getLu2004


def getBerchtold2019(self, tn=1, ta=0, tb=0):
    self.prepareData("AD20")
    atype = self.h.getSurvName('c physical activity tier')
    atypes = ['H', 'M', 'L']
    ahash = {'moderate activity':1, 'high activity':0, 'low activity':2}
    if tn == 2:
        atypes = ['H', 'L']
        ahash = {'high activity':0, 'low activity':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBerchtold2019 = getBerchtold2019


def getSrinivasan2020(self, tn=1, ta=0, tb=0):
    self.prepareData("AD21")
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['N', 'AD']
    ahash = {"Alzheimer's disease":1, 'control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSrinivasan2020 = getSrinivasan2020


def getSrinivasan2020II(self, tn=1, ta=0, tb=0):
    self.prepareData("AD22")
    atype = self.h.getSurvName('c cell type')
    ahash = {'myeloid':0, 'astrocyte':1, 'neuron':2, 'endothelial':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['N', 'AD']
    ahash = {'Control':0, 'AD':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName('c cell type')
        atypes = ['neuron', 'myeloid', 'astrocyte', 'endothelial']
        ahash = {}
    if (tn == 4):
        atype = self.h.getSurvName('c cell type')
        btype = self.h.getSurvName('c diagnosis')
        atype = [ f"{atype[k]}-{btype[k]}" for k in range(len(atype))]
        atypes = ['M-H', 'M-AD', 'A-H', 'A-AD', 'E-H', 'E-AD',
                'N-H', 'N-AD']
        ahash = {'myeloid-Control':0, 'myeloid-AD':1,
                'astrocyte-Control':2, 'astrocyte-AD':3,
                'endothelial-Control':4, 'endothelial-AD':5,
                'neuron-Control':6, 'neuron-AD':7}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSrinivasan2020II = getSrinivasan2020II


def getFriedman2018(self, tn=1, ta=0, tb=0):
    self.prepareData("AD23")
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['N', 'AD']
    ahash = {"Alzheimer's disease":1, 'control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getFriedman2018 = getFriedman2018


def getRodriguez2021(self, tn=1, ta=0, tb=0):
    self.prepareData("AD24")
    atype = self.h.getSurvName('c stimulant')
    ahash = {'':0, 'dsRNA + lipofectamine':3, 'lipopolysaccharide':4,
              'dsRNA':1, 'lipofectamine':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    dval = self.h.getSurvName('c drug_treatment')
    atype = self.h.getSurvName('c drug_concentration_um')
    atypes = ['C', '0.3', '1', '3', '10']
    ahash = {'':0, '0.3':1, '1':2, '3':3, '10':4}
    if (tn == 2):
        atype = ['' if (dval[i] == '' or dval[i] == 'dmso')
                 else atype[i] for i in range(len(atype))]            
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
        atype = [atype[i] if (dval[i] == '' or dval[i] == 'dmso'
                             or dval[i] == tb)
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRodriguez2021 = getRodriguez2021


def getNiculescu2020(self, tn=1, ta=0, tb=0):
    self.prepareData("AD25")
    atype = self.h.getSurvName('c neuropsych test results')
    atypes = ['N', 'MCI', 'AD', 'O']
    ahash = {'':0, '3 = other cognitive disorder but not MCI or ADRD':3,
             '1 = MCI':1, '2 = ADRD':2}
    if (tn == 2):
        atypes = ['N', 'AD']
        ahash = {'':0, '2 = ADRD':1}
    if (tn == 3):
        atypes = ['N', 'MCI']
        ahash = {'':0, '1 = MCI':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNiculescu2020 = getNiculescu2020


def getSood2015I(self, tn=1, ta=0, tb=0):
    self.prepareData("AD26")
    atype = self.h.getSurvName('c status')
    atypes = ['CTL', 'MCI', 'AD']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSood2015I = getSood2015I


def getSood2015II(self, tn=1, ta=0, tb=0):
    self.prepareData("AD27")
    atype = self.h.getSurvName('c status')
    atypes = ['CTL', 'MCI', 'AD', 'MCL', 'O', 'B', 'CA']
    ahash = {'MCI to CTL':3, 'OTHER':4, 'borderline MCI':5, 'CTL to AD':6}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSood2015II = getSood2015II


def getMiyashita2014(self, tn=1, ta=0, tb=0):
    self.prepareData("AD28")
    atype = self.h.getSurvName('c braak nft stage')
    atypes = ['0', 'I-II', 'III-IV', 'V-VI']
    ahash = {}
    if (tn == 2):
        atypes = ['0', 'V-VI']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMiyashita2014 = getMiyashita2014


def getSamsudin2016(self, tn=1, ta=0, tb=0):
    self.prepareData("AD29")
    atype = self.h.getSurvName('c diagnosis')
    atypes = ['C', 'AD']
    ahash = {"Probable Alzheimer's Disease":1, 'Non-demented control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSamsudin2016 = getSamsudin2016


def getMcKay2019(self, tn=1, ta=0, tb=0):
    self.prepareData("AD30")
    atype = self.h.getSurvName('c patient diagnosis')
    atypes = ['C', 'AD', 'VD']
    ahash = {'Control':0, 'Vascular dementia':2, "Alzheimer's disease":1}
    if (tn == 2):
        atypes = ['C', 'AD']
        ahash = {'Control':0, "Alzheimer's disease":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMcKay2019 = getMcKay2019


def getLow2021(self, tn=1, ta=0, tb=0):
    self.prepareData("AD31")
    atype = self.h.getSurvName('c desc')
    atypes = ['CTRL', 'DLB', 'AD', 'PDD']
    ahash = {}
    if (tn == 2):
        atypes = ['CTRL', 'AD']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLow2021 = getLow2021


def getDong2013(self, tn=1, ta=0, tb=0):
    self.prepareData("AD32")
    atype = self.h.getSurvName('c src1')
    atype = [str(k).split(" ")[2] if len(str(k).split(" ")) > 2 else k
            for k in atype]
    atypes = ['E', 'O']
    ahash = {'endogenous':0, 'overexpressed':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDong2013 = getDong2013


def getHurley2012(self, tn=1, ta=0, tb=0):
    self.prepareData("AD33")
    atype = self.h.getSurvName('c Title')
    atypes = ['C', 'KD']
    ahash = {'TCEB1_KD':0, 'ARHGAP1_KD':0, 'CDC42EP1_KD':0, 'GPR4_KD':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getHurley2012 = getHurley2012


def getPatella2020(self, tn=1, ta=0, tb=0):
    self.prepareData("AD34")
    atype = self.h.getSurvName('c treatment')
    atypes = ['C', 'EBM2', 'iVWF', 'M199', 'TNFa']
    ahash = {'nucleofected with si-VWF':2,
            'untreated, grown in EBM2_':1,
             'mock nucleofected':0,
             'untreated, grown in M199_':3,
             'treated with 50ng/ml TNFa for 3h, grown in M199_':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPatella2020 = getPatella2020


def getAvrampou2019(self, tn=1, ta=0, tb=0):
    self.prepareData("AD35")
    atype = self.h.getSurvName('c genotype')
    atypes = ['WT', 'RGS4 KO']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAvrampou2019 = getAvrampou2019


def getOlmosSerrano2016(self, tn = 1, ta = 0):
    self.prepareData("DE49", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c region');
    ahash = {'DFC':0, 'S1C':1, 'OFC':2, 'V1C':3, 'ITC':4, 'MFC':5,
             'CBC':6, 'VFC':7, 'HIP':8, 'FC':9, 'STC':10, 'IPC':11}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease status');
    atypes = ['N', 'DS']
    ahash = {'CTL':0, 'DS':1}
    if tn == 2:
        atype = self.h.getSurvName('c disease and tissue');
        atypes = ['N', 'DS']
        ahash = {'CTL ITC':0, 'CTL STC':0, 'CTL HIP':0,
                 'DS ITC':1, 'DS STC':1, 'DS HIP':1}
    if tn == 3:
        atype = self.h.getSurvName('c disease and tissue');
        atypes = ['N', 'DS']
        ahash = {'CTL V1C':0, 'DS V1C':1,}
    if tn == 4:
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getOlmosSerrano2016 = getOlmosSerrano2016


def getWilliams2009(self, tn = 1):
    self.prepareData("DE51", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c Disease State');
    atypes = ['N', 'MCI']
    ahash = {'Control':0, 'MCI':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWilliams2009 = getWilliams2009


def getAxol2021(self, tn = 1):
    self.prepareData("AD36")
    atype = self.h.getSurvName('c Cell Type');
    ahash = {'NSC':0, 'Str':1, 'iPSCs':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Disease State');
    atypes = ['C', 'AD']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if tval[i] == 0
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAxol2021 = getAxol2021


def getDunckley2006(self, tn = 1):
    self.prepareData("AD44")
    atype = self.h.getSurvName('c brain, Entorhinal Cortex')
    atype = [re.sub(".* ", "", str(k)) for k in atype]
    atypes = ['C', 'AD']
    ahash = {'tangle_e1_le1':1, 'normal_e1_le1':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getDunckley2006 = getDunckley2006


def getBlalock2011(self, tn = 1):
    self.prepareData("AD45")
    atype = self.h.getSurvName('c disease status')
    atypes = ['C', 'I', 'M', 'S']
    ahash = {'Moderate':2, 'Severe':3, 'Incipient':1, 'Control':0}
    if tn == 2:
        atypes = ['C', 'AD']
        ahash = {'Moderate':1, 'Severe':1, 'Incipient':1, 'Control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBlalock2011 = getBlalock2011


def getQian2020(self, tn = 1):
    self.prepareData("AD47")
    atype = self.h.getSurvName('c genetic manipulation')
    atypes = ['C', 'shPTB']
    ahash = {'shRNA against PTB':1, 'Mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getQian2020 = getQian2020


def getGomezCarballa2023(self, tn = 1):
    self.prepareData("MACV394")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("; 22.*", "", str(k)) for k in atype]
    atypes = ['C1', 'C2', 'A1', 'A2']
    ahash = {'ACD; Timepoint-1':2, 'ACD; Timepoint-2':3,
             'Control; Timepoint-1':0, 'Control; Timepoint-2':1}
    if (tn == 2):
        atypes = ['C', 'A']
        ahash = {'ACD; Timepoint-1':1, 'ACD; Timepoint-2':1,
             'Control; Timepoint-1':0, 'Control; Timepoint-2':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGomezCarballa2023 = getGomezCarballa2023


def getGomezCarballa2023II(self, tn = 1):
    self.prepareData("MACV394.2")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub("; 22.*", "", str(k)) for k in atype]
    atypes = ['C1', 'C2', 'A1', 'A2']
    ahash = {'ACD; Timepoint-1':2, 'ACD; Timepoint-2':3,
             'Control; Timepoint-1':0, 'Control; Timepoint-2':1}
    if (tn == 2):
        atypes = ['C', 'A']
        ahash = {'ACD; Timepoint-1':1, 'ACD; Timepoint-2':1,
             'Control; Timepoint-1':0, 'Control; Timepoint-2':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGomezCarballa2023II = getGomezCarballa2023II


def getGriggs2023ad(self, tn = 1):
    self.prepareData("MACV395")
    atype = self.h.getSurvName('c Title')
    atype = [re.sub(", sub.*region ", "_", str(k)) for k in atype]
    atypes = ['H', '9', 'AH', 'A9', 'CH', 'C9', 'ACH', 'AC9']
    ahash = {'CONTROL_ADpos_BA9':3, 'COVID_ADneg_BA9':5,
             'CONTROL_ADpos_HP':2, 'COVID_ADpos_BA9':7,
             'CONTROL_ADneg_HP':0, 'CONTROL_ADneg_BA9':1,
             'COVID_ADpos_HP':6, 'COVID_ADneg_HP':4}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGriggs2023ad = getGriggs2023ad


def getGate2024ADblood(self, tn = 1):
    self.prepareData("MACV396")
    atype = self.h.getSurvName('c disease state')
    atypes = ['H', 'AD']
    ahash = {'Healthy Control':0, 'Alzheimers Disease':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGate2024ADblood = getGate2024ADblood

def getNativio2020(self, tn=1, ta=0, tb=0):
    self.prepareData("ad59")
    atype = self.h.getSurvName('c title2')
    atypes = ['Y', 'O', "AD"]
    ahash = {'AD [RNA-seq]':2, 'Old [RNA-seq]':1, 'Young [RNA-seq]':0}
    if (tn == 2):
        atypes = ['O', "AD"]
        ahash = {'AD [RNA-seq]':1, 'Old [RNA-seq]':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNativio2020 = getNativio2020

def getPascal2011(self, tn=1, ta=0, tb=0):
    self.prepareData("ad63")
    atype = self.h.getSurvName('c gender (ch1)')
    ahash = {'F':0, 'M':1}
    gval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease (ch1)')
    ahash = {"Alzheimer's disease":0, 'Amyotrophic lateral sclerosis':1,
             "Huntington's disease":2, 'Multiple sclerosis':3,
             "Parkinson's disease":4, 'Schizophrenia':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease2 (ch1)')
    atypes = ['C', 'D']
    ahash = {'AD_Control':0, 'AD_Disease':1, 'ALS_Control':0, 'ALS_Disease':1,
             'HD_Control':0, 'HD_Disease':1, 'MS_Control':0, 'MS_Disease':1,
             'PD_Control':0, 'PD_Disease':1, 'SCHIZ_Control':0, 'SCHIZ_Disease':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta and gval[i] == tb
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getPascal2011 = getPascal2011
def getNancy2020(self, tn=1, ta=0, tb=0):
    self.prepareData("ad69")
    atype = self.h.getSurvName('c diagnosis (ch1)')
    atypes = ['C', "AD"]
    ahash = {'AD':1, 'healthy control':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getNancy2020 = getNancy2020

def getChew2019(self, tn=1, ta=0, tb=0):
    self.prepareData("ad68.3")
    atype=self.h.getSurvName('c CellType')
    ahash = {'unID':0, 'endo':1, 'oligo':2, 'astro':3, 'OPC':4, 'doublet':5,
             'neuron':6, 'mg':7}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease state (ch1)')
    atypes = ['C', "AD"]
    ahash = {'control':0, 'AD':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        atype=self.h.getSurvName('c CellType')
        atypes = ['unID', 'endo', 'oligo', 'astro', 'OPC', 'doublet',
                 'neuron', 'mg']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getChew2019 = getChew2019

def getSwarup2021(self, tn=1, ta=0, tb=0):
    self.prepareData("ad54")
    atype=self.h.getSurvName('c diagnosis (ch1)')
    atypes = ['C', "AD"]
    ahash = {'Control':0, 'AD':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSwarup2021 = getSwarup2021

def getADPooledDyn(self, tn=1, ta=0, tb=0):
    self.prepareData("AD7")
    atype = self.h.getSurvName('c AD specific');
    ahash = {'0':0, '1':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Disease State');
    atypes = ['N', 'AD']
    ahash = {'Normal':0, "Alzheimer's Disease":1, 'normal':0,
            'definite AD':0, 'Control':0}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getADPooledDyn = getADPooledDyn

def getMiller2017(self, tn=1, ta=0, tb=0):
    self.prepareData("ad6")
    atype = self.h.getSurvName('c hemisphere (ch1)');
    ahash = {'left':0, 'right':1}
    sval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c brain region (ch1)');
    ahash = {'TCx':0, 'PCx':1, 'HIP':2, 'FWM':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c new_nincds_arda_diagnosis');
    atypes = ['No Dementia', "AD" ]
    ahash = {'No Dementia':0, "Possible Alzheimer'S Disease":1,
            "Probable Alzheimer'S Disease":1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        atype = [atype[i] if tval[i] == ta and sval[i] == tb
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMiller2017 = getMiller2017

def getYang2020(self, tn=1, ta=0, tb=0):
    self.prepareData("ad103")
    atype = self.h.getSurvName('c brain region (ch1)');
    ahash = {'Superior frontal cortex':0, 'Hippocampus':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease (ch1)');
    atypes = ['Control','AD']
    ahash = {'Control':0, "AD":1}             
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYang2020 = getYang2020

def getZhang2013(self, dbid="MAC39.1"):
    self.prepareData(dbid)
    atype = self.h.getSurvName("c disease")
    atypes = ['A', 'N']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2013 = getZhang2013

def getKondo2013(self, tn = 1):
    self.prepareData("AZ27", "/Users/mgosztyl/public_html/Hegemon/explore.conf")
    atype = self.h.getSurvName('c genotype/variation');
    atypes = ['Wt','APP-IPSC']
    ahash = {'APP wild-type':0, 'APP E693delta':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getKondo2013 = getKondo2013
 
def getSun2023(self, tn = 1, ta = 0):
    self.prepareData("AD56")
    atype = self.h.getSurvName('c brain_region');
    ahash = {'Hippocampus':0, 'Prefrontal_cortex':1, 'Angular_gyrus':2,
            'Anterior_thalamus':3, 'Entorhinal_cortex':4,
            'Midtemporal_cortex':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c ADdiag2types');
    atypes = ['nonAD','AD']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSun2023 = getSun2023
 
def getSun2023II(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD56.2")
    atype = self.h.getSurvName('c celltype');
    ahash = {'Endo':0, 'Ependymal':1, 'Fib':2, 'Per':3, 'SMC':4}
    ctval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c brain_region');
    ahash = {'Hippocampus':0, 'Prefrontal_cortex':1, 'Angular_gyrus':2,
            'Anterior_thalamus':3, 'Entorhinal_cortex':4,
            'Midtemporal_cortex':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c ADdiag2types');
    atypes = ['nonAD','AD']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta and ctval[i] == tb
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSun2023II = getSun2023II
 
def getSun2023III(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD56.3")
    atype = self.h.getSurvName('c treatment');
    atypes = ['C','AB']
    ahash = {'Control':0, 'Amyloid-beta':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSun2023III = getSun2023III
 
def getSun2023IV(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD56.4")
    atype = self.h.getSurvName('c brainRegion');
    ahash = {'Hippocampus':0, 'PFC':1, 'AngularGyrus':2,
            'EntorhinalCortex':3, 'MidtemporalCortex':4, 'Thalamus':5}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c ADdiag3types');
    atypes = ['nonAD','early', 'late']
    ahash = {'lateAD':2, 'nonAD':0, 'earlyAD':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSun2023IV = getSun2023IV
 
def getGerrits2021(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD57")
    atype = self.h.getSurvName('c Type');
    ahash = {'neun':0, 'olig':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c clinical diagnosis');
    atypes = ['CTR','AD']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = self.h.getSurvName('c sample group');
        atypes = ['CTR', 'CTR+', 'AD']
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = self.h.getSurvName('c Type');
        atypes = ['neun', 'olig']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGerrits2021 = getGerrits2021
 
def getGerrits2021II(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD57.2")
    atype = self.h.getSurvName('c region');
    ahash = {'occipitotemporal cortex':0, 'occipital cortex':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c clinical diagnosis');
    atypes = ['CTR','AD']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = self.h.getSurvName('c sample group');
        atypes = ['CTR', 'CTR+', 'AD']
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = self.h.getSurvName('c sample group');
        atypes = ['CTR', 'AD']
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGerrits2021II = getGerrits2021II
 
def getGerrits2021III(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD57.3")
    atype = self.h.getSurvName('c Type');
    ahash = {'M':0, 'B':1}
    ctval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c region');
    ahash = {'occipitotemporal cortex':0, 'occipital cortex':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c clinical diagnosis');
    atypes = ['CTR','AD']
    ahash = {}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta and ctval[i] == tb
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = self.h.getSurvName('c sample group');
        atypes = ['CTR', 'CTR+', 'AD']
        atype = [atype[i] if tval[i] == ta and ctval[i] == tb
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = self.h.getSurvName('c sample group');
        atypes = ['CTR', 'AD']
        atype = [atype[i] if tval[i] == ta and ctval[i] == tb
                 else None for i in range(len(atype))]
    if (tn == 5):
        atype = [atype[i] if ctval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGerrits2021III = getGerrits2021III
 
def getGerrits2022(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD61")
    atype = self.h.getSurvName('c brain region');
    ahash = {'Frontal Cortex':0, 'Occipital Cortex':1, 'Temporal Cortex':2}
    ctval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Type');
    ahash = {'neun':0, 'olig':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c donor');
    atypes = ['C','FTD']
    ahash = {'C2':0, 'C3':0, 'C4':0, 'C5':0,
            'P1':1, 'P2':1, 'P3':1, 'P4':1, 'P5':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = [atype[i] if tval[i] == ta and ctval[i] == tb
                 else None for i in range(len(atype))]
    if tn == 4:
        atype = self.h.getSurvName('c Type');
        atypes = ['neun', 'olig']
    if tn == 5:
        atype = self.h.getSurvName('c brain region');
        atypes = ['FC', 'OC', 'TC']
        ahash = {'Frontal Cortex':0, 
                'Occipital Cortex':1, 'Temporal Cortex':2}
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getGerrits2022 = getGerrits2022
 
def getSmith2022(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD58")
    atype = self.h.getSurvName('c brain region (somatosensory or entorhinal cortex)');
    ahash = {'EC':0, 'SSC':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease state');
    atypes = ['C','AD']
    ahash = {'Non-disease control':0, 'AD':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName('c brain region (somatosensory or entorhinal cortex)');
        atypes = ['EC', 'SSC']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSmith2022 = getSmith2022
 
def getSmith2022II(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD58.2")
    atype = self.h.getSurvName('c Type');
    ahash = {'M':0, 'B':1}
    ctval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c brain region (somatosensory or entorhinal cortex)');
    ahash = {'EC':0, 'SSC':1}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c disease state');
    atypes = ['C','AD']
    ahash = {'Non-disease control':0, 'AD':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta and ctval[i] == tb
                 else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName('c brain region (somatosensory or entorhinal cortex)');
        atypes = ['EC', 'SSC']
        ahash = {}
    if (tn == 4):
        atype = self.h.getSurvName('c Type');
        atypes = ['M', 'B']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSmith2022II = getSmith2022II
 
def getAnderson2023(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD60")
    atype = self.h.getSurvName('c Diagnosis');
    atypes = ['C','AD']
    ahash = {'Unaffected':0, "Alzheimer's":1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAnderson2023 = getAnderson2023
 
def getAnderson2023II(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD60.2")
    atype = self.h.getSurvName('c predicted.id')
    ahash = {'Astrocytes':0, 'Endothelial':1, 'Excitatory':2,
            'Inhibitory':3, 'Microglia':4, 'OPCs':5,
            'Oligodendrocytes':6, 'Pericytes':7}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c Diagnosis');
    atypes = ['C','AD']
    ahash = {'Unaffected':0, "Alzheimer's":1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        atype = self.h.getSurvName('c predicted.id')
        atypes = ['Astro', 'Endo', 'Ex', 'Inh', 'Mg', 'OPCs', 'Olig', 'Peri']
        ahash = {'Astrocytes':0, 'Endothelial':1, 'Excitatory':2,
                'Inhibitory':3, 'Microglia':4, 'OPCs':5,
                'Oligodendrocytes':6, 'Pericytes':7}
    if (tn == 4):
        atype = self.h.getSurvName('c predicted.id')
        atypes = ['Astro', 'Mg']
        ahash = {'Astrocytes':0, 'Microglia':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getAnderson2023II = getAnderson2023II
 
def getZhang2022(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("ad175")
    atype = self.h.getSurvName('c title');
    atype = [re.sub("[0-9].*", "", str(k)) for k in atype]
    atypes = ['HC','AD']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZhang2022 = getZhang2022
 
def getMatarin2015Mm(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("AD62")
    atype = self.h.getSurvName('c src1');
    ahash = {'Cerebellum':0, 'Cortex':1, 'HIP':2}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c genotype');
    atypes = ['C','AD']
    ahash = {'WILD':0, 'TPM':1, 'HO_TASTPM':1, 'TAS10':1, 'TAU':1, 'HET_TASTPM':1}
    if (tn == 2):
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 3):
        ahash = {'WILD':0}
        bhash = {2:'TPM', 0:'HO_TASTPM', 3:'TAS10',
                4:'TAU', 1:'HET_TASTPM'}
        ahash[bhash[tb]] = 1
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getMatarin2015Mm = getMatarin2015Mm
 
def getJanssen2017(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("ad7","/booleanfs2/sahoo/Data/Brain/Alz-Net/sataheri.conf")
    age = self.h.getSurvName('c age (ch1)')
    ahash = {'1.5 months':0,'6 months':1,"18 months":2,'24 months':3}
    tval = [ahash[i] if i in ahash else None for i in age]
    atype = self.h.getSurvName('c genotype/variation (ch1)')
    atypes = ['WT', 'APP23 (+/-)']
    ahash = {'Wild type':0, "heterozygous APP23 (+/-)":1}
    if tn == 2:
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atype = [atype[i] if tval[i] in {2, 3}
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getJanssen2017 = getJanssen2017

def getSwarup2018(self, tn = 1, ta = 0, tb = 0):
    self.prepareData("ad17","/booleanfs2/sahoo/Data/Brain/Alz-Net/sataheri.conf")
    age = self.h.getSurvName('c age (ch1)');
    ahash = {'3 month':0,'6 month':1}
    sval = [ahash[i] if i in ahash else None for i in age]
    atype = self.h.getSurvName('c tissue (ch1)')
    ahash = {'cerebellum':0, 'cortex':1, 'hippocampus':2, 'stem':3}
    tval = [ahash[i] if i in ahash else None for i in atype]
    atype = self.h.getSurvName('c transgenic (ch1)')
    atypes = ['WT','TPR50 Tau P301S']
    ahash = {'Wt':0, 'Tg':1}
    if tn == 2:
        atype = [atype[i] if tval[i] == ta and sval[i] == tb
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSwarup2018 = getSwarup2018

def getPreuss2020(self,tn=1,ta=0, tb=0):
    self.prepareData("ad19","/booleanfs2/sahoo/Data/Brain/Alz-Net/sataheri.conf")
    age =self.h.getSurvName("c age (months) (ch1)")
    ahash = {'6', '12', '2', '4'}
    atype = self.h.getSurvName('c strain (ch1)');
    atypes = ['WT-C57BL_6J', 'APOE4 KI','Trem2-R47H','5xFAD','APOE4Trem2']
    ahash = {'Trem2-R47H':2, 'C57BL/6J':0, 'APOE4 KI':1,'5xFAD':3,'APOE4Trem2':4}
    atype = [atype[i] if age[i] == ta
             else None for i in range(len(atype))]
    if tn==2:
        atypes = ['WT_C57BL/6J','5xFAD']
        ahash = {'5xFAD':1, 'C57BL/6J':0}
    if tn==3:
        atypes = ['WT-C57BL_6J', 'APOE4 KI']
        ahash = {'APOE4 KI':1,'C57BL/6J':0}  
    if tn==4:
        atypes = ['WT-C57BL_6J', 'APOE4Trem2']
        ahash = {'APOE4Trem2':1,'C57BL/6J':0}  
    if tn==5:
        atypes = ['WT-C57BL_6J', 'APOE4-TREM2/KI']
        ahash = {'APOE4Trem2':1,'C57BL/6J':0,'APOE4 KI':1} 
    if tn==6:
        atypes = ['WT-C57BL_6J', 'Trem2-R47H']
        ahash = {'Trem2*R47H':1,'C57BL/6J':0} 
    self.initData(atype, atypes, ahash)
    return  
bone.IBDAnalysis.getPreuss2020 = getPreuss2020

def getRezaie2021(self,tn=1,ta=0, tb=0):
    self.prepareData("ad20","/booleanfs2/sahoo/Data/Brain/Alz-Net/sataheri.conf")
    tissue = self.h.getSurvName('c tissue (ch1)')
    ahash = {'hippocampus':100, 'cortex':200}
    age = self.h.getSurvName('c time (month) (ch1)');
    ahash = {'4', '8', '12', '18'}
    atype = self.h.getSurvName('c strain (ch1)');
    atypes = ['Wt','5xFAD']
    ahash = {'BL6':0, '5xFAD;BL6':1}
    if tn == 2:
        atype = [atype[i] if tissue[i] == ta and age[i] == tb
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return  
bone.IBDAnalysis.getRezaie2021 = getRezaie2021

def getNachun2020ADBlood(self,tn=1,ta=0, tb=0):
    self.prepareData("ad188","/booleanfs2/sahoo/Data/Brain/Alz-Net/sataheri.conf")
    atype = self.h.getSurvName('c apoe_status (ch1)');
    ahash = {'E2_E2':0, 'E2_E3':1, 'E2_E4':2, 'E3_E3':3, 'E3_E4':4, 'E4_E4':5}
    tval = [ahash[k] if k in ahash else None for k in atype]
    atype = self.h.getSurvName('c diagnosis (ch1)');
    atypes = ['Control','MCI','AD','bvFTD','svPPA','PSP','nfvPPA','CBS']
    ahash = {'Control':0,'MCI':1, 'AD':2,'bvFTD':3,'svPPA':4,
            'PSP':5,'nfvPPA':6, 'CBS':7}
    if tn == 2:
        atypes = [  'Control', 'AD' ]
        ahash = {'Control':0, 'AD':1}
    if tn == 3:
        atypes = [  'Control', 'AD' ]
        ahash = {'Control':0, 'AD':1}
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return  
bone.IBDAnalysis.getNachun2020ADBlood = getNachun2020ADBlood

def getNachun2020ADBloodII(self,tn=1,ta=0, tb=0):
    self.prepareData("ad96","/booleanfs2/sahoo/Data/Brain/Alz-Net/sataheri.conf")
    atype = self.h.getSurvName('c diagnosis (ch1)');
    atypes = ['Control', 'bvFTD', 'svPPA', 'PSP', 'nfvPPA', 'CBS']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return  
bone.IBDAnalysis.getNachun2020ADBloodII = getNachun2020ADBloodII

def getSamsudin2016adBlood(self, tn = 1, ta=0):
    self.prepareData("AD29")
    atype = self.h.getSurvName('c gender');
    ahash = {'male':0, 'female':1}
    tval = [ahash[k] if k in ahash else None for k in atype]
    atype = self.h.getSurvName('c diagnosis');
    atypes = ['control',"AD"]
    ahash = {"Probable Alzheimer's Disease":1, 'Non-demented control':0}
    if tn == 2:
        atypes = ['control',"AD"]
        atype = [atype[i] if tval[i] == ta
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSamsudin2016adBlood = getSamsudin2016adBlood

def getSMps19PFCmmIII(self, tn=1, ta=0):
    self.prepareData("AD55.11")
    atype = self.h.getSurvName("c gender")
    ahash = {'male':0, 'female':1}
    sval = [ahash[k] if k in ahash else None for k in atype]
    atype = self.h.getSurvName("c Type")
    atypes = ['WT', 'hTau', 'CgA-KO/hTau', 'CgA-KO']
    ahash = {'WT':0, 'CgaKO/hTau':2, 'hTau':1, 'CgaKO':3}
    if tn == 2:
        atype = [atype[i] if sval[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atypes = ['WT', 'hTau', 'CgA-KO/hTau-F', 'CgA-KO/hTau-M', 'CgA-KO']
        atype = ['CgA-KO/hTau-M' if sval[i] == 0 and atype[i] == 'CgaKO/hTau'
                 else atype[i] for i in range(len(atype))]
        atype = ['CgA-KO/hTau-F' if sval[i] == 1 and atype[i] == 'CgaKO/hTau'
                 else atype[i] for i in range(len(atype))]
        bhash = {'D5_S86', 'D2_S122', 'D3_S38', 'B1_S57', 'B2_S120', 'C1_S58'}
        bhash = {'B1_S57', 'B2_S120'}
        atype = [None if self.h.headers[i] in bhash
                 else atype[i] for i in range(len(atype))]
        ahash = {'WT':0, 'CgaKO/hTau':2, 'hTau':1, 'CgaKO':4}
    if tn == 4:
        atypes = ['hTau', 'CgA-KO/hTau-F', 'CgA-KO/hTau-M', 'CgaKO']
        atype = ['CgA-KO/hTau-M' if sval[i] == 0 and atype[i] == 'CgaKO/hTau'
                 else atype[i] for i in range(len(atype))]
        atype = ['CgA-KO/hTau-F' if sval[i] == 1 and atype[i] == 'CgaKO/hTau'
                 else atype[i] for i in range(len(atype))]
        bhash = {'D5_S86', 'D2_S122', 'D3_S38', 'B1_S57', 'B2_S120', 'C1_S58'}
        bhash = {'B1_S57', 'B2_S120', 'D5_S86'}
        atype = [None if self.h.headers[i] in bhash
                 else atype[i] for i in range(len(atype))]
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSMps19PFCmmIII = getSMps19PFCmmIII

