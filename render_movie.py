import baltic as bt
import re
import copy

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.collections import PathCollection, LineCollection
import matplotlib.transforms as transforms
from matplotlib.patches import Polygon, Circle, PathPatch
from matplotlib.path import Path
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib import animation
import matplotlib.ticker as ticker

import numpy as np
import pandas as pd

# Pulled from baltic https://github.com/blab/baltic
tree_path='data/highest_posterior.tree' ## define path to (in this case) FigTree file
print "File: %s"%(tree_path.split('/')[-1])
tipFlag=False ## used to keep track of where we are in FigTree format
tips={} ## dictionary of tip code to full tip name

for line in open(tree_path,'r'): ## iterate through FigTree lines
    l=line.strip('\n') ## strip newline characters from each line

    cerberus=re.search('dimensions ntax=([0-9]+);',l.lower()) ## check how many tips there are supposed to be
    if cerberus is not None:
        tipNum=int(cerberus.group(1))

    #####################
    cerberus=re.search('tree TREE([0-9]+) = \[&R\]',l) ## search for beginning of tree string in BEAST format
    if cerberus is not None:
        treeString_start=l.index('(') ## tree string starts where the first '(' is in the line
        ll=bt.tree() ## new instance of tree
        bt.make_tree(l[treeString_start:],ll) ## send tree string to make_tree function, provide an empty tree object
    #####################
        
    if tipFlag==True:
        cerberus=re.search('([0-9]+) ([A-Za-z\-\_\/\.\'0-9 \|?]+)',l) ## look for tip name map, where each tip is given an integer to represent it in tree
        if cerberus is not None:
            tips[cerberus.group(1)]=cerberus.group(2).strip("'") ## if you give tips an integer (in the form of a string), it will return the full name of the tip
        elif ';' not in l: ## something's wrong - nothing that matches the tip regex is being captured where it should be in the file
            print 'tip not captured by regex:',l.replace('\t','')

    if 'translate' in l.lower(): ## start looking for tips
        tipFlag=True
    if ';' in l: ## stop looking for tips
        tipFlag=False

print "Number of objects found in tree string: %d"%(len(ll.Objects))

## rename tips, find the highest tip (in absolute time) in the tree
if len(tips)==0: ## use this if tip names in the string are already the final format
    for k in ll.Objects:
        if isinstance(k,bt.leaf):
            k.name=k.numName
    highestTip=max([bt.decimalDate(x.name.strip("'").split('|')[-1].replace("_","-"),variable=True) for x in ll.Objects if isinstance(x,bt.leaf)])
else: ## there's a tip name map at the beginning, so translate the names
    ll.renameTips(tips) ## give each tip a name
    highestTip=max([bt.decimalDate(x.strip("'").split('|')[-1].replace("_","-"),variable=True) for x in tips.values()])

    
ll.sortBranches()
ll.setAbsoluteTime(highestTip)

lowestTip=min([x.absoluteTime for x in ll.Objects if x != ll.root])
timeTraversed = highestTip - lowestTip

coords = pd.read_csv("data/coordinates.csv", names=["Country", "lat", "lng", "Continent"])
coords["Country_NA_EARTH"] = coords["Country"]

# Match country names in coordinates.csv with names in shape file
coords.loc[coords["Country_NA_EARTH"]=="USA", "Country_NA_EARTH"] = "United States"
coords.loc[coords["Country_NA_EARTH"]=="Micronesia", "Country_NA_EARTH"] = "Federated States of Micronesia"
coords.loc[coords["Country_NA_EARTH"]=="Russia", "Country_NA_EARTH"] = "Russian Federation"

# Specify countries with only travel cases. 
travel_countries = ["Russia", "Canada", "United_Kingdom", "Italy", "Austria", "China", "Japan"]

groups = coords["Continent"].unique()
groups = groups[groups!="Europe"].tolist()
color_scheme = ["#4c71a5", "#47a364", "#d0684a", "#e1c62f", "#cc79a7", "#77bedb", "#7f6e85", "#476581", "#cccccc", "#f6653c", "#66ccff", "#ccc097"]
_colors = color_scheme[:len(groups)]
# Setup color scale for countries
_cmap = colors.LinearSegmentedColormap.from_list('mycmap', _colors)


z = range(1,len(groups))
hot = plt.get_cmap('hot')
cNorm  = colors.Normalize(vmin=0, vmax=len(groups))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=_cmap)

def decimalDateToMonthYear(decimalDate):
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    mIndex = int(np.floor((decimalDate%np.floor(decimalDate)) * 12))
    return months[mIndex]+" "+str(int(np.floor(decimalDate)))

f = plt.figure(figsize=(20,10))
gs = gridspec.GridSpec(5, 10, height_ratios=[0.2,1,1,1,1], width_ratios=[1,1,1,1,1,1,1,1,1,0.2])
ax = plt.subplot(gs[1:5, :5])
treeax = plt.subplot(gs[:5,5:9])
regionax =  plt.subplot(gs[:5, 9])
titleax = plt.subplot(gs[0,:5])
ax.axis("off")
titleax.axis("off")
m = Basemap(projection='eck4',lon_0=160,resolution='c', ax = ax)
m.drawcoastlines(linewidth=0, color="#FFFFFF")
m.drawmapboundary(color="aqua")
m.fillcontinents(color='#cccccc',lake_color='#FFFFFF')
# Read shapefile
m.readshapefile("data/ne_10m_admin_0_countries/ne_10m_admin_0_countries", "units", drawbounds=False)

patches = []
for info, shape in zip(m.units_info, m.units):
    poly = Polygon(np.array(shape), True)
    poly.set_facecolor('none')
    poly.set_linewidth(0)
    poly.set_edgecolor("#000000")
    poly.set_zorder(1)
    poly=ax.add_patch(poly)
    patches.append(poly)

x1,y1 = m(coords["lng"].values,coords["lat"].values)
_c = [scalarMap.to_rgba(groups.index(i)) for i in groups]
p = m.scatter(x1,y1,marker="o",alpha=1, color=_c, zorder=2, sizes = [0]*coords.shape[0])
lwidths = list(np.logspace(np.log(1), np.log(5), 300, base = np.e))

beziers = []
rounded_edges_coords = []
_colors = []
for i in ll.Objects:
    if i == ll.root or isinstance(i, bt.leaf):
        continue
    for j in i.children:
        if j.traits["Location"] != i.traits["Location"]:
            _from = coords[coords["Country"] == i.traits["Location"]]
            _to = coords[coords["Country"] == j.traits["Location"]]
            x, y = m.gcpoints(_from["lng"], _from["lat"], _from["lng"], _from["lat"], 300)
            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            _c.extend(["#000000"]*300)
            lc = LineCollection(segments, linewidths=[0]*300,colors=_c)
            lc = ax.add_collection(lc)
            beziers.append(lc)
            rounded_edges_coords.append([x[0], y[0]])
            _colors.append("#FFFFFF") # Lines drawn in individual frames. Color set in animate()
rounded_edges = m.scatter([i[0] for i in rounded_edges_coords], [i[1] for i in rounded_edges_coords], sizes = [0] * len(rounded_edges_coords), zorder = 1, color="#FFFFFF", edgecolor="#000000", linewidth= 0.4)

titleax.spines['left'].set_position('zero')
titleax.spines['bottom'].set_position('zero')
titleax.set_xlim(-1, 1)
titleax.set_xlim(-1, 1)
timeText = titleax.annotate(str(decimalDateToMonthYear(lowestTip)), (0, 0), size=25, horizontalalignment='center')

# Get longtidue where edge of map is encountered. This is to avoid ugly line across map.
eq = [m(i,90)[0] for i in range(-180, 180)]
edgeDiff = int((max(eq) - min(eq))/10**7) * 10 ** 7

# Plot tree
for i in range(int(round(lowestTip))+1, int(round(highestTip))+2):
    treeax.axvspan(i, i+.5, facecolor='#cccccc', alpha=0.5)

branchWidth=1 
nodeSize=50 
nodeAttr = []
branchSegments = []
for k in ll.Objects:
    x=k.absoluteTime
    y=k.y
    xp=k.parent.absoluteTime 
    if x==None: 
        x=0.0
    if xp==None:
        xp=x
    c='#000000'    
    if isinstance(k,bt.leaf) or k.branchType=='leaf':
        # Remove confirmed travel related cases
        if k.traits["Location"] in travel_countries:
            k.traits["Location"] = k.parent.traits["Location"]            
        _g = coords[coords["Country"]==k.traits["Location"]]["Continent"].values[0]
        print k.traits["Location"]
        _c = scalarMap.to_rgba(groups.index(_g))
        # _c = "#FFFFFF"
        nodeAttr.append({
            "size": nodeSize-30*k.height/ll.treeHeight,
            "x": x,
            "y": y,
            "color": _c
        })
    elif isinstance(k,bt.node) or k.branchType=='node': ## if node...
        branchSegments.append([[x, k.children[-1].y], [x,k.children[0].y]])
    branchSegments.append([[xp, y], [x, y]])

# Plot background 'grey' tree
greybackgroundNodeScatter = treeax.scatter([i["x"] for i in nodeAttr],[i["y"] for i in nodeAttr],sizes=[i["size"] * 1.8 for i in nodeAttr],color='#707070',zorder=10) ## plot black circle underneath
greyBranchCollection = LineCollection(branchSegments, lw=branchWidth,color="#707070",zorder=9)
treeax.add_collection(greyBranchCollection)

# Setup linecollection that plots tree with time.
backgroundNodeScatter = treeax.scatter([i["x"] for i in nodeAttr],[i["y"] for i in nodeAttr],sizes=[0 * 1.8 for i in nodeAttr],color='k',zorder=10) ## plot black circle underneath
nodeScatter = treeax.scatter([i["x"] for i in nodeAttr],[i["y"] for i in nodeAttr],sizes=[0 for i in nodeAttr],color=[i["color"] for i in nodeAttr],zorder=11) ## plot circle for every tip
branchCollection = LineCollection(branchSegments, lw=branchWidth,color="#707070",zorder=9)
treeax.add_collection(branchCollection)

treeax.axes.get_yaxis().set_visible(False)
treeax.get_xaxis().set_ticks(list(np.arange(2000, 2019, 2)))
treeax.set_ylim(-5,ll.ySpan+5)
timeline, = treeax.plot([lowestTip, lowestTip], [-5, ll.ySpan+5], lw=2, color='#000000')

cbar = mpl.colorbar.ColorbarBase(regionax, cmap=_cmap, norm=cNorm, boundaries=list(np.arange(-0.5, len(groups))))
cbar.set_ticks(ticker.MaxNLocator(integer=True))
cbar.ax.set_yticklabels(groups)
cbar.ax.get_xaxis().set_visible(False)

gs.tight_layout(f, rect=[0, 0, 0.95, 1])

def animate(i):
    print(i)
    t = timeIntervals[i]
    timeText.set_text(decimalDateToMonthYear(t))
    timeline.set_data([t, t], [-5, ll.ySpan+5])
    timeline.set_zorder(13)
    highlight = []
    bctr = 0
    _sizes = [0] * len(rounded_edges_coords)
    _coords = [[0,0]] * len(rounded_edges_coords)
    _colors = ["#FFFFFF"] * len(rounded_edges_coords)
    branchSegments = []
    nodeAttr = []
    segmentWidth = []
    branchColor = []
    for k in ll.Objects:
        nx=k.absoluteTime
        ny=k.y
        nxp=k.parent.absoluteTime        
        if nx==None: 
            nx=0.0
        if nxp==None:
            nxp=nx
        _nx = nx
        _ny = ny
        _nxp = nxp
        _width = 1
        _nx = min(nx, t)
        _nxp = min(_nxp, t)
        if k.absoluteTime > t and k.parent.absoluteTime > t:
            _width = 0            
        # Known travel cases
        if k.traits["Location"] in travel_countries:
            k.traits["Location"] = k.parent.traits["Location"]
        _g = coords[coords["Country"]==k.traits["Location"]]["Continent"].values[0]
        _c = scalarMap.to_rgba(groups.index(_g))
        if isinstance(k,bt.leaf) or k.branchType=='leaf':
            _nodeSize = nodeSize-30*k.height/ll.treeHeight
            if k.absoluteTime > t:
                _nodeSize = 0
            nodeAttr.append({
                "size": _nodeSize,
                "x": nx,
                "y": ny,
                "color": _c
            })
        elif isinstance(k,bt.node) or k.branchType=='node':
            vertical_branch_width = 1
            if k.absoluteTime > t:
                vertical_branch_width = 0    
            branchSegments.append([[_nx, k.children[-1].y], [_nx,k.children[0].y]])
            segmentWidth.append(vertical_branch_width)
            branchColor.append(_c)
        branchSegments.append([[_nxp, _ny], [_nx, _ny]])
        segmentWidth.append(_width)
        branchColor.append(_c)
    # Loop to draw great circle curves.
    for k in ll.Objects:
        if k.absoluteTime <= t:
            if k.traits["Location"] not in highlight and k.traits["Location.prob"] > 0.7:
                highlight.append(k.traits["Location"])
            if not isinstance(k, bt.leaf):
                for j in k.children:
                    if j.traits["Location"] != k.traits["Location"] and j.traits["Location"]>0.7:
                        _from = coords[coords["Country"] == k.traits["Location"]]
                        _to = coords[coords["Country"] == j.traits["Location"]]
                        x, y = m.gcpoints(_from["lng"], _from["lat"], _to["lng"], _to["lat"], 300)
                        cut_point = np.where(np.abs(np.diff(x)) > edgeDiff)[0]
                        if len(cut_point)>0:
                            cut_point = cut_point[0]
                            x = np.concatenate(
                                [x[:cut_point],
                                 [np.nan],
                                 x[cut_point+1:]]
                            )
                            y = np.concatenate(
                                [y[:cut_point], 
                                [np.nan],
                                y[cut_point+1:]]
                            )
                        curveTimeRem = (t - k.absoluteTime)/(j.absoluteTime - k.absoluteTime)
                        curveStart = min(curveTimeRem, 1)
                        curveEnd = min(curveTimeRem * 0.6, 1)
                        curveStart = int(round(curveStart * len(x)))
                        curveEnd = int(round(curveEnd * len(x)))
                        curveEnd = 0
                        _lw = lwidths[curveEnd:curveStart]
                        if len(_lw) > 0:
                            points = np.array([x[curveEnd:curveStart], y[curveEnd:curveStart]]).T.reshape(-1, 1, 2)
                            segments = np.concatenate([points[:-1], points[1:]], axis=1)
                            beziers[bctr].set_color(scalarMap.to_rgba(groups.index(coords[coords["Country"] == j.traits["Location"]]["Continent"].values[0])))
                            beziers[bctr].set_segments(segments)
                            beziers[bctr].set_linewidth(_lw)
                            _sizes[bctr] = _lw[-1] * 2
                            _coords[bctr] = [x[curveStart-1], y[curveStart-1]]
                            _colors[bctr] = scalarMap.to_rgba(groups.index(coords[coords["Country"] == j.traits["Location"]]["Continent"].values[0]))
                        if curveStart == len(x):
                            beziers[bctr].set_color(scalarMap.to_rgba(groups.index(coords[coords["Country"] == j.traits["Location"]]["Continent"].values[0])))
                            beziers[bctr].set_alpha(0.5)
                            _colors[bctr] = scalarMap.to_rgba(groups.index(coords[coords["Country"] == j.traits["Location"]]["Continent"].values[0]))
                            _sizes[bctr] = 10
                        bctr += 1
    branchCollection.set_segments(branchSegments)
    branchCollection.set_linewidth(segmentWidth)
    branchCollection.set_color(branchColor)
    branchCollection.set_alpha(0.6)
    nodeScatter.set_sizes(np.array([i["size"] for i in nodeAttr]))
    nodeScatter.set_zorder(11)
    backgroundNodeScatter.set_sizes(np.array([i["size"] * 1.8 for i in nodeAttr]))
    backgroundNodeScatter.set_zorder(10)
    rounded_edges._sizes = _sizes
    rounded_edges.set_offsets(_coords)
    # rounded_edges.set_color(_colors)
    rounded_edges.set_color("#000000")
    rounded_edges.set_zorder(20)
    _groups = coords[coords["Country"].isin(highlight)]["Continent"].tolist()
    _highlight = coords[coords["Country"].isin(highlight)]["Country_NA_EARTH"].str.lower().str.replace(" ", "_").tolist()
    pctr = 0
    for info, shape in zip(m.units_info, m.units):
        n = info["NAME_LONG"].lower().replace(" ","_")
        if n in _highlight:
            poly = patches[pctr]
            _g = _groups[_highlight.index(n)]
            poly.set_facecolor(scalarMap.to_rgba(groups.index(_g)))
            poly.set_zorder(1)
            poly.set_alpha(0.3)
        pctr += 1
    _sizes = []
    for i in coords.index:
        if coords.ix[i]["Country"] in highlight:
            _sizes.append(2)
        else:
            _sizes.append(0)
    p._sizes = _sizes
    _ = patches
    _.extend([p])
    _.extend(beziers)
    _.append(timeText)
    _.append(rounded_edges)
    _.append(branchCollection)
    _.append(backgroundNodeScatter)    
    _.append(nodeScatter)
    _.append(timeline)    
    return _

_frames = 240                   # 18 second video
_buffer = 30
# Time moves fast initially and then slows down from 2014 - 2017.
timeIntervals = highestTip - np.logspace(np.log(timeTraversed), np.log(4), 100, base=(np.e), endpoint=True)
timeIntervals = np.append(timeIntervals, np.linspace(timeIntervals[-1], highestTip+0.5, (_frames+_buffer)-100))
anim = animation.FuncAnimation(plt.gcf(), animate, frames=_frames+_buffer, interval = 20, blit=True, repeat=False)
FFwriter = animation.FFMpegWriter(fps=15, bitrate=-1)
anim.save('movie/zika_phylogeography.mp4',writer=FFwriter)
# plt.show()
plt.clf()
plt.close()
