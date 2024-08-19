import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def hideOnClick(fig=None,ax=None):
    if fig is None and ax is None:
            fig, ax = plt.subplots()
    elif fig is None and ax is not None:
        fig = ax.figure
    elif fig is not None and ax is None:
        ax = fig.gca()


    legend=ax.get_legend()
    
    lines=ax.lines[:]
    leglines=legend.get_lines()
    patches=ax.patches[:]
    legPatches=legend.get_patches()
    legElements=leglines+legPatches
    axElements=lines+patches
    graphs={}
    for count,l in enumerate(legElements):
        l.set_picker(True)
        if isinstance(l,Line2D):
            l.set_pickradius(5)
        label=l._label
        line=[x for x in axElements if x._label == label]
        if label in graphs:
            continue
        graphs[label]=line
    def on_pick(event):
        legline = event.artist
        graph=graphs[legline._label]       
        for g in graph:
            vis = not g.get_visible()
            g.set_visible(vis)
        legline.set_alpha(1.0 if vis else 0.2)
        fig.canvas.draw()
        
    fig.canvas.mpl_connect('pick_event', on_pick)