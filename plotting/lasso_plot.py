from matplotlib.widgets import Lasso
from matplotlib.nxutils import points_inside_poly
from matplotlib.colors import colorConverter
from matplotlib.collections import RegularPolyCollection

from matplotlib.pyplot import gca
from numpy import nonzero

class lasso_plot:
    # The class is designed to select the datapoints by drawing the region
    # around it. 
    # Example:
    # plot(xs, ys, ps=3 ) # first plot the data
    # las = lasso_plot.lasso_plot(xs,ys)
    # Now click on the plot and do not release the mouse button till 
    # you draw your region 
    # After that the variable las.ind will contain the indices of those 
    # points inside the region and las.verts will contain the vertices of the 
    # polygon you've just drawn 
    def __init__(self, xs, ys):
        self.axes = gca()
        self.canvas = self.axes.figure.canvas
        self.xys = [d for d in zip(xs,ys)]
        fig = self.axes.figure
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
        self.ind = None

    def callback(self, verts):
        ind = nonzero(points_inside_poly(self.xys, verts))[0]
        self.canvas.draw_idle()
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso
        self.verts = verts
        self.ind = ind

    def onpress(self, event):
        if self.canvas.widgetlock.locked(): return
        if event.inaxes is None: return
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)

