from __future__ import print_function
import scipy.spatial
import numpy as np

class MultiClicker:
    """ Class to collect the clicked points
    You initialize it as 
    >>> cl = MultiClicker(plt.gcf())
    Then all the clicks will be recorded in cl.points
    If you click the right-hand button, the recording will stop.
    You can also stop recording by calling 
    >>> cl.stop()
    """ 
    def __init__(self, fig):
        self.cid = None 
        self.points = []
        def onclick(event):
            try:
                print( 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
                        event.button, event.x, event.y, event.xdata, event.ydata))                
                if event.button==3:
                    self.stop()
                else:
                    self.points.append((event.xdata,event.ydata))
            except:
                raise
                print( 'printing failed')
        self.canvas = fig.canvas
        self.cid = self.canvas.mpl_connect('button_press_event', onclick)
    def stop(self):
        self.canvas.mpl_disconnect(self.cid)
    
    

class NearbyClicker:
    """ Class to call function on clicked points
    If you have plotted the xs,ys points
    and defined function
    def callback(i)
    Then when you create a clicker
    >>> NearbyClicker(plt.gcf(), xs, ys, callback)
    At each click a point that is closest to the clicked location is recorded 
    and the callback function is called with the integer number of this point
    
    """
    def __init__(self, fig, xs, ys, callback):
        self.cid = None 
        self.xs = xs
        self.ys = ys
        self.tree = scipy.spatial.cKDTree(np.array([xs,ys]).T)
        def onclick(event):
            try:
                if event.button==3:
                    self.stop()
                else:
                    d,pos=self.tree.query(np.r_[event.xdata,event.ydata])
                    print ('Clicked %f %f ; Selected %d: %f %f' %(event.xdata,event.ydata,pos,self.xs[pos],self.ys[pos] ))
                    callback(pos)
            except:
                print( 'callback failed')
        self.canvas = fig.canvas
        self.cid = self.canvas.mpl_connect('button_press_event', onclick)
    def stop(self):
        self.canvas.mpl_disconnect(self.cid)
    
    
def clicker(fig,xobj = None):
    cid = None 
    def onclick(event):
        try:
            print ('cid=%s button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(cid,
                event.button, event.x, event.y, event.xdata, event.ydata))
            if xobj is not None:
                if isinstance (xobj,dict):
                    xobj['x']=event.xdata
                    xobj['y']=event.ydata
        except:
            print( 'printing failed')
        event.canvas.mpl_disconnect(cid)

    cid = fig.canvas.mpl_connect('button_press_event', onclick)


def clicker_multi(fig):
    cid = None 
    res = [] 
    def onclick(event):
        try:
            print( 'cid=%s button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(cid,
                    event.button, event.x, event.y, event.xdata, event.ydata))
            if event.button==3:
                event.canvas.mpl_disconnect(cid)
            else:
                res.append((event.xdata,event.ydata))
        except:
            print( 'printing failed')

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    return res    
#cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
