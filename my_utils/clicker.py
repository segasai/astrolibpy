from __future__ import print_function

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
