> ## Plotting ##

This is mostly a set of wrappers mimicking IDL plotting commands, but also some new useful commands:

Most of the plotting commands are in idlplot package.
The set of routines includes plot(), oplot(), contour(), ploterror(), oploterror(), tvaxis(), tvhist2d().
The difference between these and the standard plot(),contour,ploterror from matplotlib is that I allow a lot of options to be specified in one single call and several options have the same name as in IDL, e.g:
```
plot(a,b,ps=3,xlog=True,ylog=True,xtitle='A',ytitle='B',yr=[1e-3,30],xr=[0.1,1e5])
```
plothist routine is mimicking astrolib plothist routine. But it also supports different kernels for histogram construction (statistics module is required), and also supports adaptive bin size selection.
There is also a useful wrapper module called idlplotInd which supports indexing in plotting, e.g.
```
plot(x,y,ind=(x+y)<3) 
```
which can be useful for interactive work

I also have two routines tvaxis and tvhist2d for plotting of 2D images and two dimensional histograms.

Another useful routine is lasso\_plot, which allows to select the points from the plot by drawing the region around them (similar to what TOPCAT allows).

## Fitting ##
  * MPFIT
This package provides the python interface to the Levenberg-Marquardt MPFIT chi-square fitter (based on the IDL version of MPFIT by Craig Markwardt and first translation to python by Mark Rivers). It is based on numpy and implements important bugfixes from the IDL version of MPFIT.

## Astrolib ##
Some astrolib utils which I converted to python :
baryvel
bprecess
convolve
cv\_coord
daycnv
euler
gal\_uvw
galage
helcorr
helio\_jd
lumdist
mwrfits
precess
precess\_xyz
premat
readcol
sphdist
xyz
zang

## Miscellaneous ##
  * Some other useful routines, such as idlsave , which gives you the opportunity to save and restore variables in the way very similar to IDL.
  * sqlutil -- Highly tuned module which allows you to query databases, upload and retrieve numpy arrays (currently tuned to works with postgresql and sqlite to some extend)
  * quick\_hist -- the N-dimensional histogram routine which is order of magnitude faster than the scipy.histogramnd
  * match\_lists -- routine to match two list of (ra,dec) by position
  * get\_dss -- DSS cutout retrieval
  * resolve -- astronomical name resolution
  * adabinner -- the routine to do adaptive 1D and 2D histogram construction (using quad-tree)
  * wav2RGB -- convert the wavelength to matplotlib color string
  * pg2hdf5 -- save the results of the big query in an hdf5 file
  * window\_func -- function which execute a user supplied function on the data which is within a given bin.
## Installation ##

There are no installation scripts for astrolibpy. So if you download the whole astrolibpy or some programs, you need to put them into your PYTHONPATH.