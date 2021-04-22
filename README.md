# astrolibpy 
### This library contains a few routines I wrote or I converted from IDL.
### Author (2008-2021) Sergey Koposov
### Email: skoposov __AT__ ed __DOT__ ac __DOT__ uk

The tools that may be still interesting to users are 
* idlsave.py -- module for saving restoring variables 
* lasso_plot -- module for selecting points by drawing a contour
* clicker -- module for extracting one or many clicks coordinates from matplotlib plot
* crossmatcher -- module for x-matching a table with the table inside postgresql database
* match_lists -- x-match lists by position on the sky using kd-tree
* mpfit -- Levenberg-Marquardt fitter (converted from IDL)
* idlplot -- set of wrappers for plotting in IDL style where your command has all 
  the options you want. Also allows efficient histogramming, subsetting etc

The astrolib part of this repo is probably obsolete with the arrival of astropy

The licensing for the programs I wrote myself is GPL3. For all other
programs (mainly converted from IDL) I guess the license is either BSD or 
they are in public domain. 
