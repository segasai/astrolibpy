from numpy import *
from numpy.fft import fft2, ifft2

def convolve(image, psf, ft_psf=None, ft_image=None, no_ft=None, correlate=None, auto_correlation=None):
   """
    NAME:
          CONVOLVE
    PURPOSE:
          Convolution of an image with a Point Spread Function (PSF)
    EXPLANATION:
          The default is to compute the convolution using a product of
          Fourier transforms (for speed).
   
    CALLING SEQUENCE:
   
          imconv = convolve( image1, psf, FT_PSF = psf_FT )
     or:
          correl = convolve( image1, image2, /CORREL )
     or:
          correl = convolve( image, /AUTO )
   
    INPUTS:
          image = 2-D array (matrix) to be convolved with psf
          psf = the Point Spread Function, (size < or = to size of image).
   
    OPTIONAL INPUT KEYWORDS:
   
          FT_PSF = passes out/in the Fourier transform of the PSF,
                  (so that it can be re-used the next time function is called).
          FT_IMAGE = passes out/in the Fourier transform of image.
   
          /CORRELATE uses the conjugate of the Fourier transform of PSF,
                  to compute the cross-correlation of image and PSF,
                  (equivalent to IDL function convol() with NO rotation of PSF)
   
          /AUTO_CORR computes the auto-correlation function of image using FFT.
   
          /NO_FT overrides the use of FFT, using IDL function convol() instead.
                  (then PSF is rotated by 180 degrees to give same result)
    METHOD:
          When using FFT, PSF is centered & expanded to size of image.
    HISTORY:
          written, Frank Varosi, NASA/GSFC 1992.
          Appropriate precision type for result depending on input image
                                  Markus Hundertmark February 2006
          Fix the bug causing the recomputation of FFT(psf) and/or FFT(image)
                                  Sergey Koposov     December 2006
   """

   n_params = 2
   psf_ft = ft_psf
   imft = ft_image
   noft = no_ft
   auto = auto_correlation
   
   sp = array(shape(psf_ft)) 
   sif = array(shape(imft))
   sim = array(shape(image))
   sc = sim / 2
   npix = array(image, copy=0).size
   
   if image.ndim!=2 or noft!=None:   
      if (auto is not None):   
         message("auto-correlation only for images with FFT", inf=True)
         return image
      else:   
         if (correlate is not None):   
            return convol(image, psf)
         else:
            return convol(image, rotate(psf, 2))
   
   if imft==None or (imft.ndim!=2) or imft.shape!=im.shape: #add the type check
      imft = ifft2(image)
   
   if (auto is not None):   
      return roll(roll(npix * real(fft2(imft * conjugate(imft))), sc[0], 0),sc[1],1)

   if (ft_psf==None or ft_psf.ndim!=2 or ft_psf.shape!=image.shape or 
            ft_psf.dtype!=image.dtype):
      sp = array(shape(psf))
      
      loc = maximum((sc - sp / 2), 0)         #center PSF in new array,
      s = maximum((sp / 2 - sc), 0)        #handle all cases: smaller or bigger
      l = minimum((s + sim - 1), (sp - 1))
      psf_ft = conjugate(image) * 0 #initialise with correct size+type according 
      #to logic of conj and set values to 0 (type of ft_psf is conserved)
      psf_ft[loc[1]:loc[1]+l[1]-s[1]+1,loc[0]:loc[0]+l[0]-s[0]+1] = \
                     psf[s[1]:(l[1])+1,s[0]:(l[0])+1]
      psf_ft = ifft2(psf_ft)
   
   if (correlate is not None):   
      conv = npix * real(fft2(imft * conjugate(psf_ft)))
   else:   
      conv = npix * real(fft2(imft * psf_ft))
   
   sc = sc + (sim % 2)   #shift correction for odd size images.
   
   return roll(roll(conv, sc[0],0), sc[1],1)

