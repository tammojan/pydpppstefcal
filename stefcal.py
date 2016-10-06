import numpy
from numpy.linalg import norm
import numpy as np
from numpy import conj, dot
from math import sqrt
from time import time
# (44, 2, 40, 1, 2, 44)

def stefcal(datacube, modeldatacube):
  ''' Solve for diagonal antenna gains using stefcal (so not full-Jones)
      datacube: [nSt, 2, nCh, nTime, 2, nSt] numpy array with visibilities
      modeldatacube: similar numpy array
      Both arrays should be in Fortran order
      Return value is a list of [nSt*2] complex gains
  '''
  oldshape = datacube.shape
  newshape = (oldshape[0]*oldshape[1], oldshape[2], oldshape[3]*oldshape[4]*oldshape[5])
  V=datacube.reshape(newshape)
  M=modeldatacube.reshape(newshape)

  G=np.ones(newshape[0])*sqrt(norm(V)/norm(M))

  nSt = newshape[0]
  nCh = newshape[1]

  iter=0
  while iter<20:
    iter+=1
    time0=time()
    Gold = G
    for st1 in range(nSt):
      w=0.
      t=0.
      for ch in range(nCh):
        z = conj(Gold) *  M[:, ch, st1] # element-wise
        w += dot(conj(z), z).real
        t += dot(conj(z), V[:, ch, st1]).real
      G[st1] = t/w;
    if iter % 2 == 0:
      dG= norm(G - Gold) / norm(G);
      G = 0.5 * G + 0.5 * Gold
      if dG<1.e-5:
        break 

  return G
