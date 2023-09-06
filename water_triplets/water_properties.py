# Shell Lab's code for water triplets
import numpy as np
import waterlib as wl

def getTripletAngs(subPos, Pos, BoxDims, lowCut=0.0, highCut=3.413):
  """This is called getCosAngs, but actually just returns the angles themselves (faster to convert
     from cos(theta) to theta in Fortran)
     Inputs:
     subPos - positions of set of atoms to measure tetrahedrality of (may be different, subset, or same as Pos)
     Pos - positions of ALL atoms that can make tetrahedral configurations (needed if subPos not same as Pos)
     BoxDims - current box dimensions to account for periodicity
     lowCut - lower cutoff for nearest-neighbor shell (default 0.0)
     highCut - higher cutoff for nearest-neighbor shell (default 3.413 - see Chaimovich, 2014, but should really
               change to reflect first peak in g(r) for the chosen water model)
     Outputs:
     angVals - all angle values for current configuration of positions supplied
     numAngs - number of angles for each central oxygen atom (i.e. number neighbors factorial)

  """

  #Set-up array to hold angle results and stack as go... list increases in size!
  angVals = np.array([])
  numAngs = np.zeros(len(subPos))

  #Find nearest neighbors for ALL atoms in subPos
  #But make sure using efficient algorithm...
  #If subPos is same as Pos, use allnearneighbors instead
  if np.array_equal(subPos, Pos):
    nearNeighbs = wl.allnearneighbors(Pos, BoxDims, lowCut, highCut).astype(bool)
  else:
    nearNeighbs = wl.nearneighbors(subPos, Pos, BoxDims, lowCut, highCut).astype(bool)

  #Loop over each position in subPos, finding angle made with all neighbor pairs
  for (i, apos) in enumerate(subPos):
    #Make sure have nearest neighbors...
    if len(Pos[nearNeighbs[i]]) > 0:
      #below returns symmetric, square array (zero diagonal)
      tempAng = wl.tetracosang(apos, Pos[nearNeighbs[i]], BoxDims) 
      #Only want half of array, flattened
      angVals = np.hstack((angVals, tempAng[np.triu_indices(len(tempAng),k=1)].tolist()))
      numAngs[i] = tempAng.shape[0]

  return angVals, numAngs
  
  