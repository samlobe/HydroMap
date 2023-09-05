#A library of code to examine properties of bulk water and near solutes
#
#Should eventually be able to handle local densities and fluctuations,
#solute-water and water-water energies, 3-body angles, hydrogen bonds,
#energy densities, and all of this as a function of space. Additionally,
#should also be able to compute interfaces, such as Willard-Chandler
#instantaneous interface, or vdW surface, SASA and volume of solute.
#
#Will work with pytraj interface for trajectory analysis, since this
#should later allow easier energy decomposition?
#If doesn't work out, will go back to sim package with netcdf plugin.
#
#Also, should have test script and some test system where know answers
#

import sys, os
import numpy as np
import scipy.optimize as optimize
from scipy.special import sph_harm
import waterlib as wl

#Define constants and unit conversions

#conversion for surface tension
kBJ = 1.38064852*(10**(-23))
temp = 300.0
tomJm2 = kBJ*temp*1000.0*(10**20) #converts kBT/Angstrom^2 to mJ/m^2

#Convert potential energy to kBT
kBTkcal = 0.0019858775*300.0

#Water density
watdens = 0.033456 # molecules or oxygens per Angstrom ^ 3 near 300 K

#Define library of useful functions

def getCosAngs(subPos, Pos, BoxDims, lowCut=0.0, highCut=3.413):
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
               This is useful for finding which angles belong to which central oxygens
               This return value was added on 07/09/2017, so any code using this function
               before then will break, unfortunately, but the fix is easy.
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
