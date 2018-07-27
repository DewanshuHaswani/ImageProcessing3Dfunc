#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 11:01:34 2018

@author: haswani
"""

import math
from skan import csr
import numpy as np
from skimage.morphology import skeletonize,skeletonize_3d
from astropy.io import fits
from skimage.morphology import binary_opening,binary_dilation,binary_erosion,closing,dilation,opening,ball,remove_small_holes
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from convert2subgraph import returnSubgraph,returnEndpointsAndJunction
import networkx as nx
from ImageProcessing3Dfunc import skeletonize3D ,pruning3D, testGraph1, testGraph2, plotAndReturnLongestFilaments 
# import imageio as iio
# from skan import draw

cube = fits.getdata('ngc6744_co21_12m+7m+tp_pbcorr_round_k_correct_mask.fits')
#cube = fits.getdata('ngc3627_co21_12m+7m+tp_mask.fits')

dskel2, pixel_graph0, coordinates0, degrees0 = skeletonize3D(cube)

parameter = 5

G = pruning3D(pixel_graph0,coordinates0,parameter)

Gtest1 = testGraph1()
Gtest2 = testGraph2()

#subcube = dskel2[125:175,300:500, 50:200]
#pixel_graph1, coordinates1, degrees1 = csr.skeleton_to_csgraph(subcube)

LongestFilamentGraph = plotAndReturnLongestFilaments(G)
