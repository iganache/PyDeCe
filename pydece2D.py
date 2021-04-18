# -*- coding: utf-8 -*-
"""
Created on Mon May 11 11:32:23 2020
A python version of the debrisflow2D code used by Luke

@author: Indujaa
"""

import sys
import os
import argparse
from pathlib import Path
import configparser
import numpy as np
# import math
from scipy import special
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
import rasterio
#import HLLCsolver



def makeGhostCells(arrays, const):
    """  Pads the input arrays with the input constant an all sides and returns the padded array  """     
    if type(arrays) == list:
        arraylist = []
        for arraynd in arrays:
            if arraynd.ndim == 2:
                arraynd = np.pad(arraynd, ((1,1), (1,1)), mode='constant', constant_values=((const,const),(const,const)))
                arraylist.append(arraynd)
            elif arraynd.ndim == 3:
                arraynd = np.pad(arraynd, ((0,0), (1,1), (1,1)), mode='constant', constant_values=((const,const),(const,const),(const,const)))
                arraylist.append(arraynd)
        return arraylist 
    
    elif type(arrays) == np.ndarray:
        arrays =  np.pad(arrays, ((1,1), (1,1)), mode='constant', constant_values=((const,const),(const,const)))
        return arrays
    
def getInput():
    """ Get nameof config file to read (and othre inputs if any) from the command line  '"""
    parser = argparse.ArgumentParser(description="Input paramters")
    
    # # mandatory arguments
    parser.add_argument('--config', dest = 'config', required=True, help = 'Path to config file')
    parser.add_argument('--type', dest = 'type', required=True, help = 'PDC initiation type; choose between collapse and fountain')
    
    
    # # arguments for collapse
    parser.add_argument('--ht', dest = 'height', type=int, default=0, nargs ='?', help = 'initial pile height')
    
    # # arguments for fountaining
    parser.add_argument('--mfps', dest = 'mass', type=int, default=0, nargs ='?', help = 'mass flux per second')
    parser.add_argument('--influx', dest = 'influx_time', type=int, default=0, nargs ='?', help = 'duration of material influx')
    
    # # arguments common to both styles
    parser.add_argument('--lambda', dest = 'pfp', type=float, default=0, nargs ='?', help = 'ratio of excess pore fluid pressure to basal stress 0 to 1')

    return parser.parse_args()

def setFlowparams(Flow, Simulation, PorePressure, args):
    """  Overrides the values input through the config file with command line input  '"""

    if args.pfp != 0:
        PorePressure['initial'] = args.pfp
    if args.mass != 0: 
        Flow['pile_ht'] = args.mass
    if args.influx_time!= 0: 
        Simulation['influx_time'] = args.influx_time
    if args.height != 0: 
        Flow['pile_ht'] = args.height
        
    return Flow, Simulation, PorePressure

def readConfig(file):
    """ Returns four dictionaries - Simulation, Propety, PorePressure, Flow - with values parssed from from the config file to the dictionaries """

    config=configparser.ConfigParser()
    config.read(file)
    
    global DEM  

    outpath = config['path']['out_path']

        
        
    # Simulation parameters
    sim_params=config['Simulation_params']

    Simulation = dict(total_time=np.int(sim_params['total_time']),
                     influx_time=np.int(sim_params['influx_time']),
                     num_frames=np.int(sim_params['num_frames']),
                     dt=np.float(sim_params['timestep']),
                     cstable = np.float(sim_params['cstable']),
                     maxstep = np.float(sim_params['maxstep']),
                     minstep = np.float(sim_params['minstep']),
                     Godunov_order = np.int(sim_params['order']))
    Simulation['hfilm'] = 1e-5
    
    
    # Physical parameters
    phys_params=config['Physical_params']

    Property = dict (gravity=np.float(phys_params['gravity']),
                          viscosity=np.float(phys_params['mu']),
                          d50=np.float(phys_params['d50']),
                          phis=np.float(phys_params['phi']),
                          rhos = np.float(phys_params['rhos']),
                          rhof = np.float(phys_params['rhof']))
    #### Additional properties ###########
    Property['phif'] = 1-Property['phis']
    Property['rho'] = Property['rhof']*Property['phif'] + Property['rhos']*Property['phis']

    # Pore pressure diffusion parameters
    ppd_params=config['PorePressure_params']
    
    PorePressure = dict(diffusivity = np.float(ppd_params['diffusivity']),
                        initial = np.float(ppd_params['lambda']),
                        minimum = np.float(ppd_params['lambda_min']))


    # Pile dimensions and location
    flow_params=config['Flow_params']
   
    Flow = dict(pile_ht=np.float(flow_params['height']),
                centerXgeo=str(flow_params['centerXgeo']),
                centerYgeo=str(flow_params['centerYgeo']),
                vmag=np.float(flow_params['velocity']),
                vdir=np.deg2rad(np.float(flow_params['direction'])))
    
    # DEM parameters
    DEM_params=config['DEM']
    DEM=DEM_params['dem_tiff']
    Flow['centerX'], Flow['centerY'], Flow['numcells'] = setDomain(DEM, Flow['centerXgeo'], Flow['centerYgeo'])                                      # set DEM, nx, ny, dx, center X and center Y

    
    return Simulation, Property, PorePressure, Flow, outpath
    
def setOutputDirectory(homedir, outputdir, Flow, Simulation, PorePressure):
    """ Creates ouput directory and output file names for the storing results """
        
    ### Check if output directory exists ###
    ht = int(Flow['pile_ht'])
    influx = int(Simulation['influx_time'])
    pfp = int(PorePressure['initial'] * 1e2)
    
    if PDCtype == 'fountain':
        vol = int(Flow['numcells'] * dx**2 * ht  * influx * 1e-9)
        dir_name = 'vol' + str(vol) + '_' + 'ht' + str(ht) + '_' + 'influx' + str(influx) + '_' + 'lambda' + str(pfp) 
    elif PDCtype == 'collapse':
        vol = int(Flow['numcells'] * dx**2 * ht  * 1e-9)
        dir_name = 'vol' + str(vol) + '_' + 'ht' + str(ht) + '_' + 'lambda' + str(pfp) 

    
    # # Create a directory path for the simulation
    # # works on both Windows and Linux
    outpath = Path(homedir) / outputdir / dir_name
    outdir = str(outpath)
    
    
    # # Check if directory exists; iff yes, delete all files inside; if not create a new directory
    # # all results from this simulation will be stored inside this directory
    if not outpath.exists():
        os.makedirs(outdir)
    elif outpath.exists():
        [f.unlink() for f in outpath.glob('*') if f.is_file()]

    outdir = outdir+'/'
    return outdir
        
    
def setDomain(file, X, Y):
    """ Reads in an elevation geotiff and assigns elevation values to global variable topo.
    Returns two arrays of pixel coordinates (p, py) of the column centers calculated from arrays of geographic coordinates (X, Y) """
    
    global nx, ny, dx
    global topo
    
    X = np.array(X.split(','), dtype=np.float32)
    Y = np.array(Y.split(','), dtype=np.float32)

    px = np.zeros((1,len(X)), dtype=int)
    py = np.zeros((1,len(Y)), dtype=int)
    
    ds = rasterio.open(file)
    nrows = ds.height
    ncols = ds.width
    nx = ncols
    ny = nrows

    topo = ds.read(1)
    topo = makeGhostCells(topo, 0)
    
    dx = ds.transform[0]
    dy = ds.transform[4]
    
    for i in range(len(X)):
        px[0,i] = int((X[i] - ds.transform[2])//dx)
        py[0,i] = int((Y[i] - ds.transform[5])//dy)

    ds.close()
    
    coords = np.transpose(np.concatenate((px, py), axis=0))
    coords_uniq = np.unique(coords, axis = 0)
    
    return coords_uniq[:,0], coords_uniq[:,1], coords_uniq.shape[0]
    
    
def setInitialConditions(Flow, Property):
    
    """ Initializes 2D arrays for variables h, hu, hv, velocity, velocity in x&y, slopes in x&y, gravitational forcing in x&y.
    Returns 3D array of (h,hu,hv), u, v, velocity, slope in x, slope in y, slope, gx, gy, gz, sx, sy"""

    ### Intitalizing flow thickness and velocity ######################
    ht = np.zeros((ny, nx), dtype = np.float64)
    vel = np.zeros((ny, nx), dtype = np.float64)
    velx = np.zeros((ny, nx), dtype = np.float64)
    vely = np.zeros((ny, nx), dtype = np.float64)
    xvelocity=Flow['vmag']*np.cos(Flow['vdir'])             # Initial flow velocity in X dir [m/s^1]
    yvelocity=Flow['vmag']*np.sin(Flow['vdir'])             # Initial flow velocity in Y dir [m/s^1]
    velocity=Flow['vmag']

    x=np.arange(0,nx,1)
    y=np.arange(0,ny,1)
    [Y,X]=np.meshgrid(y,x,indexing='ij')
    
    for i in range(len(Flow['centerX'])):
        ht[Flow['centerY'][i], Flow['centerX'][i]] = Flow['pile_ht']
        velx[Flow['centerY'][i], Flow['centerX'][i]] = xvelocity
        vely[Flow['centerY'][i], Flow['centerX'][i]] = yvelocity
        vel[Flow['centerY'][i], Flow['centerX'][i]] = velocity
    
    
    ### Initializing uh and uv ############################
    hu = ht * xvelocity
    hv = ht * yvelocity

    
    ht, hu, hv, vel, velx, vely = makeGhostCells([ht, hu, hv, vel, velx, vely], 0)
 
    
    u = np.array([ht, hu, hv])              # Conserved quantity in SWE [h hu hv]
    

    #### Set boundary conditions for topography #############
    ### Upper ###
    topo[0,:] = topo [1,:]
    
    ### Lower ###
    topo[-1,:] = topo [-2,:]
    
    ### Left ###
    topo[:,0] = topo[:,1]
    
    ### Right ###
    topo[:,-1] = topo[:,-2]
    

    slopex = np.gradient(topo, dx, axis = 1)                    # Slope in x direction 
    slopey = np.gradient(topo, dx, axis = 0)                    # Slope in y direction 
    slope = np.sqrt(slopex **2 + slopey **2 + 1)                # True slope. Adding 1 to avoid division by zero
    gy = Property['gravity'] * slopey / slope                   # Gravitational forcing in y direction
    sy = - gy * u[0,:,:]                                        # Sorce term for momentum equation in y
    gx = Property['gravity'] * slopex / slope                   # Gravitational forcing in x direction
    sx = - gx * u[0,:,:]
    gz = Property['gravity'] / slope                            # Sorce term for momentum equation in x
    
    
    # # only needed for fountaining; not ud=sed for collapse
    influx_u = u.copy()
    influx_u[1,:,:] *= 0
    influx_u[2,:,:] *= 0
    
    maxdepth = u[0,:,:].copy()
    maxvel = vel.copy()  


    return u, velx, vely, vel, slopex, slopey, slope, gx, gy, gz, sx, sy, influx_u, maxdepth, maxvel


def setTempVariables(nn, ny, nx, dtype):
    """ Initializing variables whose values are reassigned in every time step.
    The values of the temp variables are reset to 0 at the end of each time step.
    """
    #### Variables for storing shit
    
    h = np.zeros((ny, nx), dtype = np.float64)
        
    #### Sourrce terms
    betax = np.zeros((ny, nx), dtype = np.float64)
    betay = np.zeros((ny, nx), dtype = np.float64)
    
    #### Left and right quantities
    ulx = np.zeros((ny,nx-1), dtype = np.float64)
    urx = np.zeros((ny,nx-1), dtype = np.float64)
    vlx = np.zeros((ny,nx-1), dtype = np.float64)
    vrx = np.zeros((ny,nx-1), dtype = np.float64)
    
    uly = np.zeros((ny-1,nx), dtype = np.float64)
    ury = np.zeros((ny-1,nx), dtype = np.float64)
    vly = np.zeros((ny-1,nx), dtype = np.float64)
    vry = np.zeros((ny-1,nx), dtype = np.float64)
    
    #### Wave speeds
    S_x = np.zeros((nn,ny,nx-1), dtype = np.float64)
    S_y = np.zeros((nn,ny-1,nx), dtype = np.float64)
    Smaxtemp = np.zeros((ny, nx), dtype = np.float64)
    
    #### Fluxes
    tfluxx = np.zeros((nn,ny,nx-1), dtype = np.float64)
    tfluxy = np.zeros((nn,ny-1,nx), dtype = np.float64)

    
    return h, betax, betay, ulx, urx, vlx, vrx, uly, ury, vly, vry, S_x, S_y, Smaxtemp, tfluxx, tfluxy


def createOutputMovieArrays(n, ny, nx, dtype):
    """ Initializing 3D arrays to store snapshots from the model """
    
    movie_h = np.zeros((n, ny, nx), dtype=np.float64)
    movie_hu = np.zeros((n, ny, nx), dtype=np.float64)
    movie_hv = np.zeros((n, ny, nx), dtype=np.float64)
    movie_velx = np.zeros((n, ny, nx), dtype=np.float64)
    movie_vely = np.zeros((n, ny, nx), dtype=np.float64)
    movie_vel = np.zeros((n, ny, nx), dtype=np.float64)
    
    return  movie_h, movie_hu, movie_hv, movie_velx, movie_vely, movie_vel
                 
          
def updateBoundaryConditions(u, velx, vely, gx, gy, gz):

    """ Updates BCs at every time step for topo, gx, gy, gz (transmissive / solid not required)
     Updates BCs for conserved quanitites and velocity based on whether the boundary is solid or transmissive """
    
    LeftSolid = False
    RightSolid = False
    UpperSolid = False
    LowerSolid = False
    
    ## Left boundary 
    gx[:,0] = gx[:,1]
    gy[:,0] = gy[:,1]
    gz[:,0] = gz[:,1]
    vely[:,0] = vely[:,1]
    u[:,:,0] = u[:,:,1]
    if LeftSolid:
       velx[:,0] = - velx[:,1] 
       u[1,:,0] = - u[1,:,1]
    else:
       velx[:,0] = velx[:,1] 
    
    
    ## Right boundary
    gx[:,-1] = gx[:,-2]
    gy[:,-1] = gy[:,-2]
    gz[:,-1] = gz[:,-2]
    vely[:,-1] = vely[:,-2]
    u[:,:,-1] = u[:,:,-2]
    if RightSolid:
       velx[:,-1] = - velx[:,-2] 
       u[1,:,-1] = - u[1,:,-2]
    else:
       velx[:,-1] = velx[:,-2] 
       
       
    ## Upper boundary
    gx[0,:] = gx[1,:]
    gy[0,:] = gy[1,:]
    gz[0,:] = gz[1,:]
    velx[0,:] = velx[1,:] 
    u[:,0,:] = u[:,1,:]
    if UpperSolid:
       vely[0,:] = - vely[1,:]
       u[2,0,:] = - u[2,1,:]
    else:
       vely[0,:] = vely[1,:]
       
     
    ## Lower boundary
    gx[-1,:] = gx[-2,:]
    gy[-1,:] = gy[-2,:]
    gz[-1,:] = gz[-2,:]
    velx[-1,:] = velx[-2,:] 
    u[:,-1,:] = u[:,-2,:]
    if LowerSolid:
       vely[-1,:] = - vely[-2,:]
       u[2,-1,:] = - u[2,-2,:]
    else:
       vely[-1,:] = vely[-2,:]    
    
    return u, velx, vely, gx, gy, gz

def computeFluxes(u, direction, ul, ur, vl, vr, S_xy, gz, hfilm, Smaxtemp, tflux):
    """ Calculates flux across an interface using a HLLC  Reimann solver.
    Returns maximum wave speed, h/hu/hv fluxes, hstar and ustar.
    """
    
    ulr = MUSCLextrap(u, direction)
    [Smaxtemp, tflux, hstar, ustar] = hllc(ulr, ul, ur, vl, vr, S_xy, gz, hfilm, direction, Smaxtemp, tflux)
    
    # [Smaxtemp, tflux, hstar, ustar] = HLLCsolver.hllc(ulr, vl, vr, ul, ur, S_xy, gz, hfilm, direction, Smaxtemp, tflux)               ## Cython version
    
    return Smaxtemp, tflux, hstar, ustar

   

def wavespeeds(ul, hl, hl_noflow, ur, hr, hr_noflow, ustar, hstar, g_z, Sl, Sr, Sm):
    """ Returns the characteristic wave speed of left, middle and right wave """
    
    Sl = np.minimum(ul - np.sqrt(g_z*hl), ustar - np.sqrt(g_z*hstar))
    Sm = ustar.copy()
    Sr = np.maximum(ur + np.sqrt(g_z*hr), ustar + np.sqrt(g_z*hstar))
        
    ### Special cases for dry bed on left / right of cell
    Sl[hl_noflow] = ur[hl_noflow] - (2*np.sqrt(g_z[hl_noflow]*hr[hl_noflow]))
    Sm[hl_noflow] = Sl[hl_noflow]
    Sr[hl_noflow] = ur[hl_noflow] + (np.sqrt(g_z[hl_noflow]*hr[hl_noflow]))
    
    Sl[hr_noflow] = ul[hr_noflow] - (np.sqrt(g_z[hr_noflow]*hl[hr_noflow]))
    Sr[hr_noflow] = ul[hr_noflow] + (2*np.sqrt(g_z[hr_noflow]*hl[hr_noflow]))
    Sm[hr_noflow] = Sr[hr_noflow]
    
    return Sl, Sm, Sr   
    
    
def hllc(ulr, ul, ur, vl, vr, S_temp, g_z, hfilm, direction, Smaxtemp, tflux):
    """ Uses HLLC solver to calculate the fluxes at every cell interface (computation is done across the whole grid and not cell by cell).
    Returns the maximum wave speed, temporary flux, hstar and ustar"""
    
    ## Beware ulr is a 4D array. axis 0, size 2 - conserved quantity at i-1/2 and i+12; axis 1, size 3 - conserved quantities h, hu and hv;
    ## axis 3, size nrows - rows; axis 4, size ncols - values in (row, col)
    
    n = 3
    hl = ulr [0,0,:,:].copy()
    hr = ulr [1,0,:,:].copy() 
    
    if direction == 'x':
        uhl = ulr [0,1,:,:].copy() 
        uhr = ulr [1,1,:,:].copy()
        vhl = ulr [0,2,:,:].copy()
        vhr = ulr [1,2,:,:].copy()
        mx = nx-1
        my = ny
    elif direction == 'y':
        vhl = ulr [0,1,:,:].copy() 
        vhr = ulr [1,1,:,:].copy()
        uhl = ulr [0,2,:,:].copy()
        uhr = ulr [1,2,:,:].copy()
        mx = nx
        my = ny-1
    
    
    hstar = np.zeros((my,mx), dtype = np.float64)
    ustar = np.zeros((my,mx), dtype = np.float64)
    
    Sl = S_temp[0,:,:]
    Sm = S_temp[1,:,:]
    Sr = S_temp[2,:,:]


    
    hl_flow = hl > hfilm
    hl_noflow = hl <= hfilm
    
    hr_flow = hr > hfilm
    hr_noflow = hr <= hfilm
    
    ul[hl > hfilm] = uhl[hl > hfilm] / hl[hl > hfilm]
    vl[hl > hfilm] = vhl[hl > hfilm] / hl[hl > hfilm]
    
    ur[hr > hfilm] = uhr[hr > hfilm] / hr[hr > hfilm]
    vr[hr > hfilm] = vhr[hr > hfilm] / hr[hr > hfilm]
    
    
    hstar = (ul + (2 * np.sqrt(g_z*hl)) - ur - (2*np.sqrt(g_z*hr))) ** 2 / (16 * g_z)
    ustar = (0.5 * (ul + ur)) + np.sqrt(g_z*hl) - np.sqrt(g_z*hr)
         
    
    #### wave speeds ###################
    Sl, Sm, Sr = wavespeeds(ul, hl, hl_noflow, ur, hr, hr_noflow, ustar, hstar, g_z, Sl, Sm, Sr)
           
        
    ### Fluxes #############################
    cond_l = Sl >= 0
    tflux[0,:,:][cond_l] = hl[cond_l] * ul[cond_l]
    tflux[1,:,:][cond_l] = (hl[cond_l] * ul[cond_l] ** 2) + (0.5 * g_z[cond_l] * hl[cond_l] **2)

    
    cond_m = (Sl < 0) & (Sr > 0)
    tflux[0,:,:][cond_m] = (Sr[cond_m]*(hl[cond_m]*ul[cond_m]) - Sl[cond_m]*(hr[cond_m]*ur[cond_m]) + Sl[cond_m]*Sr[cond_m]*(hr[cond_m]-hl[cond_m])) / (Sr[cond_m]-Sl[cond_m])
    tflux[1,:,:][cond_m] = (Sr[cond_m]*(hl[cond_m]*ul[cond_m]**2 + 0.5*g_z[cond_m]*hl[cond_m]**2) - Sl[cond_m]*(hr[cond_m]*ur[cond_m]**2 + 0.5*g_z[cond_m]*hr[cond_m]**2) + Sl[cond_m]*Sr[cond_m]*(hr[cond_m]*ur[cond_m] - hl[cond_m]*ul[cond_m])) / (Sr[cond_m]-Sl[cond_m])

    cond_r = Sr <= 0
    tflux[0,:,:][cond_r] = hr[cond_r] * ur[cond_r]
    tflux[1,:,:][cond_r] = (hr[cond_r] * ur[cond_r] **2) + (0.5 * g_z[cond_r] * hr[cond_r] **  2)

    cond2_gt = Sm >= 0
    cond2_lt = Sm < 0
    tflux[2,:,:][cond2_gt] =  tflux[0,:,:][cond2_gt] * vl[cond2_gt]
    tflux[2,:,:][cond2_lt] = tflux[0,:,:][cond2_lt] * vr[cond2_lt]

        
    ### Maximum wave speeds for stability criteria #####
    Smaxtemp = np.maximum(np.absolute(Sl),np.absolute(Sr), np.absolute(Sm))        
    return Smaxtemp, tflux, hstar, ustar

def MUSCLextrap(u, direction):
    
    """ Monotonic Upwind Scheme blah blah something. Returns value at i-1 and i+1 for 1st order (I think) """
    
    if direction == 'x':
        val = u[:,:,:-1]
        valplus = u[:,:,1:]

        
    if direction == 'y':
        val = u[:,:-1,:]
        valplus = u[:,1:,:]

        
    return np.array([val, valplus])

def AdaptiveTimeStep(Simulation, t, Smax, dx):
    """ Checks if the CFL criterion issatisfied at each time step.
    The time step duration is adjusted accordingly to smaller values till CFL is satisfied.
    """
    Simulation['dt'] = Simulation['cstable'] * dx / Smax
    
    if Simulation['dt'] > Simulation['maxstep']:
        Simulation['dt'] = Simulation['maxstep']
    if t+Simulation['dt'] > Simulation['total_time']:
        Simulation['dt'] = Simulation['total_time'] - t
    courant = Smax * Simulation['dt'] / dx
    
    return Simulation['dt'], courant


def calculateLambda(t, h, hfilm, PorePressure):
    """ Calculates pore fluid pressure to baal normal stress ratio (lambda).
    lambda is not constant since diffusion of pore pressure is incorporated. 
    lambda needs to recalculated at eac time step.
    """
    ppd = makeGhostCells(np.zeros((ny,nx)),0)
    ppd[h>hfilm] = 1
    
    ##### Pore pressure diffusivity #################################
    if t == 0: 
        pore_lambda = PorePressure['initial'] * ppd
    else:
        pore_lambda = np.maximum(PorePressure['minimum'] * ppd, (PorePressure['initial'] * (1 - 2 * special.erfc(h/ np.sqrt(4 * PorePressure['diffusivity'] * t)))))
        
    return pore_lambda
    
def JopModel(h, velx, vely, gz, Property, pfp):
    """ Calculates coefficient of basal friction using the mu(I) rheology.
    See Jop et al. (20xx).
    """
    mu1 = np.tan(np.deg2rad(21))
    mu2 = np.tan(np.deg2rad(33))
    I0 = 0.29
    
    velx[h < 0.01] = 0
    vely[h < 0.01] = 0
    vel = np.sqrt(velx ** 2 + vely ** 2)
    
    I = makeGhostCells(np.zeros((ny,nx), dtype= np.float32),0)
    mu = makeGhostCells(np.zeros((ny,nx), dtype= np.float32),0)
    
    I[h<.01] = np.inf
    I[h>=.01] = 2 * vel[h>=.01] * Property['d50'] / h[h>=.01] / np.sqrt(Property['rho'] * gz[h>=.01] * h[h>=.01] * (1 - pfp[h>=.01]) / Property['rhos'])
    mu[I==0] = mu1
    mu[I!=0] = mu1 + ((mu2 - mu1) / (I0/I[I!=0] + 1))
    
    return mu 

def ResistanceTerms(t, u, velx, vely, gz, hfilm, PorePressure, Property, betax, betay):
    """
    Computes the total (fluid + solid) resistance terms for the x and y momentum conservation eqns.
    """
    h = u[0,:,:].copy()
    
    ## Pore pressure diffusion
    pfp = calculateLambda(t, h, hfilm, PorePressure)
    
    ##### Basal friction - Jop model #################################
    mu = JopModel(h, velx, vely, gz, Property, pfp)
    
    ### Fluid resistance
    betax[h > .01] += -1 * 2 * velx[h > .01] * Property['viscosity'] * Property['phif'] / Property['rho'] / h[h > .01]
    betay[h > .01] += -1 * 2 * vely[h > .01] * Property['viscosity'] * Property['phif'] / Property['rho'] / h[h > .01]
    
   
    ### Solid resistance
    betax[velx > 0] += -1 * (1 - pfp[velx > 0]) * gz[velx > 0] * h[velx > 0] * mu[velx > 0]
    betax[velx <=0] +=  (1 - pfp[velx <=0]) * gz[velx <=0] * h[velx <=0] * mu[velx <=0]
    
    betay[vely > 0] += -1 * (1 - pfp[vely > 0]) * gz[vely > 0] * h[vely > 0] * mu[vely > 0]
    betay[vely <= 0] += (1 - pfp[vely <= 0]) * gz[vely <= 0] * h[vely <= 0] * mu[vely <= 0]
    
    return betax, betay


    
def GravitationForcing(g, h):
    """
    Computes the gravitation source term.
    """
    return -g * h


def updateCV(CV, dt, dx, fluxx, fluxy, gravityterm=0, resistanceterm=0):
    """
    Updates the value of the conserved variables (h, hu and hv) at each time step.
    """
    oneoverdx = 1/dx
    CV_new = CV - dt * oneoverdx * np.diff(fluxx, axis=-1) - dt * oneoverdx * np.diff(fluxy, axis=-2) + dt * gravityterm + dt * resistanceterm
        
    return CV_new
          

def writeTIFF(arr, outfile, nodata, count, minvalue = 1e-6):
    """
    Writes all the model output to GeoTIFF  files.
    The geotiffs have the same spatial extent and reference as the inout DEM.

    """
    
    if nodata != None:
        arr[np.abs(arr)<minvalue] = nodata
    
    
    ds = rasterio.open(DEM)
    
    outds = rasterio.open(outfile, 'w', driver='GTiff', 
                      height = ds.height, 
                      width = ds.width, 
                      count=count, 
                      crs = ds.crs, 
                      dtype = arr.dtype,
                      transform = ds.transform,
                      nodata = nodata)
    
    if arr.ndim > 2 and count == arr.shape[0]:
        for i in range(count):
            outds.write(arr[i], i+1)
    elif arr.ndim == 2 and count == 1:      
        outds.write(arr, count)
        
    outds.close()
    ds.close()
    
def OutputVariables(u, velx, vely, maxdepth, maxvel, hfilm):
    """
    Calculates flow height, flow velocity x and flow velocity y.
    Hard-coded minimum threshold values are used to restrict thickness and velocity to meaningful values.
    """
    h = u[0,:,:].copy()
    u[0,:,:][h<=hfilm] = hfilm
    u[1,:,:][h<=hfilm] = 0
    u[2,:,:][h<=hfilm] = 0    
    
    
    velx[h>hfilm] = u[1,:,:][h>hfilm]  / u[0,:,:][h>hfilm] 
    vely[h>hfilm] = u[2,:,:][h>hfilm]  / u[0,:,:][h>hfilm]
    velx[h<=0.01] = 0
    vely[h<=0.01] = 0
    
    vel = np.sqrt(velx ** 2 + vely ** 2)
    
    maxvel = np.maximum(maxvel, vel)
    maxdepth = np.maximum(maxdepth, u[0,:,:])
    
    return u, velx, vely, vel, maxdepth, maxvel

def clearVariables(arr_list):
    """
    Reset all the temporary variables used within a timestep to 0.
    """
    new_arr_list = []
    for arr in arr_list:
        arr *= 0
        new_arr_list.append(arr)

    return new_arr_list 

    
def main():

    ## Number of conserved quantities
    nn = 3
    
    ## Get arguments from command line
    args = getInput()
    
    ## read values from config file 
    config_file = args.config
    
    ## set initiation type
    global PDCtype
    PDCtype = args.type
    
    Simulation, Property, PorePressure, Flow, outpath = readConfig(config_file)   
    
    Flow, Simulation, PorePressure = setFlowparams(Flow, Simulation, PorePressure, args)                    ## Set user-defined parameters if supplied
    
    # outdir = setOutputDirectory(outpath, Flow, Simulation)               ## Set Output directory
    outdir = setOutputDirectory(os.getcwd(), config_file.split('.')[0], Flow, Simulation, PorePressure)               ## Set Output directory
        
    
    ##### Set Initial Conditions ######
    u, velx, vely, vel, slopex, slopey, slope, gx, gy, gz, sx, sy, influx_u, maxdepth, maxvel = setInitialConditions(Flow, Property)
     
       
    ### Arrays to hold movie time series ############################
    movie_idx = 0    
    movie_h, movie_hu, movie_hv, movie_velx, movie_vely, movie_vel = createOutputMovieArrays(Simulation['num_frames'], ny, nx, np.float64)
    printinterval = Simulation['total_time'] / Simulation['num_frames']
    printat = printinterval
    
    
    h, betax, betay, ulx, urx, vlx, vrx, uly, ury, vly, vry, S_x, S_y, Smaxtemp, tfluxx, tfluxy = setTempVariables(nn, ny, nx, np.float64)
    
    ulx, urx, vlx, vrx, uly, ury, vly, vry, tfluxx, tfluxy, betax, betay = makeGhostCells([ulx, urx, vlx, vrx, uly, ury, vly, vry, tfluxx, tfluxy, betax, betay], 0)
    
    
    oneoverdx = 1/dx
    t = 0

    while t < Simulation['total_time']:

       
        u, velx, vely, gx, gy, gz = updateBoundaryConditions(u, velx, vely, gx, gy, gz)
        
        Smax = 0
        
        ## Compute fluxes #################################
        ### Flux in X-direction ###########################
        [Smaxtemp, tfluxx, hstarx, ustarx] = computeFluxes(u, 'x', ulx, urx, vlx, vrx, S_x, gz[:,:-1], Simulation['hfilm'], Smaxtemp, tfluxx)
        Smax = max(Smax, np.max(Smaxtemp))
        

        ### Flux in Y-direction ###########################
        [Smaxtemp, tfluxy, hstary, ustary] = computeFluxes(u, 'y', vly, vry, uly, ury, S_y, gz[:-1,:], Simulation['hfilm'], Smaxtemp, tfluxy)
        Smax = max(Smax, np.max(Smaxtemp))

      
        ##### SOURCE TERMS ##########################################################################################################
        betax, betay = ResistanceTerms(t, u, velx, vely, gz, Simulation['hfilm'], PorePressure, Property, betax, betay)
        
        sx = GravitationForcing(gx, u[0,:,:])
        sy = GravitationForcing(gy, u[0,:,:])
              
        
                  
        ### Adaptive time step #######################
        Simulation['dt'], courant = AdaptiveTimeStep(Simulation, t, Smax, dx)   
        if courant > 1:
            print("CFL error at t = ", t)
            t = 2 * Simulation['total_time']
            
        t += Simulation['dt']


        ### Time update ################################
        ## use np.diff with tfluxx and tfluxy ##
        if PDCtype == 'fountain' and t < Simulation['influx_time']:
            u[0,1:-1,1:-1] = updateCV(u[0,1:-1,1:-1], Simulation['dt'], dx, tfluxx[0,1:-1,:], tfluxy[0,:,1:-1], influx_u[0,1:-1,1:-1] )               
        elif PDCtype == 'fountain' and t >= Simulation['influx_time']:   
            u[0,1:-1,1:-1] = updateCV(u[0,1:-1,1:-1], Simulation['dt'], dx, tfluxx[0,1:-1,:], tfluxy[0,:,1:-1])
        elif PDCtype == 'collapse':
            u[0,1:-1,1:-1] = updateCV(u[0,1:-1,1:-1], Simulation['dt'], dx, tfluxx[0,1:-1,:], tfluxy[0,:,1:-1])                            # # same format as fountaining with t > t_influx
            
        u[1,1:-1,1:-1] = updateCV(u[1,1:-1,1:-1], Simulation['dt'], dx, tfluxx[1,1:-1,:], tfluxy[2,:,1:-1], sx[1:-1,1:-1], betax[1:-1,1:-1])
        u[2,1:-1,1:-1] = updateCV(u[2,1:-1,1:-1], Simulation['dt'], dx, tfluxx[2,1:-1,:], tfluxy[1,:,1:-1], sy[1:-1,1:-1], betay[1:-1,1:-1])
               
                   
        
        u, velx, vely, vel, maxdepth, maxvel = OutputVariables(u, velx, vely, maxdepth, maxvel, Simulation['hfilm'])
        
        #### Movie array update ####
        if t >= printat:
            print(t, courant, Simulation['dt'])
            printat += printinterval

            movie_h[movie_idx,:,:] = u[0,1:-1,1:-1]
            movie_hu[movie_idx,:,:] = u[1,1:-1,1:-1]
            movie_hv[movie_idx,:,:] = u[2,1:-1,1:-1]
            movie_vel[movie_idx,:,:] = vel[1:-1,1:-1]
            movie_velx[movie_idx,:,:] = velx[1:-1,1:-1]
            movie_vely[movie_idx,:,:] = vely[1:-1,1:-1]
            
            movie_idx += 1
            
           
            
        #### Clear all variables
        [betax, betay, ulx, urx, vlx, vrx, uly, ury, vly, vry, S_x, S_y, Smaxtemp, tfluxx, tfluxy] = clearVariables([betax, betay, ulx, urx, vlx, vrx, uly, ury, vly, vry, S_x, S_y, Smaxtemp, tfluxx, tfluxy])
        



    #### Final results as binary files #########
    writeTIFF(u[0,1:-1,1:-1], outdir+"Depth.tif", -9999, 1, .001)
    writeTIFF(u[1,1:-1,1:-1], outdir+"MomentumX.tif", 0, 1, 0)
    writeTIFF(u[2,1:-1,1:-1], outdir+"MomentumY.tif", 0, 1, 0)
    writeTIFF(vel[1:-1,1:-1], outdir+"Velocity.tif", 0, 1, .001)
    writeTIFF(velx[1:-1,1:-1], outdir+"VelocityX.tif", 0, 1, .001)
    writeTIFF(vely[1:-1,1:-1], outdir+"VelocityY.tif", 0, 1, .001)
    writeTIFF(maxdepth[1:-1,1:-1], outdir+"MaxDepth.tif", -9999, 1, 1e-5)
    writeTIFF(maxvel[1:-1,1:-1], outdir+"MaxVelocity.tif", 0, 1)
    
    
    #### Final results as movie files #########
    writeTIFF(movie_h, outdir+"DepthMovie.tif", -9999, Simulation['num_frames'], .001)
    writeTIFF(movie_hu, outdir+"MomentumXMovie.tif", 0, Simulation['num_frames'], 0)
    writeTIFF(movie_hv, outdir+"MomentumYMovie.tif", 0, Simulation['num_frames'], 0)
    writeTIFF(movie_velx, outdir+"VelocityXMovie.tif", 0, Simulation['num_frames'], .001)
    writeTIFF(movie_vely, outdir+"VelocityYMovie.tif", 0, Simulation['num_frames'], .001)
    writeTIFF(movie_vel, outdir+"VelocityMovie.tif", 0, Simulation['num_frames'], .001)
    
    ## TIFF file for classification
    u_bin = u[0,1:-1,1:-1].astype(np.int8)
    u_bin[u_bin >= 0.1] = 1                                     ## binary raster with depth>10cm = 1 and depth<10cm = 0
    u_bin[u_bin < 0.1] = 0
    writeTIFF(u_bin, outdir+"Depthbin.tif", None, 1)


        
if __name__ == '__main__':
    main()   

