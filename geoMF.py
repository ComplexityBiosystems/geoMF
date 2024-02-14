"""(c) Stefano Zapperi 2023
geoMF a python module to compute the multifractal spectrum of a geopandas dataframe.

It computes: 
 - the scaling exponents of the moments
 - the multifractal spectrum

Uses:
- pandas
- numpy
- matplotlib
- scipy
- tobler
- shapely

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from scipy.optimize import curve_fit
import tobler
import shapely

def create_grid(xmin, ymin, xmax, ymax, cell_size, crs):
    """ The function creates a geodataframe (cell) with a grid with cells of size cell_size within boundaries set by xmin, ymin, xmax, ymax and given crs """
    grid_cells = []
    for x0 in np.arange(xmin, xmax+cell_size, cell_size ):
        for y0 in np.arange(ymin, ymax+cell_size, cell_size):
        # bounds
            x1 = x0-cell_size
            y1 = y0+cell_size
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1)  )
    cell = gpd.GeoDataFrame(grid_cells, columns=['geometry'], crs=crs)
    return cell


def partition_function(gdf,blist,qlist,xmin=2.5e6,xmax=6.5e6,ymin=0.5e6,ymax=5.5e6, variable="emission", crs="EPSG:3035"):
    """ Function to compute the partition function at different scales b specified by blist, for a set of q values specified in qlist. Input a geodataframe gdf for which the specified column (variable) is treated as an extensive variable and interpolated using tobler in each grid whose cell sizes are reported in blist.

Returns a dataframe containing the q-moments for each scale b

"""
    moments=np.zeros((len(blist),len(qlist)))
    i=-1
    for b in blist:
        i=i+1 
        j=-1
# create grid of size b
        cell=create_grid(xmin, ymin, xmax, ymax, b, crs)
# merge into the grid
        area_interp_gdf = tobler.area_weighted.area_interpolate(source_df=gdf, 
                                                            target_df=cell, extensive_variables=[variable]) 
# Plot the coarse grained map
        area_interp_gdf.plot(column=variable, cmap="Greys")
# remove zeros
        area_interp_gdf=area_interp_gdf[area_interp_gdf[variable]>0]
        for q in qlist:
            j=j+1
#        print(i,j)
# compute moment            
            area_interp_gdf["eq"]=area_interp_gdf[variable]**q
            momq=area_interp_gdf["eq"].sum()
            moments[i,j]=momq
    df_q=pd.DataFrame(data=moments, columns=qlist, index=blist)
    df_q=np.log(df_q)
    df_q.index=np.log(df_q.index)
    return df_q

def get_tau_spectrum(df_q,qlist):
    """ Function to compute the list of tau(q) values associated with the dataframe df_q.

Input:
- a dataframe df_q with where the index is the log of the scale, the columns are the values of q
and the values are the log of the partition function.
- the list of q values.

Returns:
- A plot of the fits of the moments.
- the list of tau(q).

"""
    xx=df_q.index
    tau_list=[]
    for q in qlist:
        yy=df_q[q]
        popt5, pcov5 = curve_fit(lambda t, a, b: b*t+a, xx, yy)
        aa,bb= popt5
        tau_list.append(bb)
        plt.scatter(xx, yy)
        plt.plot(xx,bb*xx+aa)
    plt.xlabel(r"$\log b$")
    plt.ylabel(r"$\log Z_b(q)$")
    return tau_list

def plot_tau_spectrum(qlist,taulist):
    """ Function to plot the moment exponents $$\tau(q)$$"""
    plt.scatter(qlist,tau_list)
    plt.ylabel(r"$\tau(q)$")
    plt.xlabel(r"$q$")

    
def get_alpha_spectrum(qlist,tau_list):
    """ Function to compute the multifractal spectrum $$\alpha$$,$$f(\alpha)$$ from
    the list of moments exponents $$\tau(q)$$
    Input: 
    
    - list of q values is qlist
    - list of tau values in tau_list
    
    Returns
    
    - list of alpha.
    - list of f(alpha).
    """ 

    q0=qlist[0]
    tau0=tau_list[0]
    i=0
    alpha_list=[]
    falpha_list=[]
    for q in qlist[1:]:
        i=i+1
        dq=q-q0
        dtau=tau_list[i]-tau0
        alpha=dtau/dq
        falpha=q*alpha-tau_list[i]
        alpha_list.append(alpha)
        falpha_list.append(falpha)
        q0=q
        tau0=tau_list[i]
    return alpha_list,falpha_list