#!/bin/env python

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import math
import time
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord


class cluster:
    """
        Class hold information about the cluster and provides several methods
        required to perform the clustering

        init_center_wcs: initilizes the WCS which is used to compute the distance.
        create_cluster: Creates a new cluster and sets center co-ordinates and corresponding flux
        add_pixel: Adds a pixel to existing cluster
        compute_distnace: Computes distance of a pixel from center of the cluster


    """

    def __init__(self):
        self.center_x=0
        self.center_y=0
        self.center_flux=0
        self.numpts=0
        self.x_pixels=[]
        self.y_pixels=[]
        self.pixel_fluxs=[]

    def init_center_wcs(self):
        """
        Initializes WCS required to compute distance
        """
        self.c1 = SkyCoord(self.w.pixel_to_world(self.center_x,self.center_y))

    def create_cluster(self,x,y,flux,w):
        """
        Creates a new cluster

        Parameters:
        -----------
        x   :   int
            X position in of the center of cluster

        y   :   int
            Y position in of the center of cluster

        flux:   float or double
            flux at the center of the cluster
            
        """
        self.x_pixels.append(x)
        self.y_pixels.append(y)
        self.pixel_fluxs.append(flux)
        self.numpts=1
        self.center_x=x
        self.center_y=y
        self.center_flux=flux
        self.w=w
        self.init_center_wcs()

    def add_pixel(self,x,y,flux):
        """
        adds pixel to the existing cluster

        Parameters:
        ----------
        x   :   int
            X position of the pixel to be added to the cluster

        y   :   int
            Y position of the pixel to be added to the cluster

        flux:   float or double
            flux of the pixel to be added to the cluster

        """
        self.x_pixels.append(x)
        self.y_pixels.append(y)
        self.pixel_fluxs.append(flux)
        self.numpts=self.numpts+1

        #check how to update wcs
        if(flux>self.center_flux):
            self.center_flux=flux
            self.center_x=x
            self.center_y=y
            self.init_center_wcs()
             
    def compute_distance(self,x,y,all_pix=False):
        """
        computes distance of a pixel from the center

        Parameters:
        -----------

        x    : int
            X position of the pixel for which distance needs to computed

        y    : int
            Y position of the pixel for which distance needs to computed

        all_pix : bool, default false
                computes distance from all pixels and returns the minimum distance

        Returns:
        --------
        distance    : float
                    returns distance of the pixel from center or the shortest distance 
                    from other pixels
        """
       
        if all_pix:
            dmin=100;
            c2 = SkyCoord(self.w.pixel_to_world(x,y))
            for i in np.arange(0,len(self.x_pixels)):
                c1=SkyCoord(self.w.pixel_to_world(self.x_pixels[i],self.y_pixels[i]))
                d=self.c1.separation(c2)
            
                if d.arcmin<dmin:
                    dmin=d.arcmin
            return dmin
        else:
            c2 = SkyCoord(self.w.pixel_to_world(x,y))
            d=self.c1.separation(c2)
            return d.arcmin
       
        #check distance with all points in the cluster
        #can be run with a flag to compute distance
        #the execution time will increase
        """
        """

class utils: 
    """
    Class containing utility methods

    Methods:
    init_wcs :  Initializes the world co-ordinate system.
    """
    def init_wcs(self,img_header,dropaxis=False):
        """
        Initializes the world co-ordinate system (wcs). 
        
        Parameters
        ----------
        img_header  : fits image header 
                    The fits image header will be used to inilialize the WCS

        dropaxis    : bool, optional
                    If set to true, will delete two axis from the existing WCS. Mainly, required to process 
                    images from VLASS
                    

        Returns
        -------
        w   : WCS
            Returns object of initialized WCS

        """
        w=WCS(img_header)
        if dropaxis:
            w=w.dropaxis(2)
            w=w.dropaxis(2)    
        return w


class Sources:
    """
    The class holds informationa about the source image and performs
    source detction using clustering algorithm.

    Methods:
    -------
    read_input: Method to read the input files, image and RMS
    flag_pixels: Flags all pixels which are above the threashold*RMS.
    find_clusters: Performs clustering on the flagged pixels
    merge_clusters: Perform clustering on the clusters, this is mainly 
    to find extended sources
    plot_image: Plots the image and clusters on it.
    """
    def __init__(self):
        """
        Initilizes the cluster object
        """
        self.infile="images/VLASS1.2.se.T01t08.J033228-363000.06.2048.v1.spw6.I.iter3.image.pbcor.tt0.subim.fits"
        self.rmsfile="images/VLASS1.2.se.T01t08.J033228-363000.06.2048.v1.spw9.I.iter3.image.pbcor.tt0.rms.subim.fits"
        self.dist_threash=0.7 #in arcmin

        self.img_hdu=None
        self.img_data=None
        self.img_header=None
        self.rms_hdu=None
        self.rms_data=None
        self.cluster_list=[]
        self.numclusters=0
    
    def read_input(self,hdunum=0): 
        """
        Reads and stores the input image and RMS file. Assumption is both image and RMS
        files are following same format, Ex: Number of extensions in the file, data dimensions.

        Parameters:
        ----------
        hdunum  : int, optional, deafult is 0
                fits extension number to read data from. 
        """
        self.img_hdu=fits.open(self.infile)
        self.img_data=self.img_hdu[hdunum].data[0][0]
        self.img_header=self.img_hdu[hdunum].header
        self.rms_hdu=fits.open(self.rmsfile)
        self.rms_data=self.rms_hdu[hdunum].data[0][0]

    def flag_pixels(self,threashold=5):
        """
        Flags all the pixels which are above the threashold*rms.

        Parameters:
        ----------
        threashold  : int, optional, default is 5
                    All pixels above the threashold will be flagged
        """
        self.flags=np.zeros(np.shape(self.img_data))
        self.flags[np.where(self.img_data>threashold*self.rms_data)]=1
        self.tflags=np.where(self.flags==1)
        self.xcoor=self.tflags[0]
        self.ycoor=self.tflags[1]

    
    def find_clusters(self):
        """
        Method performs clustering on the flagged pixels and groups pixels which are 
        closer.

        Parameters:
        -----------
        dist_threash    : float, optional default 0.8 arc min
                        If distance between two pixels is less than the dist_threash,
                        then these two pixels are grouped in a single cluster.
        """
        u=utils()
        w=u.init_wcs(self.img_header,True)

        for i in np.arange(0,len(self.xcoor)):
            if(self.numclusters==0):
                #create a new cluster
                c=cluster()
                c.create_cluster(self.xcoor[i],self.ycoor[i],self.img_data[self.xcoor[i]][self.ycoor[i]],w)
                self.cluster_list.append(c)
                self.numclusters=self.numclusters+1
            else:
                isExistingCluster=False
                for c in self.cluster_list:
                    if(c.compute_distance(self.xcoor[i],self.ycoor[i])<self.dist_threash):
                        c.add_pixel(self.xcoor[i],self.ycoor[i],self.img_data[self.xcoor[i]][self.ycoor[i]])
                        isExistingCluster=True
                        break

                #check distance from center of each cluster
                #if distance <threashold: 
                #add point to the same cluster
                #else create new cluster
                if isExistingCluster==False:
                    #can be moved as a function
                    c=cluster()
                    c.create_cluster(self.xcoor[i],self.ycoor[i],self.img_data[self.xcoor[i]][self.ycoor[i]],w)
                    self.cluster_list.append(c)
                    self.numclusters=self.numclusters+1


    def merge_clusters(self):
        """
        Performs clustering on the cluster. This method tries to group all the clusters which are from
        extended source

        Parameters:
        -----------
        dist_threash    : float, optional default 0.8 arc min
                        If distance between two pixels is less than the dist_threash,
                        then these two pixels are grouped in a single cluster.

        """
        allClusterMerged = False;
        #while not allClusterMerged:
        #List all clusters and cluster merging if required. 
        for t in np.arange(0,1):
            print("Running cluster merging",t);
            pop_clusters=[]
            for i in np.arange(0,len(self.cluster_list)):
                if i not in pop_clusters:
                    for j in np.arange(i+1,len(self.cluster_list)):
                        if j not in pop_clusters:
                            d=self.cluster_list[i].compute_distance(self.cluster_list[j].center_x,self.cluster_list[j].center_y)
                            if d<self.dist_threash:
                                #create function merge clusters
                                self.cluster_list[i].x_pixels.extend(self.cluster_list[j].x_pixels)
                                self.cluster_list[i].y_pixels.extend(self.cluster_list[j].y_pixels)
                                self.cluster_list[i].pixel_fluxs.extend(self.cluster_list[j].pixel_fluxs)
                                self.cluster_list[i].numpts+=self.cluster_list[j].numpts
                                #update the center
                                if self.cluster_list[i].center_flux<self.cluster_list[j].center_flux:
                                    self.cluster_list[i].center_flux=self.cluster_list[j].center_flux
                                    self.cluster_list[i].center_x=self.cluster_list[j].center_x
                                    self.cluster_list[i].center_y=self.cluster_list[j].center_y
                                    self.cluster_list[i].init_center_wcs()    

                                pop_clusters.append(j)

            for i in  np.arange(0,len(pop_clusters)):
                self.cluster_list.pop(pop_clusters[i]-i)
    

            allClusterMerged=True
            pop_clusters.clear()




    def plot_image(self):
        """
        Plots the image and clusters on top of the image.
        """
        plt.imshow(self.img_data,interpolation="nearest",vmin=np.min(self.img_data),vmax=np.max(self.img_data),origin="lower",cmap="gray")

        #plt.imshow(flags,interpolation="nearest",vmin=0,vmax=1,origin="lower",cmap="gray")
        for c in self.cluster_list:
            plt.plot(c.y_pixels,c.x_pixels,markersize=10)
        plt.show()
        
"""


#img_data[img_data<5*rms_data]=-1
#filterd_data=np.where(img_data>5*rms_data)


"""


if __name__ == "__main__":

    start_time = time.time()
    s=Sources();
    s.read_input()
    s.flag_pixels()
    s.find_clusters()
    s.merge_clusters()
    s.plot_image()
    print("--- %s seconds ---" % (time.time() - start_time))

