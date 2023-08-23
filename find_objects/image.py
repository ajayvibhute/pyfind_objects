
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
#from . import utils as ut
#from . import cluster as cl


import utils as ut
import cluster as cl

class Image:
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
    def __init__(self,imgfile,rmsfile,dist_threash=0.7):
        """
        Initilizes the cluster object
        """
        #self.infile="/Users/avibhute/NRAO/images/VLASS1.2.se.T01t08.J033228-363000.06.2048.v1.spw6.I.iter3.image.pbcor.tt0.subim.fits"
        #self.rmsfile="/Users/avibhute/NRAO/images/VLASS1.2.se.T01t08.J033228-363000.06.2048.v1.spw9.I.iter3.image.pbcor.tt0.rms.subim.fits"
        #self.dist_threash=0.7 #in arcmin
        
        self.infile = imgfile
        self.rmsfile = rmsfile
        self.dist_threash = 40 #dist_threash
        self.img_hdu = None #hdu of fits image file
        self.img_data = None #image data
        self.img_header = None #fits image header
        self.rms_hdu = None #hdu of the rms file
        self.rms_data = None #rms data
        self.cluster_list = [] #list of clusters   
        self.numclusters = 0  #number of clusters
        self.img_freq    = None 



    def read_input(self,hdunum=0): 
        """
        Reads and stores the input image and RMS file. Assumption is both image and RMS
        files are following same format, Ex: Number of extensions in the file, data dimensions.

        Parameters:
        ----------
        hdunum  : int, optional, deafult is 0
                fits extension number to read data from. 
        """
        #reading fits hdu, data of image file
        self.img_hdu = fits.open(self.infile)
        self.img_data = self.img_hdu[hdunum].data[0][0]
        self.img_header = self.img_hdu[hdunum].header
        self.img_freq = self.img_header["RESTFRQ"] 
        #assuming key will alway present, may want to add condition to check presence of key 



        #reading data from rms file
        self.rms_hdu = fits.open(self.rmsfile)
        self.rms_data = self.rms_hdu[hdunum].data[0][0]

    def flag_pixels(self,threashold=5):
        """
        Flags all the pixels which are above the threashold*rms.

        Parameters:
        ----------
        threashold  : int, optional, default is 5
                    All pixels above the threashold will be flagged
        """

        #setting flags to zero
        self.flags=np.zeros(np.shape(self.img_data))

        #setting pixels above the threashold to 1
        self.flags[np.where( (self.img_data>threashold*self.rms_data) & (self.img_data >0)) ]=1
        
        #getting list of flagged pixels
        self.tflags=np.where(self.flags==1)

        #getting X and Y coordinates of the flagged pixels. These pixels are 
        #used to perform clustering 
        self.xcoor=self.tflags[0]
        self.ycoor=self.tflags[1]

    
    def find_clusters(self):
        """
        Method performs clustering on the flagged pixels and groups pixels which are 
        closer.

        Parameters:
        -----------
        dist_threash    : float, optional default 0.7 arc min
                        If distance between two pixels is less than the dist_threash,
                        then these two pixels are grouped in a single cluster.
        """
        #creating object of the utils to init the WCS
        u=ut.utils()
        self.w=u.init_wcs(self.img_header,True)

        #iterating through the x-pixels
        for i in np.arange(0,len(self.xcoor)):
            #if no cluster is created then create a new
            #cluster. Generally, iteration 0
            if(self.numclusters==0):
                #create a new cluster and assign
                #x, and y positions at xcoor[i], ycoor[i] to the
                #center of the cluster.
                c=cl.cluster()
                c.create_cluster(self.xcoor[i],self.ycoor[i],self.img_data[self.xcoor[i]][self.ycoor[i]],self.w)
                self.cluster_list.append(c)
                self.numclusters=self.numclusters+1
            else:
                #check if a pixel is part of existing cluster
                # A pixel is considered as a part of the existing cluster if distance of the
                #pixel from the center is less than the pre-set threashold. 
                isExistingCluster=False
                for c in self.cluster_list:
                    if(c.compute_distance(self.xcoor[i],self.ycoor[i])<self.dist_threash):
                        c.add_pixel(self.xcoor[i],self.ycoor[i],self.img_data[self.xcoor[i]][self.ycoor[i]])
                        isExistingCluster=True
                        break

                #if pixel is not part of the existing cluster, then create a new
                #cluster and assign xcoor[i],ycoor[i] (pixel under process) to the
                #center of new cluster.
                if isExistingCluster==False:
                    #can be moved as a function
                    c=cl.cluster()
                    c.create_cluster(self.xcoor[i],self.ycoor[i],self.img_data[self.xcoor[i]][self.ycoor[i]],self.w)
                    self.cluster_list.append(c)
                    self.numclusters=self.numclusters+1


        pop_clusters=[]
        #cleanup clusters
        for i in np.arange(0,len(self.cluster_list)):
            if self.cluster_list[i].numpts<5:
                pop_clusters.append(i)

        for i in  np.arange(0,len(pop_clusters)):
            self.cluster_list.pop(pop_clusters[i]-i)
            self.numclusters-=1



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
        #check if all clusters are merged if not perform cluster merging again
        allClusterMerged = False;
        #while not allClusterMerged:
        #List all clusters and cluster merging mif required. 
        #added for loop just for testing purpose, will be replaced by a while 
        for t in np.arange(0,1):
            #list of clusters merged with other cluster, and will be removed from
            #the cluster_list once merged
            pop_clusters=[]

            #iterate through all the clusters
            for i in np.arange(0,len(self.cluster_list)):
                #check if cluster is already merged with other cluster    
                if i not in pop_clusters:
                    #iterate through all clusters with which current cluster
                    #can possibly be merged
                    for j in np.arange(i+1,len(self.cluster_list)):
                        #check if cluster is already merged with other cluster
                        if j not in pop_clusters:
                            #compute distance of between two clusters
                            d=self.cluster_list[i].compute_distance(self.cluster_list[j].center_x,self.cluster_list[j].center_y)
                            
                            #if distance is less than threashold, merge the clusters 
                            #and add it's entry to the list of clusters tobe popped
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

            #removing the merged clusters
            for i in  np.arange(0,len(pop_clusters)):
                self.cluster_list.pop(pop_clusters[i]-i)
                self.numclusters-=1
    
            
            allClusterMerged=True
            #clearing the pop cluster list
            pop_clusters.clear()




    def plot_image(self):
        """
        Plots the image and clusters on top of the image.
        """


        vmax=np.max(self.img_data)
        #vmin=vmax/5
        vmin=np.min(self.img_data)/10
        
        #plot the image
        #plt.imshow(self.img_data,interpolation="nearest",norm=LogNorm(vmax=vmax/5),origin="lower",cmap="gray")
        plt.imshow(self.img_data,interpolation="nearest",origin="lower",cmap="gray")
        plt.title("Image with detected sources") 
        #plot all clusters on top of the image
        for c in self.cluster_list:
            #plt.plot(c.center_y,c.center_x,'o',markersize=1,alpha=0.5)
            plt.plot(c.y_pixels,c.x_pixels,'o',markersize=0.5,alpha=0.5)
        plt.show()
        

    def copy_clusters(self,baseimg):
        """
        Copy clusters to other images
        """
        threashold=5
        for b in baseimg.cluster_list:
            c=cl.cluster()
            if self.img_data[b.center_x][b.center_y] >0:
                c.create_cluster(b.center_x,b.center_y,self.img_data[b.center_x][b.center_y],b.w)#modify b.w to self.w
                for i in np.arange(0,len(b.x_pixels)):
                    if self.img_data[b.x_pixels[i]][b.y_pixels[i]]>threashold*self.rms_data[b.x_pixels[i]][b.y_pixels[i]] and self.img_data[b.x_pixels[i]][b.y_pixels[i]]>0 :
                        c.add_pixel(b.x_pixels[i],b.y_pixels[i],self.img_data[b.x_pixels[i]][b.y_pixels[i]])
                        c.numpts=c.numpts+1

            self.cluster_list.append(c)
       
        #update center

        for c in self.cluster_list:
            if c.numpts >0:
                maxind=np.where(c.pixel_fluxs==np.max(c.pixel_fluxs))[0][0]
                c.center_x=c.x_pixels[maxind]
                c.center_y=c.y_pixels[maxind]
                c.center_flux=c.pixel_fluxs[maxind]
                #c.center_flux=np.mean(c.pixel_fluxs)
