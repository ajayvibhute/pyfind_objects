#!/bin/env python


from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import math
import time
import os

#from . import image as image
import image as image

class Sources:
    def __init__(self):
        self.center_x=[]
        self.center_y=[]
        self.center_flux=[]
        self.freq=[]
        self.norm=None
        self.spec_ind=None
        self.ra=[]
        self.dec=[]

class ImageFile:
    """
    Class responsible for accepting inputs, initiating the processing and outputting the results
    """
    def __init__(self):
        self.test=0
        self.imgfiles=[] #path for image files
        self.rmsfiles=[] #path for rms files
        self.images =[] #store images
        self.sourcelist=[]
    """
    #accept input file while creating the source object
    def __init__(self,infile,delim=","):
        self.input(infile)

    #accept input lists while creating the source object 
    def __init__(self, imgfiles,rmsfiles):
        self.input(imgfiles,rmsfiles)

    """
    def input(self,infile,delim=","):
        """
        Reads input file containing the list of image and rms file paths.

        Parameters: 
        -----------

        infile  : str
                Name of the file containting the list of image and rms file paths

        delim   : str, optional, default is ,
                The delimiter is used to separate the values specified in the input file
        """

        #unpacking the img and rms file paths
        self.imgfiles, self.rmsfiles=np.loadtxt(infile,unpack=True,dtype=str,delimiter=delim)


    #def input(self,imgfiles,rmsfiles):
    #    """
    #    Initializes list of image and rms file paths 
    #
    #    Parameters:
    #    -----------
    #    imgfiles    : list
    #                List of input image file paths
    #
    #    rmsfiles    : list
    #                List of rms file paths
    #    """
    #    
    #    #Ensuring both the lists are of the same length
    #    """
    #    if len(imgfiles) != len(rmsfiles):
    #        print("Input lists should be of the same length");
    #        return -1
    #    """
    #    #assigning inputs 
    #    self.imgfiles=imgfiles
    #    self.rmsfiles=rmsfiles
    
    
    def init_images(self):
        """
        Creates objects for image class, assigns image, rms file and reads 
        image and rms file. 
        """
        #iterating over all the input files and reading images
        #rms and frequency
        for ind in np.arange(0,len(self.imgfiles)):
            i=image.Image(self.imgfiles[ind],self.rmsfiles[ind])
            i.read_input() 
            self.images.append(i)

        #sort the images using frequency as a key
        #self.images.sort(key=lambda x: x.img_freq)
        
    def cleanup_clusters(self):  
        #ideally check the source intensity above rms
        #implement if time permits
        pts_count=[0]*self.images[0].numclusters

        #add a flag to decide local threashold based on distance type
        local_threash=5
        for j in np.arange(0,self.images[0].numclusters):
            for i in self.images:
                if(i.cluster_list[j].numpts<5):
                    pts_count[j]=pts_count[j]+1
        removed_count=0
        for j in np.arange(0,len(pts_count)):
            if pts_count[j]>8:
                #delete cluster
                for i in self.images:
                    i.cluster_list.pop(j-removed_count)
                    i.numclusters-=1
                removed_count+=1
            

    def plot_spectra(self):
        """
        Plots spectrum (intensity as a function of frequency) of all the detected sources
        """
        for i in np.arange(0,len(self.sourcelist)):
            #print("source Id:",i,self.sourcelist[i].center_x,self.sourcelist[i].center_y,self.sourcelist[i].center_flux)
            plt.plot(self.sourcelist[i].freq, self.sourcelist[i].center_flux)
            plt.show()

    def save_source_catalog(self,delim=",",outfile="source_catalog.csv",overwrite=True):
        """
        Saves source catalog

        Parameters:
        ----------

        delim   :   str, optional, default is ","
                Delimiter to seperate the fields

        outfile : str, optional, default is source_catalog.csv
                Name of the output file

        overwrite   : bool, optional default is True
                    If true, overwrites the file otherwise appends to existing file
        """

        if overwrite:
            mode="w"
            #delete file if flag is set to overwrite
            if os.path.exists(outfile):
                os.remove(outfile)

        else:
            #open file in append mode
            mode="a"

        #create a file pointer    
        f=open(outfile,mode)
        f.write("SourceId,RA,DEC,Flux,SpectralIndex\n")
        #reformat the output
        #iterate over all the source
        for i in np.arange(0,len(self.sourcelist)):
            #write source information
            #f.write(str(i)+delim+str(self.sourcelist[i].center_x)+delim+str(self.sourcelist[i].center_y)+delim+str(self.sourcelist[i].center_flux)+"\n")
            #ideally should be the pixel repeated maximum number of times
            f.write(str(i)+delim+str(np.mean(self.sourcelist[i].ra))+delim+str(np.mean(self.sourcelist[i].dec))+delim+str(np.mean(self.sourcelist[i].center_flux))+delim+str(self.sourcelist[i].spec_ind)+"\n")
        #close the file descriptor
        f.close()
    
    
    def linear_func(self,x, a, b):
        """
        Function used to compute spectral index

        """
        return a+x*b


    def compute_spectral_index(self):
        """
        computes the spectral index for the image
        """
        
        for s in self.sourcelist:
            #rename variables, these are not positions
            ind=np.where(np.array(s.center_flux)<=0)[0].astype(int)

            freq=s.freq
            cf=s.center_flux
            for i in np.arange(0,len(ind)):
                cf.pop(ind[i]-i)
                freq.pop(ind[i]-i)
            log_xpos=np.log10(freq)
            log_ypos=np.log10(cf)
        
            if len(log_xpos)==len(log_ypos):
                norm, spec_ind = curve_fit(self.linear_func, log_xpos, log_ypos)[0]    
                s.norm=norm
                s.spec_ind=spec_ind
            else:
                print("flux and frequency length do no match")

    def process(self):
        """
        Process the images
        """
        print("Init Images")
        self.init_images()
        
        start_time = time.time()
        print("Flag Pixels")
        #flag the pixels
        self.images[0].flag_pixels()
        print("flag pix--- %s seconds ---" % (time.time()-start_time ))
        
        start_time = time.time()
        print("Perform clustering")
        #perform clsutering
        self.images[0].find_clusters()
        print("clustering--- %s seconds ---" % (time.time()-start_time ))

        start_time = time.time()
        print("Merge clusters")
        #merge clusters
        self.images[0].merge_clusters()
        print("Plot--- %s seconds ---" % (time.time()-start_time ))

        start_time = time.time()
        #copy clusters
        for i in np.arange(1,len(self.images)):
            self.images[i].copy_clusters(self.images[0])
        
        #this is required to process only this data-set
        self.cleanup_clusters()

        #create instance of source class and copy all the sources in it
        for i in np.arange(0,len(self.images[0].cluster_list)):
            s=Sources()
            self.sourcelist.append(s)
        #copying detected source information 
        #in the source clsuter
        for i in np.arange(0,len(self.images)):
            for j in np.arange(0,len(self.images[i].cluster_list)):

                if self.images[i].cluster_list[j].numpts != 0:

                    self.sourcelist[j].center_x.append(self.images[i].cluster_list[j].center_x)    
                    self.sourcelist[j].center_y.append(self.images[i].cluster_list[j].center_y)    
                    self.sourcelist[j].center_flux.append(self.images[i].cluster_list[j].center_flux)    
                    self.sourcelist[j].freq.append(self.images[i].img_freq)    
                    self.sourcelist[j].ra.append(self.images[i].cluster_list[j].ra)    
                    self.sourcelist[j].dec.append(self.images[i].cluster_list[j].dec)    


        print("Copy clusters--- %s seconds ---" % (time.time()-start_time ))

    def plot_image(self,img_index=0):
        """
        Calls plot image method in class Image to plot image and clusters

        Parameters:
        -----------

        img_index   : int, optional default is 0
                    The index of the image needs to be plotted
        """
        
        if len(self.images)<img_index:
            print("Error: Image index is greater than the number of images")
            return -1
        
        #calling plot function in image class
        self.images[img_index].plot_image()

        #plot images
        #for i in self.images:
         #   i.plot_image()
        #self.images[0].plot_image()

    
if __name__ == "__main__":

    
    infile="/Users/avibhute/NRAO/pyfind_objects/find_objects/input.txt"

    imgf=ImageFile()

    start_time = time.time()
    #set input file
    imgf.input(infile)
    print("Reading--- %s seconds ---" % (time.time()-start_time ))

    #process the image and find objects in it

    start_time = time.time()
    imgf.process()
    #imgf.plot_spectra()

    imgf.compute_spectral_index()
    imgf.save_source_catalog()
    imgf.plot_image()
    print("Process--- %s seconds ---" % (time.time()-start_time ))
    
    """
    VLASS1.2.se.T01t08.J033228-363000.06.2048.v1.spw2.I.iter3.image.pbcor.tt0.subim
    infile="/Users/avibhute/NRAO/images/VLASS1.2.se.T01t08.J033228-363000.06.2048.v1.spw2.I.iter3.image.pbcor.tt0.subim.fits"
    rmsfile="/Users/avibhute/NRAO/images/VLASS1.2.se.T01t08.J033228-363000.06.2048.v1.spw2.I.iter3.image.pbcor.tt0.rms.subim.fits"
    start_time = time.time()
    #create object of the source
    im=Image(infile,rmsfile);

    #read input
    im.read_input()

    #flag the pixels
    im.flag_pixels()
    
    #perform clsutering
    im.find_clusters()
    
    #merge clusters, 
    im.merge_clusters()
    execution_time=time.time()-start_time 
    
    #plot image and clusters
    im.plot_image()

    print("--- %s seconds ---" % (execution_time))
    """
