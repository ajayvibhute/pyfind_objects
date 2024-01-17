#!/bin/env python


from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np
import math
import time
import os

from . import image as image
# import image as image

import bdsf


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


class Host:
    """
        class is defined to store the information about the host galaxy and radio
        structures around it.
    """
    def __init__(self, center_ra, center_dec, center_x=None, center_y=None, host_id=None):
        self.center_ra = center_ra
        self.center_dec = center_dec
        self.center_x = center_x
        self.center_y = center_y
        self.host_id = host_id
        self.ellipse_list = []


class Ellipse:
    """
        Class to store information about the information about
        the ellipses fitted and fitted params.
    """
    is_assigned = False
    host_id = None

    def __init__(self, param=None, island_id=None):
        self. x0, self.y0, self.ap, self.bp, self.e, self.phi = param
        self.island_id = island_id


class Find_sources:
    def __init__(self, infile):
        self.infile = infile
        self.bdsf_img = None
        self.param = None
        self.param_temp = None
        self.hosts = []
        self.ellipse_list = []

    def search(self):
        """
        checks distance of each hosts from all the islands. if the
        distance is less than the threashold, the island will be
        assigned to the host.
        """
        for h in self.hosts:
            for e in self.ellipse_list:
                if not e.is_assigned:
                    if math.dist([h.center_x, h.center_y], [e.x0, e.y0]) < 70:
                        h.ellipse_list.append(e)
                        e.is_assigned = True

    def process(self):
        """
        Function to process image in bdsf and return the processed
        image with islands and other information
        """
        self.bdsf_img = bdsf.process_image(self.infile)

    def _check_image_boundry(self, x, y, x0, y0, x1, y1):
        """
        Checks if a given galaxy co-ordinates are within the boundry 
        of the image which is currently being processed.
        """
        # need to verify the condition with several fits images.
        # this is becuase how the co-ordinates changes and origin
        if x1 < x and x < x0 and y0 < y and y < y1:
            return True
        else:
            return False

    def init_hosts(self, host_catlog_file_name="optical_host_gal_catlog.txt"):
        """
        Function computes boundry of the fits image using 
        NAXIS1, NAXIS2 and WCS. Creates a list of host galaxies within the
        boundry of the fits image.
        The input file/catalog should contain RA, DEC of the host galaxies in degree.
        For the boundry conditions, co-ordinates look like

        (0, NAXIS2)---------(NAXIS1, NAXIS2)
        |                   |
        |                   |
        |                   |
        |                   |
        (0, 0)--------------(NAXIS1, 0) 
      """

        w = WCS(self.bdsf_img.header)
        tmp = w.pixel_to_world(0, 0, 0, 0)
        # last two arguments are dummy and not used anywhere in the process
        if tmp != 0:
            x0 = tmp[0].ra.value
            y0 = tmp[0].dec.value
        else:
            print("Error computing boundries\n Check the WCS of input image")

        tmp = w.pixel_to_world(self.bdsf_img.header["NAXIS1"], self.bdsf_img.header["NAXIS2"], 0, 0)
        # last two arguments are dummy and not used anywhere in the process
        if tmp != 0:
            x1 = tmp[0].ra.value
            y1 = tmp[0].dec.value
        else:
            print("Error computing boundries\n Check the WCS of input image")
        cat_data = np.loadtxt(host_catlog_file_name)
        gal_ra = cat_data[:, 0]
        gal_dec = cat_data[:, 1]
        for i in np.arange(0, len(gal_ra)):
            # this can be done in a smarter way
            if self._check_image_boundry(gal_ra[i], gal_dec[i], x0, y0, x1, y1):
                gal_center = w.wcs_world2pix(gal_ra[i], gal_dec[i], 0, 0, 0)
                self.hosts.append(Host(gal_ra[i], gal_dec[i], gal_center[0], gal_center[1], len(self.hosts)))
                print(gal_ra[i], gal_dec[i])

    

    def fit(self):

        for ins in self.bdsf_img.islands:
            # fitting ellipses to island
            self.param_tmp = self._fit_ellipse(ins.border[0],ins.border[1])
            if len(self.param_tmp) != 0:
                self.param = self._cart_to_pol(self.param_tmp)

                if self.param != -111:
                    el = Ellipse(self.param, ins.island_id)
                    self.ellipse_list.append(el)
        """
        if self.bdsf_img.resid_gaus_arr is None:
            low = 1.1*abs(self.bdsf_img.min_value)
        else:
            low = np.max([1.1*abs(self.bdsf_img.min_value),1.1*abs(np.nanmin(self.bdsf_img.resid_gaus_arr))])

        im=np.log10(self.bdsf_img.ch0_arr + low)
        plt.imshow(im,  cmap="gray", interpolation='nearest')

        for h in self.hosts:
            print(h.center_ra, h.center_dec, h.center_x, h.center_y)
            plt.scatter(h.center_y,h.center_x, s=6)
        for ins in self.bdsf_img.islands:
            #fitting ellipses to island
            self.param_tmp=self.fit_ellipse(ins.border[0],ins.border[1])
            if len(self.param_tmp) != 0:
                self.param = self.cart_to_pol(self.param_tmp)

                if self.param!=-111 :
                    x0, y0, ap, bp, e, phi = self.param
                    ra = 0.0
                    if  ap>=bp:
                        ra = bp/ap
                    else:
                        ra = ap/bp

                    if ra < 1:
                        xt,yt=self.get_ellipse_pts(self.param)
                        #plt.scatter(self.param[1],self.param[0],s=1)
                        #plt.scatter(yt,xt,s=1)
        """
        # plt.show()

    def plot(self):
        if self.bdsf_img.resid_gaus_arr is None:
            low = 1.1*abs(self.bdsf_img.min_value)
        else:
            low = np.max([1.1*abs(self.bdsf_img.min_value), 1.1*abs(np.nanmin(self.bdsf_img.resid_gaus_arr))])

        im = np.log10(self.bdsf_img.ch0_arr + low)
        plt.imshow(im,  cmap="gray", interpolation='nearest')
        for e in self.hosts[1].ellipse_list:
            xt, yt = self.get_ellipse_pts([e.x0, e.y0, e.ap, e.bp, e.e, e.phi])
            plt.scatter(yt, xt)
        
        plt.scatter(self.hosts[1].center_y, self.hosts[1].center_x, s=6)
        plt.show()

    # code credit: https://scipython.com/blog/direct-linear-least-squares-fitting-of-an-ellipse/
    def _fit_ellipse(self, x, y):
        """
        Fit the coefficients a,b,c,d,e,f, representing an ellipse described by
        the formula F(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0 to the provided
        arrays of data points x=[x1, x2, ..., xn] and y=[y1, y2, ..., yn].

        Based on the algorithm of Halir and Flusser, "Numerically stable direct
        least squares fitting of ellipses'.
        """

        D1 = np.vstack([x**2, x*y, y**2]).T
        D2 = np.vstack([x, y, np.ones(len(x))]).T
        S1 = D1.T @ D1
        S2 = D1.T @ D2
        S3 = D2.T @ D2
        T = -np.linalg.inv(S3) @ S2.T
        M = S1 + S2 @ T
        C = np.array(((0, 0, 2), (0, -1, 0), (2, 0, 0)), dtype=float)
        M = np.linalg.inv(C) @ M
        eigval, eigvec = np.linalg.eig(M)
        con = 4 * eigvec[0]* eigvec[2] - eigvec[1]**2
        ak = eigvec[:, np.nonzero(con > 0)[0]]
        return np.concatenate((ak, T @ ak)).ravel()

    def _cart_to_pol(self, coeffs):
        """

        Convert the cartesian conic coefficients, (a, b, c, d, e, f), to the
        ellipse parameters, where F(x, y) = ax^2 + bxy + cy^2 + dx + ey + f = 0.
        The returned parameters are x0, y0, ap, bp, e, phi, where (x0, y0) is the
        ellipse centre; (ap, bp) are the semi-major and semi-minor axes,
        respectively; e is the eccentricity; and phi is the rotation of the semi-
        major axis from the x-axis.

        """

        # We use the formulas from https://mathworld.wolfram.com/Ellipse.html
        # which assumes a cartesian form ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0.
        # Therefore, rename and scale b, d and f appropriately.
        a = coeffs[0]
        b = coeffs[1] / 2
        c = coeffs[2]
        d = coeffs[3] / 2
        f = coeffs[4] / 2
        g = coeffs[5]

        den = b**2 - a*c
        if den > 0:
            # raise ValueError('coeffs do not represent an ellipse: b^2 - 4ac must'
            #                ' be negative!')

            return -111
        # The location of the ellipse centre.
        x0, y0 = (c*d - b*f) / den, (a*f - b*d) / den

        num = 2 * (a*f**2 + c*d**2 + g*b**2 - 2*b*d*f - a*c*g)
        fac = np.sqrt((a - c)**2 + 4*b**2)
        # The semi-major and semi-minor axis lengths (these are not sorted).
        ap = np.sqrt(num / den / (fac - a - c))
        bp = np.sqrt(num / den / (-fac - a - c))

        # Sort the semi-major and semi-minor axis lengths but keep track of
        # the original relative magnitudes of width and height.
        width_gt_height = True
        if ap < bp:
            width_gt_height = False
            ap, bp = bp, ap

        # The eccentricity.
        r = (bp/ap)**2
        if r > 1:
            r = 1/r
        e = np.sqrt(1 - r)

        # The angle of anticlockwise rotation of the major-axis from x-axis.
        if b == 0:
            phi = 0 if a < c else np.pi/2
        else:
            phi = np.arctan((2.*b) / (a - c)) / 2
            if a > c:
                phi += np.pi/2
        if not width_gt_height:
            # Ensure that phi is the angle to rotate to the semi-major axis.
            phi += np.pi/2

        if isinstance(phi, complex):
            phi = phi.real % np.pi
        else:
            phi = phi % np.pi

        return x0, y0, ap, bp, e, phi

    def get_ellipse_pts(self, params, npts=10, tmin=0, tmax=2*np.pi):
        """
        Return npts points on the ellipse described by the params = x0, y0, ap,
        bp, e, phi for values of the parametric variable t between tmin and tmax.

        """
        x0, y0, ap, bp, e, phi = params
        # A grid of the parametric variable, t.
        t = np.linspace(tmin, tmax, npts)
        x = x0 + ap * np.cos(t) * np.cos(phi) - bp * np.sin(t) * np.sin(phi)
        y = y0 + ap * np.cos(t) * np.sin(phi) + bp * np.sin(t) * np.cos(phi)
        return x, y






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
        self.dist_type="euclidean"
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

        # unpacking the img and rms file paths
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
            i=image.Image(self.imgfiles[ind],self.rmsfiles[ind],self.dist_type)
            i.read_input() 
            self.images.append(i)

        #sort the images using frequency as a key
        self.images.sort(key=lambda x: x.img_freq)
        
    def cleanup_clusters(self):  
        #ideally check the source intensity above rms
        #implement if time permits
        pts_count=[0]*self.images[0].numclusters

        #add a flag to decide local threashold based on distance type
        local_threash=5
        for j in np.arange(0,self.images[0].numclusters):
            for i in self.images:
                if(i.cluster_list[j].numpts<=8):
                    pts_count[j]=pts_count[j]+1
        removed_count=0
        for j in np.arange(0,len(pts_count)):
            if pts_count[j]>=5:
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
            plt.title("Spectra for ({:.2f},{:.2f})".format(np.mean(self.sourcelist[i].ra),np.mean(self.sourcelist[i].dec)))
            plt.xlabel("Frequency")
            plt.ylabel("Flux in Jy")
            plt.plot(self.sourcelist[i].freq, self.sourcelist[i].center_flux,"ro",markersize=2)
            plt.show()

    def save_spectra(self):
        """
        saves source spectrum in the image file
        
        filename : default to spectra_<source_id>.png

        """
        #to do: Add a provision to provide file name as an input
        for i in np.arange(0,len(self.sourcelist)):
            plt.title("Spectra for ({:.2f},{:.2f})".format(np.mean(self.sourcelist[i].ra),np.mean(self.sourcelist[i].dec)))
            plt.xlabel("Frequency")
            plt.ylabel("Flux in Jy")
            plt.plot(self.sourcelist[i].freq, self.sourcelist[i].center_flux,"ro",markersize=2)
            plt.savefig("spectra_"+str(i)+".png")
            plt.clf()

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
            
            #To do: ideally should be the pixel repeated maximum number of times
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
        self.init_images()
        
        start_time = time.time()
        #flag the pixels
        self.images[0].flag_pixels()
        #print("flag pix--- %s seconds ---" % (time.time()-start_time ))
        
        start_time = time.time()
        #perform clsutering
        self.images[0].find_clusters()
        #print("clustering--- %s seconds ---" % (time.time()-start_time ))

        start_time = time.time()
        #merge clusters
        self.images[0].merge_clusters()
        #print("Plot--- %s seconds ---" % (time.time()-start_time ))

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

    #process the image and find objects in it

    start_time = time.time()
    imgf.process()


    imgf.compute_spectral_index()
    imgf.save_source_catalog()
    imgf.plot_image()
    
    #imgf.plot_spectra()
    imgf.save_spectra() 
    #print("Process--- %s seconds ---" % (time.time()-start_time ))
    
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

     
