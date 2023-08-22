#Contains defination for the cluster and other methods required to
#process the cluster
from astropy.coordinates import SkyCoord
import math
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
        self.agg_flux=0

    def init_center_wcs(self):
        """
        Initializes WCS required to compute distance
        """
        #generate instance of the Skycoord for the center of 
        #cluster. This is used to compute distance from other sky pixels
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
        #add center to the list of pixels and flux to list of flux
        self.x_pixels.append(x)
        self.y_pixels.append(y)
        self.pixel_fluxs.append(flux)
        self.numpts=1

        #Update the center position, flux
        self.center_x=x
        self.center_y=y
        self.center_flux=flux
        #assign the wcs
        self.w=w
        #initialize the wcs
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
        #check if cluster is already created or not
        #pixel can be added only if cluster exists
        if self.numpts==0:
            print("Error: cluster is not created yet\n")
            return -1
        
        #apped position, flux and update the number
        #of points in the cluster
        self.x_pixels.append(x)
        self.y_pixels.append(y)
        self.pixel_fluxs.append(flux)
        self.numpts=self.numpts+1

        #Assumption is, flux is maximum at the center. (May not be true in case of radio image
        #need to check and update if required)
        #Check flux of the newly added pixel. Update the centre position
        #if flux of the new pixel is more than the flux at the existing center
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
        #The distance can either be computed from the center or from all points 
        #in the cluster. The later is mainly required if one is looking for the extended 
        #emission of the source and will come at the cost of computing time. 
        #for the current data set (VLASS) computing distance from the center lists all the 
        #sources in the FOV
        """
        if all_pix:
            #compute distance from all points in the 
            #cluster and return minimum distance
            dmin=100;
            c2 = SkyCoord(self.w.pixel_to_world(x,y))
            for i in np.arange(0,len(self.x_pixels)):
                c1=SkyCoord(self.w.pixel_to_world(self.x_pixels[i],self.y_pixels[i]))
                d=self.c1.separation(c2)
            
                if d.arcmin<dmin:
                    dmin=d.arcmin
            return dmin
        else:
            #compute distance from centre
            c2 = SkyCoord(self.w.pixel_to_world(x,y))
            d=self.c1.separation(c2)
            return d.arcmin
        """
        return math.dist([self.center_x,self.center_y],[x,y])

        
