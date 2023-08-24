#Contains all the utility classes, methods required to perform object search

from astropy.wcs import WCS

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

        #create object of WCS with image_header
        w=WCS(img_header)

        #The VLASS data is four dimensional and two dimensions of the WCS needs to be
        #removed to perform pixel to sky conversion.
        #reference: https://github.com/aplpy/aplpy/issues/423
        
        if dropaxis:
        
            #remove dimension from WCS
            w=w.dropaxis(2)
            w=w.dropaxis(2)    

        #returning object of the wcs
        return w

