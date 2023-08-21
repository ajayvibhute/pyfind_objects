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
        w=WCS(img_header)
        if dropaxis:
            w=w.dropaxis(2)
            w=w.dropaxis(2)    
        return w

