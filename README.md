# find_objects

The package aims to find the list of the sources provided in the input image file. The package uses the clustering algorithm to find the sources in the given image. 

Developer documentation: https://github.com/ajayvibhute/pyfind_objects/wiki/Developer-Documentation

# Installation

1. Clone the git repository
```
git clone https://github.com/ajayvibhute/pyfind_objects.git
```
3. change directory to pyfind_objects
  ```
   cd pyfind_objects
  ```
4. create a .whl file
  ```
  python setup.py bdist_wheel
  ```
5. install the package from .whl file
  ```
    pip install dist/find_objects-1.0.0-py3-none-any.whl
  ```

# Usage/Example
Import the package
```
 from find_objects import find_objects as fo
```
Create an instance of ImageFile class
```
imgf=fo.ImageFile()
```
set input file, 
```
infile="input.txt"
imgf.input(infile)
```
In above example, "input.txt" is the name of the input file. Input file contains list of image and rms file paths seprated by a delimeter. The default delimeter is ','

Now, inorder to find the objects in the image with lowest frequency, run 
```
imgf.process()
```
The process function uses clustering based algorithm to find objects in the image file. The detailed algorithm along with the code description is given in the developer document.

After finding all the object in the image file, the spectral index can be computed using the function compute_spectral_index function
```
imgf.compute_spectral_index()
```
The source information can be saved using 
```
imgf.save_source_catlog()

```
By default, the catalog is stored in the "source_catalog.csv" file. However, it can be changed by passing outfile argument to the function

```
imgf.save_source_catlog(outfile="my_catalog.csv")

```
The method also takes an additional argument "overwrite", default True.

The source image and detected sources can be plotted using 

```
imgf.plot_image()

```
By default, the image for the lowest frequency will be plotted. However, one can provide the image index number for which the image needs to be plotted

```
imgf.plot_image(img_index=1)
```
The package also provides functionality to plot spectra (flux as a function of intensity) of all detected sources using 

```
imgf.plot_spectra()

```

# Software Requirements:
The package is developed and tested for the python3.11 on linux and mac operating system. The package is expected (though not tested) to work with all versions of the python3 and linux/mac operating system. In case, you face any issues, please report the error message to ajayvibhute@gmail.com


