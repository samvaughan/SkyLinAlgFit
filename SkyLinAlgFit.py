import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits

import pdb
import sys
import continuum 

import argparse




def skyfit_rows(data, a, b, r, clip_edges=True):
    d1,d2,d3= data.shape

    median=np.median(data, axis=0)


    y,x = np.ogrid[-b:d2-b, -a:d3-a]
    mask = x*x + y*y <= r*r

    if clip_edges==True:

        median[:, :4]=0
        median[:, -4:]=0


    #median[~mask]=0



    #plt.imshow(median)
    #plt.show()

    median=median.ravel()

    ones=np.ones(d3)
    continuum=np.zeros([d2, d2*d3])
    for k in range(d2):

        continuum[k, k*d3:(k+1)*d3]=ones

    model=np.column_stack((median,continuum.T))


    c_vals=np.zeros([1+continuum.shape[0], d1])
    sigma_vals=np.zeros_like(c_vals)



    for i in range(d1):
        sys.stdout.write('\r')
        
        sys.stdout.write("Fitting Slice {}".format(i))
        sys.stdout.flush()


        obj=data[i,:,:].ravel()

        c, resid, rank, sigma=np.linalg.lstsq(model, obj)

        #rint c

        """
        if i >100:

            residuals=obj-(c[0]*model.T[:,0]+c[1]*model.T[:,1])
            residuals=residuals.reshape(d2, d3)
            plt.imshow(residuals)
            plt.colorbar()
            plt.show()
        """


        c_vals[:, i]=c
        
    print "\n"
    return c_vals


def skyfit_fullimage(data, clip_edges=True):


    """
    Fit a galaxy and sky model independently to a sky-subtracted cube. The aim is to recover the sky continuum, fit a polynomial to it and subtract it off from the O-S cube.

    The sky model is an array of 1s, in the same shape as each slice of the data cube. The galaxy model is the mean of the object cube across the wavelength direction

    """


    d1,d2,d3= data.shape

    median=np.median(data, axis=0)


    

    if clip_edges==True:

        median[:, :4]=0
        median[:, -4:]=0


    #median[~mask]=0



    #plt.imshow(median)
    #plt.show()

    median=median.ravel()

    ones=np.ones(d2*d3)
    continuum=np.zeros([1, d2*d3])
    continuum[0,:]=ones

    





    model=np.column_stack((median,continuum.T))


    c_vals=np.zeros([1+continuum.shape[0], d1])
    sigma_vals=np.zeros_like(c_vals)



    for i in range(d1):
        sys.stdout.write('\r')
        
        sys.stdout.write("Fitting Slice {}".format(i))
        sys.stdout.flush()


        obj=data[i,:,:].ravel()

        c, resid, rank, sigma=np.linalg.lstsq(model, obj)

        #rint c

        
        c_vals[:, i]=c
        
    print "\n"
    return c_vals


parser = argparse.ArgumentParser(description='Fit a sky and galaxy model to an O-S cube, then fit a polynomial to the sky residuals. Subtract this continuum from each O-S cube and save')
parser.add_argument('GalName', type=str, help='Should be IC843 or NGC1277')

args = parser.parse_args()

GalName=args.GalName


if GalName=="NGC1277":
    data_list=["ms047_048.fits", "ms049_048.fits", "ms050_051.fits", "ms052_051.fits", "ms056_057.fits", "ms058_057.fits", "ms059_060.fits"]
    cubepath="/Volumes/SPV_SWIFT/SWIFT/NGC1277"


elif GalName=="IC843":
    data_list=["ms065_066.fits","ms067_066.fits","ms071_072.fits","ms073_072.fits","ms074_075.fits","ms076_075.fits","ms080_081.fits","ms082_081.fits","ms083_084.fits"]
    cubepath="/Volumes/SPV_SWIFT/SWIFT/Atlas3D/2016-03-17"

elif GalName=="GMP4928":
    data_list=["ms059_060.fits", "ms061_060.fits", "ms064_065.fits",  "ms066_065.fits"]
    cubepath="./GMP4928"
else:
    raise NameError("Please Specify a Galaxy")



final_results_rows=np.zeros([len(data_list), 45, 4112])
final_results_fullimage=np.zeros([len(data_list), 2, 4112])


continuum_fits=np.zeros([len(data_list), 4112])
lamdas=np.arange(6300, 6300+4112)


fig, axs=plt.subplots(nrows=2, ncols=2)

from matplotlib.ticker import AutoMinorLocator

for i, (name) in enumerate(data_list):
    data=fits.open("{}/{}".format(cubepath, name))[0].data
    header=fits.open("{}/{}".format(cubepath, name))[0].header

    d1,d2,d3=data.shape


    print "\n\nFitting cube {}\n\n".format(i+1)



    #final_results_rows[i, :, :]=skyfit_rows(data, a, b, r)
    final_results_fullimage[i, :, :]=skyfit_fullimage(data)

    #Find the gradient of the sky array. Skylines correspond to sharp spikes in the gradient, so we can get rid of them by ignoring pixels with a large dx.
    dx=np.gradient(final_results_fullimage[i, 1, :])
    
    fitlocs=(np.abs(dx) < 0.5)

    #Use the continuum fit function to fit polynomials, avoiding the skylines. 
    continuum_fits[i, :]=continuum.fit_continuum(lamdas, final_results_fullimage[i,1,:], np.ones_like(final_results_fullimage[i,1,:]), [2,0.5, 0.3], 7, [6300, 10412], plot=False, fitloc=fitlocs)


    import pdb
    pdb.set_trace()

    if i<2:
        j=0
        k=i
    else:
        j=1
        k=i-2


    axs[j, k].plot(lamdas, final_results_fullimage[i, 1, :], c="k")
    axs[j, k].plot(lamdas, continuum_fits[i, :], c="r", linewidth=2.0)
    axs[j, k].set_xlim([6300,10412])
    axs[j, k].set_xlabel("Wavelength")
    axs[j, k].set_ylim([-20,35])
    axs[j, k].set_ylabel("Counts")
    axs[j, k].set_title("{}".format(name))

    
    minorLocator = AutoMinorLocator()
    axs[j, k].yaxis.set_minor_locator(minorLocator)

    sky_cont_cube=np.repeat(continuum_fits[i, :], d2*d3).reshape((d1,d2,d3))
    #hdu=fits.PrimaryHDU(data-sky_cont_cube)


    cubename="{0}/{0}_sky_cont_sub_{1}".format(GalName, name)
    fits.writeto(cubename, data-sky_cont_cube, header=header, clobber=True)
    print "Saving Cube {} as {}".format(i+1, cubename)



fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
plt.savefig("{0}/{0}_sky_continuum.pdf".format(GalName), dpi=200)
plt.show()

