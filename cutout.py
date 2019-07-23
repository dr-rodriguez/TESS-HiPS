# Test script for cutouts
import os
import numpy as np
from itertools import product
from astropy import wcs
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astrocut import CutoutFactory

# try:
#     from astropy.wcs.utils import fit_wcs_from_points
# except ImportError:  # astropy version does not have the function
from astrocut.utils.wcs_fitting import fit_wcs_from_points


def my_cutout(img, xmin, xmax, ymin, ymax):
    img_cutout = img[1].data[xmin:xmax, ymin:ymax]
    # img_cutout = img[1].data[ymin:ymax, xmin:xmax]
    uncert_cutout = img[2].data[xmin:xmax, ymin:ymax]

    # Making the aperture array
    aperture = np.ones((xmax - xmin, ymax - ymin), dtype=np.int32)

    return img_cutout, uncert_cutout, aperture


def my_cutout_wcs(cutout_wcs, cutout_shape, cut):

    # Getting matched pixel, world coordinate pairs
    pix_inds = np.array(list(product(list(range(cutout_shape[1])), list(range(cutout_shape[0])))))
    world_pix = SkyCoord(cutout_wcs.all_pix2world(pix_inds, 1), unit='deg')

    # Getting the fit WCS
    linear_wcs = fit_wcs_from_points(pix_inds[:, 0], pix_inds[:, 1], world_pix, mode='wcs',
                                     # proj_point='center')
                                     proj_point=[cut.center_coord.data.lon.value, cut.center_coord.data.lat.value])

    cut.cutout_wcs = linear_wcs

    # Checking the fit
    world_pix_new = SkyCoord(linear_wcs.all_pix2world(pix_inds, 1), unit='deg')

    dists = world_pix.separation(world_pix_new).to('deg')
    sigma = np.sqrt(sum(dists.value ** 2))

    cut.cutout_wcs_fit['WCS_MSEP'][0] = dists.max().value
    cut.cutout_wcs_fit['WCS_SIG'][0] = sigma

    return (dists.max(), sigma)


def get_limits(shape=(2048, 2048), width=100, start=(44,0)):
    limit_list = []

    xrange = np.arange(start[1], shape[1], width)
    yrange = np.arange(start[0], shape[0], width)

    for x in xrange:
        for y in yrange:
            y2 = y+width
            if y2 > shape[1]:
                y2 = shape[1]
            x2 = x + width
            if x2 > shape[0]:
                x2 = shape[0]

            limit_list.append(np.array([[y, y2], [x, x2]]))

    return limit_list


def make_my_coutout(FILENAME, OUTNAME, limits, output_path='Cutouts', verbose=True, overwrite=False):
    if not overwrite:
        if os.path.isfile(os.path.join(output_path, OUTNAME)):
            # print('{} exists. Skipping'.format(OUTNAME))
            return

    if verbose: print('Working on {}'.format(FILENAME))

    cube = fits.open(FILENAME)

    cut = CutoutFactory()

    # Get the info we need from the data table (orig function fails)
    # cut._parse_table_info(cube[1].header, cube[1].data, verbose)
    # wcs_header = fits.header.Header()
    cut.cube_wcs = wcs.WCS(cube[1].header)

    cut.cutout_lims = limits
    ymin, ymax = limits[0]
    xmin, xmax = limits[1]
    if verbose:
        print("xmin,xmax: {} {}".format(xmin, xmax))
        print("ymin,ymax: {} {}".format(ymin, ymax))

    # Get center coordinates for future WCS
    w = wcs.WCS(cube[1].header)
    lon, lat = w.all_pix2world(np.mean(limits[1]), np.mean(limits[0]), 0)
    cut.center_coord = SkyCoord('{} {}'.format(lon, lat), unit='deg')
    if verbose:
        print('RA Dec: {} {}'.format(lon, lat))

    # Make the cutout (orig function fails)
    # img_cutout, uncert_cutout, aperture = cut._get_cutout(cube[1].data, verbose=verbose)
    img_cutout, uncert_cutout, aperture = my_cutout(cube, xmin, xmax, ymin, ymax)

    # Get cutout wcs info
    cutout_wcs_full = cut._get_full_cutout_wcs(cube[1].header)

    # max_dist, sigma = cut._fit_cutout_wcs(cutout_wcs_full, img_cutout.shape[1:])  # orig function fails
    try:
        max_dist, sigma = my_cutout_wcs(cutout_wcs_full, img_cutout.shape, cut)
    except Exception as e:
        print('Error when making '.format(OUTNAME))
        print('Error: {}'.format(e))
        return

    if verbose:
        print("Maximum distance between approximate and true location: {}".format(max_dist.to(u.arcsec)))
        print("Error in approximate WCS (sigma): {}".format(sigma))

    if max_dist.to(u.arcsec) > 10 * u.arcsec:
        print('ERROR: Too max distance large for {} ({})'.format(OUTNAME, max_dist.to(u.arcsec)))
        return

    # cutout_wcs_dict = cut._get_cutout_wcs_dict()

    # Build the file
    primary_hdu = cube[0]
    cut._update_primary_header(primary_hdu.header)
    cutout_hdu = fits.ImageHDU(data=img_cutout, header=cut.cutout_wcs.to_header())
    cutout_hdu_list = fits.HDUList([primary_hdu, cutout_hdu])
    cut._apply_header_inherit(cutout_hdu_list)
    # cutout_hdu_list.info()

    # Make sure the output directory exists
    if output_path and not os.path.exists(output_path):
        os.makedirs(output_path)

    # Write the file
    if verbose: print('Writting {}'.format(OUTNAME))
    cutout_hdu_list.writeto(os.path.join(output_path, OUTNAME), overwrite=True, checksum=True)

    cube.close()


# Single run
FILENAME = 'Data/tess2018206192942-s0001-1-1-0120-s_ffic.fits'
OUTNAME = 'tess_s0001-1-1_0-100_0-100.fits'
ymin, ymax = 0, 100
xmin, xmax = 0, 100
limits = np.array([[ymin, ymax], [xmin, xmax]])
make_my_coutout(FILENAME, OUTNAME, limits, output_path='Cutouts', verbose=True, overwrite=False)

# Another single run
FILENAME = 'Data/tess2018206192942-s0001-1-1-0120-s_ffic.fits'
ymin, ymax = 0, 100
xmin, xmax = 44, 144
limits = np.array([[ymin, ymax], [xmin, xmax]])
OUTNAME = 'tess_s0001-1-1_{}-{}_{}-{}.fits'.format(limits[0][0], limits[0][1], limits[1][0], limits[1][1])
make_my_coutout(FILENAME, OUTNAME, limits, output_path='Cutouts', verbose=True, overwrite=True)

# Large max dist example
"""
Working on Data/tess2018206192942-s0001-4-4-0120-s_ffic.fits
xmin,xmax: 100 200 
ymin,ymax: 0 100
RA Dec: 104.31211100371591 -56.76839717706978
Maximum distance between approximate and true location: 1356.142458902589 arcsec
Error in approximate WCS (sigma): 15.987921128646978
ERROR: Too max distance large for tess_s0001-4-4_0000-0100_0100-0200.fits (1356.142458902589 arcsec)
"""
FILENAME = 'Data/tess2018206192942-s0001-4-4-0120-s_ffic.fits'
ymin,ymax = 0, 100
xmin,xmax = 100, 200
limits = np.array([[ymin, ymax], [xmin, xmax]])
OUTNAME = 'TEST_s0001-4-4_{}-{}_{}-{}.fits'.format(limits[0][0], limits[0][1], limits[1][0], limits[1][1])
make_my_coutout(FILENAME, OUTNAME, limits, output_path='', verbose=True, overwrite=True)

# Loop over grid
# full_limits = get_limits(shape=(2078, 2136), width=100)
full_limits = get_limits(shape=(2000, 2100), width=100, start=50)
for i, limits in enumerate(full_limits):
    OUTNAME = 'tess_s0001-1-1_{:04d}-{:04d}_{:04d}-{:04d}.fits'.format(limits[0][0], limits[0][1], limits[1][0], limits[1][1])
    print(limits, OUTNAME)
    make_my_coutout(FILENAME, OUTNAME, limits, output_path='Cutouts', verbose=True, overwrite=False)
    # if i > 10: break

# Loop over everything
sector = 's0001'
for cam in (1,2,3,4):
    for ccd in (1,2,3,4):
        FILENAME = 'Data/tess2018206192942-{}-{}-{}-0120-s_ffic.fits'.format(sector, cam, ccd)
        # output_path = 'Cutouts-{}-{}'.format(cam, ccd)
        output_path = 'Cutouts'
        # full_limits = get_limits(shape=(2000, 2100), width=100, start=50)
        full_limits = get_limits()
        for i, limits in enumerate(full_limits):
            OUTNAME = 'tess_{}-{}-{}_{:04d}-{:04d}_{:04d}-{:04d}.fits'.format(sector, cam, ccd, limits[0][0],
                                                                              limits[0][1], limits[1][0], limits[1][1])
            # print(limits, OUTNAME)
            make_my_coutout(FILENAME, OUTNAME, limits, output_path=output_path, verbose=True, overwrite=False)
