import numpy as np
import healpy as hp
from astropy import units as u
import astropy.coordinates as coord
from shapely import geometry
from matplotlib.patches import Polygon
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
import statistics
from shapely.ops import linemerge, unary_union, polygonize, split
from shapely.geometry import Point
# from Teglon_Shape import *


# 2023-04-03 OLD
# from src.objects.Teglon_Shape import *
# 2023-04-03 NEW
from web.src.objects.Teglon_Shape import *


class Pixel_Element(Teglon_Shape):
    def __init__(self, index, nside, prob, pixel_id=None, mean_dist=None, stddev_dist=None):

        self.tile_references = []

        # Properties
        self.id = pixel_id
        self.__index = index
        self.__nside = nside
        self.__coord = None
        self.__prob = prob

        # Used to check if this object is within this "radius proxy" of the coordinate signularity
        self.__radius_proxy = 5 * np.sqrt(hp.nside2pixarea(self.nside, degrees=True) / np.pi)
        self.__multipolygon = None
        self.__query_polygon = None
        self.__query_polygon_string = None

        self.__epi = np.array([])

        # For all NSIDE < self.nside:
        # { NSIDE: Containing Pixel Index }
        self.parent_pixel = None

        self.mean_dist = mean_dist
        self.stddev_dist = stddev_dist

    def __str__(self):
        return str(self.__dict__)

    @property
    def index(self):
        return self.__index

    @property
    def nside(self):
        return self.__nside

    @property
    def prob(self):
        return self.__prob

    @property
    def coord(self):
        if not self.__coord:
            theta, phi = hp.pix2ang(self.nside, self.index)

            ra = np.rad2deg(phi)
            dec = np.rad2deg(0.5 * np.pi - theta)

            pix_coord = coord.SkyCoord(ra=np.rad2deg(phi),
                                       dec=np.rad2deg(0.5 * np.pi - theta),
                                       unit=(u.deg, u.deg))

            self.__coord = pix_coord

        return self.__coord

    # Telgon_Shape properties
    @property
    def radius_proxy(self):
        return self.__radius_proxy

    @property  # Returns list of polygon
    def multipolygon(self):

        if not self.__multipolygon:

            pixel_xyz_vertices = hp.boundaries(self.nside, pix=self.index)
            theta, phi = hp.vec2ang(np.transpose(pixel_xyz_vertices))

            ra_rad = phi
            dec_rad = (np.pi / 2. - theta)
            pixel_radian_vertices = [[ra, dec] for ra, dec in zip(ra_rad, dec_rad)]

            self.__multipolygon = [geometry.Polygon(pixel_radian_vertices)]
        return self.__multipolygon

    @property  # Returns list of polygon
    def projected_multipolygon(self):
        return self.multipolygon

    # NOTE -- always assuming RA coordinate of the form [0, 360]
    @property
    def query_polygon(self):

        if not self.__query_polygon:
            self.__query_polygon = self.create_query_polygon(initial_poly_in_radian=True)
        return self.__query_polygon

    # Formatted for MySQL WGS84 spatial reference system
    @property
    def query_polygon_string(self):

        if not self.__query_polygon_string:
            self.__query_polygon_string = self.create_query_polygon_string(initial_poly_in_radian=True)
        return self.__query_polygon_string

    def enclosed_pixel_indices(self, nside_out):

        # Sanity
        if nside_out < self.nside:
            raise ("Can't get enclosed pixel indices for lower resolution pixels!")

        if len(self.__epi) == 0:

            # Start with the central pixel, in case the size of the FOV is <= the pixel size
            self.__epi = np.asarray(
                [hp.ang2pix(self.nside, 0.5 * np.pi - self.__coord.dec.radian, self.__coord.ra.radian)])

            pixel_xyz_vertices = hp.boundaries(self.nside, pix=self.index)
            internal_pix = hp.query_polygon(nside_out, pixel_xyz_vertices, inclusive=False)

            if len(internal_pix) > 0:
                self.__epi = internal_pix

        return self.__epi

    def get_patch(self, bmap):

        patches = []
        query_polygon = self.query_polygon
        for p in query_polygon:
            ra_deg, dec_deg = zip(*[(coord_deg[0], coord_deg[1]) for coord_deg in p.exterior.coords])
            x2, y2 = bmap(ra_deg, dec_deg)
            lat_lons = np.vstack([x2, y2]).transpose()

            patch = Polygon(lat_lons)
            patches.append(patch)

        return patches

    def plot(self, bmap, ax_to_plot, **kwargs):
        query_polygon = self.query_polygon
        for p in query_polygon:
            ra_deg, dec_deg = zip(*[(coord_deg[0], coord_deg[1]) for coord_deg in p.exterior.coords])
            x2, y2 = bmap(ra_deg, dec_deg)
            lat_lons = np.vstack([x2, y2]).transpose()

            ax_to_plot.add_patch(Polygon(lat_lons, **kwargs))

    def plot2(self, ax_to_plot, **kwargs):

        query_polygon = self.query_polygon
        for p in query_polygon:
            ra_deg, dec_deg = zip(*[(coord_deg[0], coord_deg[1])
                                    for coord_deg in p.exterior.coords])

            list_o_coords = []
            for x, y in zip(ra_deg, dec_deg):
                list_o_coords.append((x, y))
            arr = np.asarray(list_o_coords)

            ax_to_plot.add_patch(Polygon(arr, transform=ax_to_plot.get_transform('world'), **kwargs))