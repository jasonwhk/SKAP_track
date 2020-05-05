import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
from astropy import wcs
from astropy.time import Time
from astroplan import Observer
from astroplan import download_IERS_A
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from functools import reduce
import operator
import sys
from optparse import OptionParser
from astropy.utils import iers

#download_IERS_A()

def unpackTuple(tup):

    return (reduce(operator.add, tup))


def convert_pixel_coordinate_to_equatorial(pixel_coordinates, bore_sight):
    """
    https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html
    CRVAL: coordinate system value at reference pixel
    CRPIX: coordinate system reference pixel
    CDELT: coordinate increment along axis
    CTYPE: name of the coordinate axis
    """
    step = 1 / 10000000000.

    wcs_properties = wcs.WCS(naxis=2)
    wcs_properties.wcs.crpix = [0, 0]
    wcs_properties.wcs.cdelt = [step, step]
    wcs_properties.wcs.crval = bore_sight
    wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    scaled_pixel_coordinats = np.array(pixel_coordinates) / step
    equatorial_coodinates = wcs_properties.wcs_pix2world(
        scaled_pixel_coordinats, 0)
    return equatorial_coodinates


def _append(x, y, z, x_append, y_append, z_append, rotation):
    # print(,y)
    a, b = rotate((x_append, y_append), rotation)
    x.append(a)
    #print("adding point {} {} {}".format(x_append, y_append, z_append))
    y.append(b)
    z.append(z_append)


def rotate(point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in degrees.
    """
    ox, oy = 0, 0
    px, py = point

    qx = ox + np.cos(np.deg2rad(angle)) * (px - ox) - \
        np.sin(np.deg2rad(angle)) * (py - oy)
    qy = oy + np.sin(np.deg2rad(angle)) * (px - ox) + \
        np.cos(np.deg2rad(angle)) * (py - oy)
    return qx, qy


def OTF(bore_sight, frame, time_step, start_time, step, lat, lon, alt, rotation, x_length, y_length, seperation, plot, alt_limit):
    prototype_dish = EarthLocation(
        lat=lat * u.deg, lon=lon * u.deg, height=alt * u.m)
    if frame == "altaz":
        centriod = SkyCoord(bore_sight[0], bore_sight[
                            1], unit='deg', frame="altaz", location=prototype_dish, obstime=Time.now())
        bore_sight = (centriod.icrs.ra.deg, centriod.icrs.dec.deg)
    elif frame == "galactic" or frame == "ircs":
        centriod = SkyCoord(bore_sight[0], bore_sight[
                            1], unit='deg', frame=frame)
        bore_sight = (centriod.icrs.ra.deg, centriod.icrs.dec.deg)

    step_y = int(np.ceil((y_length / seperation)))
    step_x = int(np.ceil((x_length / seperation)))
    start_x = int(step_x / 2 - step_x - 1)
    end_x = int(step_x / 2 + 2)
    start_y = int(step_y / 2 - step_y)
    end_y = int(step_y / 2 + 1)
    pixel_coordinates_x = []
    pixel_coordinates_y = []
    pixel_coordinates_z = []
    reverse = 0
    for j in range(start_y, end_y):
        x = []
        y = []
        z = []

        for i in range(start_x, end_x):
            if i == start_x and j == start_y and j % 2 == 0:
                _append(x, y, z, i - 1, j, 0, rotation)
            else:
                pass
            if i == start_x and j % 2 == 0 and j != start_y:
                #print("I am here")
                _append(x, y, z, i - 1, j - 1, 0, rotation)
                _append(x, y, z, i - 1.8661, j - 0.5, 0, rotation)
                _append(x, y, z, i - 1, j, 0, rotation)
            if i == start_x or i == end_x - 1:
                #print("i am at {} start_x + 1 = {}, start_x = {}".format(i, start_x + 1, start_x))
                _append(x, y, z, i, j, 0, rotation)
            else:
                _append(x, y, z, i, j, 1, rotation)
            if i == end_x - 1 and j % 2 != 0 and reverse != 0:
                #print("I am here")
                _append(x, y, z, i + 1, j, 0, rotation)
                _append(x, y, z, i + 1.8661, j - 0.5, 0, rotation)
                _append(x, y, z, i + 1, j - 1, 0, rotation)

        if j % 2 != 0:
            x.reverse()
            y.reverse()
            z.reverse()
            reverse = 0
            pixel_coordinates_x.append(x)
            pixel_coordinates_y.append(y)
            pixel_coordinates_z.append(z)

        else:
            reverse = 1
            pixel_coordinates_x.append(x)
            pixel_coordinates_y.append(y)
            pixel_coordinates_z.append(z)

    x = unpackTuple(pixel_coordinates_x)
    y = unpackTuple(pixel_coordinates_y)
    z = unpackTuple(pixel_coordinates_z)
    x[:] = [x * seperation for x in x]
    y[:] = [y * seperation for y in y]
    pixel_coordinates = list(zip(x, y))
    if plot == 1:
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.plot(x, y, "-o")
        plt.scatter(x, y)

    data = convert_pixel_coordinate_to_equatorial(
        pixel_coordinates, bore_sight)

    x_val = [x[0] for x in data]
    y_val = [x[1] for x in data]
    if plot == 1:
        plt.subplot(1, 2, 2)
        plt.scatter(x_val, y_val, zorder=1)
        plt.scatter(bore_sight[0], bore_sight[1], color='r', zorder=2)
        plt.plot(x_val, y_val, "-o", zorder=0)
        plt.xlabel("ra")
        plt.ylabel("dec")

    v_time = [start_time + i * u.s * time_step for i in range(0, len(x_val))]
    sc = SkyCoord(x_val, y_val, unit='deg', frame="icrs")
    prototype_dish_observer = Observer(location=prototype_dish)
    parallactic_angle = prototype_dish_observer.parallactic_angle(
        time=v_time, target=sc)
    center = SkyCoord(x_val[0], y_val[0], unit='deg', frame="icrs")
    grid_altaz = sc.transform_to(
        AltAz(obstime=start_time, location=prototype_dish))
    source_altaz = center.transform_to(
        AltAz(obstime=start_time, location=prototype_dish))
    if plot == 1:
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.scatter(grid_altaz.az.deg, grid_altaz.alt.deg, zorder=0)
        plt.plot(grid_altaz.az.deg, grid_altaz.alt.deg, "-o", zorder=1)
        plt.scatter(source_altaz.az.deg,
                    source_altaz.alt.deg, color="g", zorder=2)
        plt.xlabel("Az")
        plt.ylabel("Alt")

        sep = []
        frame = "altaz"
        for i in range(0, len(x_val)):
            c1 = SkyCoord(grid_altaz.az.deg[i] * u.deg,
                          grid_altaz.alt.deg[i] * u.deg, frame=frame)
            if i != len(x_val) - 1:
                c2 = SkyCoord(grid_altaz.az.deg[
                              i + 1] * u.deg, grid_altaz.alt.deg[i + 1] * u.deg, frame=frame)
                sep.append(c1.separation(c2).deg / time_step)

        plt.subplot(1, 2, 2)
        plt.plot(sep[:], color="y")

    alt = grid_altaz.alt.deg[:]
    az = grid_altaz.az.deg[:]
    #print("total scan time {}s".format((v_time[-1]- v_time[0])*86400))
    converted = list(zip(az, alt, v_time, z))
    for i in range(len(converted)):
        if converted[i][1] > alt_limit:
            print("{0} {1:3.8f} {2:3.8f} {3} {4}".format(converted[i][2].mjd, converted[
                  i][0], converted[i][1], converted[i][3], parallactic_angle.deg[i]), file=sys.stdout, flush=True)
        else:
            break
#print("{0} {1:3.8f} {2:3.8f}".format(time.isot, source_altaz.az.deg, source_altaz.alt.deg), file=sys.stdout, flush=True)
    if plot == 1:
        plt.show()


def cross_scan(center, width, duration, start_time, time_step, position_angle, plot, lat, lon, alt, alt_limit):
    pos_step = time_step * duration
    speed = width / duration
    start_separation = width / 2
    separation = width / pos_step
    #print("Width/2 =",start_separation)
    #print("positional sepration = ", separation)
    #print("positional step = ", pos_step)
    #print("scan speed =", speed)
    #print("scan time =", width/speed)
    position_angle = position_angle * u.deg
    x = []
    y = []
    t = []
    f = []
    c1 = SkyCoord(center[0] * u.deg, center[1] * u.deg, frame='icrs')
    #center = SkyCoord(x_val[0], y_val[0], unit='deg', frame="icrs")
    prototype_dish = EarthLocation(
        lat=lat * u.deg, lon=lon * u.deg, height=alt * u.m)
    prototype_dish_observer = Observer(location=prototype_dish)
    altaz = c1.transform_to(
        AltAz(obstime=start_time, location=prototype_dish))
    # print(altaz.alt.deg)
    c1 = SkyCoord(altaz.az.deg, altaz.alt.deg, unit='deg', frame='altaz')

    start_pos = c1.directional_offset_by(position_angle, start_separation)
    v_time = [start_time + i * u.s *
              time_step for i in range(0, int(pos_step) * 2 + 2)]
    counter = 0
    for i in range(0, int(pos_step) + 1):
        new_pos = start_pos.directional_offset_by(
            position_angle, -i * separation)
        x.append(new_pos.az.deg)
        y.append(new_pos.alt.deg)
        t.append(v_time[i])
        if i != int(pos_step):
            f.append(1)
        else:
            f.append(0)
        counter = i
    position_angle = position_angle + 90 * u.deg
    start_pos = c1.directional_offset_by(position_angle, start_separation)
    for i in range(0, int(pos_step) + 1):
        new_pos = start_pos.directional_offset_by(
            position_angle, -i * separation)
        x.append(new_pos.az.deg)
        y.append(new_pos.alt.deg)
        t.append(v_time[counter] + 10 * u.s)
        if i != 0:
            f.append(1)
        else:
            f.append(0)
        counter += 1

    if plot == 1:
        plt.figure(figsize=(5, 5))
        plt.scatter(x, y)
        plt.scatter(altaz.az.deg, altaz.alt.deg, color="r")
        plt.ylabel('Alt')
        plt.xlabel('Az')
    #prototype_dish_observer.parallactic_angle(time=time[i], target=sc).deg
    #print(len(t), len(x), len(y), len(v_time))

    for i in range(len(v_time)):
        #print(x[i], y[i])
        altaz = SkyCoord(x[i], y[i], unit='deg', frame="altaz",
                         obstime=v_time[i], location=prototype_dish)
        # altaz = AltAz(np.deg2rad(x[i])*u.rad, np.deg2rad(y[i])*u.rad:,
        # obstime=v_time[i], location=prototype_dish)
        sc = altaz.transform_to('icrs')
        print("{0} {1:3.8f} {2:3.8f} {3} {4}".format(
            t[i].mjd, x[i], y[i], f[i], prototype_dish_observer.parallactic_angle(time=v_time[i], target=sc).deg))
    if plot == 1:
        plt.show()


def simple_track(ra, dec, frame, time, input_lat, input_lon, alt, plot, alt_limit):
    prototype_dish = EarthLocation(
        lat=input_lat * u.deg, lon=input_lon * u.deg, height=alt * u.m)
    sc = SkyCoord(ra, dec, unit='deg', frame=frame, equinox="J2000",
                  obstime=time[0], location=prototype_dish)
    sc = sc.transform_to('icrs')
    prototype_dish_observer = Observer(location=prototype_dish)

    source_altaz = sc.transform_to(
        AltAz(obstime=time, location=prototype_dish))
    for i in range(len(source_altaz)):
        if source_altaz[i].alt.deg > alt_limit:
            print("{0} {1:3.8f} {2:3.8f} {3} {4}".format(time[i].mjd, source_altaz[i].az.deg,
                                                     source_altaz[i].alt.deg, 1, prototype_dish_observer.parallactic_angle(time=time[i], target=sc).deg, file=sys.stdout, flush=True))
        else:
            break
    if plot == 1:
        plt.scatter(source_altaz.az.deg, source_altaz.alt.deg)
        plt.show()


def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option('--x', dest='x', default=20, type=float,
                      help='x of the source in deg')
    parser.add_option('--y', dest='y', default=30, type=float,
                      help='y of the source in deg')
    parser.add_option('--frame', dest='frame', default='icrs', type=str,
                      help='Coordinate frame (icrs, fk5, fk4, galactic, altaz')
    parser.add_option('--start', dest='start', default=Time.now(), type=str,
                      help='Start time in ISO8601 UTC format')
    parser.add_option('--end', dest='end', default=Time.now() + 600 * u.s, type=str,
                      help='End time in ISO8601 UTC format')
    parser.add_option('--step', dest='step', default=0.5, type=float,
                      help='time step in second')
    parser.add_option('--rotation', dest='rotation', default=0, type=float,
                      help='Rotation in degree, anti-clockwise')
    parser.add_option('--xlength', dest='x_length', default=1, type=float,
                      help='x length in degrees')
    parser.add_option('--ylength', dest='y_length', default=1, type=float,
                      help='y length in degrees')
    parser.add_option('--seperation', dest='seperation', default=0.1, type=float,
                      help='Seperation in degrees in between OTF points')
    parser.add_option('--scantype', dest='type', default="track", type=str,
                      help='Scan type (track, OTF or cross_scan) for now')
    parser.add_option('--plot', dest='plot', default=0, type=int,
                      help='Plot the scan pattern or not, default = 0')
    parser.add_option('--duration', dest='duration', default=30, type=float,
                      help='Duration of the scan')
    parser.add_option('--length', dest='length', default=1, type=float,
                      help='Length of the cross scan')
    parser.add_option('--dry-run', dest='dry_run', default=None, type=str,
                      help='Dry run for SCU testing, passing "OK?"" will print out "OK!"')
    parser.add_option('--alt-limit', dest='alt_limit', default=20.0, type=float,
                      help='Elevation drive limit, default = 20.0 deg')


    # parser.add_option('-a', '--lat', dest='lat', type=float,
    #                  help='latitude of the observatory in deg')
    # parser.add_option('-o', '--lon', dest='lon', type=float,
    #                  help='lontitude of the observatory in deg')
    # parser.add_option('-l', '--alt', dest='alt', type=float,
    #                  help='height of the observatory in meter')
    (opts, args) = parser.parse_args()

    seperation = 0.1  # seperation between points
    x = opts.x
    y = opts.y
    bore_sight = (x, y)
    time_step = opts.step
    lat = -30.712
    lon = 21.411
    alt = 390
    frame = opts.frame
    rotation = opts.rotation
    x_length = opts.x_length  # in degrees
    y_length = opts.y_length  # in degrees
    start_time = opts.start
    duration = opts.duration
    length = opts.length
    plot = opts.plot
    if opts.dry_run == "OK?":
        print("OK!")
        return
    elif opts.type == "OTF":
        OTF(bore_sight, frame, time_step, start_time, seperation, lat,
            lon, alt, rotation, x_length, y_length, seperation, plot, opts.alt_limit)
    elif opts.type == "track":
        t_start, t_end = Time(opts.start), Time(opts.end) + 1 * u.s
        #duration = t_end - t_start
        duration = float(opts.duration)
       	bins = np.ceil(duration / opts.step)
        #print(bins)
        v_time = [t_start + i * u.s * opts.step for i in range(0, int(bins))]
        simple_track(x, y, frame, v_time,
                     lat, lon, alt, plot, opts.alt_limit)
    elif opts.type == "cross_scan":
        cross_scan(bore_sight, length * u.deg, duration, start_time,
                   time_step, rotation, plot, lat, lon, alt, opts.alt_limit)
    else:
        print("invaild scan mode '{}', available scan types are : 'OTF', 'track' & 'cross_scan'".format(opts.type))

if __name__ == "__main__":
    main()
