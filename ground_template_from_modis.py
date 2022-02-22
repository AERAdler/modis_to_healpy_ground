from pyhdf.SD import SD, SDC
import healpy as hp
import numpy as np
from scipy import interpolate

h = 6.62e-34
c = 3e8
k_b = 1.38e-23
T_cmb = 2.72548 

def tb2b(tb, nu):
    #Convert blackbody temperature to spectral
    x = h*nu/(k_b*tb)
    return 2*h*nu**3/c**2/(np.exp(x) - 1)

def dBdT(tb, nu):
    x = h*nu/(k_b*tb)
    slope = 2*k_b*nu**2/c**2*((x/2)/np.sinh(x/2))**2
    return slope

def rotate_to_telescope(inmap, tel_lat, tel_lon):
    nside = hp.npix2nside(inmap.shape[0])
    rlon = np.radians(tel_lon)-np.pi
    rlat = np.pi/2.-np.radians(tel_lat)
    x0, y0, z0 = hp.pix2vec(512, np.arange(hp.nside2npix(512)))
    x1 =  x0*np.cos(rlat)+z0*np.sin(rlat)
    z = -x0*np.sin(rlat)+z0*np.cos(rlat)
    x = x1*np.cos(rlon)-y0*np.sin(rlon)
    y = x1*np.sin(rlon)+y0*np.cos(rlon)
    pix_prime = hp.vec2pix(512, x,y,z)
    return inmap[pix_prime]

def telescope_view_angles(nside, telescope_altitude, surface_height=0, r_earth =6.371e6):
    """
    Calculates the coordinates of a map with a given nside projected on the telescope's view,
    given the telescope's altitude, the height of the ground, and the radius of the earth 
    at the telescope's location. Doesn't account for variations in ground height throughout 
    the fov.
    Returns the coordinates visible in the input map, and their projection in the telescope's map.
    """
    r_ground = r_earth+surface_height
    tele_h_abg = telescope_altitude - surface_height
    theta_fov = np.arccos((r_ground)/(r_ground+tele_h_abg))
    
    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    theta_visible = theta[theta<theta_fov]
    phi_visible = phi[theta<theta_fov]
    theta_from_tel = np.arctan2(r_anta*np.sin(theta_visible), 
        r_anta*(1-np.cos(theta_visible))+tele_h_abg)
    theta_from_tel = np.pi-theta_from_tel #Up is down
    phi_from_tel = np.pi-phi_visible #Therefore left is right
    return(theta_visible, phi_visible, theta_from_tel, phi_from_tel)

    def ground_template(world_map, theta_visible, phi_visible, theta_from_tel, phi_from_tel, 
                      nside_out=128, cmb=True, freq=95., frac_bwidth=.2):
    """
    Creates a ground template given a world map and sets of coordinates (see telescope_view_angles)
    Returns a filled-out ground template
    options: 
    nside_out (int): healpix NSIDE of the output map
    cmb(bool): Whether to convert the output map to CMB brightness temperature units (muK) 
    freq(float): Center of the frequency band
    frac_bwidth(float): Band width, as a fraction of the central frequency
    """
    nside_world = hp.npix2nside(world_map.shape[0])
    ground_map = np.ones(hp.nside2npix(nside_out))*hp.UNSEEN
    ground_pix = hp.ang2pix(nside_out, theta_from_tel, phi_from_tel)
    ground_map[ground_pix] = world_map_from_telescope[hp.ang2pix(nside_world, theta_visible, phi_visible)]
    if cmb:
        freq_band=np.linspace(freq*(1-frac_bwidth/2.), freq*(1+frac_bwidth/2.), 201)
        for i, tb in enumerate(ground_map):
            if tb!=hp.UNSEEN:
                bolo = np.trapz(tb2b(tb, freq_band), freq_band)
                corr = np.trapz(dBdT(T_cmb, freq_band), freq_band)
                ground_map[i] = bolo/corr*1e6
    pix_pow = int(np.log2(nside_out))
    map_low = hp.ud_grade(ground_map, 2)
    for i in range(pix_pow+1):

        nside = int(2**i)
        map_normal = hp.ud_grade(ground_map, nside)
        map_horizon = np.amin(hp.ang2pix(nside, theta_from_tel, phi_visible))
        map_normal[map_horizon:] = np.where(
                                            map_normal[map_horizon:]!=hp.UNSEEN,
                                            map_normal[map_horizon:], 
                                            hp.ud_grade(map_low, nside)[map_horizon:]
                                            )
        map_low = map_normal       
    
    return(map_normal)

def template_from_position(earth_map, lat, lon, altitude, nside_out=128, 
    cmb=True, freq=95., frac_bwidth=.2):

    nside_world = hp.npix2nside(earth_map.shape[0])
    earth_rot = rotate_to_telescope(earth_map, lat, lon)
    theta_visible, phi_visible, theta_from_tel, phi_from_tel = telescope_view_angles(
        nside_world, altitude, surface_height=0, r_earth =6.371e6)
    ground_temp = ground_template(world_map, 
                    theta_visible, phi_visible, theta_from_tel, phi_from_tel, 
                    nside_out=nside_out, cmb=cmb, freq=freq, frac_bwidth=frac_bwidth)




