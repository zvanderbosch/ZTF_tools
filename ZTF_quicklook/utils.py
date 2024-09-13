import time
import requests
import colorsys
import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
import matplotlib.colors as mc
import lmfit as lmf

from astropy.timeseries import LombScargle as ls
from astropy.coordinates import EarthLocation
from astropy.time import Time,TimeDelta
from io import StringIO
from http.client import responses
from lmfit import Model
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier


LOGIN_URL = "https://irsa.ipac.caltech.edu/account/signon/login.do"

def BaryCorr(lc_data):
    """
    Function to apply barycentric corrections to
    timestamps and convert them from beginning
    to mid-exposure.
    
    Input:
    ------
    lc_data: DataFrame
        Pandas DataFrame containing the ZTF light 
        curve downloaded from IRSA
        
    Returns:
    --------
    tbjd: array
        Mid-exposure timestamps in JD format with
        barycentric corrections applied.
    
    """
    
    loc = EarthLocation.of_site('palomar')
    c = coord.SkyCoord(
        lc_data.ra.mean(),
        lc_data.dec.mean(),
        unit='deg',
        frame='icrs'
    )
    
    t = Time(
        lc_data.mjd.values,
        scale='utc',
        format='mjd',
        location=loc
    )
    t += TimeDelta(0.5*lc_data.exptime.values * u.s)
    ltt_bary = t.light_travel_time(c)
    tbjd = t.tdb.jd + ltt_bary.jd
    
    return tbjd


def get_cookie(username, password):
    """
    Function to get a cookie from the IRSA login service
    
    Input:
    ------
    username: str
        IRSA username
    password: str
        IRSA password
        
    Returns:
    --------
    irsa_login_cookie: RequestsCookieJar
        IRSA login cookie
    """
    
    url = "{:s}?josso_cmd=login&josso_username={:s}&josso_password={:s}".format(
        LOGIN_URL, 
        username, 
        password
    )
    irsa_login_cookie = requests.get(url).cookies
    
    return irsa_login_cookie


def download_ZTF(ra,dec,radius=3.0,irsa_cookie=None,release='latest'):
    """
    Function to contruct URL and download ZTF data
    
    Input:
    ------
    ra: float
        Right ascension in decimal degrees
    dec: float
        Declination in decimal degrees
    radius: float
        Con search radius in arcseconds (default = 3.0)
    irsa_cookie: None or RequestsCookieJar
        A requests cookie based on IRSA username/password
        (default = None)
    release: str or int
        Which ZTF data release to download data from
        (default = None).
        
    Returns:
    --------
    data: DataFrame
        ZTF light curve data
    """

    # Generate URL to query the ZTF Light Curve Database
    base_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?"
    
    if release == 'latest': # Don't pass anything to COLLECTION
        pos_string = "POS=CIRCLE {:.5f} {:.5f} {:.5f}&FORMAT=csv".format(
            ra,
            dec,
            radius/3600
        )
    else:
        try:
            release_int = int(release)
        except Exception as e:
            raise(f'Invalid release input {release}...{e}' )
                  
        collection = f"ztf_dr{release_int}"
        pos_string = "POS=CIRCLE {:.5f} {:.5f} {:.5f}&COLLECTION={:s}&FORMAT=csv".format(
            ra,
            dec,
            radius/3600,
            collection
        )         
    search_url = base_url + pos_string

    # Perform request and decode response into a Pandas DataFrame
    if irsa_cookie is not None:
        response = requests.get(search_url, cookies=irsa_cookie)
    else:
        response = requests.get(search_url)
    
    # Raise error if bad status code
    response.raise_for_status()
    data = pd.read_csv(StringIO(response.content.decode('utf-8')))

    if len(data) == 0:
        data['bjd'] = []
        data['bsec'] = []
        return data
    
    # Apply barycentric corrections to time stamps
    bjd_midtimes = BaryCorr(data)
    bjdref = np.mean(bjd_midtimes)
    data['bjd'] = bjd_midtimes
    data['bsec'] = (bjd_midtimes - bjdref) * 86400.

    return data


def LC_normalize(lc_data):
    
    medmag = lc_data.mag.median()
    normflux = 10.0**(0.4*(medmag-lc_data.mag.values)) - 1.0
    normfluxerr = lc_data.magerr.values / 1.086
    
    lc_data['flux'] = normflux
    lc_data['fluxerr'] = normfluxerr

    return lc_data


def calc_lsp(t,a,f):
    """
    Function to Calculate Lomb-Scargle Periodogram (LSP) 
    with amplitude Units that match the input flux units.
    
    Input:
    ------
    t: array
        Array of light curve timestamps
    a: array
        Array of light curve flux values
    f: array
        Frequency grid for periodogram, in units that match 1/t
        
    Returns:
    --------
    norm_lsp: Numpy array
        Lomb_Scargle amplitude spectrum, normalized to 
        match the input flux units.
    """
    
    lsp = ls(t,a).power(f,normalization='psd')
    norm_lsp = np.sqrt(abs(4.0*(lsp/len(t))))
    
    return norm_lsp


def optimize_freq(time, flux, farr, freq_init, amp_init):
    """
    Function to optimize a frequency by fitting a 
    simple sineusoidal function to light curve data.
    
    Input:
    ------
    time: array
        Array of light curve timestamps
    flux: array
        Array of light curve flux values
    farr: array
        Frequency grid for periodogram, in units that match 1/time
    freq_init: float
        Initial guess for frequency
    amp_init: float
        Initial guess for amplitude
        
    Returns:
    --------
    bestfit_freq: float
        The least-squares optimized frequecy
    lsp_pw: array
        The Lomb-Scargle periodogram of the light curve
        after subtraction to best-fit sinusoidal function
    """
    
    # Sinusoidal Function
    def sine(x,freq,amp,phase):
        return amp*np.sin(2.0*np.pi*(freq*x + phase))

    # A simple constant offset
    def offset(x,offset):
        return offset
    
    # Generate an LMFIT model and parameters
    mod = Model(offset) + Model(sine)
    par = mod.make_params()

    # Set LMFIT initial values
    par['offset'].value = 0.0
    par['freq'].value = freq_init
    par['amp'].value = amp_init
    par['phase'].value = 0.0

    # Fix frequency initially for a linear least squares fit
    par['freq'].vary = False

    # Perform the linear LSQ fit
    result = mod.fit(flux, params=par, x=time)

    # Now unfix the frequenciy and perform non-linear LSQ fit
    new_par = result.params
    new_par['freq'].vary = True
    result = mod.fit(flux, params=new_par, x=time)

    # Save Fit Results
    freq_vals = []
    freq_errs = []
    for name,param in result.params.items():
        if 'freq' in name:
            freq_vals.append(param.value)
            freq_errs.append(param.stderr)
    print("  Initial Frequency: {:.6f} uHz".format(freq_init*1e6))
    print("Optimized Frequency: {:.6f} uHz".format(freq_vals[0]*1e6))
    print("  Initial Period   : {:.4f} min".format(1./freq_init/60))
    print("Optimized Period   : {:.4f} min".format(1./freq_vals[0]/60))
    bestfit_freq = freq_vals[0]

    # Subtract LMFIT model from the observed flux and recalculate the LS-Periodogram
    flux_pw = flux - result.best_fit
    lsp_pw = calc_lsp(time,flux_pw,farr)
    
    return bestfit_freq, lsp_pw


def PhaseFold(time, pfold, pdot=None, t0=None):
    """
    Function to phase fold an array of times by
    a folding period, accounting for orbital decay.
    
    This function was slightly modified from Kevin
    Burdge's LC_Tools.py code found on github.
    
    
    Input:
    ------
    time: array
        Array of light curve timestamps
    pfold: float
        Folding period, same units as time
    pdot: float or None
        Rate of period change (default = None)
    t0: float or None
        Reference epoch (default = None)
        
    Returns:
    --------
    phases: array
        Phase of each time stamp
    """


    # If no reference epoch is supplied, use earliest time
    if t0 is None:
        time -= min(time)
    else:
        time -= t0
        
    if pdot is None:
        phases = (time / pfold) % 1.0
    else:
        phases = ((time-0.5*(pdot/pfold)*time**2) / pfold) % 1.0

    return phases


def PhaseFoldBinned(t,y,dy,pfold,nbin=50,pdot=None,t0=None):
    """
    Weighted binning of a lightcurve with times t, 
    fluxes y, and flux errors dy, at period pfold.
    Also accounts for any orbital decay.
    
    This function was slightly modified from Kevin
    Burdge's LC_Tools.py code found on github.
    
    Input:
    ------
    t: array
        Array of light curve timestamps
    y: array
        Array of light curve flux/mag values
    dy: array
        Array of light curve flux/mag uncertainties
    pfold: float
        Folding period, same units as time
    nbin: int
        Number of phase bins (default = 50)
    pdot: float or None
        Rate of period change (default = None)
    t0: float or None
        Reference epoch (default = None)
        
    Returns:
    --------
    binned_lc: array
        Phase of each time stamp
    """
    
    mean_phases = np.linspace(0.0, 1.0-1./nbin, nbin)
    phases = PhaseFold(t, pfold, pdot, t0)
    lightcurve = np.array((phases,y,dy)).T
    lightcurve = lightcurve[np.argsort(lightcurve[:,0])]
    
    binned_LC = []
    for i in mean_phases:
        
        lightcurve_bin = lightcurve[
            (lightcurve[:,0]>i) &
            (lightcurve[:,0]<i+1./nbin)
        ]
        
        if len(lightcurve_bin) == 0:
            binned_LC.append(
                (i+0.5/nbin, 
                 np.nan, 
                 np.nan)
            )
            continue
            
        weights = 1./(lightcurve_bin[:,2]**2)
        weighted_mean_flux = np.sum(lightcurve_bin[:,1]*weights)/np.sum(weights)
        weighted_mean_flux_error = np.sqrt(1./np.sum(weights))
        binned_LC.append(
            (i+0.5/nbin, 
             weighted_mean_flux, 
             weighted_mean_flux_error)
        )
    binned_LC=np.array(binned_LC)

    return binned_LC


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) 
    by the given amount. Input can be matplotlib color 
    string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    
    try:
        c = mc.cnames[color]
    except:
        c = color
        
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    rgb = colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
    hexcol = mc.to_hex(rgb)
    
    return hexcol


def GaiaQuery(obj_coord):
    """
    Function to perform Gaia DR3 query
    
    Input:
    ------
    obj_coord: Astropy coord object
        Array of light curve timestamps
        
    Returns:
    --------
    gaia_data: list or None
        List of gaia DR3 parameters, including
        the Gaia designation, G-band mag,
        parallax, BP-RP color, and absolute
        G-band magnitude. None if the cone
        search query returns no results.
    """

    # Set Gaia search parameters
    Gaia.ROW_LIMIT = 5
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    
    # Query Gaia catalog
    query = True
    while query:
        try:
            gaia_job = Gaia.cone_search_async(
                obj_coord, 
                radius=3.0*u.arcsec
            )
            query = False
        except:
            print('GAIA query failed, retrying...')
            time.sleep(1.0)
            continue
    gaia_result = gaia_job.get_results().to_pandas()

    # Get results for nearest source
    if len(gaia_result) > 0:
        gaiaID = gaia_result.SOURCE_ID.iloc[0]
        designation = gaia_result.DESIGNATION.iloc[0]
        gmag = gaia_result.phot_g_mean_mag.iloc[0]
        bpmag = gaia_result.phot_bp_mean_mag.iloc[0]
        rpmag = gaia_result.phot_rp_mean_mag.iloc[0]
        parallax = gaia_result.parallax.iloc[0]
        bprp = bpmag - rpmag
        if (np.isnan(parallax)) | (parallax < 0):
            Gabs = np.nan
        else:
            Gabs = gmag + 5.0*np.log10(parallax/100)
        Gabs_lower_limit = Gabs
        Gabs_upper_limit = Gabs

        # Now query Bailer-Jones 2021 catalog for distances
        query = True
        while query:
            try:
                viz_job = Vizier.query_region(
                    obj_coord, 
                    radius=3*u.arcsec,
                    catalog='I/352/gedr3dis',
                    column_filters={
                        'Source': f'={gaiaID}'
                    }
                )
                query = False
            except:
                print('VIZIER query failed, retrying...')
                time.sleep(1.0)
                continue

        if len(viz_job) > 0:
            viz_result = viz_job[0].to_pandas()
            if len(viz_result) > 0:
                median_dist = viz_result.rgeo.iloc[0]
                lower_dist = viz_result.b_rgeo.iloc[0]
                upper_dist = viz_result.B_rgeo.iloc[0]
                Gabs = gmag + 5.0*(1.0 - np.log10(median_dist))
                Gabs_lower_limit = gmag + 5.0*(1.0 - np.log10(lower_dist))
                Gabs_upper_limit = gmag + 5.0*(1.0 - np.log10(upper_dist))

        # Generate return data
        gaia_data = {
            'designation':designation, 
            'phot_g_mean_mag':gmag, 
            'parallax':parallax, 
            'bprp':bprp, 
            'M_G':Gabs,
            'M_G_lowlim':Gabs_lower_limit,
            'M_G_upplim':Gabs_upper_limit
        }

        return gaia_data

    else:
        return None


def get_BJD_T0(filename):
    
    with open(filename) as file:
        for line in file.readlines():
            if line[0] == "#":
                if "# BJED" in line:
                    ls = line.split("=")
                    bjd = float(ls[1].split("#")[0].strip())
                    return bjd
                else:
                    continue
            else:
                break

                


def lcBin(time, flux, fluxerr, nbin=2):

    numpoints = len(time)
    med_texp = np.median(np.diff(time))

    numbins = int(np.floor(numpoints / nbin))

    binned_time = []
    binned_flux = []
    binned_fluxerr = []
    for i in range(numbins):

        # Get data points within the bin
        bin_times = time[i*nbin:i*nbin+nbin]
        bin_fluxes = flux[i*nbin:i*nbin+nbin]
        bin_fluxerrs = fluxerr[i*nbin:i*nbin+nbin]
        bin_weights = 1./(bin_fluxerrs**2)

        # Check that there aren't any major time gaps
        ttol = 0.10 # Fractional tolerance on expected time range
        bin_trange = max(bin_times) - min(bin_times)
        if bin_trange > (1+ttol) * med_texp * (nbin-1):
            continue
        else:
            tbin = np.mean(bin_times)
            fbin = np.sum(bin_weights*bin_fluxes) / np.sum(bin_weights)
            ebin = np.sqrt(np.sum(bin_fluxerrs**2))/nbin

            binned_time.append(tbin)
            binned_flux.append(fbin)
            binned_fluxerr.append(ebin)


    binned_time = np.array(binned_time)
    binned_flux = np.array(binned_flux)
    binned_fluxerr = np.array(binned_fluxerr)

    return binned_time, binned_flux, binned_fluxerr