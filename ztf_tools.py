import io
import os
import sys
import copy
import json
import time
import requests
import numpy as np
import pandas as pd
from io import StringIO
from astropy.io import fits,ascii
from astropy import wcs
from astropy.time import Time

"""
This script contains several functions useful for
retrieving a variety of ZTF data products such as:
    - Reference images (from IPAC)
    - Science images (from IPAC)
    - Light Curves (from IPAC)
    - Transient Alert Packets (from MARS/LCO)

Author: Zach Vanderbosch
"""


# Function to retrieve URLs for reference image cutouts
# centered on a given RA and DEC
def get_ref_urls(ra,dec,size=45):
    """
        ra   = Right Ascension in Decimal Degrees
        dec  = Declination in Decimal Degrees
        size = Image size in arcseonds

    Returns:
        urls = list of URLs used to download image data
    """
    
    # First get some info related to the reference image
    search_url = 'https://irsa.ipac.caltech.edu/ibe/search/ztf/products/ref'
    pos_url = search_url + '?POS={:.3f},{:.3f}'.format(ra,dec)
    r1 = requests.get(pos_url)
    im_table = ascii.read(r1.text).to_pandas()

    # Get image meta data
    urls = []
    maxbit = 33554432 # Quality cut on the INFOBITS parameter
    num_image = len(im_table[im_table.infobits < maxbit]) # Number of images in table
    for im in range(num_image):
        field = str(im_table[im_table.infobits < maxbit].field.iloc[im]).zfill(6)
        filtcode = str(im_table[im_table.infobits < maxbit].filtercode.iloc[im])
        ccdid = str(im_table[im_table.infobits < maxbit].ccdid.iloc[im]).zfill(2)
        qid = str(im_table[im_table.infobits < maxbit].qid.iloc[im])

        # Construct the Image Download URL
        data_url = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/ref'
        spec_url = '/{}/field{}/{}/ccd{}/q{}/'.format(field[0:3],field,filtcode,ccdid,qid)
        imname = 'ztf_{}_{}_c{}_q{}_refimg.fits'.format(field,filtcode,ccdid,qid)
        condis = '?center={:.5f},{:.5f}&size={}arcsec&gzip=false'.format(ra,dec,size)
        urls.append(data_url + spec_url + imname + condis)

    return urls

# Function to retrieve URLs for a science image cutouts
# centered on a given RA and DEC
def get_sci_urls(ra,dec,size=45):
    """
        ra   = Right Ascension in Decimal Degrees
        dec  = Declination in Decimal Degrees
        size = Image size in arcseconds

    Returns:
        urls = list of URLs used to download image data
    """
    
    # First get some info related to the reference image
    search_url = 'https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci'
    pos_url = search_url + '?POS={:.3f},{:.3f}'.format(ra,dec)
    r1 = requests.get(pos_url)
    im_table = ascii.read(r1.text).to_pandas()

    # Get image meta data
    urls = []
    maxbit = 33554432 # Quality cut on the INFOBITS parameter
    num_image = len(im_table[im_table.infobits < maxbit]) # Number of images in table
    for im in range(num_image):
        field = str(im_table[im_table.infobits < maxbit].field.iloc[im]).zfill(6)
        filtcode = str(im_table[im_table.infobits < maxbit].filtercode.iloc[im])
        ccdid = str(im_table[im_table.infobits < maxbit].ccdid.iloc[im]).zfill(2)
        qid = str(im_table[im_table.infobits < maxbit].qid.iloc[im])
        filefracday = str(im_table[im_table.infobits < maxbit].filefracday.iloc[im])
        imgtypecode = str(im_table[im_table.infobits < maxbit].imgtypecode.iloc[im])

        year = filefracday[0:4]
        month = filefracday[4:6]
        day = filefracday[6:8]
        fracday = filefracday[8:]

        # Construct the Image Download URL
        data_url = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci'
        date_url = '/{}/{}/{}/'.format(year,month+day,fracday)
        imname = 'ztf_{}_{}_{}_c{}_{}_q{}_sciimg.fits'.format(filefracday,field,filtcode,ccdid,imgtypecode,qid)
        condis = '?center={:.5f},{:.5f}&size={}arcsec&gzip=false'.format(ra,dec,size)
        urls.append(data_url + date_url + imname + condis)

    return urls


# Function for actually downloading an image
# using a URL generated by the get_ref_urls
# or get_sci_urls functions
def download_ztf_image(url):
    """
        url = URL used to download image data

    Returns:
        image  = Image data
        header = Image header
        wcs_solution = Image WCS solution
    """
    r = requests.get(url)
    hdu = fits.open(io.BytesIO(r.content))
    image = hdu[0].data
    header = hdu[0].header
    hdu.close()
    wcs_solution = wcs.WCS(header)
    return image,header,wcs_solution


# Function which queries IPAC for ZTF light curve
# data products using the API's cone search tool
def lightcurve_query(ra,dec,cone_radius=3.0):
    """
        ra  = RA coordinate in decimal degrees
        dec = Dec coordinate in decimal degrees

        Optional:
        cone_radius = search radius in arcseconds (default = 3)

    Returns:
        ztf_data = ZTF light curve data in a Pandas DataFrame
    """

    # Generate the light curve search URL
    cone_radius /= 3600. # Convert search radius to degrees
    ztf_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?"
    pos_string = "POS=CIRCLE {:.5f} {:.5f} {:.5f}&COLLECTION=ztf_dr3&FORMAT=csv".format(ra,dec,cone_radius)
    ztf_search_url = ztf_url + pos_string

    # Perform request and decode response into a Pandas DataFrame
    response = requests.get(ztf_search_url)
    ztf_data = pd.read_csv(StringIO(response.content.decode('utf-8')))

    return ztf_data




# Function which queries the LCO MARS alert broker
# for ZTF transient alert packets using the cone
# search functionality at given RA,DEC coords
def alert_query(RA,DEC,save=False,cone_radius=5.0,sleep_time=60.0,query_limit=1000,rb_limit=0.0):
    """
        RA  = list/array of RA coordinates in decimal degrees
        DEC = list/array of Dec coordinates in decimal degrees

        Optional:
        save        = whether to save query results to file in
                      JSON format (default = False)
        cone_radius = search radius in arcseconds (default = 5)
        sleep_time  = time in seconds between successive MARS
                      queries (default = 60)
        query_limit = number of coordinates to search per POST
                      request query (default = 1000)
        rb_limit    = Lower limit on the alert Real-Bogus score
                      (default = 0, i.e. no limit)
    """

    # Function which retrieves summary statistics
    # for an alert query
    def get_query_info(q):
        """
            q = the JSON formatted query results

        Returns:
            the total number of queries
            the number of queries with results 
            the number of objects with alerts
            the total number of alerts
        """
        objs = []
        num_q = 0
        num_a = 0
        for qi in q:
            if qi['num_alerts'] == 0:
                continue
            else:
                num_q += 1
            for qk in qi['results']:
                objs.append(qk['objectId'])
                num_a += 1
        objs_uniq = len(list(set(objs)))
        return len(q),num_q,objs_uniq,num_a


    # POST Request function formatted for MARS alert queries
    def post_request(q):
        """
            q = list of query dicts for each RA-DEC pair

        Returns:
            post = The POST request result in JSON format
        """
        base_url = 'https://mars.lco.global/'
        payload = {'queries':q}
        qcount = 0
        while True:
            try:
                post = requests.post(base_url,json=payload,timeout=50.0)
                break
            except Exception as e:
                time.sleep(0.2)
                qcount += 1
                print('   Timeout Error {}, querying again.'.format(qcount))
        return post


    # Perform check on input RA and DEC values
    Nobj = len(RA)
    if Nobj != len(DEC):
        print('ERROR: Different number of RA and DEC coordinates.')
        sys.exit(1)
    print('\nTotal Number of Queries to Perform: {}'.format(Nobj))

    # Setup a request template using LCO MARS Request Filters
    request_defs = {   'cone':'',     # Cone search values, in "ra,dec,radius" format
                    'rb__gte':'',     # Limit on Real-Bogus values
                     'format':'json'}     

    # Use RA-DEC values to generate a POST payload for the query
    print('Generating query payload......')
    queries = []
    cone_radius /= 3600 # Convert radiusfrom arcseconds to degrees
    for i in range(Nobj):
        cone_str = '{:.6f},{:.6f},{:.6f}'.format(RA[i],DEC[i],cone_radius)
        new_que = copy.deepcopy(request_defs)
        new_que['id'] = str(i)
        new_que['cone'] = cone_str
        new_que['rb__gte'] = '{:.2f}'.format(rb_limit)
        new_que['query_date'] = str(Time.now().jd)
        queries.append(new_que)

    # Perform the new query
    print('')
    numq = int(np.ceil(float(Nobj)/float(query_limit)))
    pj = []
    for i in range(numq):
        low = i*query_limit
        if i < numq-1:
            upp = low+query_limit
        else:
            upp = Nobj
        print('Querying Objects %i through %i' %(low+1,upp))
        p = post_request(queries[low:upp])

        # Add query results to full results list
        pj += p.json()['results']
        if i < numq-1:
            sleep_time = args.sleeptime
            if sleep_time > 0:
                print('Sleeping {:.2f} Minute(s)...'.format(float(sleep_time)/60.))
                time.sleep(sleep_time)

    # Get info for new query and print results
    qs,numq,num_objs,num_alerts = get_query_info(pj)
    print('\nQuery Results: {:5d} Queries, {:4d} Queries with Alerts, {:4d} Objects with Alerts, {:4d} Total Alerts'.format(
           qs,numq,num_objs,num_alerts))

    # Save combined results
    if save:
        time_str = Time.now().iso[0:10].replace("-","")
        json_fname = 'alert_results_{}.json'.format(time_str)
        if os.path.exists(json_fname):
            exists = True
            ext = 1
            while exists:
                json_fname = 'alert_results_{}_{}.json'.format(time_str,str(ext))
                if os.path.exists(json_fname):
                    ext += 1
                else:
                    exists = False
        with open(json_fname, 'w') as outfile:
            json.dump(pj, outfile)

        print('\nResults saved in Current Folder as {}\n'.format(json_fname))

    return pj
