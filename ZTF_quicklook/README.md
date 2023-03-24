## ZTF Light Curve Quick-Look Analysis

The python notebook, **ztf_quicklook**, allows you to quickly download and take a first look at ZTF light curves.

The notebook makes use of [ZTF's Light Curve API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html) to first search for objects within a cone search centered on given RA-Dec coordinates.  The default search radius used (3 arcseconds) is small to reduce the possibility of multiple sources falling within the search area, but the resulting light curve(s) may still belong to multiple objects if searching within a crowded region, so smaller cone search radii may be desired.

Following the retreival of light curves, the data is parsed into separate *g*-,*r*-, and *i*-band light curves and filtered for poor quality detections following the standard recommendations in the [ZTF data release documentation](https://irsa.ipac.caltech.edu/data/ZTF/docs/releases/ztf_release_notes_latest).

Periodograms are generated for the light curves and the highest peak's frequency is optimized using [LMFIT](https://lmfit.github.io/lmfit-py/) and then used to generate phase folded light curves.  Often the highest peak selected is consistent with being noise, but for objects with easily detectable periodic variability (such as the compact binary from [Ren et al. 2023](https://ui.adsabs.harvard.edu/abs/2023ApJS..264...39R/abstract) shown in the example notebook), the chosen peak can be significant.
