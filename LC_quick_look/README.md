## ZTF Light Curve Quick-Look Analysis

The python notebook, **view_ztf_data**, allows you to quickly download and take a first look at ZTF light curves.

The notebook makes use of [ZTF's Light Curve API](https://irsa.ipac.caltech.edu/docs/program_interface/ztf_lightcurve_api.html) to first search for objects within a cone search centered on given RA-Dec coordinates.  The search radius used (0.001 degree) is small to reduce the possibility of multiple sources falling within the search area, but the resulting light curve(s) may still belong to multiple objects if searching within a crowded region.

Following the retreival of light curves, the data is parsed into separate *g* and *r* light curves and filtered for poor quality detections following the relatively conservative recommendations of [Guidry et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv201200035G/abstract).

Periodograms are generated for the light curves and the highest peak's frequency is optimized using [LMFIT](https://lmfit.github.io/lmfit-py/) and then used to generate a phase folded light curve.  Often the highest peak selected is consistent with being noise, but for objects with easily detectable periodic variability (such as the RR Lyrae shown in the example notebook), the chosen peak can be significant.
