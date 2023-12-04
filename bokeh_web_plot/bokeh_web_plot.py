import numpy as np
import pandas as pd
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table

from bokeh.layouts import row, column, grid, layout
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.models import Range1d, Whisker, Legend, LinearAxis
from bokeh.models.glyphs import Text
from bokeh.models.widgets import Toggle
from bokeh.plotting import figure, output_file, show

""" Configure the output file """
fname = "my_bokeh_plot.html"
output_file(fname)

base_path = "./"

""" Pan-STARRS1 Magnitudes of Target """

PSmagg = 18.482
PSmagr = 18.492
PSmagi = 18.587
PSmagg_err = 0.013
PSmagr_err = 0.011
PSmagi_err = 0.012


""" Basic plot appearance params """

mjdref = 58000 # Reference MJD epoch

# Fontsizes
axis_label_fs = '10pt' # Axis labels
tick_label_fs = '9pt'  # Tick labels

# Some custom light and dark colors
lcolors = ['#97bfe0','#ffbf80','#78de78',
           '#fbc0df','#e1a685','#cda4d3',
           '#cccccc','#ed9582','#ffff6f',
           '#909090','#fbca76']
dcolors = ['#377eb8','#ff7f00','#228b22',
           '#f781bf','#a65628','#984ea3',
           '#999999','#c43a1c','#dede00',
           '#222222','#dd8f07']

# Marker colors
c_ztfg = lcolors[0]
c_ztfr = dcolors[0]
c_lcor = 'firebrick'
c_lcog = 'sandybrown'

# marker alpha values
a_ztfg = 0.75
a_ztfr = 0.75
a_lcog = 0.75
a_lcor = 0.75



""" Load in light curve data """

# Load in ZTF Data
ztf_path = base_path + "lc_dr19.csv"
ztf_data = pd.read_csv(ztf_path)
ztf_gdata = ztf_data[
    (ztf_data.filtercode == 'zg') & 
    (ztf_data.catflags == 0) &
    (ztf_data.sharp > -0.5) & 
    (ztf_data.sharp < 0.5)
].sort_values(by='mjd').reset_index(drop=True)
ztf_rdata = ztf_data[
    (ztf_data.filtercode == 'zr') & 
    (ztf_data.catflags == 0) &
    (ztf_data.sharp > -0.5) & 
    (ztf_data.sharp < 0.5)
].sort_values(by='mjd').reset_index(drop=True)

# Apply color corrections to ZTF data
ztf_gdata['mag'] = ztf_gdata.mag + ztf_gdata.clrcoeff*(PSmagg-PSmagr)
ztf_rdata['mag'] = ztf_rdata.mag + ztf_rdata.clrcoeff*(PSmagg-PSmagr)

# Load in LCO Data
lco_path = base_path + "LCO_Phot_File_binned.csv"
lco_data = pd.read_csv(lco_path)
lco_gdata = lco_data[
    (lco_data['filter'] == 'gp') & 
    (lco_data['sn'] > 15.0) & 
    (lco_data['magerr'] < 0.07) &
    (
        (lco_data['site_id'] == 'elp') |
        (lco_data['site_id'] == 'tfn')
    )
].copy(deep=True).reset_index(drop=True)
lco_rdata = lco_data[
    (lco_data['filter'] == 'rp') & 
    (lco_data['sn'] > 15.0) & 
    (lco_data['magerr'] < 0.07) &
    (
        (lco_data['site_id'] == 'elp') |
        (lco_data['site_id'] == 'tfn')
    )
].copy(deep=True).reset_index(drop=True)



""" Get initial x and y-limits for the light curve plot """

mjd_min = ztf_data.mjd.min()
mjd_max = max([lco_data.mjd.max(),ztf_data.mjd.max()])
xlow = mjd_min - 0.025*(mjd_max-mjd_min) - mjdref
xupp = mjd_max + 0.025*(mjd_max-mjd_min) - mjdref
ydiff = max(ztf_gdata.mag.max(),ztf_rdata.mag.max()) - min(min(ztf_gdata.mag),min(ztf_rdata.mag))
ylow = min(
    [
        ztf_gdata.mag.min()-ztf_gdata.mag.median(),
        ztf_rdata.mag.min()-ztf_rdata.mag.median()
    ]) - 0.15*ydiff
yupp = max(
    [
        ztf_gdata.mag.max()-ztf_gdata.mag.median(),
        ztf_rdata.mag.max()-ztf_rdata.mag.median()
    ]) + 0.15*ydiff


""" Create Column Data Source Objects """

# Pull relevant data for Column Data Sources from loaded data
def get_source_data(d):
    x = d.mjd.values - mjdref
    y = d.mag.values - d.mag.median()
    low = y - d.magerr.values
    upp = y + d.magerr.values
    return x,y,low,upp

ztf_gmjd,ztf_gmag,ztf_glow,ztf_gupp = get_source_data(ztf_gdata)
ztf_rmjd,ztf_rmag,ztf_rlow,ztf_rupp = get_source_data(ztf_rdata)
lco_gmjd,lco_gmag,lco_glow,lco_gupp = get_source_data(lco_gdata)
lco_rmjd,lco_rmag,lco_rlow,lco_rupp = get_source_data(lco_rdata)
 
# These CDS's get updated with button clicks.
source_ztf_g = ColumnDataSource(
    data=dict(
        x=ztf_gmjd, 
        y=ztf_gmag, 
        lower=ztf_glow, 
        upper=ztf_gupp
    )
)
source_ztf_r = ColumnDataSource(
    data=dict(
        x=ztf_rmjd, 
        y=ztf_rmag, 
        lower=ztf_rlow, 
        upper=ztf_rupp
    )
)
source_lco_g = ColumnDataSource(
    data=dict(
        x=lco_gmjd, 
        y=lco_gmag, 
        lower=lco_glow, 
        upper=lco_gupp
    )
)
source_lco_r = ColumnDataSource(
    data=dict(
        x=lco_rmjd, 
        y=lco_rmag, 
        lower=lco_rlow, 
        upper=lco_rupp
    )
)

# These CDS's are never changed but used to update 
# the CDS's above during button click events.
orig_ztf_g = ColumnDataSource(
    data=dict(
        x=ztf_gmjd, 
        y=ztf_gmag,
        lower=ztf_glow, 
        upper=ztf_gupp
    )
)
orig_ztf_r = ColumnDataSource(
    data=dict(
        x=ztf_rmjd, 
        y=ztf_rmag,
        lower=ztf_rlow, 
        upper=ztf_rupp
    )
)
orig_lco_g = ColumnDataSource(
    data=dict(
        x=lco_gmjd, 
        y=lco_gmag,
        lower=lco_glow, 
        upper=lco_gupp
    )
)
orig_lco_r = ColumnDataSource(
    data=dict(
        x=lco_rmjd, 
        y=lco_rmag,
        lower=lco_rlow, 
        upper=lco_rupp
    )
)


""" Initialize Light Curve Plot """

fig_lc = figure(
    height=380, 
    width=750, 
    sizing_mode='scale_both',
    x_range=[xlow,xupp], 
    y_range=[yupp,ylow],
    tools="pan,wheel_zoom,box_zoom,reset,fullscreen",
    toolbar_location="above", 
    border_fill_color="whitesmoke"
)
fig_lc.toolbar.logo = None
fig_lc.yaxis.axis_label = r'$$\Delta\mathrm{mag}$$'
fig_lc.xaxis.axis_label = r'$$\mathrm{MJD}-%.0f$$' %mjdref
fig_lc.xaxis.axis_label_text_font_size = axis_label_fs
fig_lc.yaxis.axis_label_text_font_size = axis_label_fs
fig_lc.xaxis.major_label_text_font_size = tick_label_fs
fig_lc.yaxis.major_label_text_font_size = tick_label_fs

# Add a second x-axis with fractional year coordinates
astro_mjds = Time([xlow+mjdref,xupp+mjdref],scale='utc',format='mjd')
decimal_years = astro_mjds.decimalyear
fig_lc.extra_x_ranges = {
    "fracyear": Range1d(
        start=decimal_years[0], 
        end=decimal_years[1]
    )
}
fig_lc.add_layout(LinearAxis(x_range_name="fracyear"), 'above')
fig_lc.xaxis[0].axis_label = r'$$\mathrm{Decimal\;Year}$$'
fig_lc.xaxis[0].axis_label_text_font_size = axis_label_fs
fig_lc.xaxis[0].major_label_text_font_size = tick_label_fs

# Add a second y-axis with relative flux units
tick_rellocations = [-50.,-40.,-30.,-20.,-10.,0.,10.,20.,30.]
tick_maglocations = [-2.5*np.log10(x/100+1.0) for x in tick_rellocations]
tick_labels = {}
for key,val in zip(tick_maglocations,tick_rellocations):
    tick_labels[key] = str(int(val))
fig_lc.extra_y_ranges = {"relflux": Range1d(start=yupp, end=ylow)}
fig_lc.add_layout(LinearAxis(y_range_name="relflux"), 'right')
fig_lc.yaxis[1].ticker = tick_maglocations
fig_lc.yaxis[1].major_label_overrides = tick_labels
fig_lc.yaxis[1].axis_label = r'$$\mathrm{Relative\;Flux\;Change\;(\%)}$$'
fig_lc.yaxis[1].axis_label_text_font_size = axis_label_fs
fig_lc.yaxis[1].major_label_text_font_size = tick_label_fs


# Add Error Bars
ztf_eg = Whisker(source=source_ztf_g, base="x", upper="upper", lower="lower",
                 line_color=c_ztfg,line_width=1.3,line_alpha=0.5)
ztf_er = Whisker(source=source_ztf_r, base="x", upper="upper", lower="lower",
                 line_color=c_ztfr,line_width=1.3,line_alpha=0.5)
lco_eg = Whisker(source=source_lco_g, base="x", upper="upper", lower="lower",
                 line_color=c_lcog,line_width=1.3,line_alpha=0.5)
lco_er = Whisker(source=source_lco_r, base="x", upper="upper", lower="lower",
                 line_color=c_lcor,line_width=1.3,line_alpha=0.5)
ztf_eg.upper_head.line_color = None
ztf_eg.lower_head.line_color = None
ztf_er.upper_head.line_color = None
ztf_er.lower_head.line_color = None
lco_eg.upper_head.line_color = None
lco_eg.lower_head.line_color = None
lco_er.upper_head.line_color = None
lco_er.lower_head.line_color = None
fig_lc.add_layout(ztf_eg)
fig_lc.add_layout(ztf_er)
fig_lc.add_layout(lco_eg)
fig_lc.add_layout(lco_er)

# Plot ZTF g and r light curve data
f_ztfg = fig_lc.circle('x', 'y', source=source_ztf_g, size=6, 
               # Set defaults
               fill_color=c_ztfg, line_color=c_ztfg,
               fill_alpha=a_ztfg, line_alpha=a_ztfg,
               # set visual properties for selected glyphs
               selection_color=c_ztfg,
               selection_alpha=a_ztfg,
               # set visual properties for non-selected glyphs
               nonselection_color=c_ztfg,
               nonselection_alpha=a_ztfg,
               legend_label = 'ZTF-g')
f_ztfr = fig_lc.square('x', 'y', source=source_ztf_r, size=6, 
               # Set defaults
               fill_color=c_ztfr, line_color=c_ztfr,
               fill_alpha=a_ztfr, line_alpha=a_ztfr,                    
               # set visual properties for selected glyphs
               selection_color=c_ztfr,
               selection_alpha=a_ztfr,
               # set visual properties for non-selected glyphs
               nonselection_color=c_ztfr,
               nonselection_alpha=a_ztfr,
               legend_label = 'ZTF-r')


# Plot LCO g and r light curve data
f_lcog = fig_lc.triangle('x', 'y', source=source_lco_g, size=7, 
               # Set defaults
               fill_color=c_lcog, line_color=c_lcog,
               fill_alpha=a_lcog, line_alpha=a_lcog,
               # set visual properties for selected glyphs
               selection_color=c_lcog,
               selection_alpha=a_lcog,
               # set visual properties for non-selected glyphs
               nonselection_color=c_lcog,
               nonselection_alpha=a_lcog,
               legend_label = 'LCOGT-g')
f_lcor = fig_lc.diamond('x', 'y', source=source_lco_r, size=9,
               # Set defaults
               fill_color=c_lcor, line_color=c_lcor,
               fill_alpha=a_lcor, line_alpha=a_lcor,                    
               # set visual properties for selected glyphs
               selection_color=c_lcor,
               selection_alpha=a_lcor,
               # set visual properties for non-selected glyphs
               nonselection_color=c_lcor,
               nonselection_alpha=a_lcor,
               legend_label = 'LCOGT-r')

# Modify the legend's appearance
fig_lc.legend.location = "bottom_right"
fig_lc.legend.orientation = "horizontal"
fig_lc.legend.border_line_width = 3
fig_lc.legend.border_line_color = "navy"
fig_lc.legend.border_line_alpha = 0.5
fig_lc.legend.label_text_font_style = "bold"
fig_lc.legend.label_text_font_size = '9pt'
fig_lc.legend.spacing = 7
fig_lc.legend.label_standoff = -2
fig_lc.legend.padding = 5


# Add Toggle Buttons
b1 = Toggle(label="ZTF"   ,button_type="primary",active=True,height=90,aspect_ratio=6)
b2 = Toggle(label="LCOGT" ,button_type="primary",active=True,height=90,aspect_ratio=6)
b3 = Toggle(label="g-data",button_type="primary",active=True,height=90,aspect_ratio=6)
b4 = Toggle(label="r-data",button_type="primary",active=True,height=90,aspect_ratio=6)


""" JS Callback for Buttons. Toggles the display of various data sets. """

callback = CustomJS(
    args=dict(
        b1=b1, b2=b2, b3=b3, b4=b4,
        p1=f_ztfg, p2=f_ztfr, p3=f_lcog, p4=f_lcor,
        w1=ztf_eg, w2=ztf_er, w3=lco_eg, w4=lco_er),
        code="""
            // Determine which buttons are active and change button color/type
            var blist = [b1,b2,b3,b4]
            var active = []
            for (var i = 0; i < blist.length; i++) {
                if (blist[i].active == true) {
                    active.push(i)
                    blist[i].button_type = "primary"
                } else {
                    blist[i].button_type = "default"
                }
            }

            // All the cases for different data selections
            // based on which buttons are active/inactive
            if (active.length < 2) {
                p1.visible = false
                p2.visible = false
                p3.visible = false
                p4.visible = false
                w1.visible = false
                w2.visible = false
                w3.visible = false
                w4.visible = false
            } else if (active.length == 2) {
                if (active[0] == 0 & active[1] == 2) {
                    p1.visible = true
                    p2.visible = false
                    p3.visible = false
                    p4.visible = false
                    w1.visible = true
                    w2.visible = false
                    w3.visible = false
                    w4.visible = false
                } else if (active[0] == 0 & active[1] == 3) {
                    p1.visible = false
                    p2.visible = true
                    p3.visible = false
                    p4.visible = false
                    w1.visible = false
                    w2.visible = true
                    w3.visible = false
                    w4.visible = false
                } else if (active[0] == 1 & active[1] == 2) {
                    p1.visible = false
                    p2.visible = false
                    p3.visible = true
                    p4.visible = false
                    w1.visible = false
                    w2.visible = false
                    w3.visible = true
                    w4.visible = false
                } else if (active[0] == 1 & active[1] == 3) {
                    p1.visible = false
                    p2.visible = false
                    p3.visible = false
                    p4.visible = true
                    w1.visible = false
                    w2.visible = false
                    w3.visible = false
                    w4.visible = true
                } else {
                    p1.visible = false
                    p2.visible = false
                    p3.visible = false
                    p4.visible = false
                    w1.visible = false
                    w2.visible = false
                    w3.visible = false
                    w4.visible = false
                }
            } else if (active.length == 3) {
                if (active[0] == 0 & active[1] == 1 & active[2] == 2) {
                    p1.visible = true
                    p2.visible = false
                    p3.visible = true
                    p4.visible = false
                    w1.visible = true
                    w2.visible = false
                    w3.visible = true
                    w4.visible = false
                } else if (active[0] == 0 & active[1] == 1 & active[2] == 3) {
                    p1.visible = false
                    p2.visible = true
                    p3.visible = false
                    p4.visible = true
                    w1.visible = false
                    w2.visible = true
                    w3.visible = false
                    w4.visible = true
                } else if (active[0] == 0 & active[1] == 2 & active[2] == 3) {
                    p1.visible = true
                    p2.visible = true
                    p3.visible = false
                    p4.visible = false
                    w1.visible = true
                    w2.visible = true
                    w3.visible = false
                    w4.visible = false
                } else if (active[0] == 1 & active[1] == 2 & active[2] == 3) {
                    p1.visible = false
                    p2.visible = false
                    p3.visible = true
                    p4.visible = true
                    w1.visible = false
                    w2.visible = false
                    w3.visible = true
                    w4.visible = true
                }
            } else {
                p1.visible = true
                p2.visible = true
                p3.visible = true
                p4.visible = true
                w1.visible = true
                w2.visible = true
                w3.visible = true
                w4.visible = true
            }

        """
)


# Connect buttons to callback function
b1.js_on_click(callback)
b2.js_on_click(callback)
b3.js_on_click(callback)
b4.js_on_click(callback)

# Create plot grid
l = layout(
    [
        [b1,b2,b3,b4],
        [fig_lc]
    ],
    sizing_mode='scale_width')

# Use "show" to generate HTML file and launch plot in the browser
show(l)

# Bokeh generates an annoying HTML header (<!DOCTYPE html>") 
# which displays as raw text when embedded in the website.  
# The following code removes this line of text from the html file
with open(fname, "r") as f:
    lines = f.readlines()
with open(fname, "w") as f:
    for line in lines:
        if line.strip("\n") != "<!DOCTYPE html>":
            f.write(line)





