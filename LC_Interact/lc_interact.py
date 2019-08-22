import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from astropy.io import fits
from astropy.io import ascii
from astropy.visualization import ZScaleInterval
from astropy.time import Time
from astropy.coordinates import SkyCoord

from bokeh.io import curdoc
from bokeh.layouts import column, row, layout
from bokeh.models import Slider, ColumnDataSource, Button, CustomJS
from bokeh.models import Span, Range1d, LinearColorMapper, Whisker
from bokeh.models.glyphs import Text
from bokeh.plotting import figure


# Load in the ZTF data and convert to Pandas DataFrame
data = ascii.read('lc.txt').to_pandas()

# Load in ZTF Images at this object's RA/Dec
ZS = ZScaleInterval(nsamples=10000, contrast=0.15, max_reject=0.5, 
                    min_npixels=5, krej=2.5, max_iterations=5)
fits_files = sorted(glob('ZTF_Sci_Files/*.fits'))
num_files = len(fits_files)
mjds_g,mjds_r = [],[]   # Store MJD of image
imdat_g,imdat_r = [],[] # Store Image pixel data
vmin_g,vmin_r = [],[]   # Store Z-Scale Minimums
vmax_g,vmax_r = [],[]   # Store Z-Scale Maximums
hdrs_g,hdrs_r = [],[]   # Store Image headers
for f in fits_files:
    hdr = fits.getheader(f)
    dat = fits.getdata(f)
    vmin,vmax = ZS.get_limits(dat)
    if hdr['FILTERID'] == 1:
        mjds_g.append(hdr['OBSMJD'])
        imdat_g.append(dat)
        hdrs_g.append(hdr)
        vmin_g.append(vmin)
        vmax_g.append(vmax)
    elif hdr['FILTERID'] == 2:
        mjds_r.append(hdr['OBSMJD'])
        imdat_r.append(dat)
        hdrs_r.append(hdr)
        vmin_r.append(vmin)
        vmax_r.append(vmax)

# Separate Data into g and r filters and remove poor quality epochs with CATFLAGS bit value
gdata = data[(data.filtercode == 'zg') & (data.catflags == 0)].sort_values(by='mjd').reset_index()
rdata = data[(data.filtercode == 'zr') & (data.catflags == 0)].sort_values(by='mjd').reset_index()

# Find the Image MJDs which have no associated data points in the light curves
gmjds_nodat = []
rmjds_nodat = []
for i,d in enumerate(mjds_g):
    if min(abs(gdata.mjd.values - d))*86400.0 > 1.0:
        gmjds_nodat.append(d)
for d in mjds_r:
    if min(abs(rdata.mjd.values - d))*86400.0 > 1.0:
        rmjds_nodat.append(d)

print(len(gdata), len(rdata))
print(len(gmjds_nodat), len(rmjds_nodat))

# Get data of interest from the loaded lightcurve tables
gmjds = gdata.mjd.values # MJD dates for g-data
rmjds = rdata.mjd.values # MJD dates for r-data
gmags = gdata.mag.values # Magnitudes for g-data
rmags = rdata.mag.values # Magnitudes for r-data
glower = gmags-gdata.magerr.values # Lower Magnitude Errors for g-data
gupper = gmags+gdata.magerr.values # Upper Magnitude Errors for g-data
rlower = rmags-rdata.magerr.values # Lower Magnitude Errors for r-data
rupper = rmags+rdata.magerr.values # Upper Magnitude Errors for r-data

# Get initial ylimits for the light curve plot
ydiff = max(max(gmags),max(rmags)) - min(min(gmags),min(rmags))
ylow = min(min(gmags),min(rmags)) - 0.3*ydiff
yupp = max(max(gmags),max(rmags)) + 0.3*ydiff

# Generate fixed magnitude arrays for images without any light curve data
gmags_nodat = [ylow+0.10*ydiff]*len(gmjds_nodat)
rmags_nodat = [ylow+0.10*ydiff]*len(rmjds_nodat)

# Create CDS objects for light curves
source_g = ColumnDataSource(data=dict(x=gmjds, y=gmags, lower=glower, upper=gupper))
source_r = ColumnDataSource(data=dict(x=rmjds, y=rmags, lower=rlower, upper=rupper))
source_nodat_g = ColumnDataSource(data=dict(x=gmjds_nodat, y=gmags_nodat))
source_nodat_r = ColumnDataSource(data=dict(x=rmjds_nodat, y=rmags_nodat))

# Generate an object name using the RA/Dec in the Light Curve File
ra = data['ra'].values[0]
de = data['dec'].values[0]
coord = SkyCoord(ra,de,unit="deg",frame="icrs")
coord_str = coord.to_string('hmsdms',sep="",precision=2).replace(" ","")
ztf_name = 'ZTF J{}'.format(coord_str)

# Initialize Light Curve Plot
fig_lc = figure(plot_height=320, plot_width=650,
                tools="pan,wheel_zoom,box_zoom,tap,reset",
                toolbar_location="above", border_fill_color="whitesmoke")
fig_lc.toolbar.logo = None
fig_lc.title.text = "Light Curve for {}".format(ztf_name)
fig_lc.title.offset = -10
fig_lc.yaxis.axis_label = 'Mag'
fig_lc.xaxis.axis_label = 'MJD'
fig_lc.y_range = Range1d(start=yupp, end=ylow)

# Add Error Bars
eg = Whisker(source=source_g, base="x", upper="upper", lower="lower",
             line_color="forestgreen",line_width=1.5,line_alpha=0.5)
er = Whisker(source=source_r, base="x", upper="upper", lower="lower",
             line_color="firebrick",line_width=1.5,line_alpha=0.5)
eg.upper_head.line_color = None
eg.lower_head.line_color = None
er.upper_head.line_color = None
er.lower_head.line_color = None
fig_lc.add_layout(eg)
fig_lc.add_layout(er)

# Plot g and r light curve data
fig_lc.circle('x', 'y', source=source_g, size=6, 
               # Set defaults
               fill_color="forestgreen", line_color="forestgreen",
               fill_alpha=0.5, line_alpha=0.5,
               # set visual properties for selected glyphs
               selection_color="forestgreen",
               selection_alpha=0.5,
               # set visual properties for non-selected glyphs
               nonselection_color="forestgreen",
               nonselection_alpha=0.5)
fig_lc.circle('x', 'y', source=source_r, size=6, 
               # Set defaults
               fill_color="firebrick", line_color="firebrick",
               fill_alpha=0.5, line_alpha=0.5,                    
               # set visual properties for selected glyphs
               selection_color="firebrick",
               selection_alpha=0.5,
               # set visual properties for non-selected glyphs
               nonselection_color="firebrick",
               nonselection_alpha=0.5)

# Plot diamonds for images without any corresponding light curve data
fig_lc.diamond('x', 'y', source=source_nodat_g, size=8,
               # Set defaults
               fill_color="forestgreen", line_color="forestgreen",
               fill_alpha=0.5, line_alpha=0.5,
               # set visual properties for selected glyphs
               selection_color="forestgreen",
               selection_alpha=0.5,
               # set visual properties for non-selected glyphs
               nonselection_color="forestgreen",
               nonselection_alpha=0.5)
fig_lc.diamond('x', 'y', source=source_nodat_r, size=8,
               # Set defaults
               fill_color="firebrick", line_color="firebrick",
               fill_alpha=0.5, line_alpha=0.5,                    
               # set visual properties for selected glyphs
               selection_color="firebrick",
               selection_alpha=0.5,
               # set visual properties for non-selected glyphs
               nonselection_color="firebrick",
               nonselection_alpha=0.5)

# Vertical line to indicate the cadence
vline_g = Span(location=mjds_g[0], dimension='height', 
               line_color='forestgreen',line_width=2,line_alpha=0.8)
vline_r = Span(location=mjds_r[0], dimension='height', 
               line_color='firebrick'  ,line_width=2,line_alpha=0.8)
fig_lc.add_layout(vline_g)
fig_lc.add_layout(vline_r)

# Initialize Image plots
fig_img = figure(plot_width=300, plot_height=320,
                 x_range=[0, imdat_g[0].shape[1]], 
                 y_range=[0, imdat_g[0].shape[0]],
                 title="ZTF-g Image", tools='box_zoom,wheel_zoom,reset',
                 toolbar_location="above",
                 border_fill_color="whitesmoke")
fig_imr = figure(plot_width=300, plot_height=320,
                 x_range=fig_img.x_range, 
                 y_range=fig_img.y_range,
                 title="ZTF-r Image", tools='box_zoom,wheel_zoom,reset',
                 toolbar_location="above",
                 border_fill_color="whitesmoke")
fig_img.toolbar.logo = None
fig_imr.toolbar.logo = None

# Create Z-Scale Normalized Color Maps
color_mapper_g = LinearColorMapper(palette="Greys256", low=vmin_g[0], high=vmax_g[0])
color_mapper_r = LinearColorMapper(palette="Greys256", low=vmin_r[0], high=vmax_r[0])
# Plot the images
fig_img.image(image=[imdat_g[0]], x=0, y=0, dw=imdat_g[0].shape[1], dh=imdat_g[0].shape[0], 
            dilate=True, color_mapper=color_mapper_g, name="gframe")
fig_imr.image(image=[imdat_r[0]], x=0, y=0, dw=imdat_r[0].shape[1], dh=imdat_r[0].shape[0], 
            dilate=True, color_mapper=color_mapper_r, name="rframe")


# Interactive slider widgets and buttons to select the image number
g_frame_slider = Slider(start=1,end=len(imdat_g),
                        value=1,step=1,bar_color='forestgreen',
                        title="g-Frame Slider",width=520,
                        callback_policy="throttle",
                        callback_throttle=200)
r_frame_slider = Slider(start=1,end=len(imdat_r),
                        value=1,step=1,bar_color='firebrick',
                        title="r-Frame Slider",width=520,
                        callback_policy="throttle",
                        callback_throttle=200)
rbutton_g = Button(label=">", button_type="default", width=40)
lbutton_g = Button(label="<", button_type="default", width=40)
rbutton_r = Button(label=">", button_type="default", width=40)
lbutton_r = Button(label="<", button_type="default", width=40)


# Initialize the Info Box Plots
fig_infog = figure(plot_width=325, plot_height=222,
                   x_range=[0, 1], y_range=[0, 1],
                   title="ZTF-g Image MetaData", tools='',
                   border_fill_color="whitesmoke")
fig_infor = figure(plot_width=325, plot_height=222,
                   x_range=[0, 1], y_range=[0, 1],
                   title="ZTF-r Image MetaData", tools='',
                   border_fill_color="whitesmoke")
fig_infog.toolbar.logo = None
fig_infor.toolbar.logo = None

# Configure Info Box appearance
fig_infog.xgrid.visible = False # Remove x grid
fig_infog.ygrid.visible = False # Remove y grid
fig_infor.xgrid.visible = False # Remove x grid
fig_infor.ygrid.visible = False # Remove y grid
fig_infog.xaxis.visible = False # Remove x-axis
fig_infog.yaxis.visible = False # Remove y-axis
fig_infor.xaxis.visible = False # Remove x-axis
fig_infor.yaxis.visible = False # Remove y-axis
fig_infog.title.align = "center" # Put title in center
fig_infor.title.align = "center" # Put title in center

# Function to generate text output
text_formats = ['Date-Time : {}',
                '      MJD : {:.6f}',
                '  Airmass : {:.3f}',
                '   Seeing : {:.3f}"',
                'Mag Limit : {:.2f}',
                'Moon Frac : {:.3f}',
                '     Temp : {:.1f} C',
                '     Wind : {:.1f} mph',
                ' Humidity : {:.0f} %',
                ' Infobits : {:.0f}']
hdr_keys = ['OBSMJD','OBSMJD','AIRMASS','SEEING','MAGLIM',
            'MOONILLF','TEMPTURE','WINDSPD','HUMIDITY','INFOBITS']
def gen_text(h):
    text_output = []
    for i in range(len(text_formats)):
        if i == 0:
            t = Time(h[hdr_keys[i]],scale='utc',format='mjd')
            text = text_formats[i].format(t.iso)
        else:
            text = text_formats[i].format(abs(h[hdr_keys[i]]))
        text_output.append(text)
    return text_output

# Generate text Glyphs
xlocs = [0.03]*10  # x-locations for each line of text
ylocs = [0.89 - float(i)*0.09 for i in range(10)] # y-locations for each line of text
text_init_g = gen_text(hdrs_g[0])  # Initial g-frame info to display
text_init_r = gen_text(hdrs_r[0])  # Initial r-frame info to display
text_source_g = ColumnDataSource(data=dict(xloc=xlocs,yloc=ylocs,text=text_init_g))
text_source_r = ColumnDataSource(data=dict(xloc=xlocs,yloc=ylocs,text=text_init_r))
textglyph_g = Text(x="xloc", y="yloc", text="text", 
                   text_font='courier', text_font_size='9pt')
textglyph_r = Text(x="xloc", y="yloc", text="text", 
                   text_font='courier', text_font_size='9pt')
fig_infog.add_glyph(text_source_g, textglyph_g)
fig_infor.add_glyph(text_source_r, textglyph_r)


# Callback function for g-frame slider
def update_g_frame(attr, old, new):
    newind = new['value'][0]
    fig_img.select('gframe')[0].data_source.data['image'] = [imdat_g[newind-1]]
    fig_img.select('gframe')[0].glyph.color_mapper.high = vmax_g[newind-1]
    fig_img.select('gframe')[0].glyph.color_mapper.low = vmin_g[newind-1]
    vline_g.update(location=mjds_g[newind-1])

    # Update text
    new_text = gen_text(hdrs_g[newind-1])
    text_source_g.data["text"] = new_text

    # Clear selections
    source_g.selected.indices = []
    source_nodat_g.selected.indices = []

# Callback function for r-frame slider
def update_r_frame(attr, old, new):
    newind = new['value'][0]
    fig_imr.select('rframe')[0].data_source.data['image'] = [imdat_r[newind-1]]
    fig_imr.select('rframe')[0].glyph.color_mapper.high = vmax_r[newind-1]
    fig_imr.select('rframe')[0].glyph.color_mapper.low = vmin_r[newind-1]
    vline_r.update(location=mjds_r[newind-1])

    # Update text
    new_text = gen_text(hdrs_r[newind-1])
    text_source_r.data["text"] = new_text

    source_r.selected.indices = []
    source_nodat_r.selected.indices = []

# Right button click event for g-frame
def go_right_by_one_gframe():
    existing_value = g_frame_slider.value
    if existing_value < len(imdat_g):
        g_frame_slider.value = existing_value + 1
        fake_source_g.data = dict(value=[existing_value + 1])

# Left button click event for g-frame
def go_left_by_one_gframe():
    existing_value = g_frame_slider.value
    if existing_value > 1:
        g_frame_slider.value = existing_value - 1
        fake_source_g.data = dict(value=[existing_value - 1])

# Right button click event for r-frame
def go_right_by_one_rframe():
    existing_value = r_frame_slider.value
    if existing_value < len(imdat_r):
        r_frame_slider.value = existing_value + 1
        fake_source_r.data = dict(value=[existing_value + 1])

# Left button click event for r-frame
def go_left_by_one_rframe():
    existing_value = r_frame_slider.value
    if existing_value > 1:
        r_frame_slider.value = existing_value - 1
        fake_source_r.data = dict(value=[existing_value - 1])

# Callback function which moves slider when a
# data point is clicked on in the Light Curve plot
def jump_to_lightcurve_position(event):
    if source_g.selected.indices != []:
        num_lower = np.count_nonzero(source_nodat_g.data['x'] < 
                                     source_g.data['x'][source_g.selected.indices[0]])
        g_frame_slider.value = source_g.selected.indices[0]+num_lower+1
        fake_source_g.data = dict(value=[source_g.selected.indices[0]+num_lower+1])
    if source_r.selected.indices != []:
        num_lower = np.count_nonzero(source_nodat_r.data['x'] < 
                                     source_r.data['x'][source_r.selected.indices[0]])
        r_frame_slider.value = source_r.selected.indices[0]+num_lower+1
        fake_source_r.data = dict(value=[source_r.selected.indices[0]+num_lower+1])
    if source_nodat_g.selected.indices != []:
        num_lower = np.count_nonzero(source_g.data['x'] < 
                                     source_nodat_g.data['x'][source_nodat_g.selected.indices[0]])
        g_frame_slider.value = source_nodat_g.selected.indices[0]+num_lower+1
        fake_source_g.data = dict(value=[source_nodat_g.selected.indices[0]+num_lower+1])
    if source_nodat_r.selected.indices != []:
        num_lower = np.count_nonzero(source_r.data['x'] < 
                                     source_nodat_r.data['x'][source_nodat_r.selected.indices[0]])
        r_frame_slider.value = source_nodat_r.selected.indices[0]+num_lower+1
        fake_source_r.data = dict(value=[source_nodat_r.selected.indices[0]+num_lower+1])

# Connect different objects/events to callback functions
rbutton_g.on_click(go_right_by_one_gframe)
lbutton_g.on_click(go_left_by_one_gframe)
rbutton_r.on_click(go_right_by_one_rframe)
lbutton_r.on_click(go_left_by_one_rframe)
fig_lc.on_event('tap',jump_to_lightcurve_position)

# 2-Step callback for the sliders to allow for callback throttling
fake_source_g = ColumnDataSource(data=dict(value=[]))
fake_source_r = ColumnDataSource(data=dict(value=[]))
fake_source_g.on_change('data', update_g_frame)
fake_source_r.on_change('data', update_r_frame)
g_frame_slider.callback = CustomJS(args=dict(source=fake_source_g), code="""
    source.data = { value: [cb_obj.value] }
    """)
r_frame_slider.callback = CustomJS(args=dict(source=fake_source_r), code="""
    source.data = { value: [cb_obj.value] }
    """)

# Create plot grid
l = layout([fig_lc,fig_img],
           [column(row(lbutton_g,rbutton_g,g_frame_slider),
                   row(lbutton_r,rbutton_r,r_frame_slider),
                   row(fig_infog,fig_infor)),fig_imr])

# Add everything into the Bokeh document
curdoc().add_root(l)
