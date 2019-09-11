import pandas as pd
from astropy.io import ascii

from bokeh.layouts import row, layout, widgetbox
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.models import Range1d, Whisker, Legend
from bokeh.models.glyphs import Text
from bokeh.models.widgets import Toggle
from bokeh.plotting import Figure, output_file, show

# Configure the output file
output_file("my_bokeh_plot.html")


# Load in ZTF Data
ztf_path = "/Users/zvander/data/object/ZTF/ztf_data/WDJ0139+5245/lc.txt"
ztf_data = ascii.read(ztf_path).to_pandas()
ztf_gdata = ztf_data[(ztf_data.filtercode == 'zg') & 
                     (ztf_data.catflags == 0)].copy(deep=True).reset_index(drop=True)
ztf_rdata = ztf_data[(ztf_data.filtercode == 'zr') & 
                     (ztf_data.catflags == 0)].copy(deep=True).reset_index(drop=True)

# Load in LCO Data
lco_path = "/Users/zvander/data/object/ZTF/ztf_data/WDJ0139+5245/LCO/LCO_Phot_File.csv"
lco_data = pd.read_csv(lco_path)
lco_gdata = lco_data[(lco_data['filter'] == 'gp') & 
                     (lco_data['sn'] > 20.0) & 
                     (lco_data['site_id'] == 'elp')].copy(deep=True).reset_index(drop=True)
lco_rdata = lco_data[(lco_data['filter'] == 'rp') & 
                     (lco_data['sn'] > 20.0) & 
                     (lco_data['site_id'] == 'elp')].copy(deep=True).reset_index(drop=True)

# Get initial ylimits for the light curve plot
ydiff = max(max(ztf_gdata.mag),max(ztf_rdata.mag)) - min(min(ztf_gdata.mag),min(ztf_rdata.mag))
ylow = min(min(ztf_gdata.mag),min(ztf_rdata.mag)) - 0.15*ydiff
yupp = max(max(ztf_gdata.mag),max(ztf_rdata.mag)) + 0.15*ydiff

# Pull relevant data for Column Data Sources from loaded data
def get_source_data(d):
    x = d.mjd.values
    y = d.mag.values
    low = d.mag.values - d.magerr.values
    upp = d.mag.values + d.magerr.values
    return x,y,low,upp

ztf_gmjd,ztf_gmag,ztf_glow,ztf_gupp = get_source_data(ztf_gdata)
ztf_rmjd,ztf_rmag,ztf_rlow,ztf_rupp = get_source_data(ztf_rdata)
lco_gmjd,lco_gmag,lco_glow,lco_gupp = get_source_data(lco_gdata)
lco_rmjd,lco_rmag,lco_rlow,lco_rupp = get_source_data(lco_rdata)

# Create CDS objects for ZTF light curves.  These get updated with button clicks.
source_ztf_g = ColumnDataSource(data=dict(x=ztf_gmjd, y=ztf_gmag, 
                                          lower=ztf_glow, upper=ztf_gupp))
source_ztf_r = ColumnDataSource(data=dict(x=ztf_rmjd, y=ztf_rmag, 
                                          lower=ztf_rlow, upper=ztf_rupp))
source_lco_g = ColumnDataSource(data=dict(x=lco_gmjd, y=lco_gmag, 
                                          lower=lco_glow, upper=lco_gupp))
source_lco_r = ColumnDataSource(data=dict(x=lco_rmjd, y=lco_rmag, 
                                          lower=lco_rlow, upper=lco_rupp))

# Create CDS objects for ZTF light curves.  These are never changed 
# but are used to update the CDS's above during button click events.
orig_ztf_g = ColumnDataSource(data=dict(x=ztf_gmjd, y=ztf_gmag, 
                                        lower=ztf_glow, upper=ztf_gupp))
orig_ztf_r = ColumnDataSource(data=dict(x=ztf_rmjd, y=ztf_rmag, 
                                        lower=ztf_rlow, upper=ztf_rupp))
orig_lco_g = ColumnDataSource(data=dict(x=lco_gmjd, y=lco_gmag, 
                                        lower=lco_glow, upper=lco_gupp))
orig_lco_r = ColumnDataSource(data=dict(x=lco_rmjd, y=lco_rmag, 
                                        lower=lco_rlow, upper=lco_rupp))


# Initialize Light Curve Plot
fig_lc = Figure(plot_height=320, plot_width=750,sizing_mode='scale_both',
                tools="pan,wheel_zoom,box_zoom,tap,reset",
                toolbar_location="above", border_fill_color="whitesmoke")
fig_lc.toolbar.logo = None
fig_lc.yaxis.axis_label = 'Mag'
fig_lc.xaxis.axis_label = 'MJD'
fig_lc.y_range = Range1d(start=yupp, end=ylow)


# Choose plotting colors
c_ztfg = "cornflowerblue"
c_ztfr = "firebrick"
c_lcog = "forestgreen"
c_lcor = "darkorange"

# Choose plotting alphas
a_ztfg = 1.0
a_ztfr = 1.0
a_lcog = 1.0
a_lcor = 1.0


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
               legend = 'ZTF-g')
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
               legend = 'ZTF-r')


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
               legend = 'LCOGT-g')
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
               legend = 'LCOGT-r')

# Modify the legend's appearance
fig_lc.legend.location = "bottom_right"
fig_lc.legend.border_line_width = 3
fig_lc.legend.border_line_color = "navy"
fig_lc.legend.border_line_alpha = 0.5
fig_lc.legend.label_text_font_style = "bold"
fig_lc.legend.label_text_font_size = '9pt'
fig_lc.legend.spacing = -3
fig_lc.legend.label_standoff = 0
fig_lc.legend.padding = 5


# Add Toggle Buttons
b1 = Toggle(label="ZTF"   ,button_type="primary",active=True)
b2 = Toggle(label="LCOGT" ,button_type="primary",active=True)
b3 = Toggle(label="g-data",button_type="primary",active=True)
b4 = Toggle(label="r-data",button_type="primary",active=True)


# JS Callback for Buttons. Toggles the display of various data sets.
callback = CustomJS(args=dict(b1=b1, b2=b2, b3=b3, b4=b4,
                              source1=source_ztf_g,source2=source_ztf_r,
                              source3=source_lco_g,source4=source_lco_r,
                              orig1=orig_ztf_g,orig2=orig_ztf_r,
                              orig3=orig_lco_g,orig4=orig_lco_r), code="""
    
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
        source1.data['x'] = []
        source1.data['y'] = []
        source2.data['x'] = []
        source2.data['y'] = []
        source3.data['x'] = []
        source3.data['y'] = []
        source4.data['x'] = []
        source4.data['y'] = []
    } else if (active.length == 2) {
        if (active[0] == 0 & active[1] == 2) {
            source1.data['x'] = orig1.data['x']
            source1.data['y'] = orig1.data['y']
            source2.data['x'] = []
            source2.data['y'] = []
            source3.data['x'] = []
            source3.data['y'] = []
            source4.data['x'] = []
            source4.data['y'] = []
        } else if (active[0] == 0 & active[1] == 3) {
            source1.data['x'] = []
            source1.data['y'] = []
            source2.data['x'] = orig2.data['x']
            source2.data['y'] = orig2.data['y']
            source3.data['x'] = []
            source3.data['y'] = []
            source4.data['x'] = []
            source4.data['y'] = []
        } else if (active[0] == 1 & active[1] == 2) {
            source1.data['x'] = []
            source1.data['y'] = []
            source2.data['x'] = []
            source2.data['y'] = []
            source3.data['x'] = orig3.data['x']
            source3.data['y'] = orig3.data['y']
            source4.data['x'] = []
            source4.data['y'] = []
        } else if (active[0] == 1 & active[1] == 3) {
            source1.data['x'] = []
            source1.data['y'] = []
            source2.data['x'] = []
            source2.data['y'] = []
            source3.data['x'] = []
            source3.data['y'] = []
            source4.data['x'] = orig4.data['x']
            source4.data['y'] = orig4.data['y']
        } else {
            source1.data['x'] = []
            source1.data['y'] = []
            source2.data['x'] = []
            source2.data['y'] = []
            source3.data['x'] = []
            source3.data['y'] = []
            source4.data['x'] = []
            source4.data['y'] = []
        }
    } else if (active.length == 3) {
        if (active[0] == 0 & active[1] == 1 & active[2] == 2) {
            source1.data['x'] = orig1.data['x']
            source1.data['y'] = orig1.data['y']
            source2.data['x'] = []
            source2.data['y'] = []
            source3.data['x'] = orig3.data['x']
            source3.data['y'] = orig3.data['y']
            source4.data['x'] = []
            source4.data['y'] = []
        } else if (active[0] == 0 & active[1] == 1 & active[2] == 3) {
            source1.data['x'] = []
            source1.data['y'] = []
            source2.data['x'] = orig2.data['x']
            source2.data['y'] = orig2.data['y']
            source3.data['x'] = []
            source3.data['y'] = []
            source4.data['x'] = orig4.data['x']
            source4.data['y'] = orig4.data['y']
        } else if (active[0] == 0 & active[1] == 2 & active[2] == 3) {
            source1.data['x'] = orig1.data['x']
            source1.data['y'] = orig1.data['y']
            source2.data['x'] = orig2.data['x']
            source2.data['y'] = orig2.data['y']
            source3.data['x'] = []
            source3.data['y'] = []
            source4.data['x'] = []
            source4.data['y'] = []
        } else if (active[0] == 1 & active[1] == 2 & active[2] == 3) {
            source1.data['x'] = []
            source1.data['y'] = []
            source2.data['x'] = []
            source2.data['y'] = []
            source3.data['x'] = orig3.data['x']
            source3.data['y'] = orig3.data['y']
            source4.data['x'] = orig4.data['x']
            source4.data['y'] = orig4.data['y']
        }
    } else {
        source1.data['x'] = orig1.data['x']
        source1.data['y'] = orig1.data['y']
        source2.data['x'] = orig2.data['x']
        source2.data['y'] = orig2.data['y']
        source3.data['x'] = orig3.data['x']
        source3.data['y'] = orig3.data['y']
        source4.data['x'] = orig4.data['x']
        source4.data['y'] = orig4.data['y']
    }

    // Transmit changes to the column data sources
    source1.change.emit();
    source2.change.emit();
    source3.change.emit();
    source4.change.emit();
""")


# Connect buttons to callback function
b1.js_on_click(callback)
b2.js_on_click(callback)
b3.js_on_click(callback)
b4.js_on_click(callback)

# Create plot grid
l = layout(widgetbox(row([b1,b2,b3,b4]),sizing_mode='stretch_width'),[fig_lc],
           sizing_mode='scale_both')

# Use "show" to generate HTML file and launch plot in the browser
show(l)






