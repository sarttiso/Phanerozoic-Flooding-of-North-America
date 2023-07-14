import arcpy
import numpy as np

# arcpy.env.workspace = "C:/Users/adrian/Documents/research/macrostrat/gis"
arcpy.env.workspace = "C:/Users/XPS/Documents/research/macrostrat/gis"

buffers = np.arange(100, 1800, 100)

filenames = ['Laurentia_%04dkm' % x for x in buffers]

for ii, filename in enumerate(filenames):
    arcpy.analysis.Buffer("Laurentia_no-north-slope", 
                          "Laurentia_buffered_margin_no-north-slope/"+ filename, 
                          "-%d Kilometers" % buffers[ii], 
                          "FULL", "ROUND", "LIST", "CONTINENT", "GEODESIC")