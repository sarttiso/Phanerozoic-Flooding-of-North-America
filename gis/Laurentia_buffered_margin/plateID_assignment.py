from qgis.core import *

layer_names = ['Greenland_Laurentia_%04dkm' % x for x in range(100, 1800, 100)]

for layer in layer_names:
    cur_layer = QgsProject.instance().mapLayersByName(layer)[0]
    
    with edit(cur_layer):
        for feature in cur_layer.getFeatures():
            feature['CONTINENT'] = 102
            cur_layer.updateFeature(feature)