import processing
import os

out_dir = r'/Users/emilysturdivant/PROJECTS/Haiti_biomass'
lc = 'Haiti2017_Clip'
agb = 'agb18_haiti_v6_0to310'

# Mask AGB to landcover classes 4–6
exp = f'((\"{lc}@1\" > 3)/(\"{lc}@1\" > 3)) / ((\"{lc}@1\" > 3)*1 + (\"{lc}@1\"<=3)*0) * \"{agb}@1\"'

processing.run("qgis:rastercalculator", \
{'EXPRESSION': exp,\
'LAYERS':[os.path.join(out_dir, 'R_out', 'agb18_haiti_v6_0to310.tif'),\
'CELLSIZE':None,\
'EXTENT':None,\
'CRS':QgsCoordinateReferenceSystem('EPSG:4326'),\
'OUTPUT':os.path.join(out_dir, 'LULC', 'agb_lc_over3.tif')})
