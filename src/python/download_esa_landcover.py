# Prior to running, follow directions at https://cds.climate.copernicus.eu/api-how-to
# To run using R, first try to run:
# `library('reticulate')`
# `py_run_file('src/python/download_esa_landcover.py')`
# Then, in Terminal:
# `conda activate /home/esturdivant/PROJECTS/biomass-espanola/renv/python/r-reticulate`
# `pip install cdsapi`

import cdsapi

c = cdsapi.Client()

c.retrieve(
    'satellite-land-cover',
    {
        'variable': 'all',
        'format': 'zip',
        'year': [
            '2017', '2019',
        ],
        'version': 'v2.1.1',
    },
    'download.zip')
