# -*- coding: utf-8 -*-
"""
by: Emily Sturdivant, emilysturdivant@gmail.com
requires: python 3.6+

OVERVIEW: Getting started with processing Haiti field data
"""
#%% Relate filenames to plot names
fname_plotname = {
    'Augier': {'xlsname':'FormularioDAPAugier.xls',
                            'shapefile': 'Augier.shp',
                            'plot_id':'Augier',
                            'qual_file': 'CualitativoDAugier.docx'},
    'Bombardopolis Plot 1': {'xlsname':'FormularioDAP1BO.xls',
                            'shapefile': 'Plot_1BO.shp',
                            'plot_id':'1BO',
                            'qual_file': 'CualitativoBOM1.docx'},
    'Bombardopolis Plot 2': {'xlsname':'FormatoCampo_2BO.xls',
                            'shapefile': 'Plot_2BO.shp',
                            'plot_id':'2BO'
                            'qual_file': 'CualitativoBOM2.docx'},
    'Cuyo Plot 1': {'xlsname':'FormularioDAP1Cuyo.xls',
                            'shapefile': '1Cuyo.shp',
                            # 'plot_id':'1Cuyo'
                            'qual_file': 'CualitativoCuyo.docx'},
    'Cuyo Plot 2': {'xlsname':'FormularioDAP2cuyo2.xls',
                            'shapefile': '2_cuyo.shp',
                            # 'plot_id':'2_cuyo'
                            'qual_file': 'Cualitativo2Cuyo2.docx'},
    'Decameron': {'xlsname':'FormularioDAPDeca.xls',
                            'shapefile': 'Decameron.shp',
                            # 'plot_id':'2_cuyo'
                            'qual_file': 'Cualitativo2Decameron.docx'},
    'Deluge Plot 1': {'xlsname':'FormularioDAPDeluge.xls',
                            'shapefile': 'Deluge_.shp',
                            # 'plot_id':'2_cuyo'
                            'qual_file': 'CualitativoDeluge.docx'},
    'Deluge Plot 2': {'xlsname':'FormularioDAPDeluger2.xls',
                            'shapefile': '2Delugert.shp',
                            # 'plot_id':'2_cuyo'
                            'qual_file': 'Cualitativo2deluguier.docx'},
    'Deluge Plot 3': {'xlsname':'FormularioDAPDeluger3.xls',
                            'shapefile': '3Deluger.shp',
                            # 'plot_id':'2_cuyo'
                            'qual_file': 'Cualitativo3deluguier3.docx'},
    'Gonaives Plot 1': {'xlsname':'PlotCheck_1GON.xlsx',
                            'shapefile': 'Plot_1GON.shp',
                            'plot_id':'1GON'
                            'qual_file': 'FormatoCampo_1GON.xlsx'},
    'Gonaives Plot 2': {'xlsname':'2GON.xlsx',
                            'shapefile': 'Plot_2GON.shp',
                            'plot_id':'2GON'
                            'qual_file': 'Cualitativo2GON.docx'},
    'Gonaives Plot 3': {'xlsname':'3GON.xls',
                            'shapefile': '3GON.shp',
                            'plot_id':'3GON'
                            'qual_file': 'Cualitativo3GON.docx'},
    'Campeche': {'xlsname':'FormularioDAPcampeche1.xls',
                            'shapefile': 'Campeche.shp',
                            # 'plot_id':'campeche'
                            'qual_file': 'CualitativoCampeche.docx'},
    'Casimir': {'xlsname':'FormularioDAPCasimir1.xls',
                            'shapefile': '1casimir.shp',
                            # 'plot_id':'1casimir'
                            'qual_file': 'Cualitativo2Casimir.docx'},
    'Casimir Plot 2': {'xlsname':'FormularioDAPCasimir2.xls',
                            'shapefile': '2casimir.shp',
                            # 'plot_id':'2casimir'
                            'qual_file': 'Cualitativo2Casimir2.docx'},
    'Cobana Plot 1': {'xlsname':'FormularioDAPCobana1.xls',
                            'shapefile': '1cobana.shp',
                            # 'plot_id':'1cobana'
                            'qual_file': 'CualitativoCobana.docx'},
    'Cobana Plot 2': {'xlsname':'FormularioDAPCobana2.xlsx',
                            'shapefile': '2cobana.shp',
                            # 'plot_id':'2Cobana'
                            'qual_file': 'CualitativoCobana2.docx'},
    'Lemanguier': {'xlsname':'FormularioDAPLemanguier.xlsx',
                            'shapefile': 'Lemanguier.shp',
                            # 'plot_id':'2Cobana'
                            'qual_file': 'CualitativoLemanguier.docx'},
    'Lerochekit': {'xlsname':'FormularioDAPLerochekit.xlsx',
                            'shapefile': 'Lerochekit.shp',
                            # 'plot_id':'2Cobana'
                            'qual_file': 'CualitativoLerochekit.docx'},
    'Pandiassou Plot 1': {'xlsname':'FormularioDAPPandiassou.xlsx',
                            'shapefile': 'Pandiassou.shp',
                            # 'plot_id':'2Cobana'
                            'qual_file': 'CualitativoPandiassou.docx'},
    'Pandiassou Plot 2': {'xlsname':'FormularioDAPPandiassou2.xlsx',
                            'shapefile': '2pandiassou.shp',
                            # 'plot_id':'2Cobana'
                            'qual_file': 'CualitativoPandiassou2.docx'},
    'Pewpapai': {'xlsname':'FormularioDAPewpapay.xlsx',
                            'shapefile': 'Pewpapay.shp',
                            # 'plot_id':'Pewpapay'
                            'qual_file': 'CualitativoPewpapay1.docx'},
    'PewPapay2': {'xlsname':'FormularioDAPBaspapay.xlsx',
                            'shapefile': 'Pewpapay2.shp',
                            'shapefile': '2Pewpapay2.shp',
                            # 'plot_id':'Pewpapay'
                            'qual_file': 'CualitativoPewpapay2.docx'},
    'Pewpapay': {'xlsname':'FormularioDAPBaspapay.xls',
                            'shapefile': 'baspapay.shp',
                            # 'plot_id':'Pewpapay'
                            'qual_file': ''},
