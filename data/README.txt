This file describes the content of this folder. The folder includes all data needed for the workflow developed in the Master's thesis of Nathalie Heiniger with the title:
Identifying Anthropogenic Feeding Sites from GPS Tracking Data: A Case Study for Red Kites (Milvus milvus) in Western Switzerland. (October 2020) 

Folders

Buildings
	Folder includes:
	- TLM3d data for the buildings
	- TLM_buildings: original data for Switzerland
	- tlm_buildings_size65_buff_50m_dissolved: edited data, only buildings larger 65m^2 for Switzerland
	- tlm_buildings_studyareaExtra_size65_buff_50m_dissolved.shp: edited data for study area, multipolygons
	- tlm_buildings_studyareaExtra_size65_buff_50m_sf_singlepoly.shp, edited data for study area, single polygons
	
	
Europe
	Folder includes:
	- shapefile with shape of Europe

Forest
	Folder includes:
	- vec25_forest_buff_20m_CH: shapefile with forest areas, 20m buffer around those areas for Switzerland
	- vec25_CH_NOTforest: shapefile with areas that are not forest for Switzerland (inverse version of vec25_forest_buff_20m_CH)
	- vec25_forest_buff_20m_studyarea: shapefile with forest areas, 20m buffer around those areas for study area
	- vec25_studyarea_NOTforest: shapefile with areas that are not forest for study area (inverse version of vec25_forest_buff_20m_studyarea)

Life_History
	Folder includes:
	- Excel, CSV and .rdata files with information about the life histories of the red kites. Edited by N. Heiniger with extra column for the age class in season 2018/2019.

Private_Feeders
	Folder includes:
	- CSV files with information about private feeders
	- private_feeders.csv: known feeders within Switzerland
	- private_feeders_studyarea_LV95: private feeders within study area

Study_Area
	Folder includes:
	- study_area: shapefile for extent of study area
	- studyarea_buffer_extra2km: shapefile with a extra 2km buffer around study area (necessary for creation of grid for KDE calculation)
	- studyarea_bbox_extra2km_SpatialPixels_60m.rdata: 60m grid needed for KDE calculation

Switzerland
	Folder includes:
	- shapefile for extent of Switzerland with a 5km buffer