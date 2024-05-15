import arcpy

arcpy.ImportToolbox(r"C:\Users\rzomor2\OneDrive - University of Illinois Chicago\MyDrive\Projects\CDC\RST_V2\20240503_RST_david_v2\20240503_RST_david_v2.pyt")
arcpy.toolbox.Tool(
    InputTable=r"C:\Users\rzomor2\OneDrive - University of Illinois Chicago\MyDrive\Projects\CDC\RST_V2\20240503_RST_david_v2\data\data_table_mi.csv",
    InputTableFields="GEOID Region;group Group;mortality Event;population Population",
    InputAdjacencyJSON=r"C:\Users\rzomor2\OneDrive - University of Illinois Chicago\MyDrive\Projects\CDC\RST_V2\20240503_RST_david_v2\data\miadj_cty.JSON",
    StandardPopulationYear="2000",
    OutputTable=r"C:\Users\rzomor2\OneDrive - University of Illinois Chicago\MyDrive\Projects\CDC\RST_V2\20240503_RST_david_v2\data\medians_new.csv"
)