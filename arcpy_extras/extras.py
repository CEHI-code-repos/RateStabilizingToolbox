import arcpy
import numpy as np
import pandas as pd
from collections import namedtuple

def exists(url):
    if url is not None and arcpy.Exists(url):
        return True
    return False

def get_fieldType(url, fieldName):
    if exists(url) and fieldName is not None:
        field_type = [f.type for f in arcpy.ListFields(url) if f.name == fieldName]
        if len(field_type) == 1:
            return field_type[0]
    return None

def get_fieldList(url, fieldName):
    if get_fieldType(url, fieldName) is not None:
        return [val[0] for val in arcpy.da.SearchCursor(url, fieldName)]
    return None
        
def get_pandas(url, fields):
    if exists(url) and fields is not None:
        valid_fields = [f for f in arcpy.ListFields(url) if f.name in fields]
        if len(valid_fields) == len(fields):
            return pd.DataFrame(data = arcpy.da.SearchCursor(url, fields), columns = fields)
    return None

def get_valueTableValues(valueTable):
    if valueTable.values is not None:
        return [[str(value) for value in field] for field in valueTable.values]
    else:
        return [[None for i, col in enumerate(valueTable.columns)]]

def get_fieldInfo(url, fieldName):
    fieldType = get_fieldType(url, fieldName)
    fieldList = get_fieldList(url, fieldName)
    exists = fieldType is not None
    field_info = namedtuple("field_info", ["name", "type", "list", "exists"])
    return field_info(fieldName, fieldType, fieldList, exists)
    
def set_valueTableRequired(valueTable):
    if valueTable.values is None:
        valueTable.setIDMessage('ERROR', 530)
    else:
        valueTable_str = [str(val) for row in valueTable.values for val in row]
        if "" in valueTable_str:
            valueTable.setIDMessage('ERROR', 530)

def set_parameterRequired(parameter):
    if parameter.value is None:
        parameter.setIDMessage('ERROR', 530)

def row_string(indexes):
    index_list = indexes[0:min(4, len(indexes))]
    output = ", ".join([str(x + 1) for x in index_list])
    if len(indexes) > 4:
        output += " and more"
    if len(indexes) == 1:
        output = "row " + output
    else:
        output = "rows " + output
    return output

def pandas_to_table(pd_table, out_path):
    output_np = np.rec.fromrecords(pd_table, names = list(pd_table.columns))
    arcpy.da.NumPyArrayToTable(output_np, out_path)

    arc_project = arcpy.mp.ArcGISProject("Current")
    cur_map = arc_project.activeMap
    if cur_map:
        cur_map.addDataFromPath(out_path)
        