# -*- coding: utf-8 -*-
import arcpy
import pandas as pd
import numpy as np
import helpers

import importlib
importlib.reload(helpers)

class Toolbox:
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Rate Stabilizing Toolbox"
        self.alias = "RST"

        # List of tool classes associated with this toolbox
        self.tools = [RST, IDP]


class RST:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rate Stabilizing Tool"
        self.description = "Accurate and precise estimates of local-level epidemiologic measures are critical to informing policy and program decisions, but they often require advanced statistical knowledge, programming/coding skills, and extensive computing power. In response, we developed the Rate Stabilizing Tool (RST), an ArcGIS-based tool that enables users to input their own record-level data to generate more reliable age-standardized measures of chronic disease (eg, prevalence rates, mortality rates) or other population health outcomes at the county or census tract levels. The RST uses 2 forms of empirical Bayesian modeling (nonspatial and spatial) to estimate age-standardized rates and 95% credible intervals for user-specified geographic units. The RST also provides indicators of the reliability of point estimates. In addition to reviewing the RST’s statistical techniques, we present results from a simulation study that illustrates the key benefit of smoothing. We demonstrate the dramatic reduction in root mean-squared error (rMSE), indicating a better compromise between accuracy and stability for both smoothing approaches relative to the unsmoothed estimates. Finally, we provide an example of the RST’s use. This example uses heart disease mortality data for North Carolina census tracts to map the RST output, including reliability of estimates, and demonstrates a subsequent statistical test."

    def getParameterInfo(self):
        """Define the tool parameters."""

        param_data_table = arcpy.Parameter(
            displayName="Input Table",
            name="InputTable",
            datatype="GPTableView",
            parameterType="Required",
            direction="Input"
        )

        param_data_fields = arcpy.Parameter(
            displayName="Input Table Fields",
            name="InputTableFields",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input"
        )
        param_data_fields.parameterDependencies = [param_data_table.name]
        param_data_fields.columns = [['Field', 'Region ID'], ['Field', 'Event Count'], ['Field', 'Population Count']]
        param_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'

        param_feature = arcpy.Parameter(
            displayName="Input Feature",
            name="InputFeature",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input"
        )
        param_feature.filter.list = ["Polygon"]

        param_feature_field = arcpy.Parameter(
            displayName="Input Feature Field",
            name="InputFeatureField",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input"
        )
        param_feature_field.parameterDependencies = [param_feature.name]
        param_feature_field.columns = [['Field', 'Region ID']]
        param_feature_field.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'

        param_out_table = arcpy.Parameter(
            displayName="Output Table",
            name="OutputTable",
            datatype="GPTableView",
            parameterType="Required",
            direction="Output"
        )

        param_age_grp_field = arcpy.Parameter(
            displayName="Age Group Field",
            name="AgeGroupField",
            datatype="Field",
            parameterType="Optional",
            direction="Input",
            category="Age Standardization (optional)"
        )
        param_age_grp_field.parameterDependencies = [param_data_table.name]

        param_std_pop_yr = arcpy.Parameter(
            displayName="Standard Population Year",
            name="StandardPopulationYear",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Age Standardization (optional)"
        )
        param_std_pop_yr.filter.list = ["2000", "2010"]

        param_age_std_groups = arcpy.Parameter(
            displayName="Standardized Age Groups",
            name="StandardizedAgeGroups",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input",
            category="Age Standardization (optional)"
        )
        param_age_std_groups.columns = [['String', 'Lower age value'], ['String', 'Upper age value']]
        param_age_std_groups.filters[0].type = "ValueList"
        param_age_std_groups.filters[0].list = ["0", "5", "15", "25", "35", "45", "55", "65", "75"]
        param_age_std_groups.filters[1].type = "ValueList"
        param_age_std_groups.filters[1].list = ["14", "24", "34", "44", "54", "64", "74", "84", "up"]

        param_rates_per = arcpy.Parameter(
            displayName="Rate",
            name="Rate",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input"
        )
        param_rates_per.columns = [['Long', 'Per']]
        param_rates_per.values = [[100000]]
        param_rates_per.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'

        params = [
            param_data_table,
            param_data_fields,
            param_feature,
            param_feature_field,
            param_rates_per,
            param_out_table,
            param_age_grp_field,
            param_std_pop_yr,
            param_age_std_groups
        ]
        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        data_url = parameters[0]
        data_fields = parameters[1]
        feature_url = parameters[2]
        feature_fields = parameters[3]
        rates_per = parameters[4]
        estimates_out = parameters[5]
        data_ageGrp_id = parameters[6]
        std_pop_yr = parameters[7]
        age_std_groups = parameters[8]

        # Restrict age groups to only those present within data
        data_ageGrp_info  = helpers.get_fieldInfo(data_url.valueAsText, data_ageGrp_id.valueAsText)
        if data_ageGrp_info.exists:
            ageGrp_unique = list(set(data_ageGrp_info.list))

            grp_notValid = [group for group in ageGrp_unique if group not in helpers.const_age_grps]

            if len(grp_notValid) == 0:
                lvs = sorted([int(group.split("-")[0]) for group in ageGrp_unique if group != "85up"])
                uvs = sorted([int(group.split("-")[1]) if group != "85up" else 85 for group in ageGrp_unique])
                age_std_groups.filters[0].list = [str(lv) for lv in lvs[:-1]]
                age_std_groups.filters[1].list = ["up" if uv == 85 else str(uv) for uv in uvs[1:]]

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""

        data_url = parameters[0]
        data_fields = parameters[1]
        feature_url = parameters[2]
        feature_fields = parameters[3]
        rates_per = parameters[4]
        estimates_out = parameters[5]
        data_ageGrp_id = parameters[6]
        std_pop_yr = parameters[7]
        age_std_groups = parameters[8]

        data_fields_name = helpers.get_valueTableNames(data_fields)
        data_region_info = helpers.get_fieldInfo(data_url.valueAsText, data_fields_name[0])
        data_event_info = helpers.get_fieldInfo(data_url.valueAsText, data_fields_name[1])
        data_pop_info  = helpers.get_fieldInfo(data_url.valueAsText, data_fields_name[2])
        data_ageGrp_info  = helpers.get_fieldInfo(data_url.valueAsText, data_ageGrp_id.valueAsText)
        feature_fields_name = helpers.get_valueTableNames(feature_fields)
        feature_region_info = helpers.get_fieldInfo(feature_url.valueAsText, feature_fields_name[0])

        # Check if all fields are filled in for age standardization
        if data_ageGrp_info.name is None and (std_pop_yr.valueAsText is not None or age_std_groups.valueAsText is not None):
            helpers.set_parameterRequired(data_ageGrp_id)
        if std_pop_yr.valueAsText is None and (data_ageGrp_info.name is not None or age_std_groups.valueAsText is not None):
            helpers.set_parameterRequired(std_pop_yr)
        if age_std_groups.valueAsText is None and (data_ageGrp_info.name is not None or std_pop_yr.valueAsText is not None):
            helpers.set_valueTableRequired(age_std_groups)

        # Check for Nulls
        if data_region_info.exists and None in data_region_info.list:
            data_fields.setErrorMessage("Input Table Region ID Field contains at least one Null value")
        if data_event_info.exists and None in data_event_info.list:
            data_fields.setErrorMessage("Input Table Event Count Field contains at least one Null value")
        if data_pop_info.exists and None in data_pop_info.list:
            data_fields.setErrorMessage("Input Table Population Count Field contains at least one Null value")
        if data_ageGrp_info.exists and None in data_ageGrp_info.list:
            data_ageGrp_id.setErrorMessage("Age Group Field contains at least one Null value")
        if feature_region_info.exists and None in feature_region_info.list:
            feature_fields.setErrorMessage("Input Feature Region ID Field contains at least one Null value")
        
        # Check for data types
        if data_region_info.exists and data_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            data_fields.setErrorMessage("Input Table Region ID Field is not an Integer or String")
        if data_event_info.exists and data_event_info.type not in ["SmallInteger", "Integer", "BigInteger"]:
            data_fields.setErrorMessage("Input Table Event Count Field is not an Integer")
        if data_pop_info.exists and data_pop_info.type not in ["SmallInteger", "Integer", "BigInteger"]:
            data_fields.setErrorMessage("Input Table Population Count Field is not an Integer")
        if data_ageGrp_info.exists and data_ageGrp_info.type not in ["String"]:
            data_ageGrp_id.setErrorMessage("Age Group Field is not a String")
        if feature_region_info.exists and feature_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            feature_fields.setErrorMessage("Input Feature Region ID Field is not an Integer or String")

        # Check if Output Table exists
        if helpers.exists(estimates_out.valueAsText):
            estimates_out.setErrorMessage("Output Table already exists")

        # Check if Region IDs are repeated when age is adjusted for
        if not data_ageGrp_info.exists and data_region_info.exists and len(data_region_info.list) != len(set(data_region_info.list)):
            data_fields.setWarningMessage("Repeated Input Data Region IDs are detected. Population Counts will be aggregated to totals.")

        if (data_region_info.exists and feature_region_info.exists and 
            not data_fields.hasError() and not feature_fields.hasError()):
            # Check if Region ID types are the same
            if data_region_info.type != feature_region_info.type:
                feature_fields.setErrorMessage("Input Feature Region ID Field type does not match Input Table Region ID Field type")
            else:
                # Check if Data contains Region IDs not present in Feature
                data_only_regions = set(data_region_info.list) - set(feature_region_info.list)
                if len(data_only_regions) != 0:
                    data_fields.setErrorMessage(data_only_regions)
                # Check if Feature contains Region IDs not present in Data
                feature_only_regions = set(feature_region_info.list) - set(data_region_info.list)
                if len(feature_only_regions) != 0:
                    feature_fields.setErrorMessage("Input Feature Region ID Field contains at least one value not present in Input Table Region ID Field")
            
        if data_ageGrp_info.exists and not data_ageGrp_id.hasError():
            ageGrp_unique = list(set(data_ageGrp_info.list))
            ageGrp_invalid = [group for group in ageGrp_unique if group not in helpers.const_age_grps]
            # Check if there are any invalid age groups
            if len(ageGrp_invalid) != 0:
                data_ageGrp_id.setErrorMessage('Age Group Field contains at least one invalid age group')
            elif data_region_info.exists:
                ageGrp_dict = {}
                for i, region in enumerate(data_region_info.list):
                    ageGrp_dict.setdefault(region, []).append(data_ageGrp_info.list[i])
                regions_ageGrp_nunique = [len(set(ageGrp_dict[region])) for region in ageGrp_dict]
                regions_ageGrp_n = [len(ageGrp_dict[region]) for region in ageGrp_dict]
                # Check if there duplicate age groups
                if regions_ageGrp_nunique != regions_ageGrp_n:
                    data_ageGrp_id.setErrorMessage("At least one Age Group is duplicated within a Region")
                # Check if there missing age groups
                elif len(set(regions_ageGrp_n)) != 1:
                    data_ageGrp_id.setErrorMessage("At least one Region is missing an Age Group")

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        messages.AddMessage("Preparing data...")

        data_url = parameters[0]
        data_fields = parameters[1]
        feature_url = parameters[2]
        feature_fields = parameters[3]
        rates_per = parameters[4]
        estimates_out = parameters[5]
        data_ageGrp_id = parameters[6]
        std_pop_yr = parameters[7]
        age_std_groups = parameters[8]
        
        data_fields_name = helpers.get_valueTableNames(data_fields)
        if data_ageGrp_id.valueAsText is not None:
            data_fields_name.append(data_ageGrp_id.valueAsText)
        data_region_name = data_fields_name[0]
        data_event_name = data_fields_name[1]
        data_pop_name = data_fields_name[2]
        feature_region_name = str(feature_fields.values[0][0])

        rates_per = int(rates_per.valueAsText)

        # Get the age group distribution
        age_std_groups_arr = []
        age_std_groups_names = []
        if age_std_groups.values is not None:
            for lv, uv, in age_std_groups.values:
                lv_index = [i for i, grp in enumerate(helpers.const_age_grps) if grp.startswith(lv)][0]
                uv_index = [i for i, grp in enumerate(helpers.const_age_grps) if grp.endswith(uv)][0]
                age_std_groups_arr.append(helpers.const_age_grps[lv_index:(uv_index + 1)])
                age_std_groups_names.append( lv + ("to" + uv if uv != "85up" else "up") )

        # Read in data
        data = helpers.get_pandas(data_url.valueAsText, data_fields_name)
        age_groups = [""]
        num_group = 1
        data = data.sort_values(by = [data_region_name])
        regions = data[data_region_name].unique().tolist()
        num_region = data[data_region_name].nunique()
        if data_ageGrp_id.valueAsText is not None:
            data = data.sort_values(by = [data_region_name, data_ageGrp_id.valueAsText])
            age_groups = data[data_ageGrp_id.valueAsText].unique().tolist()
            num_group = data[data_ageGrp_id.valueAsText].nunique()
        elif num_region != len(data.index):
            data = data.groupby(data_region_name).agg({data_event_name : 'sum', data_pop_name : "sum"})
            messages.addWarningMessage("Repeated Region IDs were detected. Population and Event Counts were aggregated to totals.")
        
        Y = np.array(data[data_event_name]).reshape([num_region, num_group])
        n = np.array(data[data_pop_name]).reshape([num_region, num_group])

        std_pop = 1
        # Klein, et al. 2000 standard population
        if (std_pop_yr.valueAsText == "2000"): 
            std_pop = np.array([18987, 39977, 38077, 37233, 44659, 37030, 23961, 18136, 12315, 4259])
        # 2010 standard population
        if (std_pop_yr.valueAsText == "2010"):
            std_pop = np.array([20203362, 41025851, 43626342, 41063948, 41070606, 45006716, 36482729, 21713429, 13061122, 5493433])
        if std_pop_yr.valueAsText is not None:
            std_pop = std_pop[np.isin(helpers.const_age_grps, age_groups)]

        # Create adjacency matrix
        arcpy.AddMessage("Creating adjacency matrix ...")
        arcpy.management.Sort(in_dataset = feature_url.valueAsText, out_dataset = r"in_memory\feature", sort_field = feature_region_name + " ASCENDING")
        arcpy.management.AddField(in_table=r"in_memory\feature", field_name="NUM_REGION", field_type="LONG")
        with arcpy.da.UpdateCursor(r"in_memory\feature", ["NUM_REGION"]) as cursor:
            for index, row in enumerate(cursor):
                row[0] = index + 1
                cursor.updateRow(row)
        arcpy.analysis.PolygonNeighbors(in_features=r"in_memory\feature", out_table=r"in_memory\adjTable", in_fields="NUM_REGION")
        adj_table = arcpy.da.TableToNumPyArray(r"in_memory\adjTable", ("src_NUM_REGION", "nbr_NUM_REGION"))
        adj_table = np.array([*adj_table.astype(object)]) 
        arcpy.management.Delete(r"in_memory\feature")
        arcpy.management.Delete(r"in_memory\adjTable")
        adj_dict = {}
        for source, neighbor in adj_table:
            adj_dict.setdefault(source, []).append(neighbor)
        adj = list(adj_dict.values())

        # Generate estimates
        messages.AddMessage("Generating estimates...")
        theta_out = helpers.gibbs_rucar(Y, n, adj, std_pop)
        output = helpers.expit(theta_out) * rates_per

        # If age standardized, generate age_groups
        if age_std_groups.values is not None:
            for ages in age_std_groups_arr:
                output = helpers.age_std(output, age_groups, std_pop, ages)
            age_groups.extend(age_std_groups_names)
        
        medians, ci_chart, reliable = helpers.get_medians(output, regions, age_groups)

        # If age standardized, combine output with prefixes, else just rename the output cols
        output = output_cols = []
        if age_std_groups.values is not None:
            medians = medians.add_prefix("median_")
            ci_chart = ci_chart.add_prefix("maxCI_")
            reliable = reliable.add_prefix("reliable_")
            output = pd.concat([medians, ci_chart, reliable], axis = 1)
            for i, med in enumerate(medians.columns):
                output_cols.extend([ medians.columns[i], ci_chart.columns[i], reliable.columns[i] ])
            output = output[output_cols].reset_index()
        else:
            output = pd.concat([medians, ci_chart, reliable], axis = 1).reset_index()
            output_cols = ["median", "maxCI", "reliable"]

        # Write out the final table
        output_cols = [data_region_name] + output_cols
        output_np = np.rec.fromrecords(output, names = output_cols)
        arcpy.da.NumPyArrayToTable(output_np, estimates_out.valueAsText)
    
        messages.AddMessage("Model finished!")

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class IDP:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Individual Data Processing"
        self.description = ""

    def getParameterInfo(self):
        """Define the tool parameters."""

        param_byAge = arcpy.Parameter(
            displayName="By Age",
            name="ByAge",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input"
        )
        param_byAge.value = True

        param_idv_data = arcpy.Parameter(
            displayName="Input Individual Data",
            name="InputIndividualData",
            datatype="GPTableView",
            parameterType="Required",
            direction="Input"
        )

        param_idvWAge_data_fields = arcpy.Parameter(
            displayName="Input Individual Data Fields",
            name="InputIndividualDataFieldsWAge",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input"
        )
        param_idvWAge_data_fields.parameterDependencies = [param_idv_data.name]
        param_idvWAge_data_fields.columns = [['Field', 'Region ID'], ['Field', 'Age']]
        param_idvWAge_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'
        
        param_idvWOAge_data_fields = arcpy.Parameter(
            displayName="Input Individual Data Fields",
            name="InputIndividualDataFieldsWOAge",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input"
        )
        param_idvWOAge_data_fields.parameterDependencies = [param_idv_data.name]
        param_idvWOAge_data_fields.columns = [['Field', 'Region ID']]
        param_idvWOAge_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'
        param_idvWOAge_data_fields.enabled = False

        param_pop_data = arcpy.Parameter(
            displayName="Input Population Data",
            name="InputPopulationData",
            datatype="GPTableView",
            parameterType="Required",
            direction="Input"
        )

        param_popWAge_data_fields = arcpy.Parameter(
            displayName="Input Population Data Fields",
            name="InputPopulationDataFieldsWAge",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input"
        )
        param_popWAge_data_fields.parameterDependencies = [param_pop_data.name]
        param_popWAge_data_fields.columns = [['Field', 'Region ID'], ['Field', 'Population Count'], ['Field', 'Age Group']]
        param_popWAge_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'
        
        param_popWOAge_data_fields = arcpy.Parameter(
            displayName="Input Population Data Fields",
            name="InputPopulationDataFieldsWOAge",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input"
        )
        param_popWOAge_data_fields.parameterDependencies = [param_pop_data.name]
        param_popWOAge_data_fields.columns = [['Field', 'Region ID'], ['Field', 'Population Count']]
        param_popWOAge_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'
        param_popWOAge_data_fields.enabled = False

        param_ftr = arcpy.Parameter(
            displayName="Input Feature",
            name="InputFeature",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input"
        )
        
        param_ftr_fields = arcpy.Parameter(
            displayName="Input Feature Fields",
            name="InputFeatureFields",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input"
        )
        param_ftr_fields.parameterDependencies = [param_ftr.name]
        param_ftr_fields.columns = [['Field', 'Region ID']]
        param_ftr_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'

        param_out_table = arcpy.Parameter(
            displayName="Output Table",
            name="OutputTable",
            datatype="GPTableView",
            parameterType="Required",
            direction="Output"
        )

        params = [
            param_byAge,
            param_idv_data,
            param_idvWAge_data_fields,
            param_idvWOAge_data_fields,
            param_pop_data,
            param_popWAge_data_fields,
            param_popWOAge_data_fields,
            param_ftr,
            param_ftr_fields,
            param_out_table
        ]

        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""

        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        
        byAge = parameters[0]
        idv_data_url = parameters[1]
        idvWAge_data_fields = parameters[2]
        idvWOAge_data_fields = parameters[3]
        pop_data_url = parameters[4]
        popWAge_data_fields = parameters[5]
        popWOAge_data_fields = parameters[6]
        ftr_url = parameters[7]
        ftr_fields = parameters[8]
        out_table = parameters[9]

        idvWAge_data_fields.enabled = byAge.value
        idvWOAge_data_fields.enabled = not byAge.value
        popWAge_data_fields.enabled = byAge.value
        popWOAge_data_fields.enabled = not byAge.value

        idvWAge_fields_name = helpers.get_valueTableNames(idvWAge_data_fields)
        idvWOAge_fields_name = helpers.get_valueTableNames(idvWOAge_data_fields)
        popWAge_fields_name = helpers.get_valueTableNames(popWAge_data_fields)
        popWOAge_fields_name = helpers.get_valueTableNames(popWOAge_data_fields)

        if byAge.value:
            idvWOAge_data_fields.values = [[idvWAge_fields_name[0]]]
            popWOAge_data_fields.values = [popWAge_fields_name[0:2]]
        else:
            idvWAge_data_fields.values = [[idvWOAge_fields_name[0], idvWAge_fields_name[1]]]
            popWAge_data_fields.values = [popWOAge_fields_name[0:2] + [popWAge_fields_name[2]]]
        
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        
        byAge = parameters[0]
        idv_data_url = parameters[1]
        idvWAge_data_fields = parameters[2]
        idvWOAge_data_fields = parameters[3]
        pop_data_url = parameters[4]
        popWAge_data_fields = parameters[5]
        popWOAge_data_fields = parameters[6]
        ftr_url = parameters[7]
        ftr_fields = parameters[8]
        out_table = parameters[9]

        # Get enabled fields
        if byAge.value:
            idv_data_fields = idvWAge_data_fields
            pop_data_fields = popWAge_data_fields
        else:
            idv_data_fields = idvWOAge_data_fields
            pop_data_fields = popWOAge_data_fields

        idv_fields_name = helpers.get_valueTableNames(idv_data_fields)
        idv_region_info = helpers.get_fieldInfo(idv_data_url.valueAsText, idv_fields_name[0])
        idv_age_info = helpers.get_fieldInfo(idv_data_url.valueAsText, idv_fields_name[1] if byAge.value else None)
        pop_fields_name = helpers.get_valueTableNames(pop_data_fields)
        pop_region_info = helpers.get_fieldInfo(pop_data_url.valueAsText, pop_fields_name[0])
        pop_pop_info = helpers.get_fieldInfo(pop_data_url.valueAsText, pop_fields_name[1])
        pop_ageGrp_info = helpers.get_fieldInfo(pop_data_url.valueAsText, pop_fields_name[2] if byAge.value else None)
        ftr_fields_name = helpers.get_valueTableNames(ftr_fields)
        ftr_region_info = helpers.get_fieldInfo(ftr_url.valueAsText, ftr_fields_name[0])
        
        # Make field parameters required based on off byAge
        if byAge.value:
            helpers.set_valueTableRequired(idvWAge_data_fields)
            helpers.set_valueTableRequired(popWAge_data_fields)
        else:
            helpers.set_valueTableRequired(idvWOAge_data_fields)
            helpers.set_valueTableRequired(popWOAge_data_fields)

        # Check for Nulls
        if idv_region_info.exists and None in idv_region_info.list:
            idv_data_fields.setErrorMessage("Input Individual Data Region ID Field contains at least one Null value")
        if byAge.value and idv_age_info.exists and None in idv_age_info.list:
            idv_data_fields.setErrorMessage("Input Individual Data Age Field contains at least one Null value")
        if pop_region_info.exists and None in pop_region_info.list:
            pop_data_fields.setErrorMessage("Input Population Data Region ID Field contains at least one Null value")
        if pop_pop_info.exists and None in pop_pop_info.list:
            pop_data_fields.setErrorMessage("Input Population Data Population Count Field contains at least one Null value")
        if byAge.value and pop_ageGrp_info.exists and None in pop_ageGrp_info.list:
            pop_data_fields.setErrorMessage("Input Population Data Age Field contains at least one Null value")
        if ftr_region_info.exists and None in ftr_region_info.list:
            ftr_fields.setErrorMessage("Input Feature Region ID contains at least one Null value")
        
        # Check for data types
        if idv_region_info.exists and idv_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            idv_data_fields.setErrorMessage("Input Individual Data Region ID Field is not an Integer or String")
        if byAge.value and idv_age_info.exists and idv_age_info.type not in ["SmallInteger", "Integer", "BigInteger"]:
            idv_data_fields.setErrorMessage("Input Individual Data Age Field is not an Integer")
        if pop_region_info.exists and pop_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            pop_data_fields.setErrorMessage("Input Population Data Region ID Field is not an Integer or String")
        if pop_pop_info.exists and pop_pop_info.type not in ["SmallInteger", "Integer", "BigInteger"]:
            pop_data_fields.setErrorMessage("Input Population Data Population Count Field is not an Integer")
        if byAge.value and pop_ageGrp_info.exists and pop_ageGrp_info.type not in [ "String"]:
            pop_data_fields.setErrorMessage("Input Population Data ID Field is not an String")
        if ftr_region_info.exists and ftr_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            idv_data_fields.setErrorMessage("Input Feature Region ID Field is not an Integer or String")
        
        # Check if Individual Age is negative
        if byAge.value and idv_age_info.exists and not idv_data_fields.hasError() and any(age < 0 for age in idv_age_info.list):
            idv_data_fields.setErrorMessage("Input Individual Data Age Field contains at least one negative value")

        # Check if Output Table exists
        if helpers.exists(out_table.valueAsText):
            out_table.setErrorMessage("Output Table already exists")

        # Check if 
        if not byAge.value and pop_region_info.exists and len(pop_region_info.list) != len(set(pop_region_info.list)):
            pop_data_fields.setWarningMessage("Repeated Input Population Data Region IDs are detected. Population Counts will be aggregated to totals if no age standardization is applied.")

        # Check Individual and Population data relationships
        if (idv_region_info.exists and pop_region_info.exists and 
            not idv_data_fields.hasError() and not pop_data_fields.hasError()):
            # Check if Region ID types are the same
            if idv_region_info.type != pop_region_info.type:
                pop_data_fields.setErrorMessage("Input Population Data Region ID Field type does not match Input Individual Data Region ID Field type")
            else:
                # Check if Individual Data contains Region IDs not present in Population Data
                idv_only_regions = set(idv_region_info.list) - set(pop_region_info.list)
                if len(idv_only_regions) != 0:
                    idv_data_fields.setWarningMessage("Input Individual Data Region ID Field contains at least one value not present in Input Population Data Region ID Field")

        # Check Population Data and Feature relationships
        if (pop_region_info.exists and ftr_region_info.exists and
            not pop_data_fields.hasError() and not ftr_fields.hasError()):
            # Check if Region ID types are the same
            if pop_region_info.type != ftr_region_info.type:
                ftr_fields.setErrorMessage("Input Population Data Region ID Field type does not match Input Feature Region ID Field type")
            else:
                # Check if Population Data contains Region IDs not present in Feature
                pop_only_regions = set(pop_region_info.list) - set(ftr_region_info.list)
                if len(pop_only_regions) != 0:
                    pop_data_fields.setErrorMessage("Input Population Data Region ID Field contains at least one value not present in Input Feature Region ID Field")
                # Check if Feature contains Region IDs not present in Population Data
                ftr_only_regions = set(ftr_region_info.list) - set(pop_region_info.list)
                if len(ftr_only_regions) != 0:
                    ftr_fields.setWarningMessage("Input Feature Region ID Field contains at least one value not present in Input Population Data Region ID Field\n\nPopulation within these regions will be assumed to be 0")

        # Check Individual Data and Feature relationships
        if (idv_region_info.exists and ftr_region_info.exists and 
            not idv_data_fields.hasError() and not ftr_fields.hasError()):
            if idv_region_info.type == ftr_region_info.type:
                # Check if Individual Data contains Region IDs not present in Population Data
                idv_only_regions = set(idv_region_info.list) - set(ftr_region_info.list)
                if len(idv_only_regions) != 0:
                    idv_data_fields.setErrorMessage("Input Individual Data Region ID Field contains at least one value not present in Feature Region ID Field")
        
        if pop_ageGrp_info.exists and not pop_data_fields.hasError():
            ageGrp_unique = list(set(pop_ageGrp_info.list))
            ageGrp_invalid = [group for group in ageGrp_unique if group not in helpers.const_age_grps]
            # Check if there are any invalid age groups
            if len(ageGrp_invalid) != 0:
                pop_data_fields.setErrorMessage('Age Group Field contains at least one invalid age group')
            elif pop_region_info.exists:
                ageGrp_dict = {}
                for i, region in enumerate(pop_region_info.list):
                    ageGrp_dict.setdefault(region, []).append(pop_ageGrp_info.list[i])
                regions_ageGrp_nunique = [len(set(ageGrp_dict[region])) for region in ageGrp_dict]
                regions_ageGrp_n = [len(ageGrp_dict[region]) for region in ageGrp_dict]
                # Check if there duplicate age groups
                if regions_ageGrp_nunique != regions_ageGrp_n:
                    pop_data_fields.setErrorMessage("At least one Age Group is duplicated within a Region")
                # Check if there missing age groups
                elif len(set(regions_ageGrp_n)) != 1:
                    pop_data_fields.setErrorMessage("At least one Region is missing an Age Group")

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        byAge = parameters[0]
        idv_data_url = parameters[1]
        idvWAge_data_fields = parameters[2]
        idvWOAge_data_fields = parameters[3]
        pop_data_url = parameters[4]
        popWAge_data_fields = parameters[5]
        popWOAge_data_fields = parameters[6]
        ftr_url = parameters[7]
        ftr_fields = parameters[8]
        out_table = parameters[9]

        # Get enabled fields
        idv_data_fields = pop_data_fields = None
        if byAge.value:
            idv_data_fields = idvWAge_data_fields
            pop_data_fields = popWAge_data_fields
        else:
            idv_data_fields = idvWOAge_data_fields
            pop_data_fields = popWOAge_data_fields
        
        idv_fields_name = helpers.get_valueTableNames(idv_data_fields)
        idv_region_name = idv_fields_name[0]
        idv_age_name = idv_fields_name[1] if byAge.value else None
        pop_fields_name = helpers.get_valueTableNames(pop_data_fields)
        pop_region_name = pop_fields_name[0]
        pop_pop_name = pop_fields_name[1]
        pop_ageGrp_name = pop_fields_name[2] if byAge.value else None
        ftr_fields_name = helpers.get_valueTableNames(ftr_fields)
        ftr_region_name = ftr_fields_name[0]

        idv_data = helpers.get_pandas(idv_data_url.valueAsText, idv_fields_name)
        pop_data = helpers.get_pandas(pop_data_url.valueAsText, pop_fields_name)
        ftr_data = helpers.get_pandas(ftr_url.valueAsText, ftr_fields_name)

        if byAge.value:
            idv_data["AgeGroup"] = idv_data[idv_age_name].apply(helpers.categorize_age)
            idv_data_groups = [idv_region_name, "AgeGroup"]
            pop_data_groups = [pop_region_name, pop_ageGrp_name]
            output_cols = ["RegionID", "AgeGroup", "EventCount", "PopCount"]
        else:
            idv_data_groups = [idv_region_name]
            pop_data_groups = [pop_region_name]
            output_cols = ["RegionID", "EventCount", "PopCount"]

            if pop_data[pop_region_name].nunique() != len(pop_data.index):  
                pop_data = pop_data.groupby(pop_region_name).agg({pop_pop_name : 'sum'})
                messages.addWarningMessage("Repeated Input Population Data Region IDs were detected. Population Counts were aggregated to totals.")

        # If there are regions in Feature not in Population Data, set population to 0 and warn
        if set(pop_data[pop_region_name].unique()) != set(ftr_data[ftr_region_name].unique()):
            index_names = [pop_region_name]
            if byAge.value:
                index_names.append(pop_ageGrp_name)
                index_vals = [ftr_data[ftr_region_name].unique(), pop_data[pop_ageGrp_name].unique()]
                index = pd.MultiIndex.from_product(index_vals, names = index_names)
            else:
                index = ftr_data[ftr_region_name].unique().tolist()
            pop_data = pop_data.set_index(index_names).reindex(index, fill_value=0).reset_index()
            messages.addWarningMessage("Input Feature Region ID Field contains at least one value not present in Input Population Data Region ID Field. Population within these regions were assumed to be 0.")
        
        event_data = idv_data.groupby(idv_data_groups, as_index= False).size()
        event_data = event_data.merge(pop_data, 
            right_on= pop_data_groups, 
            left_on= idv_data_groups, 
            how = "right")
        event_data["EventCount"] = event_data["size"].fillna(value=0).astype(int)
        event_data = event_data[pop_data_groups + ["EventCount", pop_pop_name]]

        # Check if Event Count is always less than Population count
        if any(event_data["EventCount"] > event_data[pop_pop_name]):
            messages.addWarningMessage("Event Count is greater than the Population Count for at least one row.")

        output_np = np.rec.fromrecords(event_data, names = output_cols)
        arcpy.da.NumPyArrayToTable(output_np, out_table.valueAsText)

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""

        return