# -*- coding: utf-8 -*-
import arcpy
import numpy as np
import pandas as pd

import arcpy_extras
import census
import model


class Toolbox:
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Rate Stabilizing Toolbox"
        self.alias = "RST"

        # List of tool classes associated with this toolbox
        self.tools = [RST, IDP, CDR]


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
        
        params_add = arcpy.Parameter(
            displayName="Additional Options",
            name="AdditionalOptions",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input"
        )
        params_add.columns = [['Double', 'Credible Level'], ['Long', 'Rate Per'], ['Long', 'Number of Years']]
        params_add.filters[0].type = "ValueList"
        params_add.filters[0].list = [0.90, 0.95, 0.99]
        params_add.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'

        params = [
            param_data_table,
            param_data_fields,
            param_feature,
            param_feature_field,
            params_add,
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
        additional_opt = parameters[4]
        estimates_out = parameters[5]
        data_ageGrp_id = parameters[6]
        std_pop_yr = parameters[7]
        age_std_groups = parameters[8]

        # Restrict age groups to only those present within data
        data_ageGrp_info  = arcpy_extras.get_fieldInfo(data_url.valueAsText, data_ageGrp_id.valueAsText)
        if data_ageGrp_info.exists:
            ageGrp_unique = list(set(data_ageGrp_info.list))

            grp_notValid = [group for group in ageGrp_unique if group not in census.constants.age_grps]

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
        additional_opt = parameters[4]
        estimates_out = parameters[5]
        data_ageGrp_id = parameters[6]
        std_pop_yr = parameters[7]
        age_std_groups = parameters[8]

        data_fields_name = arcpy_extras.get_valueTableValues(data_fields)[0]
        data_region_info = arcpy_extras.get_fieldInfo(data_url.valueAsText, data_fields_name[0])
        data_event_info = arcpy_extras.get_fieldInfo(data_url.valueAsText, data_fields_name[1])
        data_pop_info  = arcpy_extras.get_fieldInfo(data_url.valueAsText, data_fields_name[2])
        data_ageGrp_info  = arcpy_extras.get_fieldInfo(data_url.valueAsText, data_ageGrp_id.valueAsText)
        feature_fields_name = arcpy_extras.get_valueTableValues(feature_fields)[0]
        feature_region_info = arcpy_extras.get_fieldInfo(feature_url.valueAsText, feature_fields_name[0])
        age_std_groups_vals = arcpy_extras.get_valueTableValues(age_std_groups)

        # Check if all fields are filled in for age standardization
        if data_ageGrp_info.name is None and (std_pop_yr.valueAsText is not None or age_std_groups.valueAsText is not None):
            arcpy_extras.set_parameterRequired(data_ageGrp_id)
        if std_pop_yr.valueAsText is None and (data_ageGrp_info.name is not None or age_std_groups.valueAsText is not None):
            arcpy_extras.set_parameterRequired(std_pop_yr)
        if age_std_groups.valueAsText is None and (data_ageGrp_info.name is not None or std_pop_yr.valueAsText is not None):
            arcpy_extras.set_valueTableRequired(age_std_groups)

        # Check for Nulls
        if data_region_info.exists and None in data_region_info.list:
            err = "Input Table Region ID Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(data_region_info.list) if elem is None]) + "."
            data_fields.setErrorMessage(err)
        if data_event_info.exists and None in data_event_info.list:
            err = "Input Table Event Count Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(data_event_info.list) if elem is None]) + "."
            data_fields.setErrorMessage(err)
        if data_pop_info.exists and None in data_pop_info.list:
            err = "Input Table Population Count Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(data_pop_info.list) if elem is None]) + "."
            data_fields.setErrorMessage(err)
        if data_ageGrp_info.exists and None in data_ageGrp_info.list:
            err = "Age Group Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(data_ageGrp_info.list) if elem is None]) + "."
            data_fields.setErrorMessage(err)
        if feature_region_info.exists and None in feature_region_info.list:
            err = "Input Feature Region ID Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(feature_region_info.list) if elem is None]) + "."
            data_fields.setErrorMessage(err)
        
        # Check for data types
        if data_region_info.exists and data_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            data_fields.setErrorMessage("Input Table Region ID Field must be of type Integer or String.")
        if data_event_info.exists and data_event_info.type not in ["SmallInteger", "Integer", "BigInteger"]:
            data_fields.setErrorMessage("Input Table Event Count Field must be of type Integer.")
        if data_pop_info.exists and data_pop_info.type not in ["SmallInteger", "Integer", "BigInteger"]:
            data_fields.setErrorMessage("Input Table Population Count Field must be of type Integer.")
        if data_ageGrp_info.exists and data_ageGrp_info.type not in ["String"]:
            data_ageGrp_id.setErrorMessage("Age Group Field must be of type String.")
        if feature_region_info.exists and feature_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            feature_fields.setErrorMessage("Input Feature Region ID Field must be of type Integer or String.")

        # Check if Output Table exists
        if arcpy_extras.exists(estimates_out.valueAsText):
            estimates_out.setErrorMessage("Output Table already exists")

        # Check if Region IDs are repeated when age is adjusted for
        if (not data_ageGrp_info.exists and data_region_info.exists and 
            not data_fields.hasError() and 
            len(data_region_info.list) != len(set(data_region_info.list))):
            data_fields.setWarningMessage("Repeated Input Table Region IDs are detected. Population Counts will be aggregated to totals.")

        # Warn if Event Count >= Population Count
        if data_pop_info.exists and data_event_info.exists and not data_fields.hasError():
            equal100_rate_rows = [i for i, nevents in enumerate(data_event_info.list) if nevents == data_pop_info.list[i]]
            if equal100_rate_rows:
                warn = "Input Table Event Count is equal to Input Table Population Count at "
                warn += arcpy_extras.row_string(equal100_rate_rows) + "."
                data_fields.setWarningMessage(warn)
            above100_rate_rows = [i for i, nevents in enumerate(data_event_info.list) if nevents > data_pop_info.list[i]]
            if above100_rate_rows:
                err = "Input Table Event Count is greater than Input Table Population Count at "
                err += arcpy_extras.row_string(above100_rate_rows) + "."
                data_fields.setErrorMessage(err)

        if (data_region_info.exists and feature_region_info.exists and 
            not data_fields.hasError() and not feature_fields.hasError()):
            # Check if Region ID types are the same
            if data_region_info.type != feature_region_info.type:
                feature_fields.setErrorMessage("Input Feature Region ID Field type must match Input Table Region ID Field type.")
            else:
                # Check if Data contains Region IDs not present in Feature
                data_only_regions = set(data_region_info.list) - set(feature_region_info.list)
                if len(data_only_regions) != 0:
                    err = "Input Table Region ID Field must only contain values present in Input Feature Region ID Field. See "
                    err += arcpy_extras.row_string([i for i, elem in enumerate(data_region_info.list) if elem not in feature_region_info.list]) + "."
                    data_fields.setErrorMessage(err)
                # Check if Feature contains Region IDs not present in Data
                feature_only_regions = set(feature_region_info.list) - set(data_region_info.list)
                if len(feature_only_regions) != 0:
                    err = "Input Feature Region ID Field must only contain values present in Input Table Region ID Field. See "
                    err += arcpy_extras.row_string([i for i, elem in enumerate(feature_region_info.list) if elem not in data_region_info.list]) + "."
                    feature_fields.setErrorMessage(err)
            
        if data_ageGrp_info.exists and not data_ageGrp_id.hasError():
            ageGrp_unique = list(set(data_ageGrp_info.list))
            ageGrp_invalid = [group for group in ageGrp_unique if group not in census.constants.age_grps]
            # Check if there are any invalid age groups
            if len(ageGrp_invalid) != 0:
                err = "Age Group Field must not contain an invalid age group. See "
                err += arcpy_extras.row_string([i for i, elem in enumerate(data_ageGrp_info.list) if elem not in census.constants.age_grps]) + "."
                data_ageGrp_id.setErrorMessage(err)
            elif data_region_info.exists:
                ageGrp_dict = {}
                for i, region in enumerate(data_region_info.list):
                    ageGrp_dict.setdefault(region, []).append(data_ageGrp_info.list[i])
                regions_ageGrp_nunique = [len(set(ageGrp_dict[region])) for region in ageGrp_dict]
                regions_ageGrp_n = [len(ageGrp_dict[region]) for region in ageGrp_dict]
                # Check if there duplicate age groups
                if regions_ageGrp_nunique != regions_ageGrp_n:
                    data_ageGrp_id.setErrorMessage("All Regions must not have any duplicate Age Groups.")
                # Check if there missing age groups
                elif len(set(regions_ageGrp_n)) != 1:
                    data_ageGrp_id.setErrorMessage("All Regions must have the same Age Groups.")

        # Check for duplicate standardized age groups
        age_std_groups_vals = [tuple(group) for group in age_std_groups_vals]
        if len(set(age_std_groups_vals)) != len(age_std_groups_vals):
            age_std_groups.setErrorMessage("Standardized Age Groups must not be repeated.")
            
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        messages.AddMessage("Preparing data...")

        data_url = parameters[0]
        data_fields = parameters[1]
        feature_url = parameters[2]
        feature_fields = parameters[3]
        additional_opt = parameters[4]
        estimates_out = parameters[5]
        data_ageGrp_id = parameters[6]
        std_pop_yr = parameters[7]
        age_std_groups = parameters[8]
        
        data_fields_name = arcpy_extras.get_valueTableValues(data_fields)[0]
        if data_ageGrp_id.valueAsText is not None:
            data_fields_name.append(data_ageGrp_id.valueAsText)
        data_region_name = data_fields_name[0]
        data_event_name = data_fields_name[1]
        data_pop_name = data_fields_name[2]
        feature_region_name = str(feature_fields.values[0][0])

        ci_pct, rates_per, n_years = additional_opt.values
        ci_pct = float(ci_pct[1])
        rates_per = int(rates_per[1])
        n_years = int(n_years[1])

        # Get the age group distribution
        age_std_groups_arr = []
        age_std_groups_names = []
        if age_std_groups.values is not None:
            for lv, uv, in age_std_groups.values:
                lv_index = [i for i, grp in enumerate(census.constants.age_grps) if grp.startswith(lv)][0]
                uv_index = [i for i, grp in enumerate(census.constants.age_grps) if grp.endswith(uv)][0]
                age_std_groups_arr.append(census.constants.age_grps[lv_index:(uv_index + 1)])
                age_std_groups_names.append( lv + ("to" + uv if uv != "85up" else "up") )

        # Read in data
        data = arcpy_extras.get_pandas(data_url.valueAsText, data_fields_name)
        age_groups = [""]
        num_group = 1

        # Warn if Event Count == Population Count
        equal100_rate_rows = np.where(data[data_event_name] == data[data_pop_name])[0].tolist()
        if equal100_rate_rows:
            warn = "Input Table Event Count is equal to Input Table Population Count at "
            warn += arcpy_extras.row_string(equal100_rate_rows) + "."
            messages.addWarningMessage(warn)

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
            std_pop = std_pop[np.isin(census.constants.age_grps, age_groups)]

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

        # Error if a region has no adjacent regions 
        no_adjregions = [i for i, adj_regions in enumerate(adj) if len(adj_regions) == 0]
        if no_adjregions:
            err = "Each Input Feature Region must have at least one adjacent region. See "
            err += arcpy_extras.row_string(no_adjregions) + ". If using Census geographies, try using TIGER/Line Boundaries."
            messages.addErrorMessage(err)
            return

        # Generate estimates
        messages.AddMessage("Generating estimates...")
        theta_out = model.runner.gibbs_rucar(Y, n, adj, std_pop)
        output = model.param_updates.expit(theta_out) * rates_per / n_years

        # If age standardized, generate age_groups
        if age_std_groups.values is not None:
            for ages in age_std_groups_arr:
                output = model.runner.age_std(output, age_groups, std_pop, ages)
            age_groups.extend(age_std_groups_names)
        
        medians, ci_lo, ci_hi, ci_chart = model.runner.get_medians(output, regions, age_groups, ci_pct)

        # If age standardized, combine output with prefixes, else just rename the output cols
        output = output_cols = []
        if age_std_groups.values is not None:
            medians = medians.add_prefix("median_")
            ci_lo = ci_lo.add_prefix(f"CI{int(100 * ci_pct)}lower_")
            ci_hi = ci_hi.add_prefix(f"CI{int(100 * ci_pct)}upper_")
            ci_chart = ci_chart.add_prefix("maxCI_")
            output = pd.concat([medians, ci_lo, ci_hi, ci_chart], axis = 1)
            for i, med in enumerate(medians.columns):
                output_cols.extend([ medians.columns[i], ci_lo.columns[i], ci_hi.columns[i], ci_chart.columns[i] ])
            output = output[output_cols].reset_index()
        else:
            output = pd.concat([medians, ci_lo, ci_hi, ci_chart], axis = 1).reset_index()
            output_cols = ["median", f"CI{int(100 * ci_pct)}lower", f"CI{int(100 * ci_pct)}upper", "maxCI"]
        output.columns = [data_region_name] + output_cols

        arcpy_extras.pandas_to_table(output, estimates_out.valueAsText)
    
        messages.AddMessage("Model finished!")

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class IDP:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Individual Data Processor"
        self.description = ""

    def getParameterInfo(self):
        """Define the tool parameters."""

        param_byAge = arcpy.Parameter(
            displayName="Age Stratified",
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

        idvWAge_fields_name = arcpy_extras.get_valueTableValues(idvWAge_data_fields)[0]
        idvWOAge_fields_name = arcpy_extras.get_valueTableValues(idvWOAge_data_fields)[0]
        popWAge_fields_name = arcpy_extras.get_valueTableValues(popWAge_data_fields)[0]
        popWOAge_fields_name = arcpy_extras.get_valueTableValues(popWOAge_data_fields)[0]

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

        idv_fields_name = arcpy_extras.get_valueTableValues(idv_data_fields)[0]
        idv_region_info = arcpy_extras.get_fieldInfo(idv_data_url.valueAsText, idv_fields_name[0])
        idv_age_info = arcpy_extras.get_fieldInfo(idv_data_url.valueAsText, idv_fields_name[1] if byAge.value else None)
        pop_fields_name = arcpy_extras.get_valueTableValues(pop_data_fields)[0]
        pop_region_info = arcpy_extras.get_fieldInfo(pop_data_url.valueAsText, pop_fields_name[0])
        pop_pop_info = arcpy_extras.get_fieldInfo(pop_data_url.valueAsText, pop_fields_name[1])
        pop_ageGrp_info = arcpy_extras.get_fieldInfo(pop_data_url.valueAsText, pop_fields_name[2] if byAge.value else None)
        ftr_fields_name = arcpy_extras.get_valueTableValues(ftr_fields)[0]
        ftr_region_info = arcpy_extras.get_fieldInfo(ftr_url.valueAsText, ftr_fields_name[0])
        
        # Make field parameters required based on off byAge
        if byAge.value:
            arcpy_extras.set_valueTableRequired(idvWAge_data_fields)
            arcpy_extras.set_valueTableRequired(popWAge_data_fields)
        else:
            arcpy_extras.set_valueTableRequired(idvWOAge_data_fields)
            arcpy_extras.set_valueTableRequired(popWOAge_data_fields)

        # Check for Nulls
        if idv_region_info.exists and None in idv_region_info.list:
            err = "Input Individual Data Region ID Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(idv_region_info.list) if elem is None]) + "."
            idv_data_fields.setErrorMessage(err)
        if byAge.value and idv_age_info.exists and None in idv_age_info.list:
            err = "Input Individual Data Age Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(idv_age_info.list) if elem is None]) + "."
            idv_data_fields.setErrorMessage(err)
        if pop_region_info.exists and None in pop_region_info.list:
            err = "Input Population Data Region ID Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(pop_region_info.list) if elem is None]) + "."
            pop_data_fields.setErrorMessage(err)
        if pop_pop_info.exists and None in pop_pop_info.list:
            err = "Input Population Data Population Count Field must not contain Null values. See  "
            err += arcpy_extras.row_string([i for i, elem in enumerate(pop_pop_info.list) if elem is None]) + "."
            pop_data_fields.setErrorMessage(err)
        if byAge.value and pop_ageGrp_info.exists and None in pop_ageGrp_info.list:
            err = "Input Population Data Age Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(pop_ageGrp_info.list) if elem is None]) + "."
            pop_data_fields.setErrorMessage(err)
        if ftr_region_info.exists and None in ftr_region_info.list:
            err = "Input Feature Region ID Field must not contain Null values. See "
            err += arcpy_extras.row_string([i for i, elem in enumerate(ftr_region_info.list) if elem is None]) + "."
            ftr_fields.setErrorMessage(err)
        
        # Check for data types
        if idv_region_info.exists and idv_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            idv_data_fields.setErrorMessage("Input Individual Data Region ID Field must be of type Integer or String.")
        if byAge.value and idv_age_info.exists and idv_age_info.type not in ["SmallInteger", "Integer", "BigInteger", "Double"]:
            idv_data_fields.setErrorMessage("Input Individual Data Age Field must be of type Integer or Double.")
        if pop_region_info.exists and pop_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            pop_data_fields.setErrorMessage("Input Population Data Region ID Field must be of type Integer or String.")
        if pop_pop_info.exists and pop_pop_info.type not in ["SmallInteger", "Integer", "BigInteger"]:
            pop_data_fields.setErrorMessage("Input Population Data Population Count Field must be of type Integer.")
        if byAge.value and pop_ageGrp_info.exists and pop_ageGrp_info.type not in ["String"]:
            pop_data_fields.setErrorMessage("Input Population Data Age Group Field must be of type String.")
        if ftr_region_info.exists and ftr_region_info.type not in ["SmallInteger", "Integer", "BigInteger", "String"]:
            idv_data_fields.setErrorMessage("Input Feature Region ID Field must be of type Integer or String.")
        
        # Check if Individual Age is negative
        if byAge.value and idv_age_info.exists and not idv_data_fields.hasError() and any(age < 0 for age in idv_age_info.list):
            err = "Input Individual Data Age Field contains a negative value at "
            err += arcpy_extras.row_string([i for i, elem in enumerate(idv_age_info.list) if elem < 0]) + "."
            pop_data_fields.setErrorMessage(err)

        # Check if Output Table exists
        if arcpy_extras.exists(out_table.valueAsText):
            out_table.setErrorMessage("Output Table must not already exist.")

        # Check if repeated data regions
        if not byAge.value and pop_region_info.exists and len(pop_region_info.list) != len(set(pop_region_info.list)):
            pop_data_fields.setWarningMessage("Repeated Input Population Data Region IDs are detected. Population Counts will be aggregated to totals if no age standardization is applied.")

        # Check Individual and Population data relationships
        if (idv_region_info.exists and pop_region_info.exists and 
            not idv_data_fields.hasError() and not pop_data_fields.hasError()):
            # Check if Region ID types are the same
            if idv_region_info.type != pop_region_info.type:
                pop_data_fields.setErrorMessage("Input Population Data Region ID Field type must match Input Individual Data Region ID Field type.")
            else:
                # Check if Individual Data contains Region IDs not present in Population Data
                idv_only_regions = set(idv_region_info.list) - set(pop_region_info.list)
                if len(idv_only_regions) != 0:
                    err = "Input Individual Data Region ID Field must not contain a value not present in Input Population Data Region ID Field. See "
                    err += arcpy_extras.row_string([i for i, elem in enumerate(idv_region_info.list) if elem not in pop_region_info.list]) + "."
                    pop_data_fields.setErrorMessage(err)

        # Check Population Data and Feature relationships
        if (pop_region_info.exists and ftr_region_info.exists and
            not pop_data_fields.hasError() and not ftr_fields.hasError()):
            # Check if Region ID types are the same
            if pop_region_info.type != ftr_region_info.type:
                ftr_fields.setErrorMessage("Input Population Data Region ID Field must match Input Feature Region ID Field type.")
            else:
                # Check if Population Data contains Region IDs not present in Feature
                pop_only_regions = set(pop_region_info.list) - set(ftr_region_info.list)
                if len(pop_only_regions) != 0:
                    err = "Input Population Data Region ID Field must only contain values present in Input Feature Region ID Field. See "
                    err += arcpy_extras.row_string([i for i, elem in enumerate(pop_region_info.list) if elem not in ftr_region_info.list]) + "."
                    pop_data_fields.setErrorMessage(err)
                # Check if Feature contains Region IDs not present in Population Data
                ftr_only_regions = set(ftr_region_info.list) - set(pop_region_info.list)
                if len(ftr_only_regions) != 0:
                    err = "Input Feature Region ID Field contains a value not present in Input Population Data Region ID Field at "
                    err += arcpy_extras.row_string([i for i, elem in enumerate(ftr_region_info.list) if elem not in pop_region_info.list]) + "."
                    err += "\n\nPopulation within these regions will be assumed to be 0."
                    ftr_fields.setWarningMessage(err)

        # Check Individual Data and Feature relationships
        if (idv_region_info.exists and ftr_region_info.exists and 
            not idv_data_fields.hasError() and not ftr_fields.hasError()):
            if idv_region_info.type == ftr_region_info.type:
                # Check if Individual Data contains Region IDs not present in Population Data
                idv_only_regions = set(idv_region_info.list) - set(ftr_region_info.list)
                if len(idv_only_regions) != 0:
                    err = "Input Individual Data Region ID Field must only contain values present in Feature Region ID Field. See "
                    err += arcpy_extras.row_string([i for i, elem in enumerate(idv_region_info.list) if elem not in ftr_region_info.list]) + "."
                    idv_data_fields.setErrorMessage(err)
                    
        if pop_ageGrp_info.exists and not pop_data_fields.hasError():
            ageGrp_unique = list(set(pop_ageGrp_info.list))
            ageGrp_invalid = [group for group in ageGrp_unique if group not in census.constants.age_grps]
            # Check if there are any invalid age groups
            if len(ageGrp_invalid) != 0:
                err = "Age Group Field must only contain valid age groups. See "
                err += arcpy_extras.row_string([i for i, elem in enumerate(pop_ageGrp_info.list) if elem not in census.constants.age_grps]) + "."
                pop_data_fields.setErrorMessage(err)
            elif pop_region_info.exists:
                ageGrp_dict = {}
                for i, region in enumerate(pop_region_info.list):
                    ageGrp_dict.setdefault(region, []).append(pop_ageGrp_info.list[i])
                regions_ageGrp_nunique = [len(set(ageGrp_dict[region])) for region in ageGrp_dict]
                regions_ageGrp_n = [len(ageGrp_dict[region]) for region in ageGrp_dict]
                # Check if there duplicate age groups
                if regions_ageGrp_nunique != regions_ageGrp_n:
                    pop_data_fields.setErrorMessage("All Regions must not have any duplicate Age Groups.")
                # Check if there missing age groups
                elif len(set(regions_ageGrp_n)) != 1:
                    pop_data_fields.setErrorMessage("All Regions must have the same Age Groups.")

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
        
        idv_fields_name = arcpy_extras.get_valueTableValues(idv_data_fields)[0]
        idv_region_name = idv_fields_name[0]
        idv_age_name = idv_fields_name[1] if byAge.value else None
        pop_fields_name = arcpy_extras.get_valueTableValues(pop_data_fields)[0]
        pop_region_name = pop_fields_name[0]
        pop_pop_name = pop_fields_name[1]
        pop_ageGrp_name = pop_fields_name[2] if byAge.value else None
        ftr_fields_name = arcpy_extras.get_valueTableValues(ftr_fields)[0]
        ftr_region_name = ftr_fields_name[0]

        idv_data = arcpy_extras.get_pandas(idv_data_url.valueAsText, idv_fields_name)
        pop_data = arcpy_extras.get_pandas(pop_data_url.valueAsText, pop_fields_name)
        ftr_data = arcpy_extras.get_pandas(ftr_url.valueAsText, ftr_fields_name)

        idv_data = idv_data.rename(columns = {idv_region_name: "GEOID"})
        pop_data = pop_data.rename(columns = {pop_region_name: "GEOID", pop_pop_name: "PopulationCount"})
        ftr_data = ftr_data.rename(columns = {ftr_region_name: "GEOID"})

        if byAge.value:
            def classify_age(age):
                if age <= 4: return "0-4"
                elif age <= 14: return "5-14"
                elif age <= 24: return "15-24"
                elif age <= 34: return "25-34"
                elif age <= 44: return "35-44"
                elif age <= 54: return "45-54"
                elif age <= 64: return "55-64"
                elif age <= 74: return "65-74"
                elif age <= 84: return "75-84"
                else: return "85up"

            pop_data = pop_data.rename(columns = {pop_ageGrp_name: "AgeGroup"})
            idv_data = (
                idv_data
                .assign(AgeGroup = lambda x: x[idv_age_name].apply(classify_age))                
                .drop(idv_age_name, axis = 1)
            )
            data_groups = ["GEOID", "AgeGroup"]
        else:
            data_groups = ["GEOID"]

            if pop_data["GEOID"].nunique() != len(pop_data.index):  
                pop_data = pop_data.groupby("GEOID").agg({"PopulationCount" : 'sum'}).reset_index()
                messages.addWarningMessage("Repeated Input Population Data Region IDs were detected. Population Counts were aggregated to totals.")

        # If there are regions in Feature not in Population Data, set population to 0 and warn
        if set(pop_data["GEOID"].unique()) != set(ftr_data["GEOID"].unique()):
            index_names = ["GEOID"]
            if byAge.value:
                index_names.append("AgeGroup")
                index_vals = [ftr_data["GEOID"].unique(), pop_data["AgeGroup"].unique()]
                index = pd.MultiIndex.from_product(index_vals, names = index_names)
            else:
                index = ftr_data["GEOID"].unique().tolist()
            pop_data = pop_data.set_index(index_names).reindex(index, fill_value=0).reset_index()
            messages.addWarningMessage("Input Feature Region ID Field contains at least one value not present in Input Population Data Region ID Field. Population within these regions were assumed to be 0.")
        
        event_counts = idv_data.groupby(data_groups).size().reset_index(name="EventCount")
        event_data = (
            pop_data
            .merge(event_counts, on=data_groups, how="left")
            .fillna({"EventCount": 0})
            .assign(EventCount=lambda df: df["EventCount"].astype(int))
            .filter(items = data_groups + ["PopulationCount", "EventCount"])
        )

        # Warn if Event Count >= Population Count
        equal100_rate_rows = np.where(event_data["EventCount"] == event_data["PopulationCount"])[0].tolist()
        if equal100_rate_rows:
            warn = "Event Count is equal to Population Count at "
            warn += arcpy_extras.row_string(equal100_rate_rows) + "."
            messages.addWarningMessage(warn)
        above100_rate_rows = np.where(event_data["EventCount"] > event_data["PopulationCount"])[0].tolist()
        if above100_rate_rows:
            warn = "Event Count is greater than to Population Count at "
            warn += arcpy_extras.row_string(above100_rate_rows) + "."
            warn += "\n\nThe Rate Stabilizing Tool cannot produce reliable rates where the Event Count exceeds the Population Count."
            messages.addWarningMessage(warn)

        arcpy_extras.pandas_to_table(event_data, out_table.valueAsText)

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""

        return
    
class CDR:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Census Data Retriever"
        self.description = ""

    def getParameterInfo(self):
        """Define the tool parameters."""

        param_byAge = arcpy.Parameter(
            displayName="Age Stratified",
            name="ByAge",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input"
        )
        param_byAge.value = True

        param_data_fields = arcpy.Parameter(
            displayName = "Request Parameters",
            name = "InputTableFields",
            datatype = "GPValueTable",
            parameterType = "Required",
            direction = "Input"
        )
        param_data_fields.columns = [['String', 'Survey'], ['String', 'Year'], ['String', 'Geography'], ['String', 'State']]
        param_data_fields.filters[0].type = "ValueList"
        param_data_fields.filters[0].list = ["5-year ACS", "Decennial"]
        param_data_fields.filters[1].type = "ValueList"
        param_data_fields.filters[1].list = census.constants.acs_years
        param_data_fields.filters[2].type = "ValueList"
        param_data_fields.filters[2].list = ["County", "Tract"]
        param_data_fields.filters[3].type = "ValueList"
        param_data_fields.filters[3].list = list(census.constants.state_to_fips.keys())
        param_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'

        param_out_table = arcpy.Parameter(
            displayName="Output Table",
            name="OutputTable",
            datatype="GPTableView",
            parameterType="Required",
            direction="Output"
        )

        param_geom_type = arcpy.Parameter(
            displayName="Geometry Type",
            name="GeometryType",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Geography (optional)"
        )
        param_geom_type.filter.type = "ValueList"
        param_geom_type.filter.list = ["TIGER", "Cartographic"]

        param_out_feature = arcpy.Parameter(
            displayName="Output Feature",
            name="OutputFeature",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Output",
            category="Geography (optional)"
        )

        params = [
            param_byAge,
            param_data_fields,
            param_out_table,
            param_geom_type,
            param_out_feature
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
        req_param = parameters[1]
        out_table = parameters[2]
        geom_type = parameters[3]
        out_feature = parameters[4]

        req_param_val = arcpy_extras.get_valueTableValues(req_param)[0]
        req_survey = req_param_val[0]
        req_year = req_param_val[1]
        req_geography = req_param_val[2]
        req_state = req_param_val[3]
        req_geom_type = geom_type.valueAsText
        req_out_feature = out_feature.valueAsText

        if not req_survey:
            pass
        elif req_survey == "5-year ACS":
            req_param.filters[1].list = census.constants.acs_years
        elif req_survey == "Decennial":
            req_param.filters[1].list = census.constants.dec_years

        # No geometry for 2009
        if req_year == "2005-2009":
            geom_type.enabled = False
            geom_type.value = None
            out_feature.enabled = False
            out_feature.value = None
        else:
            geom_type.enabled = True
            out_feature.enabled = True

        if req_year in census.constants.acs_years and req_year in ["2007-2011", "2008-2012"]:
            geom_type.filter.list = ["TIGER"]
        else:
            geom_type.filter.list = ["TIGER", "Cartographic"]
        
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        byAge = parameters[0]
        req_param = parameters[1]
        out_table = parameters[2]
        geom_type = parameters[3]
        out_feature = parameters[4]

        req_param_val = arcpy_extras.get_valueTableValues(req_param)[0]
        req_survey = req_param_val[0]
        req_year = req_param_val[1]
        req_geography = req_param_val[2]
        req_state = req_param_val[3]
        req_geom_type = geom_type.valueAsText
        req_out_feature = out_feature.valueAsText

        if req_geom_type:
            arcpy_extras.set_parameterRequired(out_feature)
        if req_out_feature:
            arcpy_extras.set_parameterRequired(geom_type)

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        byAge = parameters[0]
        req_param = parameters[1]
        out_table = parameters[2]
        geom_type = parameters[3]
        out_feature = parameters[4]

        req_param_val = arcpy_extras.get_valueTableValues(req_param)[0]
        req_survey = req_param_val[0]
        req_year = req_param_val[1]
        req_geography = req_param_val[2]
        req_state = req_param_val[3]
        req_geom_type = geom_type.valueAsText
        req_out_feature = out_feature.valueAsText
        
        resp_df = census.data.get_census(
            req_survey, 
            req_year, 
            req_geography, 
            req_state, 
            byAge.value
        )
        arcpy_extras.pandas_to_table(resp_df, out_table.valueAsText)

        census.geometry.get_geometry(
            req_geography,
            req_geom_type, 
            req_year, 
            req_state, 
            req_out_feature
        )

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return
