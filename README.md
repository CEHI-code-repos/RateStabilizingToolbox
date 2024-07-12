<h1 align="center">
  <a href="https://github.com/CEHI-code-repos/RateStabilizingToolV2">
    Rate Stabilizing Tool
  </a>
</h1>

<p align="center">
  <strong>Generate reliable, local-level age-standardized measures of chronic disease </strong>
</p>

## Set up

1.  Click on the Green Code button and Download Zip
2.  Unzip the zip file
3.  Open ArcGIS Pro. Create a new Map Project. Within the Catalog Pane,
    right click on the Toolboxes and Add Toolbox.

## Tools

### Individual Data Processing

It is often the case that your data is collected at the individual
level. These data will be need to be aggregated to the region (and age
group, for age standardized rates) for rate calculations within the Rate
Stabilizing Tool. This tool was developed to make that process easier!

#### Inputs

- **By Age** (Required) - Specifies if output table should be stratified
  by age group.
- **Input Individual Data** - Layer or table view which contains region
  (and age, if applicable) of individuals who experienced an event. (See
  `tutorial.gdb/MI_mort_indiv` for an example)
- **Input Individual Data Fields**
  - **Region ID** (Required) - String or Integer identifier that
    uniquely identifies each geographic region in which an event
    occured.
  - **Age** (Required if By Age checked, Optional if By Age unchecked) -
    Integer age of each individual who experienced an event.
- **Input Population Data** - Layer or table view which contains
  population data aggregated by region and if applicable, age group.
  (See `tutorial.gdb/MI_pop` or `tutorial.gdb/MI_pop_grouped` for
  an example)
- **Input Population Data Fields**
  - **Region ID** (Required) - String or Integer identifier that
    uniquely identifies each geographic region.
  - **Population Count** (Required) - Integer population within region
    (among valid age groups, if applicable).
  - **Age Group** (Required if By Age checked, Optional if By Age
    unchecked) - String designating regional age group.
- **Input Feature** - Polygon feature of geographic areas
  within analysis. (See `tutorial.gdb/MI_carto` for an example) 
- **Input Feature Fields**
  - **Region ID** (Required) - String or Integer identifier that
    uniquely identifies each geographic region.
- **Output Table** (Required) - Output table location and name.

> [!NOTE]
>
> The following are valid age groups: \[“0-4”, “5-14”, “15-24”, “25-34”,
> “35-44”, “45-54”, “55-64”, “65-74”, “75-84”, “85up”\]. When provided age
> stratified data age groups must be spelled **EXACTLY** as written above.

### Rate Stabilizing Tool

The Rate Stabilizing Tool offers an easy-to-use interface for the
restricted CAR model to generate estimates of event rates for geographic
areas with small population sizes or small counts.

#### Inputs

- **Input Table** (Required) - Layer or table view which contains
  incidence and population data aggregated by region and if applicable,
  age group. Typically, this table is the output of the Individual Data
  Processing tool. (See `tutorial.gdb/MI_mort_pop` or
  `tutorial.gdb/MI_mort_pop_grouped` for an example)
- **Input Table Fields**
  - **Region ID** (Required) - String or integer uniquely identifying
    each geographic area.
  - **Event Count** (Required) - Integer number of events in region
    among total population or if applicable, age group
  - **Population Count** (Required) - Integer population within region
    among the total population or if applicable, among an age group
- **Input Feature** (Required) - Polygon feature of geographic areas
  within the input table. (See `tutorial.gdb/MI_carto` for an example)
- **Input Feature Field**
  - **Region ID** (Required) - String or integer field uniquely
    identifying each geographic area.
- **Rate** (Required) - Integer specifying desired denominator for
  output rates.
- **Output Table** (Required) - Output table location and name.
- **Age Standardization Parameters** (Optional, but must all be filled
  for Age Standardization)
  - **Age Group Field** - String designating regional age group within
    Input Table.
  - **Standardized Population Year** - Standard population year used to
    adjust rate
  - **Standardized Age Groups** - Designates lower and upper bounds of
    age groups that will be standardized.