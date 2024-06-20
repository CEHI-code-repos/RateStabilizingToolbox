# Rate Stabilizing Tool

## Install

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
  `TestData.gdb/MI_indiv_data` for an example)
- **Input Individual Data Fields**
  - **Region ID** (Required) - String or Integer identifier that
    uniquely identifies each geographic region in which an event
    occured.
  - **Age** (Required if By Age checked, Optional if By Age unchecked) -
    Integer age of each individual who experienced an event.
- **Input Population Data** - Layer or table view which contains
  population data aggregated by region and if applicable, age group.
  (See `TestData.gdb/MI_popByAge_data` or `TestData.gdb/MI_pop_data` for
  an example)
- **Input Population Data Fields**
  - **Region ID** (Required) - String or Integer identifier that
    uniquely identifies each geographic region in which an event
    occured.
  - **Population Count** (Required) - Integer population within region
    (among valid age groups, if applicable).
  - **Age Group** (Required if By Age checked, Optional if By Age
    unchecked) - String designating regional age group.
- **Output Table** (Required) - Output table location and name.

The following are valid age groups: \[“0-4”, “5-14”, “15-24”, “25-34”,
“35-44”, “45-54”, “55-64”, “65-74”, “75-84”, “85up”\]. When provided age
stratified data age groups must be spelled EXACTLY as written above.

> [!NOTE]
>
> If the Individual Data Processing tool detects that the Input
> Population Data contains multiple rows with the same region when By
> Age is not checked, it will sum the values of Population Count for
> those regions with multiple rows. This allows for age stratified
> population data to be provided to the Individual Data Processing tool
> even when By Age is not checked.
>
> If you encounter the following warning without the intention of using
> age stratified population data to generate a non age stratified output
> table, check your population dataset for duplicate regions.
>
> <img src="warning_IDP.png" data-fig-align="center" />

#### Examples

Non Age Stratified Example:  
<img src="crude_IDP.png" data-fig-align="center" />

Age Stratified Example:  
<img src="ageStd_IDP.png" data-fig-align="center" />

### Rate Stabilizing Tool

The Rate Stabilizing Tool offers an easy-to-use interface for the
restricted CAR model to generate estimates of event rates for geographic
areas with small population sizes or small counts.

#### Inputs

- **Input Table** (Required) - Layer or table view which contains
  incidence and population data aggregated by region and if applicable,
  age group. Typically, this table is the output of the Individual Data
  Processing tool. (See `TestData.gdb/MI_mort_byAgeCounty` or
  `TestData.gdb/MI_mort_SingleGroupByCounty` for an example)
- **Input Table Fields**
  - **Region ID** (Required) - String or integer uniquely identifying
    each geographic area.
  - **Event Count** (Required) - Integer number of events in region
    among total population or if applicable, age group
  - **Population Count** (Required) - Integer population within region
    among the total population or if applicable, among an age group
- **Input Feature** (Required) - Polygon feature of geographic areas
  within the input table.
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

    To generate age standardized rates for individuals 35 to 64, 65 and
    up, and 35 and up:

    | Lower age value | Upper age value |
    |-----------------|-----------------|
    | 35              | 64              |
    | 65              | up              |
    | 35              | up              |

> [!NOTE]
>
> If the Rate Stabilizing Tool detects that the Input Table contains
> multiple rows with the same region when age standardized parameters
> are not provided, it will sum the values of Event Count and Population
> Count for those regions with multiple rows. This allows for age
> stratified data to be provided to the Rate Stabilizing Tool even when
> calculating crude rates.
>
> If you encounter the following warning without the intention of using
> age stratified data to generate crude rates, check your population
> dataset for duplicate regions.
>
> <img src="warning_RST.png" data-fig-align="center" />

#### Examples

Crude Rate Example:  
<img src="crude_RST.png" data-fig-align="center" />

Age Standardized Rates Example:  
<img src="ageStd_RST.png" data-fig-align="center" />
