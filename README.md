# data-pipeline-nwis

A data pipeline used to retrieve, process and visualize time series data from the National Water Information System (NWIS). 

<img src="./images/example_choropleth_plot.gif?raw=true" width="1920px">

## How to use

This data analysis workflow uses Snakemake (installation instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)) as a pipelining tool to retrieve, process and visualize environmenal time series data from NWIS. The project files are organized following the [conventions](https://github.com/DOI-USGS/ds-pipelines-targets-1-course/blob/main/course-instructions.md) outlined in the USGS data science branch.

First, create a Conda environment with all the required packages by running the following command: `
conda env create -f environment.yaml
`

Once in the new environment, we can execute the snakemake pipeline with this command: `snakemake --cores 1 -s Snakefile.txt`

When the jobs are done, a choropleth plot displaying processed timeseries data and associated metadata will be in a newly created out folder in [3_plot/](3_plot/).

Modifications to the data query can be made via the [Snakefile](Snakefile.txt) by changing the params key values in the `get_data` rule, and by changing the hydrologic unit code (huc) inputs found in [1_fetch/src/in](1_fetch/src/in).
