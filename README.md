# Wastewater genomic surveillance captures early detection of Omicron in Utah

This repository contains the R code used in the analysis and visualization of genomic surveillance data of the Omicron variant of SARS-CoV-2 in Utah. The manuscript associated with this code is titled "Wastewater genomic surveillance captures early detection of Omicron in Utah" and is available at https://www.medrxiv.org/content/10.1101/2022.11.24.22282643v2

## Data preprocessing

The R code reads a raw data file, formats the dates, filters the data based on the study period, and groups the lineages for further analysis. It also reads a corrected lineage abundance data file and applies additional filtering and grouping.

## Analysis

The analysis focuses on summarizing the data, grouping lineages, and creating visualizations. The code includes the creation of a detailed lineage plot for the supplementary materials and a lineage plot for the main text.

## Visualizations

Visualizations created using ggplot2 include:

1. Detailed lineage plot (Supplementary)
2. Bar plot for Delta and Omicron lineages
3. Bar plot for Omicron lineages

## Dependencies

The following R packages are required to run the code:

- dplyr
- tidyverse
- ggplot2
- ggpubr
- ggthemes
- gridExtra
- lubridate
- tibbletime
- scales
- viridis
- patchwork
- broom.mixed

## File Structure

- R code: The main R code file is provided in this repository.

## Usage

To run the code, clone the repository and ensure that the required R packages are installed. Open the R script file and update the file paths to match the location of the data files on your system. Run the script to generate the visualizations and output data file.

## Citation

If you use this code or data, please cite the associated manuscript:

Wastewater genomic surveillance captures early detection of Omicron in Utah
