# LRP clustering and spatial visualization Tool

## Step 1: Load the LRP Dataset and Perform Preprocessing Steps

- Load the LRP dataset
- Create a Seurat object
- Perform preprocessing steps
- Perform dimensionality reduction and clustering
- Identify marker genes
- Annotate the clusters
- Predict the pseudotime

## Step 2: Generation of Realistic Layout

- Load lateral root image as a spatial object using the sf package
- Convert to binary image
- Define the cells as connected components
- Label the connected components
- Run tests on coloring the layout with random datasets

## Step 3: Realistic Visualization

- Load your LRP data as a data frame
- Identify marker genes for LRP using Seurat 
- Merge the gene expression data with the spatial object using the cluster names as a common identifier
- Set the pseudotime values to map different developmental stages 
- Join the gene data to the lateral root spatial object
- Plot the lateral root image with marker genes overlaid
- Plot clustering visualization        

## Experimental Steps to Confirm Data

- Confocal analysis of existing reporter lines
- Design and generation of new reporter marker lines

## Optional / Additional Steps

1. Alter project parameters based on findings

2. Use shiny to create a UI for users to interact with the data and visualizations by either maintaining existing interface or updating and improving upon it. (This is a feasible temporal solution)

3. Create API to interact with this UI in order to update datasets and visualize other Seurat accessible plots on existing UI. Potentially allow user to upload own datasets using the UI? This is a better way of doing things, but you might not be in time with it.
