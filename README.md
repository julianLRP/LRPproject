
# Load the LRP dataset
# Create a Seurat object
# Perform preprocessing steps
# Perform dimensionality reduction and clustering
# Identify marker genes
# Annotate the clusters
# Predict the pseudotime. using
Step 2. Generation of realistic layout
# Load lateral root image as a spatial object using the sf package
#Convert to binary image
#Define the cells as connected components
#Label the connected components
# Run tests on coloring the layout with random datasets.

Step 3. Realistic vizualizaton
# Load your LRP data as a data frame
# Identify marker genes for LRP using Seurat - (THIS has been done on step 1)
# Merge the gene expression data with the spatial object using the cluster names as a common identifier
#Set the pseudotime values to map different developmental stages (write the rules explicitly) 
# Join the gene data to the lateral root spatial object
# Plot the lateral root image with marker genes overlaid
# Plot clustering visualization (Umap tSNE or both?)        

Experimental steps to confirm data(Immunostaining/ In situ hybridization of LRP marker genes?, or other experiments?)
Confocal analysis of existing reporter lines
Design and generation(?) of new reporter marker lines.

Optional /additional 1: alter project parameters based on findings

(Optional 2)
#using shiny to create a UI for users to interact with the data and visualizations by either maintaining existing interface or updating and improving upon it.
That is a feasible temporal solution.

(Optional 3)
#Create API to interact with this UI in order to update datasets and visualize other Seurat accessible plots on existing UI.( potentially allow user to upload own datasets using the UI? Or would this be too complicated/ unnecessary?)
This is a better way of doing things, but you might not be in time with it.


