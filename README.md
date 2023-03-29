Setup gitversion control for sharing scripts, potentially also possible to share datasets using gitversion control and separate data storage facility not needed? ( Remote repository on cloud, or on local server? If Cloud which cloud service will be used for this? Paid service/ free service?)
For scripts – a high priority task, for datasets – a low priotity.

Find and create image of LPR + all stages to be used included to be processed.
Take one of those in published papers shared by Maria (propose a few examples this week)

Make the R project:
Packages that will be used:
SF, Seurat, GGplot2, package for reading and processing of png’s( choices: PNG/Imager/magick ) choose the least amount of packages in order to achieve project goal.
In the end of your project, we might develop a package prototype for single cell data spatial visualization

Note that some of these packages may require additional system-level dependencies to be installed (Find this out and document the findings obtained from R studios while writing)
You usually find it when the code does not work 😊
Step 1. Single cell data preprocessing
# Load the LRP dataset
# Create a Seurat object
# Perform preprocessing steps
# Perform dimensionality reduction and clustering
# Identify marker genes
# Annotate the clusters
# Predict the pseudotime. It is using Slingshot:
https://bioconductor.org/packages/release/bioc/html/slingshot.html

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

Make presentation explaining the project:
-	General explanation of transcriptomics
Starts with single cell transcriptomics right away
-	Very brief Explanation of R, programming basics and the use of these visualization techniques to the researchers.
No need for that, but there is need in introduction of Seurat and standard analysis pipeline.

-	Explanation of LRP, benefit of integrating this in the Root Cell Atlas.
-	Project layout and Experimental design + explanation of verification techniques.
-	Future research regarding project scalability/optimization/ increased functionality of tool.
-	Implications for the field and future research uses if proven to be successful.

