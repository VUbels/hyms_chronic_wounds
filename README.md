# Information to use the hyms_chronic_wounds project
**First set your working directory to the appropriate project location in CLI:**

cd "your-directory"\
(For example this project would be: cd hyms_chronic_wounds)

**To run proseg pipeline to resegment cells simply run the following**
in CLI and change directory to your base working directory through:

-i: input folder containing the spatial xenium output (Function is recursive so only point to base location)
-o: output-folder will recursively output data to folder within working directory

(No need for function prefix as script is executable through chmod):\
./scripts/run_proseg.sh -i /mnt/d/HYMS/chronic_wounds -o ./proseg_results

**Building scRNA reference dataset to project annotation labels**
If you need to rebuild the complete scRNA reference data run the following:

Rscript ./scripts/prep_CW.R\
Rscript ./scripts/prep_HSCA.R\
Rscript ./scripts/prep_reference.R\
Rscript ./scripts/reference_build.R

##Then to build the seurat objects with appropriate labeling run:

Rscript ./scripts/prepare_objects.R



################################################################################

You can opt to look at the segmentation and troubleshoot yourselves by running 

gene_diagnostic.R or compare_segmentation.R in any IDE (Rstudio for example) and
run through each line individually so you can see/compare plots and data easily.
