# InTiCAR

![Figure_1_for_Github](https://github.com/user-attachments/assets/5964d5aa-7665-40fc-be30-f144db5db6f2)
</br>
We developed a network analysis method to find Inter-Tissue Communicators for Autoimmune diseases by Random walk with restart (InTiCAR).
This repository shares a executable tool to run InTiCAR for the user's choice of network, inter-tissue communicators, and genes of interest (e.g. disease genes).

## Systems Requirements

Python version

- 'python' == 3.8.18

Required packages

- 'pandas' == 1.5.3
- 'scipy' == 1.10.1
- 'numpy' == 1.24.3
- 'networkx' == 3.1

## Input Files Descriptions

The following files are required to run InTiCAR for the interested researchers: </br>
The delimeters for all files are '|'. </br>
The example files are provided in the directory "prep_files". </br>
You can simply prepare your own version of files as needed, and add their file names as arguments when running InTiCAR.

The input arguments are the following:

- <strong>'background_network' (-b)</strong> </br> The network to run the analysis from. All the gene pairs need to be in ENSG IDs. </br> Please note that the name of this file will be used as the name of the network throughout the analysis. </br> Refer to 'example_network.csv' for the required format. </br> Default set as '../prep_files/example_network.csv'

- <strong>disease_genes_of_interest (-g)</strong> </br> A list of genes that you would like to search the related ITCs for. </br> Please note that the name of this file will be used as the name of the genes-of-interest throughout the analysis. </br> Refer to 'example_GOI.csv' for the required format. </br> Default set as '../prep_files/example_GOI.csv'

- <strong>disease_genes_full_collection (-d)</strong> </br> A dataframe with the prior knowledge-based disease genes. </br> We recommend using the file provided (i.e. "DiseaseGenes_parsed.csv") with the prior knowledge-based genes for 265 diseases. </br> This table will provide the reference to compare when calculating the modified Z score for your genes-of-interest. </br> Refer to 'DiseaseGenes_parsed.csv' for the required format. </br> Default set as '../prep_files/DiseaseGenes_parsed.csv'

- <strong>modified_z_threshold (-t)</strong> </br> The threshold to use for the modified z-score to find ITCs specific to your genes of interest. </br> Default set as 5.

- <strong>parallel_num (-p)</strong> </br> A number of cores to use for the parallel processing. </br> Parallel processing is highly recommended. </br> Default set as 50.
  </br>

## The example command to run InTiCAR: </br>

(At the ./src directory...)
./run_inticar.py -b {dir_for_network}/my_network.csv -g {dir_for_gene_of_interest}/my_GOI.csv -p 20

Currently, the whole process can take days depending on the size of the given network. There may be updates to optimize the tool further. In the mean time, we made sure to report the time for each step while the code is running, so you can check which point of the process you are in.
