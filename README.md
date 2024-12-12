# fRank
This R script is designed to process and annotate transcription start site (TSS) data obtained from TSSpredator pipeline. 
TSS Annotation and Ranking Script
This repository contains an R script for annotating transcription start sites (TSS) and ranking them based on enrichment factors and step heights. The script processes data from TSSpredator and gene annotation files to identify and rank primary and secondary TSS across different genomic conditions.
Prerequisites
Before running the script, ensure you have the following R packages installed:
	•	`dplyr`
	•	`tibble`
You can install them using:
install.packages(c("dplyr", "tibble"))

Script Overview
1.	Setup and Data Loading:
•	Clears the R environment.
•	Loads necessary libraries.
•	Sets the working directory.
•	Reads in a CSV file containing cleaned TSS data (`ChlData_Cleaned.csv`).
•	Reads in a gene annotation file (`NC_003028.v3.17.ncrna.genes`).

2.	Data Preparation:
•	Splits the data into top and complementary strands based on `SuperStrand`.
•	Converts relevant columns to numeric types for accurate processing.

3.	Annotation Process:
•	Annotates 5’ UTRs, genic regions, and antisense regions for both strands.
•	Handles edge cases for positions before the first gene.

4.	Primary TSS Ranking:
•	Filters and ranks primary TSS based on enrichment factors and step heights.
•	Separates data into different genomic conditions (`NDCt60`, `Chl1q_t60`, `Chl3q_t60`).
•	Annotates primary TSS with updated information.

5.	Secondary TSS Analysis:
•	Removes primary TSS positions from the dataset.
•	Repeats the annotation process for secondary TSS.
•	Ranks secondary TSS similarly to primary ones.

6.	Output:
•	Combines primary and secondary ranked data.
•	Writes the final annotated and ranked data to a CSV file (`ChlData_fRanked.csv`).

Usage
To run the script, ensure your working directory is set correctly, and that the input files (`ChlData_Cleaned.csv` and `NC_003028.v3.17.ncrna.genes`) are available in the specified directory. You can change these names according to your corresponding TSSpredator file obtained and annotation file, respectively.
Execute the script in R or RStudio to perform the annotation and ranking process. The output will be saved as `ChlData_fRanked.csv`.
Important Notes
	•	The script assumes specific column names in the input files; ensure your files match these expectations.
	•	Adjust paths and filenames as necessary to suit your local setup.
	•	Review comments within the script for detailed explanations of each step.
