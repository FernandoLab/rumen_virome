Analysis to recreate the manuscript "Digestible energy content drives the dynamic response of bovine rumen viral communities to dietary change".

The analysis is separated into several RMarkdown files.

You can setup the same environment used to analyze the data and render the markdown files using the instructions below.

Generally, the markdown files rely on previously generate outputs, so rendering in this order likely matters.

	1. data_curation.Rmd
	2. qc_viral.Rmd
	3. qc_microbial.Rmd
	4. qc_deep_viral.Rmd
	5. assembly_and_orfs.Rmd
	6. map_reads.Rmd
	7. annotation.Rmd (due to KEGG licensing, not fully reproducible but outputs available in intermediate_results directory)
	8. abundance_tables.Rmd (due to KEGG licensing, not fully reproducible but outputs available in intermediate_results directory)
	9. beta_div_viral.Rmd
	10. beta_div_microbial.Rmd
	11. circos_plots_viral.Rmd
	12. circos_plots_microbial.Rmd
	13. plsr.Rmd
	14. enzyme_enrichment.Rmd (due to KEGG licensing, not fully reproducible but outputs available in intermediate_results directory)

Due to licensing issues, USEARCH can not be included in the setup. To obtain a download link, go to the USEARCH download page and select version USEARCH v8.0.1623 for linux. A link (expires after 30 days) will be sent to the provided email. Use the link as an argument for shell script below.

Clone the github repository and run the setup.sh script (provide the link to download your licensed USEARCH v8.0.1623 as an argument):

- git clone https://github.com/chrisLanderson/rumen_virome.git
- cd rumen_virome
- ./setup.sh usearch_link

Anaconda is downloaded and prompts you during installataion of the packages. The prompts are roughly as follows:

- Press enter to view the license agreement
- Press enter to read the license and q to exit
- Accept the terms
- Prompts you where to install anaconda. Simply type anaconda to create a directory within the current directory. Should be: [/Users/user/anaconda] >>> anaconda
- No to prepend anaconda to your path. Choosing yes should not impact the installation though.
- Will be asked a few times if you wish to proceed with installing the packages...agree to it.
- After installation, enter 'source anaconda/bin/activate rumenVirome' on the command line to activate the virtual enviornment

To render one of the RMarkdown files:
R CMD BATCH --no-restore --no-save '--args file1.Rmd' render.R
