Analysis to recreate the rumen virome manuscript.

The analysis is separated into several RMarkdown files.
You can render them using the instruction below, in the order below.
Generally, later RMarkdown files rely on outputs generated earlier.

	1. setup.sh
	2. qc.Rmd
	3. assembly_and_orfs.Rmd
	4. beta_div_viral_mg.Rmd
	5. beta_div_total_mg.Rmd
	6. circos_plots_viral_mg.Rmd
	7. circos_plots_total_mg.Rmd
	8. annotation.Rmd
	9. enzyme_analysis.Rmd
	
- git clone https://github.com/chrisLanderson/rumen_virome.git
- R
- library(knitr)
- library(rmarkdown)
- render("Rmd_file", "all")

Versions already rendered can be viewed in this repository with extension .md or .html.

