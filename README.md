Analysis to recreate the rumen virome manuscript.

The analysis is separated into several RMarkdown files.
You can render them using the instruction below, in the order below.
Generally, later RMarkdown files rely on outputs generated earlier.

	1. setup.Rmd
	2. get_data.Rmd
	3. qc.Rmd
	4. assembly.Rmd
	5. beta_div.Rmd
	6. circos_plots.Rmd
	7. annotation.Rmd
	8. enzyme_analysis.Rmd
	
- git clone
- R
- library(knitr)
- library(rmarkdown)
- render("Rmd_file", "all")

Versions already rendered can be viewed in this repository with suffix .md or .html.

