subworkflow fluWorkFlow:
	workdir:
		"seasonal-flu-nextstrain"
subworkflow zikaWorkFlow:
	workdir:
		"zika-nextstrain"	
rule all:
	input:
		"docs/FullLinkedChartClickableZika.html",
		"docs/FullLinkedChartClickableFlu.html",
		"docs/cartography.pdf",
		"docs/cartography.html"
rule run_json:
	input:
		fluWorkFlow(expand("auspice/flu_seasonal_h3n2_ha_2y_tree.json), 
		zikaWorkFlow("auspice/zika-cartography_tree.json")
	output:
		"zika-nextstrain/aligned.fasta",
		"zika-nextstrain/auspice/zika-cartography_tree.json",
		"seasonal-flu/aligned_cdc_h3n2_ha_2y_cell_hi.fasta",
		"seasonal-flu/auspice/flu_seasonal_h3n2_ha_2y.json",
		touch("mytask.done")
rule run_zika:
	input:
		"mytask.done",
		"zika-nextstrain/aligned.fasta",
		"notebooks/dropped_Strains_zika.txt",
		"zika-nextstrain/auspice/zika-cartography_tree.json",
		"notebooks/Data/name_of_disease_zika.txt"
	output:
		"docs/FullViolinPlotZika.png",
		"docs/FullScatterplotZika.png",
		"docs/FullLinkedChartClickableZika.pdf",
		"docs/FullLinkedChartClickableZika.html",
	notebook:											   
		"notebooks/2019-08-08FinalNotebookFlu.ipynb"
rule run_flu:
	input:
		"mytask.done",
		"seasonal-flu/aligned_cdc_h3n2_ha_2y_cell_hi.fasta",
		"notebooks/dropped_Strains_flu.txt",
		"seasonal-flu/auspice/flu_seasonal_h3n2_ha_2y.json",
		"notebooks/Data/clade_names.txt",
		"notebooks/Data/name_of_disease_flu.txt"
	output:
		"docs/FullViolinPlotFlu.png",
		"docs/FullScatterplotFlu.png",
		"docs/FullLinkedChartClickableFlu.pdf",
		"docs/FullLinkedChartClickableFlu.html"
	notebook:											   
		"notebooks/2019-08-08FinalNotebookFlu.ipynb"
		
def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
	    return "cartography.html"
    elif suffix == "pdf":
	    return "cartography.pdf"
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")
	
rule report:
	params: options = pandoc_options
	output: "docs/cartography.{suffix}"
	wildcard_constraints: suffix = "((html)|(pdf))"
	shell:
		"""
		pandoc --filter pandoc-citeproc --bibliography=docs/cartography.bib -s docs/index.md -o docs/cartography.{wildcards.suffix}
		"""