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
        fluWorkFlow("auspice/flu_seasonal_h3n2_ha_2y_tree.json"), 
        zikaWorkFlow("auspice/zika-cartography_tree.json")
    output:
        "zika-nextstrain/aligned.fasta",
        "zika-nextstrain/auspice/zika-cartography_tree.json",
        "seasonal-flu-nextstrain/aligned.fasta",
        "seasonal-flu/auspice/flu_seasonal_h3n2_ha_2y_tree.json",
rule run_zika:
    input:
        "zika-nextstrain/aligned.fasta",
        "notebooks/dropped_Strains_zika.txt",
        "zika-nextstrain/auspice/zika-cartography_tree.json",
        "notebooks/Data/name_of_disease_zika.txt"
    output:
        "docs/FullViolinPlotZika.png",
        "docs/FullScatterplotZika.png",
        "docs/FullLinkedChartClickableZika.pdf",
        "docs/FullLinkedChartClickableZika.html",
    conda: "cartography.yml"
    notebook:											   
        "notebooks/2019-08-08FinalNotebookFlu.ipynb"
rule run_flu:
    input:
        "seasonal-flu-nextstrain/aligned.fasta",
        "notebooks/dropped_Strains_flu.txt",
        "seasonal-flu-nextstrain/auspice/flu_seasonal_h3n2_ha_2y_tree.json",
        "notebooks/Data/clade_names.txt",
        "notebooks/Data/name_of_disease_flu.txt"
    output:
        "docs/FullViolinPlotFlu.png",
        "docs/FullScatterplotFlu.png",
        "docs/FullLinkedChartClickableFlu.pdf",
        "docs/FullLinkedChartClickableFlu.html"
    conda: "cartography.yml"
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
    conda: "cartography.yml"
    shell:
        """
        pandoc --filter pandoc-citeproc --bibliography=docs/cartography.bib -s docs/index.md -o docs/cartography.{wildcards.suffix}
        """