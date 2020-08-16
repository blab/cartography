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
		
rule run_zika:
    input:
        tree = zikaWorkFlow("auspice/zika-cartography_tree.json"),
        disease_name = "notebooks/Data/name_of_disease_zika.txt",
        pca = zikaWorkFlow("results/embed_pca.csv"),
        mds = zikaWorkFlow("results/embed_mds.csv"),
        tsne = zikaWorkFlow("results/embed_tsne.csv"),
        umap = zikaWorkFlow("results/embed_umap.csv"),
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
		tree = fluWorkFlow("auspice/flu_seasonal_h3n2_ha_2y_tree.json"),
        clade_names = "notebooks/Data/clade_names.txt",
        disease_name = "notebooks/Data/name_of_disease_flu.txt"
    output:
        "docs/FullViolinPlotFlu.png",
        "docs/FullScatterplotFlu.png",
        "docs/FullLinkedChartClickableFlu.pdf",
        "docs/FullLinkedChartClickableFlu.html"
    conda: "cartography.yml"
    notebook:											   
        "notebooks/2019-08-08FinalNotebookFlu.ipynb"
	shell:
	"""
		snp-sites -o flu-nextstrain/results-20182020/variable_sites.fasta flu-nextstrain/results-20182020/aligned.fasta \
		snp-sites -o flu-nextstrain/results/variable_sites.fasta flu-nextstrain/results/aligned.fasta \
		
        {python:q} scripts/bases_missing.py \
            --strains {input.strains} \
            --dropped_strains {input.dropped_strains} \
            --disease_name {input.disease_name} \
			
		{python:q} scripts/hamming_distance.py \
	"""
		
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