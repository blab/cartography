rule all:
    input:
        "docs/cartography.pdf",
        "docs/cartography.html"

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