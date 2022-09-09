"""Functions for manipulating and plotting pairwise distances and reduced
dimensionality embeddings of distance matrices.
"""
from augur.io import read_sequences
import Bio.SeqIO
from collections import OrderedDict
import numpy as np
import altair as alt
import pandas as pd
import re
from scipy.spatial.distance import squareform, pdist
from scipy.stats import linregress
import statsmodels

def get_embedding_columns_by_method(method):
    if method in ("pca"):
        return list(f"{method}1 {method}2 {method}3 {method}4 {method}5 {method}6 {method}7 {method}8 {method}9 {method}10".split())
    if method in ("mds"):
        return list(f"{method}1 {method}2 {method}3 {method}4".split())
    if method in ("t-sne"):
        return list("tsne_x tsne_y".split())
    else:
        return list(f"{method}_x {method}_y".split())

def get_PCA_feature_matrix(alignment, sequence_names):
    """Calculate PCA feature matrix from the alignment and the strains that should be kept
    within the analysis
    Parameters
    ----------
    alignment : string
        string corresponding to the address of the alignment file (local or global)
    Returns
    -------
    numpy array feature map of PCA
    """
    sequences_by_name = OrderedDict()

    for sequence in read_sequences(alignment):
        if sequence.id in sequence_names:
            sequences_by_name[sequence.id] = str(sequence.seq)

    sequence_names_val = list(sequences_by_name.keys())
    assert(len(sequence_names_val) == len(sequence_names))

    # Create a matrix representation of sequences with one row per sequence and
    # one integer value per character state in the sequence. Order rows to match
    # the order of strain names in the `sequence_names` input.
    numbers = [sequences_by_name[strain] for strain in sequence_names]
    for i, strain in enumerate(sequence_names):
        numbers[i] = re.sub(r'[^AGCT]', '5', numbers[i])
        numbers[i] = list(numbers[i].replace('A','1').replace('G','2').replace('C', '3').replace('T','4'))
        numbers[i] = [int(j) for j in numbers[i]]

    numbers = np.array(numbers)
    return numbers


def get_hamming_distances(genomes):
    """Calculate pairwise Hamming distances between the given list of genomes
    and return the nonredundant array of values for use with scipy's squareform function.
    Bases other than standard nucleotides (A, T, C, G) are ignored. Treat indels as a single event.
    Parameters
    ----------
    genomes : list
        a list of strings corresponding to genomes that should be compared
    Returns
    -------
    list
        a list of distinct Hamming distances as a vector-form distance vector
    >>> genomes = ["ATGCT", "ATGCT", "ACGCT"]
    >>> get_hamming_distances(genomes)
    [0, 1, 1]
    >>> genomes = ["AT-GCT", "AT--CT", "AC--CT"]
    >>> get_hamming_distances(genomes)
    [1, 2, 2]
    >>> genomes = ["ACTGG", "A--GN", "A-NGG"]
    >>> get_hamming_distances(genomes)
    [1, 1, 1]
    >>> genomes = ["ACTGTA", "A--CCA", "A--GT-"]
    >>> get_hamming_distances(genomes)
    [3, 2, 4]
    """

    # Define an array of valid nucleotides to use in pairwise distance calculations.
    # Using a numpy array of byte strings allows us to apply numpy.isin later.
    nucleotides = np.array([b'A', b'T', b'C', b'G'])

    # Convert genome strings into numpy arrays to enable vectorized comparisons.
    genome_arrays = [
        np.frombuffer(genome.encode(), dtype="S1")
        for genome in genomes
    ]

    # Precalculate positions of valid bases (A, T, C, and G) in each genome to speed up later comparisons.
    valid_bases = [
        np.isin(genome_array, nucleotides)
        for genome_array in genome_arrays
    ]

    # Calculate Hamming distance between all distinct pairs of genomes at valid bases.
    # The resulting list is a reduced representation of a symmetric matrix that can be
    # converted to a square matrix with scipy's squareform function:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html
    hamming_distances = []
    for i in range(len(genomes)):
        # Only compare the current genome, i, with all later genomes.
        # This avoids repeating comparisons or comparing each genome to itself.
        np_genomes = np.array(list(re.sub("-", "0", genomes[i])))

        for j in range(i + 1, len(genomes)):
            np_genomes_b = np.array(list(re.sub("-", "0", genomes[j])))
            result = np.where(np_genomes_b == "0")
            np_genomes[result] = 0
            a = (np_genomes=="0")
            num_indel = (a&~np.r_[[False],a[:-1]]).sum()
            # Find all mismatches between these two genomes.
            mismatches = genome_arrays[i] != genome_arrays[j]

            # Count the number of mismatches where both genomes have valid bases.
            hamming_distances.append(((mismatches & valid_bases[i] & valid_bases[j]).sum()) + num_indel)

    return hamming_distances


def get_y_positions(tree):
    """Create a mapping of each clade to its vertical position. Dict of {clade:
    y-coord}. Coordinates are negative, and integers for tips.
    We use the y position layout function from BioPython [1]. This function is
    hidden inside the top-level draw function, so we cannot reuse it.
    [1] https://github.com/biopython/biopython/blob/d1d3c0d6ab33de12057201e09eb48bdb1964521a/Bio/Phylo/_utils.py#L471-L495
    Parameters
    ----------
    tree : Bio.Phylo.BaseTree
        a tree from BioPython
    Returns
    -------
    dict
        mapping of BioPython Clade instances to y-axis coordinates
    """
    maxheight = tree.count_terminals()
    # Rows are defined by the tips
    heights = {
        tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
    }

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (
            heights[clade.clades[0]] + heights[clade.clades[-1]]
        ) / 2.0

    if tree.root.clades:
        calc_row(tree.root)
    return heights


def concatenate_results_with_strain_data(principal_Df, result_metadata, fields):
    """Takes the data from data reductions (T-SNE, MDS, etc) and pairs up each strain's euclidean plotpoints with its metadata

    Parameters
    ----------
    principal_Df: Pandas Dataframe
        Data from data reduction
    result_metadata: Pandas Dataframe
        the metadata that is being read in (Pandas DataFrame)
    fields: list
        the parts of metadata that should be concatenated with princiapl_Df (eg. "strain", "region", "country")

    Returns
    --------
    finalDf: Pandas Dataframe
        the concatenated Pandas Dataframe
    """
    finalDf = pd.concat([principal_Df, result_metadata[fields]], axis=1)
    return finalDf

def scatterplot_with_tooltip(finalDf, x, y, Titlex, Titley, ToolTip, color):
    """Creates an interactive scatterplot in altair

    Parameters
    -----------
    finalDf: Pandas Dataframe
        the data that is used to generate the scatter plot
    x: string
        the data for the x axis
    y: string
        the data for the y axis
    Titlex: string
        the name for the x axis
    Titley: string
        the name for the y axis
    Tooltip: list
        the data available when scanning over a plot
    Color: string
        what the scatterplot is colored by

    Returns
    --------
    an Altair chart
    """
    brush = alt.selection(type='interval', resolve='global')
    chart = alt.Chart(finalDf).mark_circle(size=60).encode(
        x=alt.X(x, title=Titlex),
        y=alt.X(y, title=Titley),
        color=color,
        tooltip=ToolTip
    ).interactive()
    # chart.display()
    return chart

def scatterplot_with_tooltip_interactive(finalDf, x, y, Titlex, Titley, ToolTip, color, domain=None, range_=None):
    """Creates an interactive scatterplot in altair

    Parameters
    -----------
    finalDf: Pandas Dataframe
        the data that is used to generate the scatter plot
    x: string
        the data for the x axis
    y: string
        the data for the y axis
    Titlex: string
        the name for the x axis
    Titley: string
        the name for the y axis
    Tooltip: list
        the data available when scanning over a plot
    Color: string
        what the scatterplot is colored by

    Returns
    --------
    an Altair chart
    """
    brush = alt.selection(type='interval', resolve='global')
    chart = alt.Chart(finalDf).mark_circle(size=60).encode(
        x=alt.X(x, title=Titlex),
        y=alt.X(y, title=Titley),
        color=alt.Color(color, scale=alt.Scale(domain=domain, range=range_)),
        tooltip=ToolTip
    ).interactive()
    # chart.display()
    return chart


def linking_tree_with_plots_brush(dataFrame, list_of_data, list_of_titles, color, ToolTip, domain=None, range_=None):
    """Creates a linked brushable altair plot with the tree and the charts appended
    Parameters
    -----------
    dataframe: Pandas Dataframe
        dataframe including node data and dimensionality reduction data
    list_of_data: list
        list of all the names of the columns in the dataframe for which you want graphs: goes in the order of [x1,y1,x2,y2,x3,y3] etc.
    list_of_titles: list
        list of all the TITLES you want for each axis: goes in order of[x1,y1,x2,y2,x3,y3] etc.
    color: string
        what the data should be colored by (ex. by clade, by region)
    ToolTip: list
        when hovering over the data, what data should be shown

    Returns
    ---------
    A brushable altair plot combining the tree with the plots of columns passed in
    """

    list_of_chart = []
    if(len(list_of_data) % 2 != 0 or len(list_of_titles) % 2 != 0):
        raise Exception(
            'The length of list_of_data and the length of list_of_titles should not be odd.')
    else:
        base = alt.Chart(dataFrame)
        brush = alt.selection(type='interval', resolve='global')
        tree_name = base.mark_circle().encode(
            x=alt.X(
                "date:Q",
                scale=alt.Scale(
                    domain=(dataFrame["date"].min() - 0.2, dataFrame["date"].max() + 0.2)),
                title="Date",
                axis=alt.Axis(labels=False, ticks=False)
            ),
            y=alt.Y(
                "y:Q",
                title="",
                axis=alt.Axis(labels=False, ticks=False)
            ),
            color=alt.condition(brush, if_false=alt.ColorValue('gray'), if_true=alt.Color(color, scale=alt.Scale(domain=domain, range=range_))),
            tooltip=ToolTip
        ).add_selection(brush).properties(width=560, height=250)
        list_of_chart.append(tree_name)

        for i in range(0, len(list_of_data) - 1, 2):
            if(i == len(list_of_data)):
                break
            chart = base.mark_circle(size=60).encode(
                x=alt.X(list_of_data[i], title=list_of_titles[i], axis=alt.Axis(labels=False, ticks=False)),
                y=alt.X(list_of_data[i + 1], title=list_of_titles[i + 1], axis=alt.Axis(labels=False, ticks=False)),
                color=alt.condition(
                    brush,
                    if_false=alt.ColorValue('gray'),
                    if_true=alt.Color(color, scale=alt.Scale(domain=domain, range=range_), legend=alt.Legend(symbolLimit=len(domain)))
                ),
                tooltip=ToolTip
            ).add_selection(
                brush
            ).properties(
                width=250,
                height=250
            )
            list_of_chart.append(chart)
        return list_of_chart


def linking_tree_with_plots_clickable(dataFrame, list_of_data, list_of_titles, colors, fields, ToolTip, domain=None, range_=None):
    """
    Parameters
    -----------
    dataframe: Pandas Dataframe
        dataframe including node data and dimensionality reduction data
    list_of_data: list
        list of all the names of the columns in the dataframe for which you want graphs: goes in the order of [x1,y1,x2,y2,x3,y3] etc.
    list_of_titles: list
        list of all the TITLES you want for each axis: goes in order of[x1,y1,x2,y2,x3,y3] etc.
    color: string
        what the data should be colored by (ex. by clade, by region)
    ToolTip: list
        when hovering over the data, what data should be shown

    Returns
    ---------
    A clickable altair plot combining the tree with the plots of columns passed in
    """
    list_of_chart = []
    if(len(list_of_data) % 2 != 0 or len(list_of_titles) % 2 != 0):
        raise Exception(
            'The length of list_of_data and the length of list_of_titles should not be odd.')
    else:
        base = alt.Chart(dataFrame)
        selection = alt.selection_multi(fields=fields)

        color = alt.condition(selection,
                              alt.Color(colors, legend=None),
                              alt.value('lightgray'))
        tree_name = base.mark_circle().encode(
            x=alt.X(
                "date:Q",
                scale=alt.Scale(
                    domain=(dataFrame["date"].min() - 0.2, dataFrame["date"].max() + 0.2)),
                title="Date"
            ),
            y=alt.Y(
                "y:Q",
                title=""
            ),
            color=color,
            tooltip=ToolTip
        ).add_selection(selection).properties(width=400, height=250)

        list_of_chart.append(tree_name)
        for i in range(0, len(list_of_data) - 1, 2):
            if(i == len(list_of_data)):
                break
            chart = base.mark_circle(size=60).encode(
                x=alt.X(list_of_data[i], title=list_of_titles[i]),
                y=alt.X(list_of_data[i + 1], title=list_of_titles[i + 1]),
                color=color,
                tooltip=ToolTip
            ).add_selection(
                selection
            ).properties(
                width=250,
                height=250
            )
            list_of_chart.append(chart)
        legend = base.mark_point().encode(
            y=alt.Y(colors, axis=alt.Axis(orient='right')),
            color=colors
        ).add_selection(
            selection
        )
        list_of_chart.append(legend)

        return list_of_chart


def scatterplot_xyvalues(strains, similarity_matrix, embedding_df, column_list, type_of_embedding):
    """Returns a unraveled similarity matrix Pandas dataframe of pairwise and euclidean distances for each strain pair
     Parameters
    -----------
    strains: list
        list of strains for the build (ex. A/Oman/5263/2017)
    similarity_matrix: Pandas Dataframe
        a similarity matrix using hamming distance to compare pairwise distance between each strain
    df_merged: Pandas Dataframe
        dataframe containing the euclidean coordinates of each strain in an embedding (PCA, MDS, t-SNE, UMAP)
    column_list: list
        string list contaning the names of the columns to find the distance between
    type_of_embedding: string
        "MDS", "PCA", "TSNE", or "UMAP"
    Returns
    ---------
    A Pandas Dataframe of pairwise and euclidean distances for every strain pair
    """
    pairwise_distance_array = np.array(similarity_matrix)[
        np.triu_indices(len(embedding_df), k=0)]

    indexes_tuple = np.triu_indices(len(embedding_df), k=0)
    row_number = indexes_tuple[0]
    column_number = indexes_tuple[1]
    row_strain = pd.DataFrame([strains[x] for x in row_number])
    column_strain = pd.DataFrame([strains[x] for x in column_number])

    pairwise_df = pd.DataFrame(pairwise_distance_array)
    row_column = row_strain.merge(
        column_strain, how='outer', left_index=True, right_index=True)
    row_column.columns = ["row", "column"]

    distances = pdist(embedding_df[column_list])
    euclidean_df = pd.DataFrame({"distance": distances})
    euclidean_df["embedding"] = type_of_embedding
    euclidean_df.columns = ["euclidean", "embedding"]

    row_column_pairwise = row_column.merge(
        pairwise_df, how='outer', left_index=True, right_index=True)
    row_column_pairwise = row_column_pairwise.where(
        row_column_pairwise["row"] != row_column_pairwise["column"]).dropna().set_index(euclidean_df.index)

    total_df = row_column_pairwise.merge(
        euclidean_df, how='inner', left_index=True, right_index=True).dropna()
    total_df.columns = ["row", "column", "genetic",
                        "euclidean", "embedding"]

    return total_df


def scatterplot_tooltips(strains, similarity_matrix, df_merged, column1, column2, type_of_embedding, n_sample):
    """uses scatterplot_xyvalues, returns a pairwise vs euclidean scatterplot

    Parameters
    -----------
    strains: list
        list of strains for the build (ex. A/Oman/5263/2017)
    similarity_matrix: Pandas Dataframe
        a similarity matrix using hamming distance to compare pairwise distance between each strain
    df_merged: Pandas Dataframe
        dataframe containing the clade information and euclidean coordinates of each strain in an embedding (PCA, MDS, t-SNE, UMAP)
    column1: string
        the name of the first column in df_merged to compare distance between
    column2: string
        the name of the second column in df_merged to compare distance between
    type_of_embedding: string
        "MDS", "PCA", "TSNE", or "UMAP"
    n_sample:
        how many strains to sample (altair cannot currently run with that many strains, sample around 1000 - 2000)

    Returns
    ----------
    An altair pairwise vs Euclidean scatterplot with tooltips
    """
    alt.data_transformers.disable_max_rows()

    total_df = scatterplot_xyvalues(strains, similarity_matrix,
                                    df_merged, column1, column2, type_of_embedding)
    regression = linregress(total_df["genetic"], total_df["euclidean"])
    slope, intercept, r_value, p_value, std_err = regression

    chart = alt.Chart(total_df.sample(n=n_sample)).mark_circle(size=60).encode(
        x=alt.X('genetic', title="genetic distance"),
        y=alt.X('euclidean', title="Euclidean distance"),
        tooltip=['row', 'column', 'genetic', 'euclidean']
    ).properties(title="Genetic vs. Euclidean scatterplot: " + type_of_embedding + "  (R^2 = " + str((r_value ** 2).round(3)) + ")", height=200, width=300)
    return chart

def get_euclidean_data_frame(sampled_df, column_for_analysis, embedding, column_list=None):
    """Gives a dataframe of euclidean distances for embedding columns to use in plotting and analysis
    Parameters
    -----------
    sampled_df: pandas DataFrame
        a dataframe of euclidean coordinate points containing the two columns passed in
    column1: string
        the name of the first column in sampled_df
    column2: string
        the name of the second column in sampled_df
    column_for_analysis: string
        the name of the column which the dataframe will construct "between" and "within" from (ex. Host, Clade_membership, etc)
    Returns
    ----------
    A data frame of Euclidean distances for the requested embedding columns.
    """
    # Traverse pairs of samples from left-to-right, top-to-bottom
    # along the upper triangle of the pairwise matrix and collect
    # the clade status of each pair as either within- or between-clades.
    # This traversal excludes self-self comparisons along the diagonal.
    clade_status = []
    clade_memberships = sampled_df[column_for_analysis].values
    for i in range(sampled_df.shape[0] - 1):
        for j in range(i + 1, sampled_df.shape[0]):
            if clade_memberships[i] == clade_memberships[j]:
                clade_status.append(1)
            else:
                clade_status.append(0)

    # Calculate pairwise distances between samples for the requested columns.
    # The resulting array is in the same left-to-right, top-to-bottom order
    # as the clade statuses above.
    if (column_list is not None):
        sampled_distances = pdist(sampled_df[column_list])

    else:
        sampled_distances = squareform(sampled_df.drop([column_for_analysis, "strain"], axis=1))

    # Align clade status with pairwise distance for each pairwise comparison.
    sampled_distances_df = pd.DataFrame(
        {"distance": sampled_distances, "clade_status": clade_status})

    # Annotate the requested embedding.
    sampled_distances_df["embedding"] = embedding

    return sampled_distances_df
