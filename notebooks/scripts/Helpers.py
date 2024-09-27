"""Functions for manipulating and plotting pairwise distances and reduced
dimensionality embeddings of distance matrices.
"""
from augur.io import read_sequences
from augur.utils import read_tree
import Bio.SeqIO
from collections import OrderedDict
import numpy as np
import altair as alt
import pandas as pd
import re
from scipy.spatial.distance import squareform, pdist
from scipy.stats import linregress, entropy
import statsmodels

def joint_entropy(X, Y):
    n = sum(len(x) for x in X)
    H_X_Y = 0.0
    for x in X:
        for y in Y:
            p_x_y = len(set(x) & set(y)) / n
            if p_x_y > 0:
                H_X_Y -= (p_x_y * np.log2(p_x_y))

    return H_X_Y

def variation_of_information(X, Y, normalized=False):
    """Calculate variation of information (VI) score between the ground truth
    clustering and the proposed clustering. The score is 0 for an exact match,
    and the log of the total samples in the case of a complete difference.

    Parameters
    ----------
    X : list of list
        a list of total elements n (where n is number of total samples), partitioned into lists defining
        separated clusters
    Y : list of list
        a list of total elements n (where n is number of total samples), partitioned into lists defining
        separated clusters
    normalized : boolean, default = False
        determines if the VI score is normalized or not

    Returns
    -------
    float :
        the variation of information score between two separate clusterings.


    For maximally separated clusters, VI should be log_2(n) = log_2(4) = 2.

    >>> X = [[1], [2], [3], [4]]
    >>> Y = [[1, 2, 3, 4]]
    >>> variation_of_information(X, Y) == np.log2(len(X))
    True

    For the same maximally separated clusters, the normalized VI should be 1.0.

    >>> variation_of_information(X, Y, normalized=True) == 1.0
    True

    For identical clusters, VI should be 0.

    >>> X = [[1], [2], [3], [4]]
    >>> Y = [[1], [2], [3], [4]]
    >>> variation_of_information(X, Y) == 0.0
    True

    """
    H_X = entropy([len(k) for k in X], base=2)
    H_Y = entropy([len(k) for k in Y], base=2)
    H_X_Y = joint_entropy(X, Y)

    # From Equation 22 of Meila 2007
    VI = (2 * H_X_Y) - H_X - H_Y

    if normalized:
        n = sum(len(x) for x in X)
        return VI / np.log2(n)
    else:
        return VI

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


def get_clade_label_chart(
    df, x_column, y_column, text_column, drop_labels=None,
    xoffset_by_label=None,
    yoffset_by_label=None,
):
    """Build an Altair chart of text labels from the given data frame, placing
    labels at the mean x and y position of the given x and y columns and using
    text from the given text column.

    """
    if drop_labels is None:
        drop_labels = {"other"}
    else:
        drop_labels = set(drop_labels) | {"other"}

    clade_label_positions = df.loc[
        (~df["is_internal_node"]) & (~df[text_column].isin(drop_labels)),
        [text_column, x_column, y_column]
    ].groupby(
        text_column
    ).aggregate({
        x_column: "mean",
        y_column: "mean",
    }).reset_index()

    if xoffset_by_label is not None:
        clade_label_positions[x_column] = clade_label_positions.apply(
            lambda record: record[x_column] + xoffset_by_label.get(record[text_column], 0),
            axis=1
        ).values

    if yoffset_by_label is not None:
        clade_label_positions[y_column] = clade_label_positions.apply(
            lambda record: record[y_column] + yoffset_by_label.get(record[text_column], 0),
            axis=1
        ).values

    clade_labels_chart = alt.Chart(clade_label_positions).mark_text().encode(
        x=f"{x_column}:Q",
        y=f"{y_column}:Q",
        text=f"{text_column}:N",
        size=alt.SizeValue(10),
    )

    return clade_labels_chart


def make_node_branch_widths(tree_path, min_width=1.0, max_width=4.0):
    """Build a data frame of branch widths for a given Newick tree path bounded
    by the given min and max width. Branch widths reflect the log of the number
    of leaves descended from each node normalized to range from 0 to 1. Nodes
    receive a branch width of the minimum width plus the remaining range between
    the max and min values scaled by the normalized log leaves.
    """
    tree = read_tree(tree_path)

    node_branch_widths = []
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            node.branch_width = 1
        else:
            node.branch_width = sum((child.branch_width for child in node.clades))

        node_branch_widths.append({
            "node": node.name,
            "leaves": node.branch_width
        })

    node_branch_widths = pd.DataFrame.from_records(node_branch_widths)
    node_branch_widths["branch_weight"] = np.sqrt(node_branch_widths["leaves"] / node_branch_widths["leaves"].max())

    width_range = max_width - min_width
    node_branch_widths["branch_width"] = min_width + width_range * node_branch_widths["branch_weight"]

    return node_branch_widths


def make_branch_lines_for_columns(embedding_segments, column1, column2, color_domain=None, color_range=None, color_spec="clade_membership_color:N"):
    """Return an Altair chart representing branch line segments defined in the
    given data frame. Each segment spans from `column1` to `{column1}_parent`
    and `column2` to `{column2}_parent`. Color branches by clade membership with
    an optional domain and range. Otherwise, color branches gray.

    """
    if color_domain and color_range:
        color = alt.Color(color_spec).scale(domain=color_domain, range=color_range)
    else:
        color = alt.ColorValue("#cccccc")

    branch_lines = alt.Chart(embedding_segments).mark_rule().encode(
        x=f"{column1}_parent:Q",
        x2=f"{column1}:Q",
        y=f"{column2}_parent:Q",
        y2=f"{column2}:Q",
        color=color,
        strokeWidth=alt.StrokeWidth("branch_width", legend=None),
    )

    return branch_lines


def linking_tree_with_plots_brush(dataFrame, list_of_data, list_of_titles, color, legend_title, ToolTip, domain=None, range_=None, legend_columns=1, plot_legend=True, color_branches=True):
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
        what the data should be colored by (ex. by clade, by region) including the Vega type specification (e.g., "clade:N")
    legend_title: string
        title to use for the color legend
    ToolTip: list
        when hovering over the data, what data should be shown

    Returns
    ---------
    A brushable altair plot combining the tree with the plots of columns passed in
    """
    x_nice = 15.0
    y_nice = True

    list_of_chart = []
    if(len(list_of_data) % 2 != 0 or len(list_of_titles) % 2 != 0):
        raise Exception(
            'The length of list_of_data and the length of list_of_titles should not be odd.')
    else:
        # Split data frame into internal nodes and tips.
        internal_nodes_df = dataFrame[dataFrame["is_internal_node"]].copy()
        tips_df = dataFrame[~dataFrame["is_internal_node"]].copy()

        # Optionally, split tips data frame into gray tips and color tips. Gray
        # tips represent "unassigned" values and get plotted below the color
        # tips. We have to check for the gray color in the color scale to find
        # the value that corresponds to the unassigned state and then check the
        # tips data frame for those values to make sure we don't have an empty
        # gray tips data frame.
        color_map = dict(zip(range_, domain))
        unassigned_color = "#999999"
        unassigned_value = color_map.get(unassigned_color)

        # Extract the color field name from the Vega field specification which
        # includes type (e.g., "clade:N").
        color_field = color.split(":")[0]
        color_values = tips_df[color_field].drop_duplicates().values

        if unassigned_color in color_map and unassigned_value in color_values:
            tips_gray_df = tips_df[tips_df[color_field] == unassigned_value].copy()
            tips_color_df = tips_df[tips_df[color_field] != unassigned_value].copy()
        else:
            tips_gray_df = None
            tips_color_df = tips_df

        # Prepare tree plot which is circles for tips and lines for branches
        # (internal nodes).
        brush = alt.selection(type='interval', resolve='global')
        tips_color = alt.Chart(tips_color_df).mark_circle().encode(
            x=alt.X(
                "divergence:Q",
                scale=alt.Scale(
                    domain=(dataFrame["parent_mutation"].min(), dataFrame["divergence"].max()), nice=x_nice),
                title="Divergence",
                axis=alt.Axis(labels=True, ticks=True)
            ),
            y=alt.Y(
                "y_value:Q",
                title="",
                axis=alt.Axis(labels=False, ticks=False)
            ),
            color=alt.condition(brush, if_true=alt.Color(color, scale=alt.Scale(domain=domain, range=range_)), if_false=alt.ColorValue('gray')),
            tooltip=ToolTip
        ).add_selection(brush)

        if tips_gray_df is not None:
            tips_gray = alt.Chart(tips_gray_df).mark_circle().encode(
                x=alt.X(
                    "divergence:Q",
                    scale=alt.Scale(
                        domain=(dataFrame["parent_mutation"].min(), dataFrame["divergence"].max()), nice=x_nice),
                    title="Divergence",
                    axis=alt.Axis(labels=True, ticks=True)
                ),
                y=alt.Y(
                    "y_value:Q",
                    title="",
                    axis=alt.Axis(labels=False, ticks=False)
                ),
                color=alt.condition(brush, if_true=alt.Color(color, scale=alt.Scale(domain=domain, range=range_)), if_false=alt.ColorValue('gray')),
                tooltip=ToolTip
            )
            tips = tips_gray + tips_color
        else:
            tips = tips_color

        # Look for parent color field (e.g., "parent_clade_membership_color") in
        # data frame to determine whether to color tree branches by that field.
        # This field allows us to color branches by the value of the node from
        # which the branch descends.
        if color_branches and color_field in dataFrame.columns:
            tree_branch_color_spec = f"{color_field}:N"
            tree_branch_color = alt.Color(tree_branch_color_spec).scale(domain=domain, range=range_)
        else:
            tree_branch_color = alt.ColorValue("#cccccc")

        if "branch_width" in dataFrame.columns:
            tree_branch_width = alt.StrokeWidth("branch_width", legend=None)
        else:
            tree_branch_width = alt.StrokeValue(1.25)

        horizontal_lines = alt.Chart(dataFrame).mark_rule().encode(
            x=alt.X("parent_mutation:Q", scale=alt.Scale(domain=(dataFrame["parent_mutation"].min(), dataFrame["divergence"].max()), nice=x_nice)),
            x2="divergence:Q",
            y=alt.Y("y_value:Q", scale=alt.Scale(domain=(dataFrame["y_value"].min(), dataFrame["y_value"].max()), nice=y_nice)),
            color=tree_branch_color,
            strokeWidth=tree_branch_width,
        )

        vertical_line_aggregates = {
            "y_min": ("y_value", "min"),
            "y_max": ("y_value", "max"),
            "x": ("parent_mutation", "min"),
        }

        if color_branches and f"parent_{color_field}" in dataFrame.columns:
            parent_tree_branch_color_spec = f"parent_{color_field}:N"
            parent_tree_branch_color = alt.Color(parent_tree_branch_color_spec).scale(domain=domain, range=range_)
            vertical_line_aggregates[f"parent_{color_field}"] = (f"parent_{color_field}", "first")
        else:
            parent_tree_branch_color = alt.ColorValue("#cccccc")

        if "branch_width" in dataFrame.columns:
            vertical_line_aggregates["branch_width"] = ("branch_width", "max")

        vertical_line_positions = dataFrame.groupby("parent_name").aggregate(**vertical_line_aggregates)
        vertical_lines = alt.Chart(vertical_line_positions).mark_rule().encode(
            x="x:Q",
            y="y_min:Q",
            y2="y_max:Q",
            color=parent_tree_branch_color,
            strokeWidth=tree_branch_width,
        )

        tree = (horizontal_lines + vertical_lines + tips).properties(width=560, height=250)
        list_of_chart.append(tree)

        # Prepare embedding plots which are tips only.
        for i in range(0, len(list_of_data) - 1, 2):
            if(i == len(list_of_data)):
                break

            if plot_legend:
                legend = alt.Legend(
                        symbolLimit=len(domain),
                        columns=legend_columns,
                        title=legend_title,
                    )
            else:
                legend = None

            chart_color = alt.Chart(tips_color_df).mark_circle(size=60).encode(
                x=alt.X(list_of_data[i], title=list_of_titles[i], axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y(list_of_data[i + 1], title=list_of_titles[i + 1], axis=alt.Axis(labels=False, ticks=False)),
                color=alt.condition(
                    brush,
                    if_false=alt.ColorValue('gray'),
                    if_true=alt.Color(color, scale=alt.Scale(domain=domain, range=range_), legend=legend)
                ),
                tooltip=ToolTip
            ).add_selection(
                brush
            )

            # As with the tree plot, plot gray tips first, if they exist, and
            # then plot color tips second to place them on top.
            if tips_gray_df is not None:
                if plot_legend:
                    legend = alt.Legend(
                            symbolLimit=len(domain),
                            columns=legend_columns,
                            title=legend_title,
                        )
                else:
                    legend = None

                chart_gray = alt.Chart(tips_gray_df).mark_circle(size=60).encode(
                    x=alt.X(list_of_data[i], title=list_of_titles[i], axis=alt.Axis(labels=False, ticks=False)),
                    y=alt.Y(list_of_data[i + 1], title=list_of_titles[i + 1], axis=alt.Axis(labels=False, ticks=False)),
                    color=alt.condition(
                        brush,
                        if_false=alt.ColorValue('gray'),
                        if_true=alt.Color(color, scale=alt.Scale(domain=domain, range=range_), legend=legend)
                    ),
                    tooltip=ToolTip
                )
                chart = chart_gray + chart_color
            else:
                chart = chart_color

            chart = chart.properties(
                width=250,
                height=250,
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
