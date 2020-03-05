"""Functions for manipulating and plotting pairwise distances and reduced
dimensionality embeddings of distance matrices.
"""
import numpy as np


def get_hamming_distances(genomes):
    """Calculate pairwise Hamming distances between the given list of genomes
    and return the nonredundant array of values for use with scipy's squareform function.
    Bases other than standard nucleotides (A, T, C, G) are ignored.

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
    [0, 1, 1]
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
        for j in range(i + 1, len(genomes)):
            # Find all mismatches between these two genomes.
            mismatches = genome_arrays[i] != genome_arrays[j]

            # Count the number of mismatches where both genomes have valid bases.
            hamming_distances.append((mismatches & valid_bases[i] & valid_bases[j]).sum())

    return hamming_distances


"""
principal_Df -- Data from data reduction (-T-SNE, MDS, etc) (pandas DataFrame)
result_metadata -- the metadata that is being read in (Pandas DataFrame)
fields, the parts of metadata that should be concatenated with princiapl_Df (list)
"""


def concatenate_results_with_strain_data(principal_Df, result_metadata, fields):
    import pandas as pd
    finalDf = pd.concat([principal_Df, result_metadata[fields]], axis=1)
    return finalDf


"""
Defining Fields:
finalDf: The data that is used to generate the scatter plot (pandas DataFrame)
x, the data you want on the x axis (string)
y, the data you want on the y axis (string)
Titlex,the name you want on the x axis (string)
Titley, the name you want on the y axis (string)
Tooltip, when scanning over a point, the data you want avaiable (list)
Color, what the scatterplot is colored by (String)
"""


def scatterplot_with_tooltip_interactive(finalDf, x, y, Titlex, Titley, ToolTip, color):
    import altair as alt
    brush = alt.selection(type='interval', resolve='global')
    chart = alt.Chart(finalDf).mark_circle(size=60).encode(
        x=alt.X(x, title=Titlex),
        y=alt.X(y, title=Titley),
        color=color,
        tooltip=ToolTip
    ).interactive()
    # chart.display()
    return chart


"""
dataframe: dataframe including node data and dimensionality reduction data (Pandas Dataframe)
list_of_data: list of all the names of the columns in the dataframe for which you want graphs: goes in the order of [x1,y1,x2,y2,x3,y3] etc.(list)
list_of_titles: list of all the TITLES you want for each axis: goes in order of[x1,y1,x2,y2,x3,y3] etc.(list)
color: what the data should be colored by (string)
ToolTip: when hovering over the data, what data should be shown (list)
"""


def linking_tree_with_plots_brush(dataFrame, list_of_data, list_of_titles, color, ToolTip):
    import altair as alt
    import pandas as pd
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
                title="Date"
            ),
            y=alt.Y(
                "y:Q",
                title=""
            ),
            color=alt.condition(brush, color, alt.ColorValue('gray')),
            tooltip=ToolTip
        ).add_selection(brush).properties(width=400, height=250)
        list_of_chart.append(tree_name)

        for i in range(0, len(list_of_data) - 1, 2):
            if(i == len(list_of_data)):
                break
            chart = base.mark_circle(size=60).encode(
                x=alt.X(list_of_data[i], title=list_of_titles[i]),
                y=alt.X(list_of_data[i + 1], title=list_of_titles[i + 1]),
                color=alt.condition(brush, color, alt.ColorValue('gray')),
                tooltip=ToolTip
            ).add_selection(
                brush
            ).properties(
                width=250,
                height=250
            )
            list_of_chart.append(chart)
        return list_of_chart


def linking_tree_with_plots_clickable(dataFrame, list_of_data, list_of_titles, colors, fields, ToolTip):
    import altair as alt
    import pandas as pd
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


def scatterplot_xyvalues(strains, similarity_matrix, df_merged, column1, column2, type_of_embedding):
    import pandas as pd
    import altair as alt
    import numpy as np
    from scipy.spatial.distance import squareform, pdist
    from scipy.stats import linregress
    import pandas as pd
    import numpy as np
    embedding_df = df_merged[[column1, column2, 'strain']]
    pairwise_distance_array = np.array(similarity_matrix)[
        np.triu_indices(len(df_merged), k=0)]

    indexes_tuple = np.triu_indices(len(df_merged), k=0)
    row_number = indexes_tuple[0]
    column_number = indexes_tuple[1]
    row_strain = pd.DataFrame([strains[x] for x in row_number])
    column_strain = pd.DataFrame([strains[x] for x in column_number])

    pairwise_df = pd.DataFrame(pairwise_distance_array)
    row_column = row_strain.merge(
        column_strain, how='outer', left_index=True, right_index=True)
    row_column.columns = ["row", "column"]

    euclidean_df = get_euclidean_data_frame(
        df_merged, column1, column2, "MDS").dropna()
    euclidean_df.columns = ["euclidean", "clade_status", "embedding"]

    row_column_pairwise = row_column.merge(
        pairwise_df, how='outer', left_index=True, right_index=True)
    row_column_pairwise = row_column_pairwise.where(
        row_column_pairwise["row"] != row_column_pairwise["column"]).dropna().set_index(euclidean_df.index)

    total_df = row_column_pairwise.merge(
        euclidean_df, how='inner', left_index=True, right_index=True).dropna()
    total_df.columns = ["row", "column", "genetic",
                        "euclidean", "clade_status", "embedding"]

    return total_df


def scatterplot_tooltips(strains, similarity_matrix, df_merged, column1, column2, type_of_embedding, n_sample):
    import pandas as pd
    import altair as alt
    import numpy as np
    from scipy.spatial.distance import squareform, pdist
    from scipy.stats import linregress
    import pandas as pd
    import numpy as np
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


def LOESS_scatterplot(strains, similarity_matrix, df_merged, column1, column2, type_of_embedding, n_sample):
    import pandas as pd
    import altair as alt
    import numpy as np
    from scipy.spatial.distance import squareform, pdist
    from scipy.stats import linregress
    import pandas as pd
    import numpy as np
    import statsmodels
    total_df = scatterplot_xyvalues(strains, similarity_matrix,
        df_merged, column1, column2, type_of_embedding)
    y_values = statsmodels.nonparametric.smoothers_lowess.lowess(
        total_df["euclidean"], total_df["genetic"], frac=0.6666666666666666, it=3, delta=0.0, is_sorted=False, missing='drop', return_sorted=True)

    PD_Y_values = pd.DataFrame(y_values)
    PD_Y_values.columns = ["LOWESS_x", "LOWESS_y"]

    line = alt.Chart(PD_Y_values).mark_line(color='red').encode(
        alt.X('LOWESS_x'),
        alt.Y('LOWESS_y')
    )

    chart = scatterplot_tooltips(strains, similarity_matrix,
        df_merged, column1, column2, type_of_embedding, n_sample)
    return chart + line


def get_euclidean_data_frame(sampled_df, column1, column2, embedding):
    import pandas as pd
    import altair as alt
    import numpy as np
    from scipy.spatial.distance import squareform, pdist
    from scipy.stats import linregress
    import pandas as pd
    import numpy as np
    import statsmodels
    """
    Returns a data frame of Euclidean distances for the requested embedding columns.
    
    The given `sampled_df` MUST include a "clade_membership" column.
    """
    # Traverse pairs of samples from left-to-right, top-to-bottom
    # along the upper triangle of the pairwise matrix and collect
    # the clade status of each pair as either within- or between-clades.
    # This traversal excludes self-self comparisons along the diagonal.
    clade_status = []
    clade_memberships = sampled_df["clade_membership"].values
    for i in range(sampled_df.shape[0] - 1):
        for j in range(i + 1, sampled_df.shape[0]):
            if clade_memberships[i] != clade_memberships[j]:
                clade_status.append("between")
            else:
                clade_status.append("within")

    # Calculate pairwise distances between samples for the requested columns.
    # The resulting array is in the same left-to-right, top-to-bottom order
    # as the clade statuses above.
    sampled_distances = pdist(sampled_df[[column1, column2]])

    # Align clade status with pairwise distance for each pairwise comparison.
    sampled_distances_df = pd.DataFrame(
        {"distance": sampled_distances, "clade_status": clade_status})

    # Annotate the requested embedding.
    sampled_distances_df["embedding"] = embedding

    return sampled_distances_df
