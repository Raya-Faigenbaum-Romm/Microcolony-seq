import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd

### Create_bins returns an equal-width (distance) partitioning. It returns an ascending list of tuples, representing the intervals.
def create_bins(lower_bound, width, quantity):
    bins = []
    for low in range(lower_bound,
                     lower_bound + quantity * width + 1, width):
        bins.append((low, low + width))
    return bins

### find_bin() gets a number in value and bins. bins is a list of tuples, like [(0,20), (20, 40), (40, 60)], binning returns the smallest index i of bins so that bin[i][0] <= value < bin[i][1]
def find_bin(value, bins):
    for i in range(0, len(bins)):
        if bins[i][0] <= value < bins[i][1]:
            return i
    return -1

### fromDataToBins gets a df and bins and assigns a bin for each gene based on its mean expression (log10(basemean))
def fromDataToBins(df, bins):
    for i, row in df.iterrows():  # Initialize for loop
        bin_index = find_bin(float(row['log10_basemean']) * 100, bins)
        df.loc[i,'Bin'] = bin_index

### Main function that reads from a csv file normalized reads of counts, each biological sample in a separate columns and calculates for each gene if it is a variable gene or not based on the dispersion and mean expression of each gene
### The normalized reads are calculated by DESeq2 before running defineOutliersFromNormCounts() function
### The functions saves a plot of dispersion vs. mean expression and an excel table with the infromation for each gene if it is a dispersion outliers, i.e., a variable gene
def defineOutliersFromNormCounts():
    df = pd.read_csv(r"I:\My_Microcolony_seq_project\My_normalized_counts.txt", sep="\t")

    #Part 1: calculate basemean, variance, dispersion.
    #Speficy the names of the columns containing the normalized reads counts for each sample in the df
    #Calcualte mean expression for each gene and save it in the basemean column
    df['basemean'] = df[['UTI_1', 'UTI_2', 'UTI_3', 'UTI_4', 'UTI_5', 'UTI_7', 'UTI_9', 'UTI_11', 'UTI_13',  'UTI_14', 'UTI_15', 'UTI_16', 'UTI_17', 'UTI_20', 'UTI_21', 'UTI_24', 'UTI_25', 'UTI_27', 'UTI_28', 'UTI_29']].mean(axis=1)

    #Take out genes with basemean <10 which probably include only noise
    df = df[df['basemean'] >= 10]

    # Calcualte variance for each gene and save it in the var column
    df['var'] = df[['UTI_1', 'UTI_2', 'UTI_3', 'UTI_4', 'UTI_5', 'UTI_7', 'UTI_9', 'UTI_11', 'UTI_13',  'UTI_14', 'UTI_15', 'UTI_16', 'UTI_17', 'UTI_20', 'UTI_21', 'UTI_24', 'UTI_25', 'UTI_27', 'UTI_28', 'UTI_29']].var(axis=1)

    #Calcualte dispersion for each gene
    df['dispersion'] = (df['var'] - df['basemean'])/np.square(df['basemean'])
    df['log10_basemean'] = np.log10(df['basemean'])

    #Take out genes with negative dipersion
    df = df[df['dispersion'] > 0]
    df['log10_dispersion'] = np.log10(df['dispersion'])

    # Part 2: Separation to binds according to mean expression
    bins_1 = create_bins(lower_bound=100, width=25,  quantity=6)
    bins_2 = create_bins(lower_bound=275, width=800, quantity=1) # A single bin for genes with high expression
    bins = bins_1 + bins_2
    print(bins)
    #Get the bin for each gene
    fromDataToBins(df, bins)
    numOfGenesInBins = df.groupby('Bin')['log10_dispersion'].count()
    print (numOfGenesInBins)

    #Part 3: Find dispersion outliers, i.e., variable genes
    df['log10_dispersion'] = df['log10_dispersion'].astype(float)

    stdOfBins = df.groupby('Bin')['log10_dispersion'].std()
    meanOfBins = df.groupby('Bin')['log10_dispersion'].mean()

    stdOfBins = {i: v for i, v in enumerate(stdOfBins)}
    meanOfBins = {i: v for i, v in enumerate(meanOfBins)}

    df['meanOfBin'] = df['Bin'].map(meanOfBins)
    df['stdOfBin'] = df['Bin'].map(stdOfBins)
    df['disp_minus_mean'] = df['log10_dispersion'] - df['meanOfBin']

    df['outlier_by_bins'] = 'False'
    df.loc[df['disp_minus_mean'] >= df['stdOfBin'] * 2  , 'outlier_by_bins'] = 'True'

    #Part 4: Scatter plot
    fig = px.scatter(df,x='basemean',y='dispersion', color='outlier_by_bins' ,color_discrete_sequence=['teal','red'] ,
                     log_x=True, log_y=True,
                     labels={
                     "basemean": "Mean expression",
                     "dispersion": "Dispersion",
                     "outlier_by_bins": "Dispersion outlier"
                 },)
    fig.update_xaxes(dtick=1)
    fig.update_yaxes(dtick=1)
    fig.update_traces(marker_size=14)
    fig.update_layout(yaxis=dict(tickfont=dict(size=22)), xaxis=dict(tickfont=dict(size=22)), font=dict(size=22,color="black")) #family="Ariel"
    fig.update_layout(template="plotly_white")
    fig.update_layout(legend=dict(yanchor="top", y=1.5, xanchor="center", x=0.9))
    #Save the plot under the desired destination
    fig.write_image(r"I:\My_Microcolony_seq_project\Dispersion_vs_mean_expression_scatter_plot.jpeg", width=6, height=3, scale=2)

    #Save table with variable genes under the desired destination
    df.to_excel(r"I:\My_Microcolony_seq_project\UTIs_unicycler_contigs_tech_std_2_210224.xlsx")

### Main
defineOutliersFromNormCounts()
