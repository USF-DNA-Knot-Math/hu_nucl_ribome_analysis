# Human Nuclear Ribome analysis

This folder contains the code used to analyze genomic DNA samples treated with <i>Escherichia coli</i> RNase HII, an enzyme that introduces nicks at sites containing embedded ribonucleoside monophosphates (rNMPs). The treated DNA samples were subsequently denatured and separated by agarose gel electrophoresis.

Gel electrophoresis images were processed to quantify DNA fragment distributions for each treatment condition. These fragment-length distributions were then analyzed using mathematical simulations to estimate the number of embedded rNMPs present in each experiment.

The simulation code used for fragment-size estimation is available in the 
<a href="./fragment size analysis/">fragment size analysis</a> folder.

## Input Data
The gel electrophoresis pictures are included in the folder <a href="./gel_regular_figures/">gel_regular_figures</a>. These were inverted the colors using regular image editing software. The results of this transformation are included in the folder <a href="./gel_inverted_images/">gel_inverted_images</a>.

## Main Analysis
The gel images were analysed with the Python Jupyter Notebook <a href="main_figure_analysis.ipynb".main_figure_analysis.ipynb</a>. The uploaded notebook has the results for the gel figure "sample cropped gel.jpg", which were repeated for all the gels in <a href="./gel_regular_figures/">gel_regular_figures</a>. We summarize this analysis here:

1. We estimated the total number of DNA nucleotides loaded per lane, using the average molecular weight of a nucleotide (330 g/mol) and Avogadro's constant:

$$ N = \frac{W}{330 \times 10^9} \times 6.023 \times 10^{23},$$

where \(W\) is the DNA mass in nanograms (\(W=250\) for the example).

2. We inverted the colors in the gel image to produce the file "gel_grey.png".
3. We used Gelpy to set up gel images and lane pixel ranges. These were adjusted manually to coincide with the image lanes.
4. 



