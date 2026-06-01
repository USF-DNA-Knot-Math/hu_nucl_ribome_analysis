# Human Nuclear Ribome analysis

This folder contains the code used to analyze genomic DNA samples treated with <i>Escherichia coli</i> RNase HII, an enzyme that introduces nicks at sites containing embedded ribonucleoside monophosphates (rNMPs). The treated DNA samples were subsequently denatured and separated by agarose gel electrophoresis.

Gel electrophoresis images were processed to quantify DNA fragment distributions for each condition, both before and after treatment with RNase HII. These fragment-length distributions were then analyzed using mathematical simulations to estimate the number of embedded rNMPs present in each experiment. The simulation code used for fragment-size estimation is available in the 
<a href="./fragment size analysis/">fragment size analysis</a> folder.

## Input Data

The original gel electrophoresis images are provided in the folder <a href="./gel_regular_figures/">gel_regular_figures</a>. For analysis with the Python library GelPy, the images were converted to inverted grayscale versions, which are available in the folder <a href="./gel_inverted_images/">gel_inverted_images</a>.

## Main Analysis
The gel images were analyzed using the Python Jupyter notebook <a href="./main_figure_analysis.ipynb">main_figure_analysis.ipynb</a>. The notebook included in this repository contains the analysis of the gel image [sample cropped gel.jpg](./sample%20cropped%20gel.jpg), and the same workflow was applied to each gel image in the folder <a href="./gel_regular_figures/">gel_regular_figures</a>. We summarize the analysis pipeline below.

1. We estimated the total number of DNA nucleotides loaded per lane, using the average molecular weight of a nucleotide (330 g/mol) and Avogadro's constant:

$$ N = \frac{W}{330 \times 10^9} \times 6.023 \times 10^{23},$$

where \(W\) is the DNA mass in nanograms (\W=250\ for the example).

2. We converted the gel image to grayscale and inverted its intensity values to produce [`gel_grey.png`](./gel_grey.png), which was used as the input image for subsequent analysis in GelPy.

3. We imported [`gel_grey.png`](./gel_grey.png) into GelPy and created separate gel objects for the untreated and treated samples. We assigned and manually adjusted lane labels and pixel ranges to align with the corresponding gel lanes, using a shared ladder lane as a reference in both analyses. We then used GelPy to extract intensity profiles for the selected lanes. The resulting lane-profile plots and gel overview figures are displayed directly in the notebook output.

4. We extracted the normalized lane intensity profiles from GelPy for both samples and normalized each profile to represent the percentage of the total signal intensity. Using the ladder lane, we identified the pixel positions corresponding to DNA markers of known molecular weight and fit an exponential decay model relating pixel position to fragment size. We then used the fitted model to estimate the molecular weight of every pixel row in the gel image. The file [`ladder_values.csv`](./ladder_values.csv) contains a smoothed intensity profile of the ladder lane. The fitted curve and ladder marker positions are shown in [`alkaline_pixels_MW.png`](./alkaline_pixels_MW.png).
  
5. We estimated the number of DNA fragments represented at each pixel position in every lane by assuming that SYBR Gold fluorescence intensity is proportional to the amount of DNA present. After removing low-intensity noise using a 5% quantile threshold, we converted normalized intensity values into fragment counts by scaling with the estimated number of nucleotides loaded per lane and dividing by the molecular weight assigned to each pixel position. To focus on the signal-containing portion of each lane, we truncated the profiles using a common cutoff derived from the intensity distributions. We then compared untreated and RNase HII-treated samples using Mood’s median test and the Mann–Whitney U test on both the intensity and fragment-count distributions.

This analysis produced [`ladder_fragments.csv`](./ladder_fragments.csv), containing the estimated fragment distribution for the ladder lane; [`gel24_main_fragments_nontreated.csv`](./gel24_main_fragments_nontreated.csv) and [`gel24_main_fragments_treated.csv`](./gel24_main_fragments_treated.csv), containing the estimated fragment counts as a function of molecular weight for each sample, and [`stats.csv`](./stats.csv), containing the statistical comparisons between untreated and treated samples.

6. We visualized the normalized lane intensity profiles for both untreated and treated samples by plotting smoothed intensity distributions across each lane. These plots provided a qualitative comparison of DNA fragment distributions between samples and experimental conditions. We also computed the total estimated number of fragments for each lane by summing the fragment-count distributions obtained in the previous step and displayed these values in the notebook output.

This analysis produced the figures [`intensity_no_treatment_full.png`](./intensity_no_treatment_full.png) and [`intensity_treatment_full.png`](./intensity_treatment_full.png), showing the smoothed intensity profiles. The total estimated fragment counts for each sample lane were printed directly on the notebook.

7. Finally, we computed summary statistics for the estimated fragment-count distributions in each lane. 

This analysis produced the figure [`fragmentsmedians.png`](./fragmentsmedians.png), which displays the median fragment positions for each sample before and after RNase HII treatment. In addition, in the notebook we printed the estimated median fragment molecular weights for each lane.


