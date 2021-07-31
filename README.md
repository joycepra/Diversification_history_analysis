Species Delimitation and Biogeography 

Scripts and input files used in the paper published in Molecular Phylogenetics and Evolution:

Prado JR, Knowles LL, Percequillo AR. 2021. New species boundaries and the diversification history of marsh rat taxa clarify historical connections among ecologically and geographically distinct wetlands of South America. Molecular Phylogenetics and Evolution 155: 106992.

Input files of all analyzes presented in the manuscript

(1)Ancestral Range Estimation
- Holochilus_tree.nwk: Phylogentic tree in newick format
- geo_data: Geographic data from each species

Biogeobears scripts can be found on the software's website

(2)Divergence Time 
Input Files from both approaches.
- Concatenated.xml
- snapp.xml

(3) Morphological _Variation
- morphometric.csv: file containing all the morphological data used in the manuscript

(4) Processing_Sequence_Data
- Processing_Seq_Data.R: 
Visualize and filter loci to trim sites with high variation suggestive of sequencing 
and/or alignment errors and remove suspicious clusters of possible paralogs based
on a maximum pairwise sequence divergence

(5) Species_Delimitation 
- BPP > RAxML_InputTree > 200L
Input files 
200_3-2-3-2.ctl
200_3-2-3-2.ctl
200_3-002-3-2.ctl
200_3-002-3-002.ctl
L200.txx

- BPP > RAxML_InputTree > 500L
Input files 
500_3-2-3-2.ctl
500_3-2-3-2.ctl
500_3-002-3-2.ctl
500_3-002-3-002.ctl
L500.txx


- BPP > SVDQuartets_InputTree > 200L
Input files 
200_3-2-3-2.ctl
200_3-2-3-2.ctl
200_3-002-3-2.ctl
200_3-002-3-002.ctl
L200.txx

- BPP > SVDQuartets_InputTree > 500L
Input files 
500_3-2-3-2.ctl
500_3-2-3-2.ctl
500_3-002-3-2.ctl
500_3-002-3-002.ctl
L500.txx

- iBPP input files
morpho_ibpp : input file with the morphological data

-iBPP > RAxML_InputTree 
Input files 
200_1-10-1-10.ctl
200_1-10-2-2000.ctl
200_2-2000-1-10.ctl
200_2-2000-2-2000.ctl

-iBPP > SVDQuartets_InputTree 
Input files 
200_1-10-1-10.ctl
200_1-10-2-2000.ctl
200_2-2000-1-10.ctl
200_2-2000-2-2000.ctl

(6) SVDQuartets 
Nexus input file used in each step
- SVD_ind_tree.nex 
- SVD_species_tree
