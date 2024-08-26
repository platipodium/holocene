<!--
# SPDX-FileContributor: Kai Wirtz <kai.wirtz@hereon.de>
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
-->

# MATLAB and R-scripts for the analysis of archaeological C14 and paleoclimate proxy data

Code and data availability on this repository:  https://github.com/kaiwirtz3/holocene

The MATLAB and R code enables standalone or parallel processing.

> On High Performance Computers (HPC) with the SLURM scheduler, parallel processing of batch jobs can be invoked on HPC by `sbatch --array=1-64 slurm.sh`.

## Requirenments

### MATLAB

 * needs packages: 'M_MAP' and 'pca_toolbox'
 * insert/correct installation path in `addpath('~/tools/m_map')` found in the header of some MATLAB-scripts (foremost `load_par.m`), or `addpath('~/tools/pca_toolbox')` in `proxy_dtw.m`

###  R

 * needs libraries: rworldmap, rcarbon, R.matlab, RColorBrewer, cluster, sp/sf
   for rcarbon version <1.5 see comments in grid_growth.r ("spweights")

## Input data
  (in directories c14mat/ paleoclim/ and data/)

  Note that the input data are already available in pre-processed binary (MATLAB) format.

 The  C14 dates from the p3k14c dataset were downloaded in Mar 2023 from
  https://www.nature.com/articles/s41597-022-01118-7/tables/4.
  Updates can be found at https://www.p3k14c.org
  The p3k14c csv/xls file can be converted by "read_dASIS.m" (see outcommented settings) to
  "c14mat/p3k14c.mat" (needed by "merge_C14data.m")
  The new compilation refers to the database dASIS which in turn collects previous databases
  and adds own dates compiled by Gronenborn/Antunes at LEIZA Mainz
  The xls file "c14mat/C14_dASIS.xlsx" was converted by runing "read_dASIS.m"
  to "c14mat/C14_dASIS.mat"
  Using the script "merge_C14data.m", all binary files for pooled and gridded were produced

   c14mat/C14_europe_[1-64].mat:  radiocarbon dates for Europe on a coarse grid
   c14mat/C14_europe0.mat:        radiocarbon dates for Europe pooled
   c14mat/C14_[cont].mat:         radiocarbon dates for continents
   c14mat/C14_NIreland.mat:       radiocarbon dates for North Ireland
   c14mat/C14_dASIS.xlsx:         compilation of (dASIS) C14v dates in XLS format

  collected paleoclimate proxy time-series and meta-info in
    paleoclim/InEur_27_920.mat:   index ('InEur') to proxies in EU domain
    paleoclim/proxydescription_488_0_11.00_1.70.mat: MATLAB structure ("evinfo") with paleoclimate info
      for 209 global records. the 98 used indices (for Europe) are given by 'InEur'
      meta-info:
       e.g., "Proxy" for proxy type (e.g. d18O), "Plotname" for site name, "Latitude", "Longitude"
        or "Source" for bibtex-key of original reference
      evinfo.time:  cells with vectors containing the time values of a record
      evinfo.value: cells with vectors containing the normalized values
    Note that the MATLAB structure 'evinfo' and its fields can be converted into an R dataframe
     using the R.matlab package

  misc time-series or fields
    SA_spd_rgr.mat: RGR South America from digitized and calculated from Riris et al 2019
    Steinhilber2012_Solar.dat: Total Solar Irradiance (TSI) from Steinhilber et al 2012
    bog_std.mat: Northern Irish bog data by Rowan McLaughlin
    seamask_norm_0.05.mat: landmask Europe at 0.05 degree resolution

## Execution:

all scripts can be invoked by a single master shell command: `./master.sh`
 this serial processing may take some time depending on computational resources;
time consuming script processing such as of 'spd_growth.r' or 'cluster.r' can be run at a HPC in parallel mode (see also slurm.sh)

for running R-scripts (e.g."r-example.r"):  start 'R' terminal and type 'source("r-example.r")'
                                or from (bash) terminal, type  'Rscript r-example.r'
for running MATLAB-scripts (e.g."m-example.m"): start MATLAB terminal and type 'm-example'
                                or from (bash) terminal, type  'matlab -nodesktop -r "try; m-example; catch; end; quit" '
--------------------------------------------------------------------
sequence of script execution and dependencies (see also master.sh)

1. grid_growth.r   # calculates spatial statistics of SPDs on a 4°x4° grid

2. collect.r       # collect/merge grid cells and prepares clustering

3. cluster.r       # create region clusters of C14 sites

4. make_grid.m     # spatial kriging on a grid based on cluster points
  (uses cl_distance.m, make_grid_regions.m)

5. spd_growth.r    # calculates Summed Probability Density (SPD)
                #   and related growth (RGR) for each region and time slice
  (uses movavg.r) # reads PrePop_ clusti writes AllPop for time slices

6. spd_pooled.r    # calculates SPD and RGR for pooled method
                #reads C14_europe0 (or C14_EA, C14_NIreland), writes AllPop_all

7. plot_RGR.m      # process and plot calculated growth rates (RGR)
  (uses movavg.m, calc_aravg_rgr.m)
                # reads AllPop_ writes AllPop_tag_all avg_rgr_  RGR_Comp.png

### Process climate proxy data
proxy_dtw.m     # applies DTW writes dtw_proxydata
  (uses dtw.m by T. Felty, movweighavg.m)

collect_ts.m    # integrate and smooth time-series, merges DTW time segments data
                # and includes other times series (e.g. RGR, solar forcing)
  (uses movavg.m, movweighavg.m) # reads AllPop_EA_all avg_rgr_ data/* writes target_ts_0

### collects RGR from time slices
### plot_varmap_slice %reads AllPop_i

glmloop.r       # run GLM model for different sets of input variables
   (uses do_lgm.r)  # reads target_ts_0, writes glmres*
   
[ add_logitres.m  # add logit model results to matrix ]
overlap_ts.m    # calc and plot overlap between time-series
   (uses add_logitres.m and calc_overlap.m)

## License

If not stated otherwise, the entire analysis software is licensed under
the GNU Public License version 3 or later. See <http://www.gnu.org/licenses/gpl-3.0.txt> for the complete terms.

## Documentation

see text and equations in Material and Methods of Wirtz et al, Multicentennial cycles in continental demography synchronous with solar activity and climate stability, subm.
--------------


## Date preprocessing

### Archaeological data

The original dASIS dataset was filtered end enriched with the P3K14C dataset.  The combined dataset was spatially fuzzified and overlaps were removed.  If you use this dataset for Europe, please cite our dataset, and - either/or the dASIS and P314K dataset.

#### Radon
Open Access 

Rinne et al. 2024: Christoph Rinne/Jutta Kneisel/Martin Hinz/Martin Furholt/Nina Krischke/Johannes Müller/Dirk Raetzel-Fabian/Marcel Rodens/Karl-Göran Sjögren/Helle Vandkilde and Hans-Peter Wotzka, Rado.NB. In: https://radonb.ufg.uni-kiel.de. 


Martin Hinz, Martin Furholt, Johannes Müller, Dirk Raetzel-Fabian, Christoph Rinne, Karl-Göran Sjögren, Hans-Peter Wotzka, RADON - Radiocarbon dates online 2012. Central European database of 14C dates for the Neolithic and Early Bronze Age. www.jungsteinsite.de, 2012, 1-4.

#### P3K14C

CC0-1.0

Bird, D., Miranda, L., Vander Linden, M. et al. p3k14c, a synthetic global database of archaeological radiocarbon dates. Sci Data 9, 27 (2022). https://doi.org/10.1038/s41597-022-01118-7

#### AIDA

CC-BY-4.0

Alessio Palmisano, Andrew Bevan, Alexander Kabelindde, Neil Roberts, & Stephen Shennan. (2022). AIDA (Archive of Italian radiocarbon DAtes) (v3.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5846835

#### MedAfriCarbon

CC-BY-4.0

Lucarini, G., Wilkinson, T., Crema, E. R., Palombini, A., Bevan, A., & Broodbank, C. (2020). The MedAfriCarbon radiocarbon database and web application. Archaeological dynamics in Mediterranean Africa, ca. 9600-700 BC [Data set]. In Journal of Open Archaeology Data (1.0.3). Zenodo. https://doi.org/10.5281/zenodo.3689716

#### Aegean History

CC-BY-4.0

Katsianis, Markos; Bevan, Andrew; Styliaras, Giorgos; Maniatis, Yannis (2020). Dataset for: An Aegean history and archaeology written through radiocarbon dates. University College London. Dataset. https://doi.org/10.5522/04/12489137.v1

#### Chapple

CC-BY-4.0

Chapple, R. M. (2018) ‘Catalogue of radiocarbon determinations & dendrochronology dates’. Zenodo. doi: 10.5281/zenodo.3367518.

#### EUBAR

CAPUZZO G. 2014. Space-temporal analysis of radiocarbon evidence and associated archaeological record: from Danube to Ebro rivers and from Bronze to Iron ages, Ph.D. Thesis discussed at the Autonomous University of Barcelona.

#### NAKALA 

CC-BY-4.0

Perrin, Thomas (2021) «BDA - Base de Données Archéologiques» [Dataset] NAKALA. https://doi.org/10.34847/nkl.dde9fnm8


### Paleoclimate data

The paleoclimate data was obtained from the PANGAEA and NOAA archives.  Many of the data not available in these archives was digitized from the graphics in the original publication.  All data sources are cited in the supplementary material. 

If you use the combined dataset provided here, we ask you to cite this publication, and -- if highlighting a particular source -- also the original dataset.