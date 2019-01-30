# Biological-Interactions-Analysis

### Introductions 
------
Biological interactions network analysis which related to a list of specific (interested) human genes.<br/>

### Data and Code
------
In this project Biological Interactions Network analysis related to specific human [seed genes](https://github.com/AAbasinejad/Biological-Interactions-Analysis/blob/master/seed_genes.txt) has been carried out by using Human [Integrated Interactions Database](http://iid.ophid.utoronto.ca/static/download/human_annotated_PPIs.txt.gz) (a.k.a IID) which is an on-line database of detected and predicted protein-protein interactions (PPIs) and [BioGRID](https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.168/BIOGRID-ORGANISM-3.5.168.tab2.zip) which is an interaction repository with data compiled through comprehensive curation efforts.<br />
The BioGRID and IID datasets mentioned above is exactly the ones that has been used in this project but you can find other formats of same biogrid dataset [here](https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.5.168/) or other PPI networks databases provided by IID [here](http://iid.ophid.utoronto.ca/search_by_proteins/).<br />
Furthermore, all basic and essential needed data has been fetched (`Basic_info.py`) from National Center for Biotechnology Information [NCBI](https://www.ncbi.nlm.nih.gov/) website which was approved also by HUGO Gene Nomenclature Committee (a.k.a [HGNC](https://www.genenames.org/)) website.<br/>

In order to run this code you have to put all *needed files* in a same directory and run this command in terminal:<br />


`python main.py <seed_genes.txt> <BIOGRID-ORGANISM-Homo_sapiens dataset> <IID dataset>`

**Note**: Needed files are `main.py`, `Basic_info.py`, `Interactions.py`, `Network_Analysis.py` plus the two mentioned datasets. <br/>
**Note**: This project has been done in `Python 3.x`, so ... .<br/>
**P.S.**: running of this code in **lunch break** is strongly recommended. :D <br/>

**Specific Libraries**

```python
import requests
import pandas as pd
from bs4 import BeautifulSoup
import networkx as nx
import markov_clustering
import community
from scipy.stats import hypergeom
```


### Modules
------
TBD soon...<br/>
A brief explaination of each module(both technically and conceptually) will be here ASAP. ;) <br/>

### Terms
------
For better understanding of this project, some terms and abbreviations will be defined in follow:<br/>

**PPI:** Protein-Protein Interaction
**Uniprot AC:** Uniprot AC ‘Accession Number’ of each gene (a.k.a Uniprot entry, e.g. P01344, P15502, etc.).<br/> 
**GeneSymbol:** Scienctific symbol of each gene (e.g. IGF2, ELN, PTPRC, etc.).<br/>
*Note:* In general GeneSymbol is more important for both practical and scientific purposes since it's more understandable.<br/>
**SG:** in code it stands for Seed Gene (in fact it refers to seed_genes list).<br/>
**SGI:** Seed Gene Interactions which refers to interactions that involves seed genes only. (from both DBs)<br/>
**Union_Interactions:** It represents all interactions that involves at least one seed gene. (from both DBs)<br/>
**Intersection_Interactions:** all interactions that involves at least one seed gene, confirmed by both DBs.<br/>
*Note:* In general in both Union and Intersection interactions, the interactions between interactomes which has direct interaction with a [seed_gene](https://github.com/AAbasinejad/Biological-Interactions-Analysis/blob/master/seed_genes.txt) has been considered.<br/>
*Note:* In code, variables with sgi, u an I signs refers to SGI, Union_Interactions and Intersection_Interactions respectively.<br/>
**p-value:** p-value is to measure under- or over-enrichment based on the cumulative distribution function (CDF) of the hypergeometric distribution.<br/> 
*Note:* for computing p-value you can use [Hypergeometric p-value calculator](http://systems.crump.ucla.edu/hypergeometric/index.php) or `hypergeom` library from `scipy.stats`.
**Note:** in this project we will often use the terms ‘gene’ and ‘protein’ as synonyms, even if they are not, from the purely biological point of view.<br/>


### putative disease proteins using the DIAMOnD tool
------
Using the tool DIAMOnD, compute the putative disease protein list using as reference interactome (“network_file”) the latest BioGrid interactome already used to collect [PPIs](https://github.com/AAbasinejad/Biological-Interactions-Analysis/blob/master/PPI.txt).
The [`DIAMOnD.py`](https://github.com/AAbasinejad/Biological-Interactions-Analysis/blob/master/DIAMOnD.py) file up here is compatible with `Python 3.x` while you can find comatible version with `Python 2.x` on [DIAMOnD source page](https://github.com/barabasilab/DIAMOnD.git).<br/)
In order to find how to use it, check the source page as well. ;)
