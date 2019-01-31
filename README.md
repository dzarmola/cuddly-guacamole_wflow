# cuddly-guacamole_wflow
Workflow for aligning protein family profiles.

Uses hhsuite for profile creation and pairwise alignments, cd-hit for representative selection,
mcl for clustering [optional], and clustalo for adding single sequences to the family alignments
[subject to change].

```
usage: full_workflow.py [-h] [--run_name RUN_NAME] [--evalue EVALUE]
                        [--num_representatives NUM_REPRESENTATIVES]
                        [--mparam MPARAM] [--inflation INFLATION]
                        directories [directories ...]

Workflow for alignment of family profiles [with extracting just part ofthe
profile based on structure]. Needs hhmake, hhsearch, clustalo, cd-hit and mcl
(if clustering) installed and available on system path.

positional arguments:
  directories           Directories with groups of families. Each directory
                        must contain an 'org_fastas' directory with fasta
                        formatted alignments (.fa or .fasta) of desired
                        families. Additionally, if a 'ranges.json' files and
                        'struct_hhms' directory are presented all family
                        profiles in this group will be aligned to profiles
                        from struct_hhms, and just columns correspnding to
                        ranges specified in the ranges.json will be used in
                        further analysis.

optional arguments:
  -h, --help            show this help message and exit
  --run_name RUN_NAME, -n RUN_NAME
                        Name to be given to this run. Directory will be
                        created in the current working dir.
  --evalue EVALUE, -e EVALUE
                        E-value cutoff for significant hhsearch hits
  --num_representatives NUM_REPRESENTATIVES, -r NUM_REPRESENTATIVES
                        Number of representative sequences for each family
  --mparam MPARAM, -m MPARAM
                        Cutoff parameter for columns in hhmake
  --inflation INFLATION, -i INFLATION
                        Inflation value for mcl clustering
  ```
