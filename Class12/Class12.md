Class 12: Structural Bioinformatics and Drug Discvoery
================

Clean up our protein target structure
-------------------------------------

First we download a target (i.e. protein receptor) structure from the main PDB database. We will pick PDB ID "1hsg"

``` r
library(bio3d)

pdb.code <- "1hsg"
file.name <- get.pdb(pdb.code)
```

    ## Warning in get.pdb(pdb.code): ./1hsg.pdb exists. Skipping download

Extract the protein only segments of this PDB entry and write out a new PDB format file. We will also do the same for the bound ligand.

``` r
hiv <- read.pdb(file.name)
hiv
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

Protein extraction first

``` r
prot <- trim.pdb(hiv, "protein")
prot
```

    ## 
    ##  Call:  trim.pdb(pdb = hiv, "protein")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
prot.filename <- paste(pdb.code, "_protein.pdb", sep="")
#write.pdb(prot, file = "1hsg_protein.pdb") will hardcode the protein into the script. Using the paste function allows us to have a more friendly code to use in the future. 
write.pdb(prot, file =prot.filename)
```

Do the same thing for the ligand...

Ligand extraction

``` r
lig <- trim.pdb(hiv, "ligand")
lig
```

    ## 
    ##  Call:  trim.pdb(pdb = hiv, "ligand")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
lig.filename <- paste(pdb.code, "_ligand.pdb", sep="")
write.pdb(lig, file =lig.filename)
```

Convert our docking results for viewing in VMD
----------------------------------------------

``` r
res <- read.pdb("all.pdbqt", multi = TRUE)
res
```

    ## 
    ##  Call:  read.pdb(file = "all.pdbqt", multi = TRUE)
    ## 
    ##    Total Models#: 13
    ##      Total Atoms#: 50,  XYZs#: 1950  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 50  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, xyz, calpha, call

``` r
write.pdb(res, file = "results.pdb")
```

``` r
# res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1] 10.446 11.590 10.773 11.857 11.228 13.357 11.716  9.146 10.033 11.926
    ## [11]  9.103  9.361 13.066

``` r
rmsd(res)
```

    ## Warning in rmsd(res): No indices provided, using the 50 non NA positions

    ##         [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]
    ##  [1,]  0.000  4.593  2.459  8.482  4.214  9.814  8.493 12.286 12.697
    ##  [2,]  4.593  0.000  5.495  9.180  2.384 10.665  9.775 11.409 14.497
    ##  [3,]  2.459  5.495  0.000  7.422  4.804  8.652  7.598 12.767 11.930
    ##  [4,]  8.482  9.180  7.422  0.000  8.317  3.955  5.685 14.975  9.841
    ##  [5,]  4.214  2.384  4.804  8.317  0.000  9.619  8.766 11.722 13.681
    ##  [6,]  9.814 10.665  8.652  3.955  9.619  0.000  5.699 16.849  9.151
    ##  [7,]  8.493  9.775  7.598  5.685  8.766  5.699  0.000 16.151  8.409
    ##  [8,] 12.286 11.409 12.767 14.975 11.722 16.849 16.151  0.000 16.552
    ##  [9,] 12.697 14.497 11.930  9.841 13.681  9.151  8.409 16.552  0.000
    ## [10,]  8.418  9.743  7.477  5.439  8.745  5.028  2.471 15.899  8.557
    ## [11,] 12.570 11.545 13.129 15.542 11.901 17.388 16.336  2.507 16.643
    ## [12,] 12.730 11.840 13.399 15.405 12.143 17.267 16.466  2.911 16.817
    ## [13,]  9.523 11.407  8.320  6.201 10.542  4.907  6.735 17.705  8.384
    ##        [,10]  [,11]  [,12]  [,13]
    ##  [1,]  8.418 12.570 12.730  9.523
    ##  [2,]  9.743 11.545 11.840 11.407
    ##  [3,]  7.477 13.129 13.399  8.320
    ##  [4,]  5.439 15.542 15.405  6.201
    ##  [5,]  8.745 11.901 12.143 10.542
    ##  [6,]  5.028 17.388 17.267  4.907
    ##  [7,]  2.471 16.336 16.466  6.735
    ##  [8,] 15.899  2.507  2.911 17.705
    ##  [9,]  8.557 16.643 16.817  8.384
    ## [10,]  0.000 16.211 16.227  6.551
    ## [11,] 16.211  0.000  3.360 18.088
    ## [12,] 16.227  3.360  0.000 18.131
    ## [13,]  6.551 18.088 18.131  0.000
