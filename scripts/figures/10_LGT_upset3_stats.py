import pandas as pd

og_lgt_genes_file = "data/orthogroups/og_lgt_genes.csv"

og_lgt_genes = pd.read_csv(og_lgt_genes_file, header="infer", sep="\t").drop_duplicates("id")


"Stats"

og_lgt_genes = og_lgt_genes.rename(columns={"carpe": "C. membranifera",
                            "kbiala": "K. bialata",
                            "HIN": "H. inflata",
                            "trepo": "Trepomonas pc1-LGT",
                            "spiro": "S. salmonicida",
                            "wb": "G. intestinalis",
                            "muris": "G. muris"})


# Get the subset of genes for blue category"
blue_subset = og_lgt_genes[
    (og_lgt_genes['K. bialata'] == False) &
    (og_lgt_genes['C. membranifera'] == False) &
    (og_lgt_genes['H. inflata'] == True) &
    (og_lgt_genes['Trepomonas pc1-LGT'] == True) &
    (og_lgt_genes['S. salmonicida'] == True) &
    (og_lgt_genes['G. intestinalis'] == True) &
    (og_lgt_genes['G. muris'] == True)]

# Get the count of genes in this filtered subset
number_of_blue_genes = len(blue_subset)
print("")
print(f"Number of genes highlighted in blue: {number_of_blue_genes}")

# Get the subset of genes for red categroy"
red_subset = og_lgt_genes[
    (og_lgt_genes['K. bialata'] == False) &
    (og_lgt_genes['C. membranifera'] == False) &
    (og_lgt_genes['H. inflata'] == True) &
    (og_lgt_genes['Trepomonas pc1-LGT'] == True) &
    (og_lgt_genes['S. salmonicida'] == False) &
    (og_lgt_genes['G. intestinalis'] == False) &
    (og_lgt_genes['G. muris'] == False)]

# Get the count of genes in this filtered subset
number_of_red_genes = len(red_subset)
print(f"Number of genes highlighted in red: {number_of_red_genes}")

# Get the subset of genes for green categroy"
green_subset = og_lgt_genes[
    (og_lgt_genes['K. bialata'] == True) &
    (og_lgt_genes['C. membranifera'] == True) &
    (og_lgt_genes['H. inflata'] == True) &
    (og_lgt_genes['Trepomonas pc1-LGT'] == True) &
    (og_lgt_genes['S. salmonicida'] == False) &
    (og_lgt_genes['G. intestinalis'] == False) &
    (og_lgt_genes['G. muris'] == False)]

# Get the count of genes in this filtered subset
number_of_green_genes = len(green_subset)
print(f"Number of genes highlighted in green: {number_of_green_genes}")

trepo= og_lgt_genes[(og_lgt_genes['Trepomonas pc1-LGT'] == True)].drop_duplicates(subset= "id")
hin= og_lgt_genes[(og_lgt_genes['H. inflata'] == True)].drop_duplicates(subset= "id")
spiro= og_lgt_genes[(og_lgt_genes['S. salmonicida'] == True)].drop_duplicates(subset= "id")
wb= og_lgt_genes[(og_lgt_genes['G. intestinalis'] == True)].drop_duplicates(subset= "id")
muris= og_lgt_genes[(og_lgt_genes['G. muris'] == True)].drop_duplicates(subset= "id")


""" this is only printing the case where there is one single gene"""
print(" number of genes in OG from trepo= ", len(trepo))
print(" number of genes in OG from hin= ", len(hin))
print(" number of genes in OG from spiro= ", len(spiro))
print(" number of genes in OG from wb= ", len(wb))
print(" number of genes in OG from muris= ", len(muris))