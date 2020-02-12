#!/usr/bin/python2

# This script has to be run using python 2 or the mapping of TF names to motifs 
# throw an error because of how lists + dictionaries work in python3. I need
# to find a workaround so we can use python3 for everything. 

import os
import argparse
import pandas as pd
import numpy as np
import multiprocessing
from joblib import Parallel, delayed

###-Environment setup-###
# Set up the argument parser with all of the options and defaults set up.
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--MOODS", 
                dest='IN_MOODS', 
                help="Input MOODS file: sep = |", 
                required=True)

parser.add_argument("-T", "--TF_META", 
                dest='TF_META', 
                help="Input TF meta file: sep = \t", 
                required=True)

parser.add_argument("-K", "--Keep_list", 
                dest='Keep_LIST', 
                help="Input TF list to keep: 1 column : 1 TF per row", 
                required=True)

parser.add_argument("-o", "--OUT_DIR", 
                dest='OUT_DIR', 
                help="Output dir", 
                required=True)

args = parser.parse_args()

# create output directory if it does not exist
if not os.path.exists(args.OUT_DIR):
        os.makedirs(args.OUT_DIR)

basename = os.path.splitext(os.path.basename(args.IN_MOODS))[0]

num_cores = multiprocessing.cpu_count()

class MOODS:
    """
    This class object will parse the MOODS output: 
    MOODS Object: parse MOODS output so every TF has the correct coordinates
    MAP TF names to motifs: Map the names associated with each motif
    Remove TFs not in list: Will remove TFs that were expanded with the motif meta data
    write bed files: write a bed file for each TF
    """
    def __init__(self, moods_input, tf_meta):
        column_names = ["ID", "Motif", "TF_pos", 'Match_Score', "Motif_sequence"]
        column_info = {
            'ID':             "category",
            'Motif':          "category",
            'TF_pos':         "uint16",
            'Match_Score':    "float",
            'Motif_sequence': "category"
            }
       
        #self.MOODS_DF holds the moods data
        self.MOODS_DF = pd.read_csv(moods_input, 
                         usecols=[0,1,2,4,5], 
                         names=column_names, 
                         dtype=column_info, 
                         header=None, 
                         sep="|",
                         low_memory=False)

        #Get chromosome and coordinates from ID column
        self.MOODS_DF["Chr"], self.MOODS_DF["Coord"] = self.MOODS_DF["ID"].str.split(':').str

        #Get the TF start position
        self.MOODS_DF["TF_start"] = self.MOODS_DF["Coord"].str.split("-").str[0].apply(int) + self.MOODS_DF["TF_pos"]

        #Get the TF stop position
        self.MOODS_DF["TF_end"] = self.MOODS_DF["TF_start"] + self.MOODS_DF["Motif_sequence"].str.len()
                
        #Drop all columns except these listed
        self.MOODS_DF = self.MOODS_DF[["Chr", "TF_start", "TF_end", "Motif", "Match_Score", "ID"]]

        #Create the df to hold the TF-motif meta data file
        df_TF_meta = pd.read_csv(tf_meta, 
                                 sep="\t", 
                                 header=0, 
                                 low_memory=False)

        #Group the motifs and TF names together
        df_gb = df_TF_meta.groupby(["Motif_ID", "TF_Name"]).count()

        #Create a dictionary of motif-TF associations, key: motif, value: list of TF names and then map
        self.MOODS_DF["Motif"] = self.MOODS_DF["Motif"].map({k:list(df_gb.loc[k].index) for k in df_gb.index.levels[0]})

        #Drop rows that have na in the Motif columns
        self.MOODS_DF.dropna(subset=["Motif"], inplace=True)

    def expand_list(self, x, fill_value='', preserve_index=False):
        # make sure `x` is list-alike
        if (x is not None
            and len(x) > 0
            and not isinstance(x, (list, tuple, np.ndarray, pd.Series))):
            x = [x]
        
        # all columns except `x`
        idx_cols = self.MOODS_DF.columns.difference(x)
        
        # calculate lengths of lists
        lens = self.MOODS_DF[x[0]].str.len()
        
        # preserve original index values    
        idx = np.repeat(self.MOODS_DF.index.values, lens)
        
        # create "exploded" DF
        self.MOODS_DF = (pd.DataFrame({
                    col:np.repeat(self.MOODS_DF[col].values, lens)
                    for col in idx_cols},
                    index=idx)
                .assign(**{col:np.concatenate(self.MOODS_DF.loc[lens>0, col].values)
                                for col in x}))
        
        # append those rows that have empty lists
        if (lens == 0).any():
            # at least one list in cells is empty
            self.MOODS_DF = (self.MOODS_DF.append(self.MOODS_DF.loc[lens==0, idx_cols], sort=False)
                    .fillna(fill_value))
        
        # revert the original index order
        self.MOODS_DF = self.MOODS_DF.sort_index()
        
        # reset index if requested
        if not preserve_index:        
            self.MOODS_DF = self.MOODS_DF.reset_index(drop=True)

        #Clean up dataframe        
        self.MOODS_DF = self.MOODS_DF[["Chr", "TF_start", "TF_end", "Motif", "Match_Score", "ID"]]
        
        #Create a list of unique TFs in the dataframe
        self.TF_list = self.MOODS_DF["Motif"].unique().tolist()

    def removeTFs(self, x):   
        #read file of TFs to keep into a python list     
        keep_list = list(line.rstrip('\n') for line in open(x))

        #remove TFs not in list
        self.MOODS_DF = self.MOODS_DF.loc[self.MOODS_DF['Motif'].isin(keep_list)]

        #Update TF_list to only include those passing filter
        self.TF_list = self.MOODS_DF["Motif"].unique().tolist()

    def write_bed(self, basename, TF):
        #Create a dataframe of only the TF in loop, sort, drop duplicates, and write
        temp_MOODS_DF = self.MOODS_DF[self.MOODS_DF.Motif == TF].copy()
                    
        temp_MOODS_DF.sort_values(by=['Chr', "TF_start", "TF_end", "Match_Score"], inplace=True)
    
        temp_MOODS_DF.drop_duplicates(subset=['Chr', "TF_start", "TF_end"], keep='first', inplace=True)
        #temp_MOODS_DF.groupby(['Chr', "TF_start", "TF_end", "Motif", "ID"], as_index=False).max()

        temp_MOODS_DF.to_csv(basename + "_" + TF + ".bed.gz", 
                                sep="\t", 
                                index=False, 
                                header=False, 
                                compression='gzip')

###-Start-###
print("Creating MOODS Object")
MOODS_obj = MOODS(args.IN_MOODS, args.TF_META)

print("Expanding Lists")
MOODS_obj.expand_list("Motif")

print("Removing TFs not in list")
MOODS_obj.removeTFs(args.Keep_LIST)

os.chdir(args.OUT_DIR)

#For every unique TF in the dataframe, write a bed file for each of its motif positions
if __name__ == "__main__":
    print("Write BED for each TF")

    Parallel(n_jobs=num_cores)(delayed(MOODS_obj.write_bed)(basename, i) for i in MOODS_obj.TF_list)
