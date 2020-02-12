#!/usr/bin/python3

import pandas as pd
import numpy as np
import os
import argparse
import sys

###-Environment setup-###
# Set up the argument parser with all of the options and defaults set up.
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--INPUT", 
                dest='IN_NPY', 
                help="Input NPY array", 
                required=True)

args = parser.parse_args()

def Leopard_predictions_to_BED(x):    
    basename_file = os.path.basename(x).split(".npy")[0]

    print (basename_file)
    
    predictions = np.load(x)

    # Round the prediction scores to 3 decimal places
    predictions = np.round(predictions, decimals=3)
       
    # Create a dataframe from the NP array
    preds_DF = pd.DataFrame(predictions)
    
    # add a column with the chromosome number, in this case we are only using chr1
    preds_DF["chr"] = "chr1"
    
    # add a column with the start position based on the index
    preds_DF["Start"] = preds_DF[0].index
    
    # add a column with the stop position 1 bp from the start
    preds_DF["Stop"] = preds_DF[0].index + 1
    
    # rename the columns
    preds_DF.columns = ["Score", "chr", "Start", "Stop"]
    
    # Reorganize the columns
    preds_DF = preds_DF[["chr", "Start", "Stop", "Score"]]
    
    # FIll na values with a 0; might need to change later to be more accurate
    preds_DF.fillna(0, inplace=True)

    # write to csv file. 
    preds_DF.to_csv(basename_file + "_1bp.bed", sep="\t", index=False, header=False)

Leopard_predictions_to_BED(args.IN_NPY)