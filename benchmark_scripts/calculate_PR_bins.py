#!/usr/local/bin/python3

import pandas as pd
import argparse
import numpy as np
import os
from sklearn import metrics

########################## Settings ###################################################

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--IN_GS", 
                dest='IN_GS', 
                help="Input GS file intersected with the bins", 
                required=True)

parser.add_argument("-p", "--IN_PREDS", 
                dest='IN_PREDS', 
                help="Input PREDS file intersected with the bins", 
                required=True)

parser.add_argument("-bins", "--IN_BINS", 
                dest='IN_BINS', 
                help="Input BINS", 
                required=True)

parser.add_argument("-TF", "--TF_NAME", 
                dest='TF_NAME', 
                help="Input BINS", 
                required=True)

parser.add_argument("-o", "--output_dir", 
                dest="OUT_DIR",
                help="The output directory",
                required=False)

args = parser.parse_args()

if not os.path.exists(args.OUT_DIR):
        os.makedirs(args.OUT_DIR)

########################## Functions ###################################################

class PREDICTIONS_PR:
    """
    This class object will perform the precision and recall analysis of binned ChIP-seq and binned TF
    binding predictions. 
    """
    def __init__(self, Predictions, GoldStandard):
        # This is the basename of the predictions. Could also be the project name instead        
        self.basename_PREDS = os.path.splitext(os.path.basename(Predictions))[0]
        
        print(self.basename_PREDS)

        self.predictions = Predictions
        
        self.GoldStandard = GoldStandard
        
        self.basename_GS = os.path.splitext(os.path.basename(GoldStandard))[0]
        
        self.PR_curve = []

        self.TF_NAME = args.TF_NAME

        self.Resolution = "w200"

        self.output_filename = (self.basename_PREDS + "_PR.tsv")                                                           
        
    def parse_GS(self, TOTAL_BINS_GENOME):
        GS_column_names = ["BIN_chr", "BIN_start", "BIN_stop", "Score"]

        GS_df = pd.read_csv(self.GoldStandard, sep="\t", header=None, usecols=[0, 1, 2, 3], names=GS_column_names, low_memory=False)

        GS_df["BIN"] = GS_df["BIN_chr"] + ":" + GS_df["BIN_start"].map(str) + "-" + GS_df["BIN_stop"].map(str)
        
        self.TOTAL_BINS_GS = len(GS_df["BIN"].unique())

        self.random_Precision = (self.TOTAL_BINS_GS/int(args.IN_BINS))

        self.GS_df= GS_df[["BIN", "Score"]]

        self.GS_gb = self.GS_df.groupby(["BIN"]).count()
        
        self.GS_gb_index = self.GS_gb.index

    def parse_PREDS(self):
        PREDS_df_columns = ["BIN_chr", "BIN_start", "BIN_stop", "Score"]

        PREDS_df = pd.read_csv(self.predictions, sep="\t", header=None, usecols=[0,1,2,3], names=PREDS_df_columns, low_memory=False)

        PREDS_df["BIN"] = PREDS_df["BIN_chr"] + ":" + PREDS_df["BIN_start"].map(str) + "-" + PREDS_df["BIN_stop"].map(str)

        self.PREDS_df = PREDS_df[["BIN", "Score"]]
        
        self.PREDS_df["Score"] = self.PREDS_df["Score"].round(3)

        self.PREDS_UNIQUE_RANKS = self.PREDS_df["Score"].sort_values(ascending=False).unique()
        
        self.PREDS_df.set_index("BIN", inplace=True)

        self.PREDS_UNIQUE_RANKS_len = len(self.PREDS_df["Score"].unique())

    def calculate_PR(self, rankings_list, rank_number):
        """This function will take two groupby objects from count_TF_bin_overlap and output a row of 
        statistics for the sample"""
        Temp_PR_DF = self.PREDS_df.loc[self.PREDS_df['Score'].isin(rankings_list)].index
        
        Number_of_predictions = len(Temp_PR_DF)

        GS_recovered = len(Temp_PR_DF.intersection(self.GS_gb_index))

        Precision = GS_recovered/Number_of_predictions

        Recall = GS_recovered/self.TOTAL_BINS_GS

        self.PR_curve.append(pd.DataFrame.from_dict({"Precision": [Precision],
                        "Recall": [Recall],
                        "TF_NAME": [self.TF_NAME],
                        "Rank": [rank_number],
                        "Random_Precision": [self.random_Precision],
                        "Number_Predictions": [Number_of_predictions],
                        "GS_Recovered": [GS_recovered],
                        "Unique_ranks": [self.PREDS_UNIQUE_RANKS_len],
                        "Unique_bins_GS": [self.TOTAL_BINS_GS],
                        "Resolution": [self.Resolution]}))

rankings_list = []

rank_number = 1

PR_obj = PREDICTIONS_PR(args.IN_PREDS, args.IN_GS)

PR_obj.parse_GS(args.IN_BINS)

PR_obj.parse_PREDS()

print ("Unique ranks: " + str(PR_obj.PREDS_UNIQUE_RANKS_len))

for rank in PR_obj.PREDS_UNIQUE_RANKS:

    rankings_list.append(rank)

    PR_obj.calculate_PR(rankings_list, rank_number)

    rank_number = rank_number + 1

PR_CURVE_DF = pd.concat(PR_obj.PR_curve)

if 0 not in PR_CURVE_DF["Recall"]:
    first_precision_point=PR_CURVE_DF["Precision"][0]

    modDfObj = dfObj.append({'Precision' : first_precision_point, 
                             'Recall' : 0, 
                             "TF_NAME": [self.TF_NAME], 
                             "Rank": [rank_number], 
                             "Random_Precision": [self.random_Precision],
                             "Number_Predictions": [Number_of_predictions],
                             "GS_Recovered": [GS_recovered],
                             "Unique_ranks": [self.PREDS_UNIQUE_RANKS_len],
                             "Unique_bins_GS": [self.TOTAL_BINS_GS],
                             "Resolution": [self.Resolution]}, ignore_index=True)

else:
    print("Pass")

try:
    PR_CURVE_DF["AUPR"] = metrics.auc(PR_CURVE_DF["Recall"], PR_CURVE_DF["Precision"])

except:
    PR_CURVE_DF["AUPR"] = 0

os.chdir(args.OUT_DIR)

PR_CURVE_DF.to_csv(PR_obj.output_filename, sep="\t", header=True, index=False)
