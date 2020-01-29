import pyBigWig
import argparse
import os
import numpy as np
import gzip

class peakclass:
    size=10240
    num_channel=6
    peak_features=[]
    peakvec = []
    def __init__(self,peakline):
        self.peakline = peakline.split('\t')
        self.chr = self.peakline[0]
        self.start = self.peakline[1]
        self.end = self.peakline[2]
        self.label=self.chr+'_'+str(self.start)+'_'+str(self.end)
        self.lenth=self.end-self.start
        self.peak_center=int((self.end+self.start)/2)
        self.padded_start=self.peak_center-0.5*peakclass.size
        self.padded_end=self.padded_start+peakclass.size

        # features
        self.feature=[]   # 10240 * 6

    def derive_feature(self,list_dna,dict_dna,feature_avg,feature_test):
        image = np.zeros((peakclass.num_channel, peakclass.size))
        num=0
        for j in np.arange(len(list_dna)):
            the_id = list_dna[j]
            image[num, :] = dict_dna[the_id].values(self.chr, self.padded_start, self.padded_end)
            num += 1
        # feature & diff
        image[num, :] = np.array(feature_test.values(self.chr, self.padded_start, self.padded_end))
        avg = np.array(feature_avg.values(self.chr, self.padded_start, self.padded_end))
        image[num + 1, :] = image[num, :] - avg

        self.feature=image.T   # size * channel

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-datapath', '--datapath', default='../data', type=str,
                        help='data folder')
    parser.add_argument('-tepath', '--tepath', default='dnase_bigwig', type=str,
                        help='test cell type folder')
    parser.add_argument('-te', '--test', default='K562', type=str,
                        help='test cell type')
    parser.add_argument('-peak', '--peak', default='', type=str,
                        help='peak file name')
    parser.add_argument('-peakpath', '--peakpath', default='', type=str,
                        help='peak file folder')
    parser.add_argument('-out', '--out', default='', type=str,
                        help='output file name')
    parser.add_argument('-outpath', '--outpath', default='', type=str,
                        help='output file folder')
    args=parser.parse_args()
    return args

args = get_args()
path1 = os.path.join(args.datapath, 'dna_bigwig')  # dna path
path2 = os.path.join(args.datapath,args.tepath)    # dnase or atac signal path
path3= os.path.join(args.datapath,args.peakpath)   # peak file path
path4=os.path.join(args.datapath,args.outpath)     # output path


# open bigwig files
list_dna = ['A', 'C', 'G', 'T']
dict_dna = {}
for the_id in list_dna:
    dict_dna[the_id] = pyBigWig.open(path1 + the_id + '.bigwig')
feature_avg = pyBigWig.open(os.path.join(path2,'avg.bigwig'))        # load bigwig of average data
feature_test = pyBigWig.open(os.path.join(path2,args.te+'.bigwig'))  # load bigwig of test data

# loading peaks:
if args.peak.find('gz')>=0:
    with gzip.open(os.path.join(path3,args.peak),'rb') as f:
        data=f.readlines()
        for i in data:
            peak=peakclass(i)
            peakclass.peakvec.append(peak)
else:
    with open(os.path.join(path3,args.peak),'r') as f:
        data=f.readlines()
        for i in data:
            peak=peakclass(i)
            peakclass.peakvec.append(peak)

# get features
for i in peakclass.peakvec:
    i.derive_feature(list_dna,dict_dna,feature_avg,feature_test)
    peakclass.peak_features.append(i.feature)

# close bigwig files
for the_id in list_dna:
    dict_dna[the_id].close()
feature_avg.close()
feature_test.close()

# output
peakclass.peak_features=np.array(peakclass.peak_features)
os.system('mkdir -p '+path4)
np.save(os.path.join(path4,args.out),peakclass.peak_features)   # peak number * size * channel

