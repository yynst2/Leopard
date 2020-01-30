import pyBigWig
import argparse
import os
import numpy as np
import gzip


class peakclass:
    # constant
    size=10240
    num_channel=6

    def __init__(self,peakline,chr_map):
        self.peakline = peakline.split('\t')
        self.chr = self.peakline[0]
        self.start = int(self.peakline[1])
        self.end = int(self.peakline[2])
        self.label=self.chr+'_'+str(self.start)+'_'+str(self.end)
        self.length=self.end-self.start
        self.peak_center=int((self.end+self.start)/2)

        if int(self.peak_center-0.5*peakclass.size)<0:
            self.padded_start= 0
            self.padded_end = int(self.padded_start + peakclass.size)
        elif int(self.peak_center+0.5*peakclass.size)>=chr_map[self.chr]:
            self.padded_end = chr_map[self.chr]
            self.padded_start = self.padded_end  -peakclass.size
        else:
            self.padded_start= int(self.peak_center - 0.5 * peakclass.size)
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
    parser.add_argument('-datapath', '--datapath', default='./data', type=str,
                        help='data folder')
    parser.add_argument('-te', '--test', default='K562', type=str,
                        help='test cell type')
    parser.add_argument('-tepath', '--testpath', default='dnase_bigwig', type=str,
                        help='test cell type folder')
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
path2 = os.path.join(args.datapath,args.testpath)    # dnase or atac signal path
path3= os.path.join(args.datapath,args.peakpath)   # peak file path
path4=os.path.join(args.datapath,args.outpath)     # output path

peak_features = []
peakvec = []
# open bigwig files
list_dna = ['A', 'C', 'G', 'T']
dict_dna = {}
for the_id in list_dna:
    dict_dna[the_id] = pyBigWig.open(os.path.join(path1, the_id + '.bigwig'))
feature_avg = pyBigWig.open(os.path.join(path2,'avg.bigwig'))        # load bigwig of average data
feature_test = pyBigWig.open(os.path.join(path2,args.test+'.bigwig'))  # load bigwig of test data

# define boundary
chr_all = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
           'chr22', 'chrX']
num_bp = np.array(
    [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210,
     78077248, 59128983, 63025520, 48129895, 51304566, 155270560])
chr_map = {}
for idx, i in enumerate(chr_all):
    chr_map[i]=num_bp[idx]

# loading peaks:
if args.peak.find('gz')>=0:
    with gzip.open(os.path.join(path3,args.peak),'r') as f:
        data=f.readlines()
        for i in data:
            peak=peakclass(i.decode().rstrip(),chr_map)
            peakvec.append(peak)
else:
    with open(os.path.join(path3,args.peak),'r') as f:
        data=f.readlines()
        for i in data:
            peak=peakclass(i.rstrip(),chr_map)
            peakvec.append(peak)

print('total peak number: {}'.format(len(peakvec)))
# get features
for idx,i in enumerate(peakvec):
    # print('handling : {} {} {} {}                                   '.format(idx,i.chr,i.start,i.end), end="\r")
    print('handling : {}'.format(idx), end="\r")
    i.derive_feature(list_dna,dict_dna,feature_avg,feature_test)
    peak_features.append(i.feature)

# close bigwig files
for the_id in list_dna:
    dict_dna[the_id].close()
feature_avg.close()
feature_test.close()

# output
peak_features=np.array(peak_features)
os.system('mkdir -p '+path4)
np.save(os.path.join(path4,args.out+'.npy'),peak_features)   # peak number * size * channel

