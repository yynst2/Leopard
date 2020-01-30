#!/home/hyangl/anaconda3/bin/python
import pyBigWig
import argparse
import os
import sys
import numpy as np
import tensorflow as tf

#tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)  ## disable tf warning!
stderr = sys.stderr  # disable printing "Using TensorFlow backend."
sys.stderr = open(os.devnull, 'w')
import keras

sys.stderr = stderr
from keras import backend as K
import unet
from util.auc import score_record, calculate_auc

K.set_image_data_format('channels_last')  # TF dimension ordering in this code

###### PARAMETER ###############

size = 2 ** 11 * 5  # 10240
num_channel = 6
write_pred = True  # whether generate .vec prediction file
size_edge = int(100)  # chunk edges to be excluded
batch = 100


# argv
def get_args():
    parser = argparse.ArgumentParser(description="predict")
    parser.add_argument('-m', '--model', default='weights_1.h5', type=str,
                        help='model name')
    parser.add_argument('-tf', '--transcription_factor', default='CTCF', type=str,
                        help='transcript factor')
    parser.add_argument('-inputpath', '--inputpath', default='preloaded_features', type=str,
                        help='preloaded input data path')
    parser.add_argument('-input', '--input', default='', type=str,
                        help='preloaded input data for prediction')
    parser.add_argument('-te', '--test', default='K562', type=str,
                        help='test cell type')
    parser.add_argument('-para', '--parallel', default=1, type=int,
                        help='control GPU memory usage when running multiple models in parallel')
    args = parser.parse_args()
    return args


args = get_args()

path_computer = 'data/'
path1 = path_computer + 'dna_bigwig/'                   # dna
# path2 = path_computer + 'dnase_bigwig/'               # dnase
path4=  os.path.join(path_computer,args.inputpath)      # input data path

from keras.backend.tensorflow_backend import set_session

config = tf.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 1 / float(args.parallel) - 0.05
set_session(tf.Session(config=config))

# print(sys.argv)
name_model = args.model
the_tf = args.transcription_factor
cell_test = args.test

model = unet.get_unet(the_lr=1e-3, num_class=1, num_channel=num_channel, size=size)
model.load_weights(name_model)
# model.summary()

test_data=np.load(os.path.join(path4,args.input))    # peak number * size * channel
shape = np.shape(test_data)
print('Number of predictions: {}'.format(shape[0]))
if __name__ == '__main__':

    the_name = name_model.split('/')[-1].split('.')[0]


    ## make predictions ################

    output = model.predict(test_data)                            # output: peak number * size * channel (=1)
    output = np.reshape(output,(shape[0],shape[1]))              # output:  peak number * size

    if write_pred:
        os.system('mkdir -p output')
        np.save('./output/pred_' + the_tf + '_' + cell_test + '_Full_' + the_name, output)
