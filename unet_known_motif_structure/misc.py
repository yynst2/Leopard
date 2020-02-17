import os
import gzip
import pickle
import numpy as np
import pandas as pd
import datetime
import string
import random
import argparse
import copy

# AUROC AUPRC
from sklearn.preprocessing import normalize
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import auc

# visualization
from IPython.display import SVG, display
from keras.utils.vis_utils import plot_model
from keras.utils.vis_utils import model_to_dot
import matplotlib
import matplotlib.pylab as plt

#
matplotlib.use('agg')

# functions
def saveObj(obj, filename):
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)
def loadObj(filename):
    with open(filename, 'rb') as f:
        _ = pickle.load(f)
    return _
def gen_roc_curve(true_label,
                  predicted_label,
                  outputfilename,
                  ):
    # Compute ROC curve and area the curve
    plt.figure()
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k',
             label='random guess', alpha=.8)
    for _class in range(2):
        fpr, tpr, thresholds = roc_curve(true_label[:, _class], predicted_label[:, _class])
        auc_value = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=2, alpha=0.3,
                 label='(class {0:d} AUC = {1:0.2f})'.format(_class, auc_value))
    plt.legend(loc='best')
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.savefig(outputfilename)
def gen_auprc_curve(true_label,
                    predicted_label,
                    outputfilename,
                    ):
    # Compute ROC curve and area the curve
    plt.figure()
    # plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',label='Chance', alpha=.8)
    for _class in range(2):
        precision, recall, thresholds = precision_recall_curve(true_label[:, _class], predicted_label[:, _class])
        auc_value = auc(recall, precision, True)
        plt.plot(recall, precision, lw=2, alpha=0.3,
                 label='(class {0:d} AUPRC = {1:0.2f})'.format(_class, auc_value))
    plt.legend(loc='best')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.savefig(outputfilename)

# keras visual class
class keras_visual():
    @staticmethod
    def plot_n_display(keras_model, output_file_path):
        plot_model(model=keras_model,
                   to_file=os.path.join(output_file_path, 'model.pdf'),
                   show_shapes=True)
        # visualize in notebook view
        display(SVG(model_to_dot(keras_model, show_shapes=True).create(prog='dot', format='svg')))
    @staticmethod
    def plot(keras_model, output_file_path, output_file_name):
        plot_model(model=keras_model,
                   to_file=os.path.join(output_file_path, output_file_name),
                   show_shapes=True)
    @staticmethod
    def display(keras_model):
        display(SVG(model_to_dot(keras_model, show_shapes=True).create(prog='dot', format='svg')))
def get_model_memory_usage(batch_size, model):
    import numpy as np
    from keras import backend as K

    shapes_mem_count = 0
    internal_model_mem_count = 0
    for l in model.layers:
        layer_type = l.__class__.__name__
        if layer_type == 'Model':
            internal_model_mem_count += get_model_memory_usage(batch_size, l)
        single_layer_mem = 1
        for s in l.output_shape:
            if s is None:
                continue
            single_layer_mem *= s
        shapes_mem_count += single_layer_mem

    trainable_count = np.sum([K.count_params(p) for p in set(model.trainable_weights)])
    non_trainable_count = np.sum([K.count_params(p) for p in set(model.non_trainable_weights)])

    number_size = 4.0
    if K.floatx() == 'float16':
        number_size = 2.0
    if K.floatx() == 'float64':
        number_size = 8.0

    total_memory = number_size * (batch_size * shapes_mem_count + trainable_count + non_trainable_count)
    gbytes = np.round(total_memory / (1024.0 ** 3), 3) + internal_model_mem_count
    return gbytes