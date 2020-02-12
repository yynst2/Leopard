import os
import sys
import numpy as np
from keras.models import Model
from keras.layers import Input, concatenate, Conv1D, MaxPooling1D, Conv2DTranspose,Lambda,BatchNormalization
from keras.optimizers import Adam
import tensorflow as tf
from keras import backend as K
import keras
K.set_image_data_format('channels_last')  # TF dimension ordering in this code

ss=10


def gelu(x):
    """Gaussian Error Linear Unit.
    This is a smoother version of the RELU.
    Original paper: https://arxiv.org/abs/1606.08415
    Args:
        x: float Tensor to perform activation.
    Returns:
        `x` with the GELU activation applied.

    from google BERT model

    this could be useful as negative correction can be potentially found between PWMs
    """
    cdf = 0.5 * (1.0 + tf.tanh(
        (np.sqrt(2 / np.pi) * (x + 0.044715 * tf.pow(x, 3)))))
    return x * cdf


def crossentropy_cut(y_true,y_pred):
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    y_pred_f= tf.clip_by_value(y_pred_f, 1e-7, (1. - 1e-7))
    mask=K.greater_equal(y_true_f,-0.5)
    losses = -(y_true_f * K.log(y_pred_f) + (1.0 - y_true_f) * K.log(1.0 - y_pred_f))
    losses = tf.boolean_mask(losses, mask)
    masked_loss = tf.reduce_mean(losses)
    return masked_loss

def Conv1DTranspose(input_tensor, filters, kernel_size, strides=2, padding='same'):
    x = Lambda(lambda x: K.expand_dims(x, axis=2))(input_tensor)
    x = Conv2DTranspose(filters=filters, kernel_size=(kernel_size, 1), strides=(strides, 1), padding=padding)(x)
    x = Lambda(lambda x: K.squeeze(x, axis=2))(x)
    return x

def dice_coef(y_true, y_pred):
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    mask=K.cast(K.greater_equal(y_true_f,-0.5),dtype='float32')
    intersection = K.sum(y_true_f * y_pred_f * mask)
    return (2. * intersection + ss) / (K.sum(y_true_f * mask) + K.sum(y_pred_f * mask) + ss)

def dice_coef_loss(y_true, y_pred):
    return -dice_coef(y_true, y_pred)

def pcc(layer_in, num_filter, size_kernel, activation='relu', padding='same'):
    x = MaxPooling1D(pool_size=2)(layer_in)
    x = BatchNormalization()(Conv1D(num_filter,size_kernel,activation=activation,padding=padding)(x))
    x = BatchNormalization()(Conv1D(num_filter,size_kernel,activation=activation,padding=padding)(x))
    return x

def ucc(layer_in1,layer_in2, num_filter, size_kernel, activation='relu', padding='same'):
    x = concatenate([Conv1DTranspose(layer_in1,num_filter,2,strides=2,padding=padding), layer_in2], axis=2)
    x = BatchNormalization()(Conv1D(num_filter,size_kernel,activation=activation,padding=padding)(x))
    x = BatchNormalization()(Conv1D(num_filter,size_kernel,activation=activation,padding=padding)(x))
    return x

def get_unet(the_lr=1e-1, num_class=2, num_channel=10, size=2048*25,known_motif_array_shape=[]):
    inputs = Input((size, num_channel))
    #inputs2=Input((size,num_channel-2))
#    print(inputs.shape)

    num_blocks=5 
    initial_filter=15
    scale_filter=1.5
    size_kernel=7
    activation='relu'
    padding='same'    

    layer_down=[]
    layer_up=[]


    inputs2=Lambda(lambda x: K.slice(x,[0,0,0],[-1,-1,4]), name='DNA_slice')(inputs)

    # first conv+bn
    conv0 = BatchNormalization()(Conv1D(initial_filter, size_kernel,
                                        activation=activation, padding=padding)(inputs))

    # # known motif conv+bn
    conv1 =Conv1D(filters=known_motif_array_shape[-1], kernel_size=known_motif_array_shape[1],
           padding=padding,name='known_motif_scan',use_bias=False,activation=gelu)(inputs2)

    # concat known motif into original layer
    conv0 = keras.layers.concatenate(axis=-1, inputs=[conv0, conv1],name='concat_known_motif_results')

    # second conv+bn
    layer_down.append(conv0)
    num=initial_filter

    for i in range(num_blocks):
        num=int(num * scale_filter)
        the_layer=pcc(layer_down[i], num, size_kernel, activation=activation, padding=padding)
        layer_down.append(the_layer)

    layer_up.append(the_layer)
    for i in range(num_blocks):
        num=int(num / scale_filter)
        the_layer=ucc(layer_up[i],layer_down[-(i+2)],num, size_kernel, activation=activation, padding=padding)
        layer_up.append(the_layer)

    convn = Conv1D(num_class, 1, activation='sigmoid', padding=padding)(layer_up[-1])

    model = Model(inputs=[inputs], outputs=[convn])

    model.compile(optimizer=Adam(lr=the_lr,beta_1=0.9, beta_2=0.999,decay=1e-5), loss=crossentropy_cut, metrics=[dice_coef])



    return model


