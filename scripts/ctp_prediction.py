# Copyright 2024, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"
## =========================================================================
## Import of modules
## =========================================================================
import sys
import numpy as np
import h5py
import pickle

#suppress warnings
import warnings
warnings.filterwarnings('ignore')

## =========================================================================
## Define the logistic sigmoid function
## =========================================================================


def _sigmoid(x):

    return(1./(1. + np.exp(-x)))

## =========================================================================
## Define the softmax function
## =========================================================================


def _softmax(x):

    exp_x = np.exp(x - np.max(x))

    return(exp_x / np.sum(exp_x))

## =========================================================================
## Define the rectified linear unit activation function
## =========================================================================


def _relu(x):

    output = x
    output[output < 0] = 0.

    return(output)

## =========================================================================
## Define the hyperbolic tangent activation function
## =========================================================================


def _tanh(x):

    return ((np.exp(x) - np.exp(-x))/(np.exp(x) + np.exp(-x)))

## =========================================================================
## Define the artifical neural network fitting routine
## =========================================================================


def _ann_fit(features=None,
             n_hl=None,
             weights_hl1=None,
             intercepts_hl1=None,
             weights_hl2=None,
             intercepts_hl2=None,
             weights_hl3=None,
             intercepts_hl3=None,
             weights_hl4=None,
             intercepts_hl4=None,
             weights_ll=None,
             intercepts_ll=None,
             activation=None):

    n_samples, n_vars = np.shape(features)
    n_neurons = len(intercepts_hl1)
    n_surfs = len(intercepts_ll)

    result = np.zeros((n_samples, n_surfs))
    for iSamples in range(0, n_samples):

        gamma = 0.
        for iVars in range(0, n_vars):
            gamma = gamma + (features[iSamples, iVars]*weights_hl1[iVars])

        if activation == 'sigmoid':
            Gamma = _sigmoid(gamma+intercepts_hl1)
        if activation == 'relu':
            Gamma = _relu(gamma+intercepts_hl1)
        if activation == 'tanh':
            Gamma = _tanh(gamma+intercepts_hl1)

        if n_hl == 1:
            for iSurfs in range(0, n_surfs):
                result[iSamples, iSurfs] = np.sum(Gamma[:] *
                                                  weights_ll[:, iSurfs]) + \
                    intercepts_ll[iSurfs]
        else:

            gamma2 = 0.
            for iNeurons in range(0, n_neurons):
                gamma2 = gamma2 + Gamma[iNeurons]*weights_hl2[iNeurons]

            if activation == 'sigmoid':
                Gamma2 = _sigmoid(gamma2+intercepts_hl2)
            if activation == 'relu':
                Gamma2 = _relu(gamma2+intercepts_hl2)
            if activation == 'tanh':
                Gamma2 = _tanh(gamma2+intercepts_hl2)

            if n_hl == 2:
                for iSurfs in range(0, n_surfs):
                    result[iSamples, iSurfs] = np.sum(Gamma2[:] *
                                                      weights_ll[:, iSurfs]) + \
                        intercepts_ll[iSurfs]
            else:

                gamma3 = 0.
                for iNeurons in range(0, n_neurons):
                    gamma3 = gamma3 + Gamma2[iNeurons]*weights_hl3[iNeurons]

                if activation == 'sigmoid':
                    Gamma3 = _sigmoid(gamma3+intercepts_hl3)
                if activation == 'relu':
                    Gamma3 = _relu(gamma3+intercepts_hl3)
                if activation == 'tanh':
                    Gamma3 = _tanh(gamma3+intercepts_hl3)

                if n_hl == 3:
                    for iSurfs in range(0, n_surfs):
                        result[iSamples, iSurfs] = np.sum(Gamma3[:] *
                                                          weights_ll[:, iSurfs]) + \
                            intercepts_ll[iSurfs]
                else:

                    gamma4 = 0.
                    for iNeurons in range(0, n_neurons):
                        gamma4 = gamma4 + \
                            Gamma3[iNeurons]*weights_hl4[iNeurons]

                    if activation == 'sigmoid':
                        Gamma4 = _sigmoid(gamma4+intercepts_hl4)
                    if activation == 'relu':
                        Gamma4 = _relu(gamma4+intercepts_hl4)
                    if activation == 'tanh':
                        Gamma4 = _tanh(gamma4+intercepts_hl4)

                    for iSurfs in range(0, n_surfs):
                        result[iSamples, iSurfs] = np.sum(Gamma4[:] *
                                                          weights_ll[:, iSurfs]) + \
                            intercepts_ll[iSurfs]

    output = np.zeros((len(result), 2))
    for iSamples in range(0, n_samples):
        if np.isfinite(result[iSamples, 0]):
            output[iSamples] = _softmax(result[iSamples])
        else:
            output[iSamples] = [np.nan,np.nan]
    return(output)

## =========================================================================
## Create a class that contains weights data
## =========================================================================


class _Weights_object():

    def __init__(self,
                 Standardization_Brightness_Temperatures_Mean_CM,
                 Standardization_Brightness_Temperatures_Std_CM,
                 Normalization_Brightness_Temperatures_Min_CTP,
                 Normalization_Brightness_Temperatures_Max_CTP,
                 Normalization_Labels_Min_CTP,
                 Normalization_Labels_Max_CTP,
                 Weights_Hidden_Layer_1,
                 Intercepts_Hidden_Layer_1,
                 Weights_Hidden_Layer_2,
                 Intercepts_Hidden_Layer_2,
                 Weights_Hidden_Labels_Layer,
                 Intercepts_Hidden_Labels_Layer,
                 Bands,
                 MIFs,
                 Channels_Band_1,
                 Channels_Band_2,
                 Channels_Band_3,
                 Channels_Band_4,
                 Channels_Band_5,
                 Channels_Band_6,
                 Channels_Band_7,
                 Channels_Band_8,
                 Channels_Band_9,
                 Channels_Band_10,
                 Activation_Function):

        self.bt_mean_cm = Standardization_Brightness_Temperatures_Mean_CM
        self.bt_std_cm = Standardization_Brightness_Temperatures_Std_CM
        self.bt_min_ctp = Normalization_Brightness_Temperatures_Min_CTP
        self.bt_max_ctp = Normalization_Brightness_Temperatures_Max_CTP
        self.lb_min_ctp = Normalization_Labels_Min_CTP
        self.lb_max_ctp = Normalization_Labels_Max_CTP
        self.w_l1 = Weights_Hidden_Layer_1
        self.i_l1 = Intercepts_Hidden_Layer_1
        self.w_l2 = Weights_Hidden_Layer_2
        self.i_l2 = Intercepts_Hidden_Layer_2
        self.w_lo = Weights_Hidden_Labels_Layer
        self.i_lo = Intercepts_Hidden_Labels_Layer
        band_dummy = []
        for i_band in range(0, len(Bands)):
            band_dummy.append(Bands[i_band].decode('utf-8'))
        self.bands = np.array(band_dummy)
        self.mifs = MIFs
        self.ch_b1 = Channels_Band_1
        self.ch_b2 = Channels_Band_2
        self.ch_b3 = Channels_Band_3
        self.ch_b4 = Channels_Band_4
        self.ch_b5 = Channels_Band_5
        self.ch_b6 = Channels_Band_6
        self.ch_b7 = Channels_Band_7
        self.ch_b8 = Channels_Band_8
        self.ch_b9 = Channels_Band_9
        self.ch_b10 = Channels_Band_10
        self.actf = Activation_Function[0].decode('utf-8')

## =========================================================================
## Read weights
## =========================================================================


def _read_weights(weights_file=None):

    file = h5py.File(weights_file, 'r')

    # Features CM
    Standardization_Brightness_Temperatures_Mean_CM = file[
        'Standardization_Brightness_Temperatures_Mean_CM'][:]
    Standardization_Brightness_Temperatures_Std_CM = file[
        'Standardization_Brightness_Temperatures_Std_CM'][:]

    # Features CTP
    Normalization_Brightness_Temperatures_Min_CTP = file[
        'Normalization_Brightness_Temperatures_Min_CTP'][:]
    Normalization_Brightness_Temperatures_Max_CTP = file[
        'Normalization_Brightness_Temperatures_Max_CTP'][:]

    # Labels
    Normalization_Labels_Min_CTP = file['Normalization_Labels_Min_CTP'][:]
    Normalization_Labels_Max_CTP = file['Normalization_Labels_Max_CTP'][:]

    # Weights
    Weights_Hidden_Layer_1 = file['Weights_Hidden_Layer_1'][:]
    Intercepts_Hidden_Layer_1 = file['Intercepts_Hidden_Layer_1'][:]
    Weights_Hidden_Layer_2 = file['Weights_Hidden_Layer_2'][:]
    Intercepts_Hidden_Layer_2 = file['Intercepts_Hidden_Layer_2'][:]
    Weights_Hidden_Labels_Layer = file['Weights_Hidden_Labels_Layer'][:]
    Intercepts_Hidden_Labels_Layer = file['Intercepts_Hidden_Labels_Layer'][:]

    # Radiance information
    Bands = file['Bands'][:]
    MIFs = file['MIFs'][:]
    Channels_Band_1 = file['Channels_Band_#1'][:]
    Channels_Band_2 = file['Channels_Band_#2'][:]
    Channels_Band_3 = file['Channels_Band_#3'][:]
    Channels_Band_4 = file['Channels_Band_#4'][:]
    Channels_Band_5 = file['Channels_Band_#5'][:]
    Channels_Band_6 = file['Channels_Band_#6'][:]
    Channels_Band_7 = file['Channels_Band_#7'][:]
    Channels_Band_8 = file['Channels_Band_#8'][:]
    Channels_Band_9 = file['Channels_Band_#9'][:]
    Channels_Band_10 = file['Channels_Band_#10'][:]

    # Other information
    Activation_Function = file['Activation_Function'][:]

    file.close()

    data_out = _Weights_object(Standardization_Brightness_Temperatures_Mean_CM,
                               Standardization_Brightness_Temperatures_Std_CM,
                               Normalization_Brightness_Temperatures_Min_CTP,
                               Normalization_Brightness_Temperatures_Max_CTP,
                               Normalization_Labels_Min_CTP,
                               Normalization_Labels_Max_CTP,
                               Weights_Hidden_Layer_1,
                               Intercepts_Hidden_Layer_1,
                               Weights_Hidden_Layer_2,
                               Intercepts_Hidden_Layer_2,
                               Weights_Hidden_Labels_Layer,
                               Intercepts_Hidden_Labels_Layer,
                               Bands,
                               MIFs,
                               Channels_Band_1,
                               Channels_Band_2,
                               Channels_Band_3,
                               Channels_Band_4,
                               Channels_Band_5,
                               Channels_Band_6,
                               Channels_Band_7,
                               Channels_Band_8,
                               Channels_Band_9,
                               Channels_Band_10,
                               Activation_Function)

    return(data_out)

## =========================================================================
## Read L1BRAD data
## =========================================================================


def _read_l1brad(l1brad_file=None,
                 l1b_bands=None):

    # Read file
    R = {}
    file = h5py.File(l1brad_file, 'r')
    for i_bands in range(len(l1b_bands)):
        R_D = file[l1b_bands[i_bands]][:, :, :]
        R_D2 = file[l1b_bands[i_bands]+' Baseline'][:, :]
        R_DD = np.zeros_like(R_D)
        for i_mifs in range(0, len(R_DD[0, :, 0])):
            R_DD[:, i_mifs, :] = np.array(R_D[:, i_mifs, :] + R_D2[:, :])
            R[i_bands] = R_DD
    file.close()

    return(R)

## =========================================================================
## Read L1BOA data
## =========================================================================


def _read_l1boa(l1boa_file=None):

    # Read file.
    # We ar eonly interested in GeodLat of the GHz receiver
    file = h5py.File(l1boa_file, 'r')
    l1boa_lat = file['/GHz/GeodLat'][:, :]
    l1boa_lon = file['/GHz/Lon'][:, :]
    l1boa_alt = file['/GHz/GeodAlt'][:, :]
    file.close()

    return(l1boa_lat,
           l1boa_lon,
           l1boa_alt)

## =========================================================================
## Precision
## =========================================================================


def _define_precision():

    # The GBDT precision is based on the RMSD between ANN and colocated MODIS
    # retrievals from the test data set.
    gbdt_prec = np.array([70.73858429])

    return(gbdt_prec)

## =========================================================================
## Prediction wrapper for H2O
## =========================================================================


def _ctp_prediction(l1bradg_file=None,
                    l1boa_file=None,
                    weights_file=None,
                    decision_file=None,
                    out_file=None):

    # Read weights
    weights = _read_weights(weights_file=weights_file)

    # Read L1b GHz data
    ghz_rad = _read_l1brad(l1brad_file=l1bradg_file,
                           l1b_bands=weights.bands)
    # Read L1b GHz data
    l1boa_lat, l1boa_lon, l1boa_alt = _read_l1boa(l1boa_file=l1boa_file)

    # We average l1boa_lat and l1boa_lon over all MIFs
    lat = np.mean(l1boa_lat[0:len(ghz_rad[0]), :], axis=1)
    lon = np.mean(l1boa_lon[0:len(ghz_rad[0]), :], axis=1)

    # Construct input matrix
    n_obs = len(ghz_rad[0])
    n_mifs = len(weights.mifs)
    n_ch1 = len(weights.ch_b1)
    n_ch2 = len(weights.ch_b2)
    n_ch3 = len(weights.ch_b3)
    n_ch4 = len(weights.ch_b4)
    n_ch5 = len(weights.ch_b5)
    n_ch6 = len(weights.ch_b6)
    n_ch7 = len(weights.ch_b7)
    n_ch8 = len(weights.ch_b8)
    n_ch9 = len(weights.ch_b9)
    n_ch10 = len(weights.ch_b10)
    n_feat = n_ch1+n_ch2+n_ch3+n_ch4+n_ch5+n_ch6+n_ch7+n_ch8+n_ch9+n_ch10
    n_feat *= n_mifs

    bt_dummy = 1
    bt1 = np.zeros((n_obs, n_feat))
    bt2 = np.zeros((n_obs, n_feat))
    bt3 = np.zeros((n_obs, n_feat))
    ind_l1b = np.arange(0, n_obs, 1)
    count_b = 0

    for i_bands in range(1, 11):
        n_ch = np.copy(n_ch1)
        if i_bands == 5:
            n_ch = np.copy(n_ch5)
        if i_bands == 10:
            n_ch = np.copy(n_ch10)
        for i_mifs in range(0, len(weights.mifs)):
            for i_ch in range(0, n_ch):
                if i_bands == 1:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b1[i_ch]-1]
                if i_bands == 2:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b2[i_ch]-1]
                if i_bands == 3:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b3[i_ch]-1]
                if i_bands == 4:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b4[i_ch]-1]
                if i_bands == 5:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b5[i_ch]-1]
                if i_bands == 6:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b6[i_ch]-1]
                if i_bands == 7:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b7[i_ch]-1]
                if i_bands == 8:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b8[i_ch]-1]
                if i_bands == 9:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b9[i_ch]-1]
                if i_bands == 10:
                    bt_dummy = ghz_rad[i_bands-1][:,
                                                  weights.mifs[i_mifs]-1, weights.ch_b10[i_ch]-1]
                bt1[:, count_b] = bt_dummy[np.roll(ind_l1b, 1)]
                bt2[:, count_b] = bt_dummy[ind_l1b]
                bt3[:, count_b] = bt_dummy[np.roll(ind_l1b, -1)]
                count_b += 1

    features = np.zeros((n_obs, 3*n_feat))
    features[:, 0:n_feat] = bt1
    features[:, n_feat:2*n_feat] = bt2
    features[:, 2*n_feat:3*n_feat] = bt3
    features = np.append(features, np.reshape(
        lat, (n_obs, 1)), axis=1)
    features = np.append(features, np.reshape(
        lon, (n_obs, 1)), axis=1)
    features = np.append(features, np.reshape(
        l1boa_alt[:, weights.mifs-1], (n_obs, n_mifs)), axis=1)

    # Predict cloud flag
    pred_cm = np.zeros((len(features), 1))
    pred_cm[:, 0] = -999.99

    x = (features[:, :3*n_feat] - weights.bt_mean_cm[:]) / weights.bt_std_cm[:]

    dummy = _ann_fit(features=x,
                     n_hl=2,
                     weights_hl1=weights.w_l1[:, :],
                     intercepts_hl1=weights.i_l1[:],
                     weights_hl2=weights.w_l2[:, :],
                     intercepts_hl2=weights.i_l2[:],
                     weights_ll=weights.w_lo[:, :],
                     intercepts_ll=weights.i_lo[:],
                     activation=weights.actf)
    
    for i_pred in range(0, n_obs):
        pred_cm[i_pred, 0] = np.argmax(dummy[i_pred])
        if np.isnan(dummy[i_pred, 0]):
            pred_cm[i_pred, 0] = -999.99

    # Predict cloud top pressure
    ctp_model = pickle.load(open(decision_file, 'rb'))

    pred = np.zeros((len(features), 1))
    pred[:, 0] = -999.99

    x = (features[:, :] - weights.bt_min_ctp[:]) / \
        (weights.bt_max_ctp[:] - weights.bt_min_ctp[:])

    dummy = ctp_model.predict(x)
    ind_cm = np.where(pred_cm[:, 0] == 1)[0]
    if len(ind_cm) > 0:
        pred[ind_cm, 0] = 10**dummy[ind_cm]

    # Precision
    prec = _define_precision()
    prec = np.tile(prec, (len(features), 1))

    # Set first and last value to -999.99
    prec[0, 0] = -999.99
    prec[-1, 0] = -999.99

    # Check for faulty latitudes and set precision to -999.99.
    ind_lat = np.where((lat < -90))[0]
    if len(ind_lat) > 0:
        prec[ind_lat, :] = -999.99

    # Check for unrealistic predictions outside the expected range
    ind_valid = np.where(pred[:, 0] > 0)[0]
    ind_thresh_lb = np.where((pred[ind_valid, 0] < 0.9*10**weights.lb_min_ctp) |
                          (pred[ind_valid, 0] > 1.1*10**weights.lb_max_ctp))[0]
    if len(ind_thresh_lb) > 0:
        prec[ind_valid[ind_thresh_lb], 0] = -999.99

    # We see whether for each observation there is an outlier feature
    delta = weights.bt_max_ctp[:]-weights.bt_min_ctp[:]
    thresh_min = weights.bt_min_ctp[:]-0.25*delta
    thresh_max = weights.bt_max_ctp[:]+0.25*delta
    for i_obs in range(0, n_obs):
        ind_thresh_bt = np.where(
            (features[i_obs, :]-thresh_min < 0) | (features[i_obs, :]-thresh_max > 0) )[0]
        if len(ind_thresh_bt) > 0:
            prec[i_obs, 0] = -999.99

    # Output file
    file = h5py.File(out_file, 'w')
    dset = file.create_dataset('ANN_Prediction', data=pred, shape=[
                               len(features), 1], dtype=float)
    dset = file.create_dataset('ANN_Precision', data=prec, shape=[
                               len(features), 1], dtype=float)
    file.close()

    return (None)


## =========================================================================
## Read the terminal input from the command line
## =========================================================================
args = {}
args['python_file'] = sys.argv[0]
args['l1bradg_file'] = sys.argv[1]
args['l1boa_file'] = sys.argv[2]
args['weights_file'] = sys.argv[3]
args['decision_file'] = sys.argv[4]
args['out_file'] = sys.argv[5]
l1bradg_file = args['l1bradg_file']
l1boa_file = args['l1boa_file']
weights_file = args['weights_file']
decision_file = args['decision_file']
out_file = args['out_file']

## =========================================================================
## Run main
# Call: python ctp_prediction.py '/data/emls/l1b/v05.00/2019/239/MLS-Aura_L1BRADG_v05-00-c01_2019d239.h5' '/data/emls/l1b/v05.00/2019/239/MLS-Aura_L1BOA_v05-00-c01_2019d239.h5' '/users/fwerner/Documents/database/neural_network_weights/weights/v05-0x/MLS-Aura_ANN-CM_v05-0x-02_20240308.h5' '/users/fwerner/Documents/database/neural_network_weights/weights/v05-0x/MLS-Aura_ANN-CTP_v05-0x-02_20240308.pkl' '/users/fwerner/Documents/software/python/ann/test_data/MLS-Aura_CTP_v05-00-c01_2019d239.h5'
## =========================================================================
result = _ctp_prediction(l1bradg_file=l1bradg_file,
                         l1boa_file=l1boa_file,
                         weights_file=weights_file,
                         decision_file=decision_file,
                         out_file=out_file)

## =========================================================================
## Revisions
## =========================================================================
## Revision 1.1  2024/03/06 fwerner
## Tweaked thresh_min and thresh_max from +/-0.1*delta to +/-0.25*delta.
## Those deltas are still <10% of the absolute values. We still
## eliminate noticeable outliers, while allowing for some reasonable
## extrapolation. Note that this affects obs where <=3 features exceed thresh.
##
## Revision 1.2  2024/03/15 fwerner
## Added a check for NAN before returning results in _ann_fit. This is 
## necessary because _softmax sets all outputs to NAN if there is even a
## single NAN in the input. pred_cm[i_pred, 0] = -999.99 for NAN predictions.
##
