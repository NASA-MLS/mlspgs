# Copyright 2023, by the California Institute of Technology. ALL
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

#suppress warnings
import warnings
warnings.filterwarnings('ignore')

## =========================================================================
## Define the logistic sigmoid function
## =========================================================================
def _sigmoid(x):

  return(1./(1. + np.exp(-x)))

## =========================================================================
## Define the rectified linear unit activation function
## =========================================================================
def _relu(x):

  output =  x
  output[output<0] = 0.

  return(output)

## =========================================================================
## Define the hyperbolic tangent activation function
## =========================================================================
def _tanh(x):

    return ((np.exp(x) - np.exp(-x))/(np.exp(x) + np.exp(-x)))

## =========================================================================
## Define the artifical neural network fitting routine
## =========================================================================
def _ann_fit(features=None, \
             n_hl=None, \
             weights_hl1=None, \
             intercepts_hl1=None, \
             weights_hl2=None, \
             intercepts_hl2=None, \
             weights_hl3=None, \
             intercepts_hl3=None, \
             weights_hl4=None, \
             intercepts_hl4=None, \
             weights_ll=None, \
             intercepts_ll=None, \
             activation=None):

    n_samples, n_vars = np.shape(features)
    n_neurons = len(intercepts_hl1)
    n_surfs = len(intercepts_ll)
    
    result = np.zeros ((n_samples,n_surfs))
    for iSamples in range ( 0, n_samples ):
        
        gamma = 0.
        for iVars in range ( 0, n_vars ):
            gamma = gamma + (features[iSamples,iVars]*weights_hl1[iVars])
        
        if activation=='sigmoid':
            Gamma = _sigmoid(gamma+intercepts_hl1)
        if activation=='relu':
            Gamma = _relu(gamma+intercepts_hl1)
        if activation=='tanh':
            Gamma = _tanh(gamma+intercepts_hl1)
        
        if n_hl == 1:
            for iSurfs in range(0, n_surfs):
                result[iSamples,iSurfs] = np.sum( Gamma[:] * \
                                          weights_ll[:,iSurfs] ) + \
                                          intercepts_ll[iSurfs]
        else:
           
            gamma2 = 0.
            for iNeurons in range ( 0, n_neurons ):
                gamma2 = gamma2 + Gamma[iNeurons]*weights_hl2[iNeurons]
            
            if activation=='sigmoid':
                Gamma2 = _sigmoid(gamma2+intercepts_hl2)
            if activation=='relu':
                Gamma2 = _relu(gamma2+intercepts_hl2)
            if activation=='tanh':
                Gamma2 = _tanh(gamma2+intercepts_hl2)

            if n_hl == 2:
                for iSurfs in range(0, n_surfs):
                    result[iSamples,iSurfs] = np.sum( Gamma2[:] * \
                                              weights_ll[:,iSurfs] ) + \
                                              intercepts_ll[iSurfs]
            else:

                gamma3 = 0.
                for iNeurons in range ( 0, n_neurons ):
                    gamma3 = gamma3 + Gamma2[iNeurons]*weights_hl3[iNeurons]
                
                if activation=='sigmoid':
                    Gamma3 = _sigmoid(gamma3+intercepts_hl3)
                if activation=='relu':
                    Gamma3 = _relu(gamma3+intercepts_hl3)
                if activation=='tanh':
                    Gamma3 = _tanh(gamma3+intercepts_hl3)
        
                if n_hl == 3:
                    for iSurfs in range(0, n_surfs):
                        result[iSamples,iSurfs] = np.sum( Gamma3[:] * \
                                                  weights_ll[:,iSurfs] ) + \
                                                  intercepts_ll[iSurfs]
                else:
           
                    gamma4 = 0.
                    for iNeurons in range ( 0, n_neurons ):
                        gamma4 = gamma4 + Gamma3[iNeurons]*weights_hl4[iNeurons]
                    
                    if activation=='sigmoid':
                        Gamma4 = _sigmoid(gamma4+intercepts_hl4)
                    if activation=='relu':
                        Gamma4 = _relu(gamma4+intercepts_hl4)
                    if activation=='tanh':
                        Gamma4 = _tanh(gamma4+intercepts_hl4)

                    for iSurfs in range(0, n_surfs):
                        result[iSamples,iSurfs] = np.sum( Gamma4[:] * \
                                                  weights_ll[:,iSurfs] ) + \
                                                  intercepts_ll[iSurfs]

    return(result)

## =========================================================================
## Create a class that contains weights data
## =========================================================================
class _Weights_object():

    def __init__ ( self, \
                   Standardization_Brightness_Temperatures_Mean, \
                   Standardization_Brightness_Temperatures_Std, \
                   Standardization_Labels_Mean, \
                   Standardization_Labels_Std, \
                   Normalization_Labels_Min, \
                   Normalization_Labels_Max, \
                   Weights_Hidden_Layer_1, \
                   Intercepts_Hidden_Layer_1, \
                   Weights_Hidden_Layer_2, \
                   Intercepts_Hidden_Layer_2, \
                   Weights_Hidden_Labels_Layer, \
                   Intercepts_Hidden_Labels_Layer, \
                   Bands, \
                   MIFs, \
                   Channels_Band_1, \
                   Channels_Band_2, \
                   Channels_Band_3, \
                   Channels_Band_4, \
                   Latitude_Bins, \
                   Output_Pressure_Levels_Indices, \
                   Activation_Function ):

        self.bt_mean = Standardization_Brightness_Temperatures_Mean
        self.bt_std = Standardization_Brightness_Temperatures_Std
        self.labels_mean = Standardization_Labels_Mean
        self.labels_std = Standardization_Labels_Std
        self.labels_min = Normalization_Labels_Min
        self.labels_max = Normalization_Labels_Max
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
        self.lats = Latitude_Bins
        self.surfs = Output_Pressure_Levels_Indices
        self.actf = Activation_Function[0].decode('utf-8')

## =========================================================================
## Read weights
## =========================================================================
def _read_weights(weights_file=None):

    file = h5py.File(weights_file, 'r')

    # Features
    Standardization_Brightness_Temperatures_Mean = file['Standardization_Brightness_Temperatures_Mean'][:]
    Standardization_Brightness_Temperatures_Std = file['Standardization_Brightness_Temperatures_Std'][:]

    # Labels
    Standardization_Labels_Mean = file['Standardization_Labels_Mean'][:]
    Standardization_Labels_Std = file['Standardization_Labels_Std'][:]
    Normalization_Labels_Min = file['Normalization_Labels_Min'][:]
    Normalization_Labels_Max = file['Normalization_Labels_Max'][:]

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

    # Other information
    Latitude_Bins = file['Latitude_Bins'][:]
    Output_Pressure_Levels_Indices = file['Output_Pressure_Levels_Indices'][:]
    Activation_Function = file['Activation_Function'][:]
        
    file.close()

    data_out = _Weights_object(Standardization_Brightness_Temperatures_Mean, \
                               Standardization_Brightness_Temperatures_Std, \
                               Standardization_Labels_Mean, \
                               Standardization_Labels_Std, \
                               Normalization_Labels_Min, \
                               Normalization_Labels_Max, \
                               Weights_Hidden_Layer_1, \
                               Intercepts_Hidden_Layer_1, \
                               Weights_Hidden_Layer_2, \
                               Intercepts_Hidden_Layer_2, \
                               Weights_Hidden_Labels_Layer, \
                               Intercepts_Hidden_Labels_Layer, \
                               Bands, \
                               MIFs, \
                               Channels_Band_1, \
                               Channels_Band_2, \
                               Channels_Band_3, \
                               Channels_Band_4, \
                               Latitude_Bins, \
                               Output_Pressure_Levels_Indices, \
                               Activation_Function)

    return(data_out)

## =========================================================================
## Read L1BRAD data
## =========================================================================
def _read_l1brad(l1brad_file=None, \
                 l1b_bands=None):

    # Read file
    R = {}
    file = h5py.File(l1brad_file, 'r')
    for i_bands in range(len(l1b_bands)):
        R_D = file[l1b_bands[i_bands]][:,:,:]
        R_D2 = file[l1b_bands[i_bands]+' Baseline'][:,:]
        R_DD = np.zeros_like ( R_D )
        for i_mifs in range(0, len(R_DD[0,:,0])):
            R_DD[:,i_mifs,:] = np.array(R_D[:,i_mifs,:] + R_D2[:,:])
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
    l1boa_lat = file['/GHz/GeodLat'][:,:]
    file.close()

    return(l1boa_lat)

## =========================================================================
## Precision
## =========================================================================
def _define_precision(surfs=None):

    # The ANN precision is based on the RMSD between ANN and L2
    # from the test data set (averaged over all latitude bins).
    rmsd_ann = np.array([-9.99990000e+02, -9.99990000e+02, -9.99990000e+02, -9.99990000e+02,
                          4.20425664e-01,  4.35085096e-01,  4.27523999e-01,  4.66817550e-01,
                          4.97767777e-01,  5.43183935e-01,  5.86066063e-01,  4.69996715e-01,
                          4.68247813e-01,  4.92537926e-01,  5.45789214e-01,  6.25295421e-01,
                          7.75882944e-01,  1.08322854e+00, -9.99990000e+02, -9.99990000e+02,
                         -9.99990000e+02, -9.99990000e+02, -9.99990000e+02, -9.99990000e+02,
                         -9.99990000e+02, -9.99990000e+02, -9.99990000e+02, -9.99990000e+02,
                         -9.99990000e+02, -9.99990000e+02, -9.99990000e+02, -9.99990000e+02,
                         -9.99990000e+02, -9.99990000e+02, -9.99990000e+02, -9.99990000e+02,
                         -9.99990000e+02])

    # Here we use the annual median of precision per level, 
    # averaged over all years.
    rmsd_l2 = np.array([-9.99000000e+11, -9.99000000e+11, -9.99000000e+11, -9.99000000e+11,
                         4.52555556e-01,  4.87722222e-01,  4.99555556e-01,  5.31944444e-01,
                         5.97000000e-01,  6.96833333e-01,  8.34722222e-01,  6.83666667e-01,
                         7.55055556e-01,  8.40500000e-01,  9.44277778e-01,  1.08450000e+00,
                         1.34322222e+00,  1.79127778e+00, -9.99000000e+11, -9.99000000e+11,
                        -9.99000000e+11, -9.99000000e+11, -9.99000000e+11, -9.99000000e+11,
                        -9.99000000e+11, -9.99000000e+11, -9.99000000e+11, -9.99000000e+11,
                        -9.99000000e+11, -9.99000000e+11, -9.99000000e+11, -9.99000000e+11,
                        -9.99000000e+11, -9.99000000e+11, -9.99000000e+11, -9.99000000e+11,
                        -9.99000000e+11])

    # Total precision via RSS
    nn_prec = np.zeros ( (37) )
    nn_prec[:] = np.sqrt ( rmsd_l2**2 + rmsd_ann**2 )/1e9

    # Set other values to -999.99
    ind = np.where(rmsd_ann<-999)[0]
    nn_prec[ind] = -999.99

    return(nn_prec)

## =========================================================================
## Prediction wrapper for H2O
## =========================================================================
def _hno3_prediction(l1bradg_file=None, \
                   l1bradd_file=None, \
                   l1boa_file=None, \
                   weights_file=None, \
                   out_file=None):

    # Read weights
    weights = _read_weights(weights_file=weights_file)

    # Read L1b GHz data
    ghz_rad = _read_l1brad(l1brad_file=l1bradg_file, \
                           l1b_bands=weights.bands[0:4])

    # Read L1b GHz data
    l1boa_lat = _read_l1boa(l1boa_file=l1boa_file)
    
    # Slice radiance arrays for valid MIFs and channels.
    # For some reason I have to do dimensions individually
    b1 = ghz_rad[0][:,weights.mifs-1,:]
    b1 = b1[:,:,weights.ch_b1-1]

    b2 = ghz_rad[1][:,weights.mifs-1,:]
    b2 = b2[:,:,weights.ch_b2-1]

    b3 = ghz_rad[2][:,weights.mifs-1,:]
    b3 = b3[:,:,weights.ch_b3-1]

    b4 = ghz_rad[3][:,weights.mifs-1,:]
    b4 = b4[:,:,weights.ch_b4-1]
    
    # We average l1boa_lat over all valid MIFs
    lat = np.mean(l1boa_lat[0:len(b1),weights.mifs-1],axis=1)

    # Construct input matrix
    n_obs = len(b1)
    n_mifs = len(weights.mifs)
    n_ch1 = len(weights.ch_b1)
    n_ch2 = len(weights.ch_b2)
    n_ch3 = len(weights.ch_b3)
    n_ch4 = len(weights.ch_b4)

    features = np.reshape(b1,(n_obs,n_mifs*n_ch1))
    features = np.append(features,np.reshape(b2,(n_obs,n_mifs*n_ch2)),axis=1)
    features = np.append(features,np.reshape(b3,(n_obs,n_mifs*n_ch3)),axis=1)
    features = np.append(features,np.reshape(b4,(n_obs,n_mifs*n_ch4)),axis=1)

    # Predict
    pred = np.zeros ( (len(features),37) )
    pred[:,:] = -999.99
    for i_lat in range ( 0, len(weights.lats) ):
        lat_bin = weights.lats[i_lat]
        ind_lat = np.where((np.min(features,axis=1)>-100) & (lat>=lat_bin[0]) & (lat<lat_bin[1]))[0]
        if len(ind_lat)>0:

            x = (features[ind_lat,:] - weights.bt_mean[:,i_lat]) / \
                weights.bt_std[:,i_lat]

            dummy = _ann_fit(features=x, \
                             n_hl=2, \
                             weights_hl1=weights.w_l1[:,:,i_lat], \
                             intercepts_hl1=weights.i_l1[:,i_lat], \
                             weights_hl2=weights.w_l2[:,:,i_lat], \
                             intercepts_hl2=weights.i_l2[:,i_lat], \
                             weights_ll=weights.w_lo[:,:,i_lat], \
                             intercepts_ll=weights.i_lo[:,i_lat], \
                             activation=weights.actf)
            dummy = dummy * \
                weights.labels_std[:,i_lat] + \
                    weights.labels_mean[:,i_lat]
            pred[ind_lat,weights.surfs[0]-1:weights.surfs[-1]] = dummy

    # Precision
    prec = _define_precision(surfs=weights.surfs)
    prec = np.tile(prec,(len(features),1))

    # Check for faulty latitudes and set precision to -999.99.
    # Predictions are already -999.99.
    ind_lat = np.where((lat<-90))[0]
    if len(ind_lat)>0:
        prec[ind_lat,:] = -999.99

    # Check for unrealistic predictions outside the expected range
    for i_lat in range ( 0, len(weights.lats) ):
        lat_bin = weights.lats[i_lat]
        ind_lat = np.where((lat>=lat_bin[0]) & (lat<lat_bin[1]))[0]
        if len(ind_lat)>0:
            for i_surfs in range(weights.surfs[0]-1, weights.surfs[-1]):
                # We need to be careful: weights.labels_min and weights.labels_max
                # are only provided for valid surfs (4:18)
                min_clim = weights.labels_min[i_surfs-weights.surfs[0]+1,i_lat]
                max_clim = weights.labels_max[i_surfs-weights.surfs[0]+1,i_lat]
                delta = np.zeros_like(max_clim) #max_clim - min_clim
                thresh_min = min_clim - delta
                thresh_max = max_clim + delta

                ind_thresh = np.where( (pred[ind_lat,i_surfs]<thresh_min) | \
                                       (pred[ind_lat,i_surfs]>thresh_max) )[0]
                if len(ind_thresh)>0:
                    prec[ind_lat[ind_thresh],:] = -999.99

    # Output file
    file = h5py.File(out_file, 'w')
    dset = file.create_dataset('ANN_Prediction', data=pred, shape=[len(features),37], dtype=float)
    dset = file.create_dataset('ANN_Precision', data=prec, shape=[len(features),37], dtype=float)
    file.close()

    return ( None )

## =========================================================================
## Read the terminal input from the command line
## =========================================================================
args = {}
args['python_file'] =  sys.argv[0]
args['l1bradg_file'] =  sys.argv[1]
args['l1bradd_file'] =  sys.argv[2]
args['l1boa_file'] =  sys.argv[3]
args['weights_file'] =  sys.argv[4]
args['out_file'] =  sys.argv[5]
l1bradg_file = args['l1bradg_file']
l1bradd_file = args['l1bradd_file']
l1boa_file = args['l1boa_file']
weights_file = args['weights_file']
out_file = args['out_file']

## =========================================================================
## Run main
# Call: python hno3_prediction.py '/data/emls/nrt/v05.01.NRT.16/2022/272/MLS-Aura_L1BRADG_v05-01-NRT-16-c01_2022d272t2110.h5' '/data/emls/nrt/v05.01.NRT.16/2022/272/MLS-Aura_L1BRADD_v05-01-NRT-16-c01_2022d272t2110.h5' '/data/emls/nrt/v05.01.NRT.16/2022/272/MLS-Aura_L1BOA_v05-01-NRT-16-c01_2022d272t2110.h5' '/users/fwerner/Documents/database/neural_network_weights/weights/v05-0x/MLS-Aura_ANN-HNO3_v05-0x-01_20220930.h5' '/users/fwerner/Documents/software/python/ann/test_data/MLS-Aura_HNO3_v05-01-NRT-16-c01_2022d142t1410.h5'
## =========================================================================
result = _hno3_prediction(l1bradg_file=l1bradg_file, \
                          l1bradd_file=l1bradd_file, \
                          l1boa_file=l1boa_file, \
                          weights_file=weights_file, \
                          out_file=out_file)

## =========================================================================
## Revisions
## =========================================================================
## n/a
