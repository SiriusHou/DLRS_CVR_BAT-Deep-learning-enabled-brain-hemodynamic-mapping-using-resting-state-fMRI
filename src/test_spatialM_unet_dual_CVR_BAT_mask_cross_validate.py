# %%
import time
import os
import glob
import fnmatch
import numpy as np
import torch
import torch.nn as nn
from torch.utils import data
from torch.utils.data.sampler import SubsetRandomSampler
import matplotlib.pyplot as plt
import nibabel as nib
import math
from model_spatialM_unet_dual_CVR_BAT import UNet_dual

class CVRDataset(data.Dataset):
    'Characterizes a dataset for PyTorch'

    def __init__(self, data_root):
        'initialization'
        self.input_IDs = []
        self.mask_IDs = []
        'Find the file location'
        for subname in os.listdir(data_root):
            if os.path.isdir(os.path.join(data_root, subname)):
                for subfile in os.listdir(os.path.join(data_root, subname)):
                    if fnmatch.fnmatch(subfile, 'bold_spatialM_ce_nm_0[0-9][0-9].nii'):
                        sub_filepath = os.path.join(data_root, subname, subfile)
                        self.input_IDs.append(sub_filepath)

    def __len__(self):
        'Denotes the total number of samples'
        return len(self.input_IDs)

    def __getitem__(self, index):
        'Generates one sample of data'

        # Load input
        input_ID = self.input_IDs[index]

        parts = input_ID.split('.')
        ID = parts[0]

        img = nib.load(input_ID)
        train_image = img.get_fdata()
        train_image = train_image.astype(np.float32)
        img.uncache()

        train_image = np.transpose(train_image, (2, 0, 1))

        A_affine = img.affine
        A = torch.Tensor(train_image).type(torch.FloatTensor)

        return [ID, A, A_affine]

if __name__ == "__main__":

    for pod in range(5):
        # Creating data indices for training and validation splits:
        pod_test_name = 'pod_'+ str(pod+1)
        print(pod_test_name)

        test_dir = '/data1/xhou/ML_CVR/smooth8/ML_data/cross-validate/' + pod_test_name

        os.chdir(test_dir)

        train_dir = '/data1/xhou/ML_CVR/smooth8/ML_data/cross-validate'
        img_mean = nib.load(os.path.join(train_dir, 'CVR_mean.nii'))

        # Input test data
        test_set = CVRDataset(test_dir)
        test_set_size = len(test_set)

        test_loader = data.DataLoader(test_set, batch_size=1, shuffle=False, num_workers=10)

        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

        model_path = '/data1/xhou/ML_CVR/smooth8/ML_data/cross-validate/' + pod_test_name + '/Unet_dual_corr_beta0_spatialM_ce_nm_Epoch_100_AdaB1_1e12_final2_mask_soft_4.pt'
        model = UNet_dual(in_channels=136, n_classes=1, depth=5, batch_norm=True, padding=True, up_mode='upconv')
        state_dict = torch.load(model_path, map_location='cpu')

        model.load_state_dict(state_dict)
        model.to(device)

        model.eval()
        with torch.no_grad():
            for X_ID, X, X_affine in test_loader:
                X_dir = X_ID[0]
                # X = X[:, :-1, :, :]
                X = X.to(device)
                prediction = model(X)  # [N, H, W]

                cvr = np.squeeze(prediction[0, 0].cpu().clone().numpy())
                X_affine = np.squeeze(X_affine.cpu().clone().numpy())
                img_output = nib.Nifti1Image(cvr, X_affine)
                img_output.get_data_dtype() == np.dtype(np.int16)
                nib.save(img_output, '_'.join([X_dir, 'ML_Unet_dual_result_CVR.nii']))

                bat = np.squeeze(prediction[0, 1].cpu().clone().numpy())
                bat_output = nib.Nifti1Image(bat, X_affine)
                bat_output.get_data_dtype() == np.dtype(np.int16)
                nib.save(bat_output, '_'.join([X_dir, 'ML_Unet_dual_result_BAT.nii']))

        img_size = (91, 109, 91)
        for subname in glob.glob('*/'):

            img_cvr_data_3D = np.zeros(img_size)
            img_bat_data_3D = np.zeros(img_size)
            for subfile in os.listdir(os.path.join(test_dir, subname)):
                if fnmatch.fnmatch(subfile, 'bold_spatialM_ce_nm_[0-9][0-9][0-9]_ML_Unet_dual_result_CVR.nii'):
                    parts = subfile.split('_')
                    ID = int(parts[-6]) - 1

                    sub_filepath = os.path.join(test_dir, subname, subfile)
                    img = nib.load(sub_filepath)
                    img_data = img.get_fdata()
                    img.uncache()

                    img_cvr_data_3D[:, :, ID] = img_data[2:-3, 1:-2]

                if fnmatch.fnmatch(subfile, 'bold_spatialM_ce_nm_[0-9][0-9][0-9]_ML_Unet_dual_result_BAT.nii'):
                    parts = subfile.split('_')
                    ID = int(parts[-6]) - 1

                    sub_filepath = os.path.join(test_dir, subname, subfile)
                    img = nib.load(sub_filepath)
                    img_bat_data = img.get_fdata()
                    img.uncache()

                    img_bat_data_3D[:, :, ID] = img_bat_data[2:-3, 1:-2]

            img_output = nib.Nifti1Image(img_cvr_data_3D, img_mean.affine)
            img_output.get_data_dtype() == np.dtype(np.int16)
            nib.save(img_output, os.path.join(test_dir, subname, 'ML_result_spatialM_ce_nm_Unet_dual_CVR4.nii'))

            img_bat_output = nib.Nifti1Image(img_bat_data_3D, img_mean.affine)
            img_bat_output.get_data_dtype() == np.dtype(np.int16)
            nib.save(img_bat_output, os.path.join(test_dir, subname, 'ML_result_spatialM_ce_nm_Unet_dual_BAT4.nii'))

