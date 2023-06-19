# %%
import os
import fnmatch
import numpy as np
import torch
from torch.utils import data
import nibabel as nib
from model_spatialM_unet_dual_CVR_BAT import UNet_dual


class CVRDataset(data.Dataset):
    'Characterizes a dataset for PyTorch'

    def __init__(self, data_root):
        'initialization'
        self.input_IDs = []

        for dir in os.listdir(data_root):
            if os.path.isdir(os.path.join(data_root, dir)) and dir != "output":
                dlrs_input_path = os.path.join(data_root, dir, "DLRS_input")
                # Get all files in the DLRS_input folder
                files_in_dlrs_input = [
                    os.path.join(dlrs_input_path, file) for file in os.listdir(dlrs_input_path)
                    if fnmatch.fnmatch(file, 'DLRS_input_layer_[0-9][0-9][0-9].nii')
                ]
                # Add the files to the input_IDs list
                self.input_IDs.extend(files_in_dlrs_input)

    def __len__(self):
        'Denotes the total number of samples'
        return len(self.input_IDs)

    def __getitem__(self, index):
        'Generates one sample of data'

        # Load input
        input_ID = self.input_IDs[index]

        filename, extension = os.path.splitext(input_ID)

        img = nib.load(input_ID)
        train_image = img.get_fdata()
        train_image = train_image.astype(np.float32)
        img.uncache()

        train_image = np.transpose(train_image, (2, 0, 1))

        A_affine = img.affine
        A = torch.Tensor(train_image).type(torch.FloatTensor)

        return [filename, A, A_affine]


def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


if __name__ == "__main__":

    test_dir = './../data'
    os.chdir(test_dir)

    # Input test data
    test_set = CVRDataset(test_dir)
    test_set_size = len(test_set)

    test_loader = data.DataLoader(test_set, batch_size=1, shuffle=False, num_workers=10)

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    device = 'cpu'

    model_path = './../model/pretrained_DLRS_CVR_BAT_model.pt'
    model = UNet_dual(in_channels=136, n_classes=1, depth=5, batch_norm=True, padding=True, up_mode='upconv')
    state_dict = torch.load(model_path, map_location='cpu')

    model.load_state_dict(state_dict)
    model.to(device)

    model.eval()
    cvr = []
    bat = []
    with torch.no_grad():
        for X_ID, X, X_affine in test_loader:
            X_dir = X_ID[0]
            X = X.to(device)
            prediction = model(X)  # [N, H, W]

            cvr = np.squeeze(prediction[0, 0, 2:-3, 1:-2].cpu().clone().numpy())
            X_affine = np.squeeze(X_affine.cpu().clone().numpy())
            img_output = nib.Nifti1Image(cvr, X_affine)
            nib.save(img_output, '_'.join([X_dir, 'CVR.nii']))

            bat = np.squeeze(prediction[0, 1, 2:-3, 1:-2].cpu().clone().numpy())
            bat_output = nib.Nifti1Image(bat, X_affine)
            nib.save(bat_output, '_'.join([X_dir, 'BAT.nii']))

    output_dir = "./../data/output"
    create_folder(output_dir)

    for dir in os.listdir(test_dir):
        if os.path.isdir(os.path.join(test_dir, dir)) and dir != "output":
            create_folder('/'.join([output_dir, dir]))

            dlrs_input_path = os.path.join(test_dir, dir, "DLRS_input")

            # Get all 2D CVR files in the DLRS_input folder
            cvr_files = sorted([
                os.path.join(dlrs_input_path, file) for file in os.listdir(dlrs_input_path)
                if fnmatch.fnmatch(file, 'DLRS_input_layer_[0-9][0-9][0-9]_CVR.nii')
            ])

            cvr_3D = np.dstack([nib.load(cvr_file).get_fdata() for cvr_file in cvr_files])
            cvr_3D = nib.Nifti1Image(cvr_3D, X_affine)
            nib.save(cvr_3D, '/'.join([output_dir, dir, "DLRS_CVR.nii"]))

            # Get all 2D BAT files in the DLRS_input folder
            bat_files = sorted([
                os.path.join(dlrs_input_path, file) for file in os.listdir(dlrs_input_path)
                if fnmatch.fnmatch(file, 'DLRS_input_layer_[0-9][0-9][0-9]_BAT.nii')
            ])

            bat_3D = np.dstack([nib.load(bat_file).get_fdata() for bat_file in bat_files])
            bat_3D = nib.Nifti1Image(bat_3D, X_affine)
            nib.save(bat_3D, '/'.join([output_dir, dir, "DLRS_BAT.nii"]))
