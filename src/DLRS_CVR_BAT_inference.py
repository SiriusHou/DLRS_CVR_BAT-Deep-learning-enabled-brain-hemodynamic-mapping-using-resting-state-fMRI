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
        filename, extension = os.path.splitext(self.input_IDs[index])
        img = nib.load(self.input_IDs[index])

        image, img_affine = img.get_fdata(), img.affine
        image = np.transpose(image, (2, 0, 1))
        image = torch.Tensor(image).type(torch.FloatTensor)
        img.uncache()

        return [filename, image, img_affine]


def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def main():

    wkdir = './../data'
    model_path = './../model/pretrained_DLRS_CVR_BAT_model.pt'
    os.chdir(wkdir)

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # Load data
    inference_set = CVRDataset(wkdir)
    data_loader = data.DataLoader(inference_set, batch_size=1, shuffle=False, num_workers=4)

    # Load model
    model = UNet_dual(in_channels=136, n_classes=1, depth=5, batch_norm=True, padding=True, up_mode='upconv')
    state_dict = torch.load(model_path, map_location='cpu')
    model.load_state_dict(state_dict)
    model.to(device)

    model.eval()
    with torch.no_grad():
        for X_ID, X, X_affine in data_loader:
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

    output_dir = "/".join([wkdir, "output"])
    create_folder(output_dir)

    for subdir in os.listdir(wkdir):
        if os.path.isdir(os.path.join(wkdir, subdir)) and subdir != "output":
            create_folder('/'.join([output_dir, subdir]))

            dlrs_input_path = os.path.join(wkdir, subdir, "DLRS_input")

            # Get all 2D CVR files in the DLRS_input folder
            cvr_files = sorted([
                os.path.join(dlrs_input_path, file) for file in os.listdir(dlrs_input_path)
                if fnmatch.fnmatch(file, 'DLRS_input_layer_[0-9][0-9][0-9]_CVR.nii')
            ])

            cvr_3D = np.dstack([nib.load(cvr_file).get_fdata() for cvr_file in cvr_files])
            cvr_3D = nib.Nifti1Image(cvr_3D, X_affine)
            nib.save(cvr_3D, '/'.join([output_dir, subdir, "DLRS_CVR.nii"]))

            # Get all 2D BAT files in the DLRS_input folder
            bat_files = sorted([
                os.path.join(dlrs_input_path, file) for file in os.listdir(dlrs_input_path)
                if fnmatch.fnmatch(file, 'DLRS_input_layer_[0-9][0-9][0-9]_BAT.nii')
            ])

            bat_3D = np.dstack([nib.load(bat_file).get_fdata() for bat_file in bat_files])
            bat_3D = nib.Nifti1Image(bat_3D, X_affine)
            nib.save(bat_3D, '/'.join([output_dir, subdir, "DLRS_BAT.nii"]))


if __name__ == "__main__":
    main()
