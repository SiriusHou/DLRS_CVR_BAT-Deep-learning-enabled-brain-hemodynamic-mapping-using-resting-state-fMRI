import os
import fnmatch
import numpy as np
import torch
import torch.nn as nn
from torch.utils import data
from torch.utils.data.sampler import SubsetRandomSampler
import matplotlib.pyplot as plt
import nibabel as nib
import math
from DLRS_CVR_BAT_model import UNet_dual
from adabelief_pytorch import AdaBelief


class CVRDataset(data.Dataset):
    'Characterizes a dataset for PyTorch'

    def __init__(self, data_root):
        'initialization'
        self.input_IDs = []
        self.cvr_IDs = []
        self.bat_IDs = []
        'Find the file location'
        for dir in os.listdir(data_root):
            if os.path.isdir(os.path.join(data_root, dir)) and dir != "output":
                dlrs_input_path = os.path.join(data_root, dir, "DLRS_input")
                # Get all files in the DLRS_input folder
                for file in os.listdir(dlrs_input_path):
                    # Add the DL input files to the input_IDs list
                    if fnmatch.fnmatch(file, 'DLRS_input_layer_[0-9][0-9][0-9].nii'):
                        self.input_IDs.append(os.path.join(data_root, dir, "DLRS_input", file))

                        # Add the corresponding HC CVR/BAT files to the cvr_IDs/bat_IDs list
                        index_with_ext = file.split('_')[-1]
                        self.cvr_IDs.append(os.path.join(data_root, dir, 'CVR_HC', 'CVR_HC_clean_' + index_with_ext))
                        self.bat_IDs.append(os.path.join(data_root, dir, 'BAT_HC', 'BAT_HC_clean_' + index_with_ext))


    def __len__(self):
        'Denotes the total number of samples'
        return len(self.input_IDs)

    def __getitem__(self, index):
        'Generates one sample of data'
        input_ID = self.input_IDs[index]

        # Load input
        input_file = nib.load(input_ID)
        input = input_file.get_fdata()
        input = np.transpose(input, (2, 0, 1))
        input = torch.Tensor(input).type(torch.FloatTensor)
        input_file.uncache()

        # Load HC CVR
        cvr_ID = self.cvr_IDs[index]
        cvr = nib.load(cvr_ID).get_fdata()
        cvr = torch.Tensor(cvr).type(torch.FloatTensor).unsqueeze(0)

        # Load HC BAT
        bat_ID = self.bat_IDs[index]
        bat = nib.load(bat_ID).get_fdata()
        bat = torch.Tensor(bat).type(torch.FloatTensor).unsqueeze(0)

        label = torch.cat((cvr, bat), dim=0)
        return [input, label]


def main():

    validation_split = 0.10
    random_seed = 53
    shuffle_dataset = True

    training_set = CVRDataset('./../data')
    training_set_size = len(training_set)
    indices = list(range(training_set_size))
    split = int(np.floor(validation_split * training_set_size))
    if shuffle_dataset:
        np.random.seed(random_seed)
        np.random.shuffle(indices)

    train_indices, val_indices = indices[split:], indices[:split]

    train_sampler = SubsetRandomSampler(train_indices)
    valid_sampler = SubsetRandomSampler(val_indices)

    # Generators
    train_loader = data.DataLoader(training_set, batch_size=64, sampler=train_sampler, num_workers=12)
    validation_loader = data.DataLoader(training_set, batch_size=64, sampler=valid_sampler, num_workers=12)
    model = UNet_dual(in_channels=136, n_classes=1, depth=5, batch_norm=True, padding=True, up_mode='upconv')

    # single gpu
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    gpu_index = 1
    model = model.to(device)

    optim = AdaBelief(model.parameters(), lr=0.00002, eps=1e-12)

    criterion = nn.L1Loss()
    epochs = 100
    print_every = 10

    train_loss_list = []
    test_loss_list = []

    cc_cost_list = []
    for t in range(epochs):

        train_loss = 0.0
        running_loss = 0.0
        step_count = 0.0

        test_loss = 0.0
        test_total_loss = 0.0
        test_step_count = 0.0

        for i, (X, y) in enumerate(train_loader):
            X = X.to(device)
            y = y.to(device)
            prediction = model(X)
            loss = criterion(prediction, y)

            optim.zero_grad()
            loss.backward()
            optim.step()

            if math.isnan(loss.item()) == False:
                train_loss += loss.item()
                step_count += 1

            if (i + 1) % print_every == 0:
                running_loss += train_loss
                print('Epoch number {}, Step [{}/{}], Training Loss: {:.4f}'
                      .format(t + 1, i + 1, len(train_loader), loss.item()))
                train_loss = 0.0

                # validation
                model.eval()
                with torch.no_grad():
                    for X, y in validation_loader:
                        X = X.to(device)
                        y = y.to(device)
                        prediction = model(X)
                        loss = criterion(prediction, y)

                        if math.isnan(loss.item()) == False:
                            test_loss = loss.item()
                            test_step_count += 1

                        test_total_loss += test_loss
                        test_loss = 0

                print('Epoch number {}, Step [{}/{}], Validation Loss: {:.4f}'
                      .format(t + 1, i + 1, len(train_loader), loss.item()))

        test_loss_list.append(test_total_loss / test_step_count)
        train_loss_list.append(running_loss / step_count)

        if (t + 1) % 10 == 0:

            model_parts = ('./../model/Unet_dual_Epoch_', str(t + 1), '.pt')

            model_name = ''.join(model_parts)
            if gpu_index == 2:
                torch.save(model.module.state_dict(), model_name)
            else:
                torch.save(model.state_dict(), model_name)

    ## Plotting batch-wise train loss curve:
    plt.plot(train_loss_list, '-o', label='train_loss', color='blue')
    plt.plot(test_loss_list, '-o', label='valid_loss', color='orange')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.savefig('./../model/Unet_dual_loss_CVR_BAT.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
