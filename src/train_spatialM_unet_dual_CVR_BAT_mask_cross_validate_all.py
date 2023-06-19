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
from adabelief_pytorch import AdaBelief

# cwd = os.getcwd()
# print(cwd)
class CVRDataset(data.Dataset):
    'Characterizes a dataset for PyTorch'

    def __init__(self, data_root):
        'initialization'
        self.input_IDs = []
        self.label_IDs = []
        self.label_re_IDs = []
        self.bat_IDs = []
        self.bat_re_IDs = []
        'Find the file location'
        for pod in os.listdir(data_root):
            if os.path.isdir(os.path.join(data_root, pod)):
                for subname in os.listdir(os.path.join(data_root, pod)):
                    if os.path.isdir(os.path.join(data_root, pod, subname)):
                        for subfile in os.listdir(os.path.join(data_root, pod, subname)):
                            if fnmatch.fnmatch(subfile, 'bold_spatialM_ce_nm_[0-2][0-9][0-9].nii') or fnmatch.fnmatch(subfile, 'bold_spatialM_ce_nm_[0-2][0-9][0-9][xyz].nii'):
                                sub_filepath = os.path.join(data_root, pod, subname, subfile)
                                self.input_IDs.append(sub_filepath)

                                parts = sub_filepath.split('/')
                                parts_prefix = '/'.join(parts[:-1])
                                parts_end = parts[-1]
                                parts_end = parts_end.split('_')
                                parts_end = parts_end[-1]

                                self.label_IDs.append(os.path.join(parts_prefix, 'HC', 'CVR_HC_clean_' + parts_end))
                                self.label_re_IDs.append(os.path.join(parts_prefix, 'HC', 'CVR_HC_clean_re_' + parts_end))

                                self.bat_IDs.append(os.path.join(parts_prefix, 'bat_v3', 'bat_HC_clean_' + parts_end))
                                self.bat_re_IDs.append(os.path.join(parts_prefix, 'bat_v3', 'bat_HC_clean_re_' + parts_end))

    def __len__(self):
        'Denotes the total number of samples'
        return len(self.input_IDs)

    def __getitem__(self, index):
        'Generates one sample of data'

        'initialization'
        input_x = []
        label_y = []

        # Select input
        input_ID = self.input_IDs[index]
        # print(input_ID)

        # Load input
        img = nib.load(input_ID)
        train_image = img.get_fdata()
        train_image = train_image.astype(np.float32)
        #         img.uncache()

        # # For single channel
        # train_image = np.expand_dims(train_image, axis=2)

        # Transpose
        train_image = np.transpose(train_image, (2, 0, 1))

        A = torch.Tensor(train_image).type(torch.FloatTensor)

        # Select label
        label_ID = self.label_IDs[index]
        #         print(label_ID)

        # Load label
        taimg = nib.load(label_ID)
        label_image = taimg.get_fdata()
        label_image = label_image.astype(np.float32)
        #         taimg.uncache()
        B = torch.Tensor(label_image).type(torch.FloatTensor).unsqueeze(0)

        label_re_ID = self.label_re_IDs[index]

        taimg_re = nib.load(label_re_ID)
        label_re_image = taimg_re.get_fdata()
        label_re_image = label_re_image.astype(np.float32)
        B_re = torch.Tensor(label_re_image).type(torch.FloatTensor).unsqueeze(0)

        # Select bat
        bat_ID = self.bat_IDs[index]
        #         print(label_ID)

        # Load label
        batimg = nib.load(bat_ID)
        bat_image = batimg.get_fdata()
        bat_image = bat_image.astype(np.float32)
        C = torch.Tensor(bat_image).type(torch.FloatTensor).unsqueeze(0)

        bat_re_ID = self.bat_re_IDs[index]

        batimg_re = nib.load(bat_re_ID)
        bat_re_image = batimg_re.get_fdata()
        bat_re_image = bat_re_image.astype(np.float32)
        C_re = torch.Tensor(bat_re_image).type(torch.FloatTensor).unsqueeze(0)

        D = torch.cat((B, C, B_re, C_re), dim=0)
        return [A, D]


# %%

if __name__ == "__main__":

    validation_split = 0.10
    random_seed = 53
    shuffle_dataset = True

    training_set = CVRDataset('/data1/xhou/ML_CVR/smooth8/ML_data/cross-validate')
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
    train_loader = data.DataLoader(training_set, batch_size=4, sampler=train_sampler, num_workers=12)
    validation_loader = data.DataLoader(training_set, batch_size=4, sampler=valid_sampler, num_workers=12)
    model = UNet_dual(in_channels=136, n_classes=1, depth=5, batch_norm=True, padding=True, up_mode='upconv')

    # single gpu
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    gpu_index = 1
    model = model.to(device)

    # #multiple gpu
    # device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    # model = model.to(device)
    # if torch.cuda.device_count() > 1:
    #     gpu_index = 2
    #     print("Let's use", torch.cuda.device_count(), "GPUs!")
    #     model = nn.DataParallel(model)  # enabling data parallelism

    # optim = torch.optim.Adam(model.parameters(), lr=0.0001)
    optim = AdaBelief(model.parameters(), lr=0.00002, eps=1e-12)

    # optim = torch.optim.SGD(model.parameters(), lr=0.0001, momentum=0.9)
    criterion = nn.L1Loss()
    epochs = 1
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
            # X = X[:, :-1, :, :]
            X = X.to(device)  # [N, 6, H, W]
            # y = y[:, :-1, :, :]
            y = y.to(device)  # [N, H, W] with class indices (0, 1)
            prediction = model(X)  # [N, H, W]
            # loss = criterion(prediction, y[:, :-2, :, :]) + 4*criterion(prediction * y[:, 2:, :, :], y[:, :-2, :, :] * y[:, 2:, :, :])
            mask_index = torch.sum(y[:, 2:, :, :])
            if mask_index == 0:
                print('Warning!')
            loss = criterion(prediction, y[:, :-2, :, :]) + 4*torch.sum(torch.abs( (prediction- y[:, :-2, :, :]) * y[:, 2:, :, :]))/(torch.sum(y[:, 2:, :, :])+0.00001)

            #cross-correlation
            cvr_pred = prediction[:, 0, :, :]
            cvr_pred_mean = torch.mean(torch.mean(cvr_pred, 2), 1).unsqueeze(1).unsqueeze(2)
            cvr_pred_mean = cvr_pred_mean.repeat(1, 96, 112)
            vx = cvr_pred - cvr_pred_mean

            cvr_y = y[:, 0, :, :]
            cvr_y_mean = torch.mean(torch.mean(cvr_y, 2), 1).unsqueeze(1).unsqueeze(2)
            cvr_y_mean = cvr_y_mean.repeat(1, 96, 112)
            vy = cvr_y - cvr_y_mean

            cc_cost = 0
            cc_index = 0
            for jj in range(cvr_pred.shape[0]):
                vx_1 = vx[jj, :, :]
                vy_1 = vy[jj, :, :]
                cc_cost += torch.sum(vx_1 * vy_1) / (torch.sqrt(torch.sum(vx_1 ** 2)) * torch.sqrt(torch.sum(vy_1 ** 2)) +0.00001)
                cc_index += 1

            cc_cost = cc_cost/cc_index

            optim.zero_grad()
            loss.backward()
            optim.step()

            if math.isnan(loss.item()) == False:
                train_loss += loss.item()
                step_count += 1

            if (i + 1) % print_every == 0:
                running_loss += train_loss
                print('Epoch number {}, Step [{}/{}], Training Loss: {:.4f}, CC Loss: {:.4f}'
                      .format(t + 1, i + 1, len(train_loader), loss.item(), cc_cost.item()))
                train_loss = 0.0

                # validation
                model.eval()
                with torch.no_grad():
                    for X, y in validation_loader:
                        # X = X[:, :-1, :, :]
                        X = X.to(device)  # [N, 6, H, W]
                        # y = y[:, :-1, :, :]
                        y = y.to(device)  # [N, H, W] with class indices (0, 1)
                        prediction = model(X)  # [N, H, W]
                        # loss = criterion(prediction, y[:, :-2, :, :]) + 4*criterion(prediction * y[:, 2:, :, :],
                        #                                                           y[:, :-2, :, :] * y[:, 2:, :, :])
                        loss = criterion(prediction, y[:, :-2, :, :]) + 4*torch.sum(
                            torch.abs( (prediction - y[:, :-2, :, :]) * y[:, 2:, :, :]) ) / (torch.sum(y[:, 2:, :, :])+0.00001)

                        # cross-correlation
                        cvr_pred = prediction[:, 0, :, :]
                        cvr_pred_mean = torch.mean(torch.mean(cvr_pred, 2), 1).unsqueeze(1).unsqueeze(2)
                        cvr_pred_mean = cvr_pred_mean.repeat(1, 96, 112)
                        vx = cvr_pred - cvr_pred_mean

                        cvr_y = y[:, 0, :, :]
                        cvr_y_mean = torch.mean(torch.mean(cvr_y, 2), 1).unsqueeze(1).unsqueeze(2)
                        cvr_y_mean = cvr_y_mean.repeat(1, 96, 112)
                        vy = cvr_y - cvr_y_mean

                        cc_cost = 0
                        cc_index = 0
                        for jj in range(cvr_pred.shape[0]):
                            vx_1 = vx[jj, :, :]
                            vy_1 = vy[jj, :, :]
                            cc_cost += torch.sum(vx_1 * vy_1) / (
                                        torch.sqrt(torch.sum(vx_1 ** 2)) * torch.sqrt(torch.sum(vy_1 ** 2)) + 0.00001)
                            cc_index += 1

                        cc_cost = cc_cost / cc_index

                        if math.isnan(loss.item()) == False:
                            test_loss = loss.item()
                            test_step_count += 1

                        test_total_loss += test_loss
                        test_loss = 0

                # cc_cost_list.append(cc_cost)
                print('Epoch number {}, Step [{}/{}], Validation Loss: {:.4f}, CC Loss: {:.4f}'
                      .format(t + 1, i + 1, len(train_loader), loss.item(), cc_cost.item()))

        test_loss_list.append(test_total_loss / test_step_count)
        train_loss_list.append(running_loss / step_count)

        if (t + 1) % 10 == 0:
            # model_parts = ('/data1/xhou/ML_CVR/smooth8/ML_data/train/Unet_dual_corr_beta0_spatialM_ce_aal_Epoch_', str(t+1) , '_AdaB1_1e12_new_nomask.pt')

            model_parts = ('/data1/xhou/ML_CVR/smooth8/ML_data/cross-validate/Unet_dual_corr_beta0_spatialM_ce_nm_Epoch_', str(t+1) , '_AdaB1_1e12_final3_mask_soft_4.pt')
            # model_parts = ('/data1/xhou/ML_CVR/smooth8/ML_data/train/Unet_dual_corr_beta0_spatialM_ce_nm_Epoch_', str(t+1) , '_AdaB1_1e12_final_mask_soft_4.pt')
            # model_parts = ('/data1/xhou/ML_CVR/smooth8/ML_data/train/Unet_dual_corr_CVR_spatialM_ce_nm_Epoch_', str(t+1) , '_AdaB1_1e12_final_mask_soft_2.pt')

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
    plt.savefig('/data1/xhou/ML_CVR/smooth8/ML_data/train/Unet_dual_corr_bold0_spatialM_ce_nm_Epoch_100_valid_loss_CVR_BAT.png', dpi=300)
    plt.show()
