function calc_FD(subjdir, rp_file_name, thresh_FD)

subjdir = deblank(subjdir);
rp = [subjdir, filesep, rp_file_name];
disp(rp)

%%Framewise Displacement (FD, Power et al. 2012)
filerp = fopen(rp,'r');
rp_data = fscanf(filerp,'%f %f %f %f %f %f',[6,inf])';
dynnum = length(rp_data);
fclose(filerp);

alpha = rp_data(:,4); %change radian to mm (Power et al)
beta = rp_data(:,5);
gamma = rp_data(:,6);

for r = 1:dynnum 
    if r == 1
        dx(r,1) = 0;
        dy(r,1) = 0;
        dz(r,1) = 0;
        dalpha(r,1) = 0;
        dbeta(r,1) = 0;
        dgamma(r,1) = 0;
    else
        dx(r,1) = abs(rp_data(r,1) - rp_data(r-1,1)); %absolute relative displacement
        dy(r,1) = abs(rp_data(r,2) - rp_data(r-1,2));
        dz(r,1) = abs(rp_data(r,3) - rp_data(r-1,3));
        dalpha(r,1) = abs(alpha(r,1) - alpha(r-1,1));
        dbeta(r,1) = abs(beta(r,1) - beta(r-1,1));
        dgamma(r,1) = abs(gamma(r,1) - gamma(r-1,1));
    end
end


FD = (dx+dy+dz+dalpha+dbeta+dgamma); % FD (Power)

FD_index = find(FD>thresh_FD);

flag_index = FD_index;
dummy_index=ones(size(rp_data,1),1);

dummy_index(flag_index)=0;
dummy_index(flag_index+1)=0;
dummy_index(flag_index-1)=0;

flag(1,1:2) = [sum(FD>thresh_FD) length(find(dummy_index==0))];

%% plot FD and DVARS
p = figure('visible', 'off');
plot(FD, 'r');
xlabel('Dynnum')
ylabel('FD (mm)');
ylim([0 2]);
hold on
plot([0 dynnum], [thresh_FD thresh_FD], 'r-.');
xlim([0 dynnum])
hold off

print(p, '-djpeg', fullfile(subjdir, 'FD.jpg'));

csvwrite(fullfile(subjdir, 'flag.csv'), flag);  %% 1 = number of DVAR>5, 2 = FD > 0.5
csvwrite(fullfile(subjdir, 'censor.csv'), dummy_index);  %% 1 = number of DVAR>5, 2 = FD > 0.5
csvwrite(fullfile(subjdir, 'flag_index.csv'), flag_index);

end
