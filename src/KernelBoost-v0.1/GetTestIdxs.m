function GetTestIdxs

rng(2014);

DATASET_NAME = 'DRIVE';

params = setup_config_L(DATASET_NAME);

params.codename = date;

params = setup_lists(params);
params = setup_directories_L(params);
save(fullfile(params.results_dir,'params.mat'),'params','-v7.3');

%% Data loading
fprintf('  Loading data...\n');
[params,data] = load_data(params);

%% Sample acquisition
fprintf('  Getting good sampling positions...\n');
[good_sampling_idx,tot_idx] = get_sampling_pos(params,data);
params.pos_samples_no = min(params.pos_samples_no,tot_idx.pos);
params.neg_samples_no = min(params.neg_samples_no,tot_idx.neg);

[samples_idx] = get_samples_idx(params,good_sampling_idx,data.train.gts.X,tot_idx);