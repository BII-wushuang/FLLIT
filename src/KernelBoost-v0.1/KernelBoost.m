%  Source code for the paper:
%  "Supervised Feature Learning for Curvilinear Structure Segmentation",
%  C. Becker, R. Rigamonti, V. Lepetit, P. Fua, MICCAI 2013

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

clear;
rng(2014);

DATASET_NAME = 'DRIVE';

addpath(genpath('mex'));
addpath('/cvlabdata1/home/rigamont/libs/dip/common/dipimage');

%% Configuration
fprintf('== KernelBoost - v0.1 ==\n');
fprintf('  Setting up configuration and directories...\n');
params = setup_config(DATASET_NAME);
params = setup_lists(params);
params = setup_directories(params);
save(fullfile(params.results_dir,'params.mat'),'params','-v7.3');

%% Data loading
fprintf('  Loading data...\n');
[params,data] = load_data(params);

%% Sample acquisition
fprintf('  Getting good sampling positions...\n');
[good_sampling_idx,tot_idx] = get_sampling_pos(params,data);
params.pos_samples_no = min(params.pos_samples_no,tot_idx.pos);
params.neg_samples_no = min(params.neg_samples_no,tot_idx.neg);

fprintf('  Getting sample indices...\n');
smpls_fname = fullfile(params.results_dir,'samples.mat');
if (~exist(smpls_fname,'file'))
    [samples_idx] = get_samples_idx(params,good_sampling_idx,data.train.gts.X,tot_idx);
    save(smpls_fname,'samples_idx','-v7.3');
else
    load(smpls_fname);
end

% Recompute samples numbers according to the collected samples numbers
params.pos_samples_no = sum(samples_idx(:,end)==1);
params.neg_samples_no = sum(samples_idx(:,end)==-1);

% number to sample for filter search (T1) and tree learning (T2)
params.T1_size = round((params.pos_samples_no+params.neg_samples_no)/3);
params.T2_size = round((params.pos_samples_no+params.neg_samples_no)/3);


params.pos_to_sample_no = params.T1_size;
params.neg_to_sample_no = params.T1_size;

%% KB training
wl_fname = fullfile(params.results_dir,'weak_learners.mat');
if (~exist(wl_fname,'file'))
    weak_learners = train_boost_general(params,data,samples_idx);
    save(wl_fname,'weak_learners','-v7.3');
else
    load(wl_fname);
end

%% KB testing
fprintf('Final classifier stored on disk, now classifying test images...\n');
test_classifiers(params,data,weak_learners);
