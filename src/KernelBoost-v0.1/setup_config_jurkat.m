function [params] = setup_config_jurkat(dataset_name)

%% Dataset parameters
params.dataset_name = dataset_name;

%% Boosting parameters
params.loss_type = 'squared';
params.wl_no = 3;
params.tree_depth = 3;
params.pool_filters = false;
params.sp_pooling = true;
params.sp_p_image_only = true;
params.use_zoned_sp_pool = true;
params.use_sp_diff_ops = true;
params.use_ilastik = false;
params.use_clusters = false;
params.cluster_centers_no = [12;5];
params.use_downscaling = true;
params.use_colors_luv = false;
params.shrinkage_factor = 0.1;
params.zoned_pool_radius = [-2;0;4;6];

%% Sampling parameters
params.sample_size = [51;51];
params.border_size = 2*(params.sample_size(1)-1);
params.border_skip_size = 4*params.sample_size(1);

params.pos_samples_no = 5000;
params.neg_samples_no = 5000;
params.T1_size = round((params.pos_samples_no+params.neg_samples_no)/3);
params.T2_size = round((params.pos_samples_no+params.neg_samples_no)/3);
params.pos_to_sample_no = params.T1_size;
params.neg_to_sample_no = params.T1_size;
params.rand_samples_no = 400;

params.use_qws = false;
params.use_2nd_deriv = false;
params.use_uniform_random_sampling = true;

%% Channels
params.ch_list = {'imgs'};
if (params.use_ilastik)
    params.ch_list{end+1} = 'ilastik';
end
params.ch_no = length(params.ch_list);

%% Image channel
params.imgs.learner_type = 'conv_reg';
params.imgs.min_f_size = 3;
params.imgs.max_f_size = 17;
params.imgs.filters_no = 200;
params.imgs.conv_reg_values = [500;1000;2000;3000];
params.imgs.perp_filters_no = 0;
params.imgs.zero_mean_flag = false;

%% ilastik features
params.ilastik.learner_type = 'conv_reg';
params.ilastik.min_f_size = 3;
params.ilastik.max_f_size = 17;
params.ilastik.filters_no = 200;
params.ilastik.conv_reg_values = [500;1000;2000;3000];
params.ilastik.perp_filters_no = 0;
params.ilastik.zero_mean_flag = false;

%% Filter pooling
params.pooling.filters_no = 100;
params.pooling.type = 'max';
params.pooling.nonlinearity = 'posneg';
params.pooling.pool_size = 5;
params.pooling.gauss_sigma = 3;

if (strcmp(params.pooling.type,'gauss'))
    h = fspecial('gaussian',params.pooling.pool_size,params.pooling.gauss_sigma);
    params.pooling.gauss_weights = h/norm(h);
end

%% Performance plotting parameters
params.PR_thresholds_no = 200;
params.cham_threshold = 0;
params.validate_set_flag = false;

%% Operations
params.ops = {'train','test'};
params.ops_no = length(params.ops);

%% Define codename
params.codename = input('    Please enter a valid codename for the simulation: ','s');
while (~all(isstrprop(params.codename,'alphanum')|params.codename=='_')||isempty(params.codename))
    fprintf('      INVALID INPUT! Only alphanumeric strings are accepted!\n');
    params.codename = input('      Please enter again a valid codename for the simulation: ','s');
end

pool_str = '';
if (params.pool_filters)
    pool_str = '__POOL__';
end
sppool_str = '';
if (params.sp_pooling)
    if (params.sp_p_image_only)
        sppool_str = '__SPPIX_IMG__';
    else
        sppool_str = '__SPPIX_ALL__';
    end
end
zoned_sppool_str = '';
if (params.use_zoned_sp_pool)
    zoned_sppool_str = ['__ZONEDSP_',sprintf('%d_',params.zoned_pool_radius),'_'];
end
ilastik_str = '';
if (params.use_ilastik)
    ilastik_str = '__ILASTIK__';
end
clusters_str = '';
if (params.use_clusters)
    clusters_str = sprintf('__CLUSTERS_%d__',params.cluster_centers_no);
end
dwn_str = '';
if (params.use_downscaling)
    dwn_str = '__DWN__';
end
params.codename = [params.codename,sprintf('____LOSS_%s__WL_%d__PS_%d_NS_%d__TreeD_%d__',params.loss_type,params.wl_no,params.pos_samples_no,params.neg_samples_no,params.tree_depth),pool_str,sppool_str,zoned_sppool_str,ilastik_str,clusters_str,dwn_str,'__Z'];

end
