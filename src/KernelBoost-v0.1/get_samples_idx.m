function [samples_idx] = get_samples_idx(params,good_sampling_idx,imgs,tot_idx)

% Extract the sampling positions and their labels
%
%  authors: Carlos Becker, Roberto Rigamonti, CVLab EPFL
%  e-mail: name <dot> surname <at> epfl <dot> ch
%  web: http://cvlab.epfl.ch/
%  date: February 2014

fprintf('+++++++++++++++++++ SAMPLES +++++++++++++++++++\n');
fprintf('There are %d positives and %d negatives\n',tot_idx.pos,tot_idx.neg);
fprintf('===============================================\n');

% Structure of samples_idx:
% 1st index = image number
% 2nd index = row
% 3rd index = column
% 4th index = label

if (tot_idx.pos-params.pos_samples_no<0)
    error('Not enough positive training samples!');
end
if (tot_idx.neg-params.neg_samples_no<0)
    error('Not enough negative training samples!');
end

imgs_no = length(imgs);

tmp_pos = zeros(tot_idx.pos,3);
tmp_neg = zeros(tot_idx.neg,3);
pos_collected_no = 0;
neg_collected_no = 0;

for i_img = 1:imgs_no
    ns_pos = length(good_sampling_idx{i_img}.pos);
    ns_neg = length(good_sampling_idx{i_img}.neg);

    tmp_pos(pos_collected_no+1:pos_collected_no+ns_pos,1) = i_img*ones(ns_pos,1);
    tmp_neg(neg_collected_no+1:neg_collected_no+ns_neg,1) = i_img*ones(ns_neg,1);
    [tmp_pos(pos_collected_no+1:pos_collected_no+ns_pos,2),tmp_pos(pos_collected_no+1:pos_collected_no+ns_pos,3)] = ind2sub(size(imgs{i_img}),good_sampling_idx{i_img}.pos);
    [tmp_neg(neg_collected_no+1:neg_collected_no+ns_neg,2),tmp_neg(neg_collected_no+1:neg_collected_no+ns_neg,3)] = ind2sub(size(imgs{i_img}),good_sampling_idx{i_img}.neg);            

    pos_collected_no = pos_collected_no+ns_pos;
    neg_collected_no = neg_collected_no+ns_neg;
end

rnd_idx_pos = randperm(size(tmp_pos,1));
rnd_idx_neg = randperm(size(tmp_neg,1));

sampled_idx.pos = tmp_pos(rnd_idx_pos(1:min(length(rnd_idx_pos),params.pos_samples_no)),:);
sampled_idx.neg = tmp_neg(rnd_idx_neg(1:min(length(rnd_idx_neg),params.neg_samples_no)),:);

labels = [ones(size(sampled_idx.pos,1),1);-ones(size(sampled_idx.neg,1),1)];
samples_idx = [cat(1,sampled_idx.pos,sampled_idx.neg),labels];

end
