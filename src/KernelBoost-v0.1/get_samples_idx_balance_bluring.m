function [samples_idx] = get_samples_idx_balance_bluring(params,good_sampling_idx,tot_idx)

% Extract the sampling positions and their label proportions

% Unlike the previous methods, here attempts to balance the bluring and
% unbluring region


fprintf('+++++++++++++++++++ SAMPLES +++++++++++++++++++\n');
fprintf('There are %d positives and %d negatives\n',tot_idx.pos,tot_idx.neg);
fprintf('===============================================\n');

% Structure of samples_idx:
% 1st index = image number
% 2nd index = index
% 3th index = label alpha

if (tot_idx.pos-params.pos_samples_no<0)
    error('Not enough positive training samples!');
end
if (tot_idx.neg-params.neg_samples_no<0)
    error('Not enough negative training samples!');
end

imgs_no = length(good_sampling_idx);

tmp_pos = zeros(tot_idx.pos,3);
tmp_neg = zeros(tot_idx.neg,3);
pos_collected_no = 0;
neg_collected_no = 0;

for i_img = 1:imgs_no
    ns_pos = length(good_sampling_idx{i_img}.pos);
    ns_neg = length(good_sampling_idx{i_img}.neg);

    tmp_pos(pos_collected_no+1:pos_collected_no+ns_pos,1) = i_img*ones(ns_pos,1);
    tmp_neg(neg_collected_no+1:neg_collected_no+ns_neg,1) = i_img*ones(ns_neg,1);
    
    
    tmp_pos(pos_collected_no+1:pos_collected_no+ns_pos,2) = good_sampling_idx{i_img}.pos;
    tmp_neg(neg_collected_no+1:neg_collected_no+ns_neg,2) = good_sampling_idx{i_img}.neg;            

    tmp_pos(pos_collected_no+1:pos_collected_no+ns_pos,3) = good_sampling_idx{i_img}.pos_val;    
    tmp_neg(neg_collected_no+1:neg_collected_no+ns_neg,3) = good_sampling_idx{i_img}.neg_val;
    
    pos_collected_no = pos_collected_no+ns_pos;
    neg_collected_no = neg_collected_no+ns_neg;
end

% b refers to the blur part

tmp_pos_b = tmp_pos(tmp_pos(:,3) < 0.99,:);

% b refers to the determined part

tmp_pos_d = tmp_pos(tmp_pos(:,3) > 0.99,:);

rnd_idx_pos_b = randperm(size(tmp_pos_b,1));

rnd_idx_pos_d = randperm(size(tmp_pos_d,1));

pos_samples_no = params.pos_samples_no;

pos_samples_no = round(pos_samples_no / 2);



sampled_idx.pos = tmp_pos_b(rnd_idx_pos_b(1:min(length(rnd_idx_pos_b),pos_samples_no)),:);

sampled_idx.pos = [sampled_idx.pos ;...
    tmp_pos_d(rnd_idx_pos_d(1:min(length(rnd_idx_pos_d),pos_samples_no)),:)];


tmp_neg_b = tmp_neg(tmp_neg(:,3) > 0.01,:);

tmp_neg_d = tmp_neg(tmp_neg(:,3) < 0.01,:);

rnd_idx_neg_b = randperm(size(tmp_neg_b,1));

rnd_idx_neg_d = randperm(size(tmp_neg_d,1));

neg_samples_no = params.neg_samples_no;

neg_samples_no = round(neg_samples_no / 2);


sampled_idx.neg = tmp_neg_b(rnd_idx_neg_b(1:min(length(rnd_idx_neg_b),neg_samples_no)),:);

sampled_idx.neg = [sampled_idx.neg ;...
tmp_neg_d(rnd_idx_neg_d(1:min(length(rnd_idx_neg_d),neg_samples_no)),:)];


% rnd_idx_neg = randperm(size(tmp_neg,1));

% sampled_idx.pos = tmp_pos(rnd_idx_pos(1:min(length(rnd_idx_pos),params.pos_samples_no)),:);
% sampled_idx.neg = tmp_neg(rnd_idx_neg(1:min(length(rnd_idx_neg),params.neg_samples_no)),:);

samples_idx = cat(1,sampled_idx.pos,sampled_idx.neg); 

% 
% labels = [ones(size(sampled_idx.pos,1),1);-ones(size(sampled_idx.neg,1),1)];
% samples_idx = [cat(1,sampled_idx.pos,sampled_idx.neg),labels];

end
