% MAIN FUNCTION TO RUN MinSeqFind ALGORITHM  

function minseqfind_1(File1,LibFile,Llib,Rlib,param)

nmer_cut=200; % THRESHOLD TO REMOVE SEQUENCES WHICH APPEAR IN EXCESS WITHOUT THEIR REVERSE COMPLEMENT

% MINSEQ RELATED PARAMETERS
min_half_nmer=param.min_half_nmer;
max_half_nmer=param.max_half_nmer;
min_nmer=param.min_nmer;
max_nmer=param.max_nmer;
max_gap=param.max_gap;
pow_factor=param.pow_factor;
max_cut_arr=param.max_cut_arr;
flip3to5=param.flip3to5;
seeds_to_cont=param.seeds_to_cont;

% LIBRARY RELATED PARAMETERS
lib_nmer=param.lib_nmer;
gap_mod=param.gap_mod;
interpol=param.interpol;
lib_cut=param.lib_cut;

% PWM RELATED PARAMETERS
max_cut_pwm=param.max_cut_pwm;
min_ct=param.min_ct;
seeds_to_cont2=param.seeds_to_cont2;
num_motifs_max=param.num_motifs_max;
num_Ns=param.num_Ns;

% THRESHOLDING TO REMOVE SOME MINSEQS
enrich_cut=param.enrich_cut;
enrich_perc_cut=param.enrich_perc_cut;
weighted_enrich_perc_cut=param.weighted_enrich_perc_cut;

max_cut=min(max_cut_arr);
model_order_max=lib_nmer;
if max_half_nmer~=-1
    model_order_max=min(lib_nmer,2*max_half_nmer-1);
end
model_order_min=min(model_order_max,2*min_half_nmer-1);

File2x=LibFile(1:end-4);
if interpol==1
    File2x=[File2x '_i1_' num2str(lib_cut)];
end
File2x=[File2x '_4.'];

% LIBRARY MODEL CREATION AND NAMING LIBRARY FILE
if gap_mod
    fname_txt=[File2x num2str(model_order_min) '.order.model.1_1_1.txt'];
else
    fname_txt=[File2x num2str(model_order_min) '.order.model.1.txt'];
end

if exist(fname_txt,'file') % ASSUMING IF 1 FILE EXISTS ALL EXISTS AND SAME PARAMETERS 
    [l,l_seqs_lib,l_seqs_lib2,cont_seqs1l,lib_model_dat1]=...
        lib_model_load_var_v2(LibFile,model_order_min,model_order_max,lib_cut,...
        interpol,gap_mod,flip3to5);
else
    [l,l_seqs_lib,l_seqs_lib2,cont_seqs1l,lib_model_dat1]=...
        lib_model_v2_2(LibFile,model_order_min,model_order_max,lib_cut,...
        interpol,gap_mod,flip3to5);
end

% INITIAL DATA PROCESSING
[l_seqs_raw,l_seqs2,cont_seqs1,cont_seqs2,seq_colp5n,count,seqs_uq,to_contf]=...
    init_data_process_v1(File1,Llib,Rlib,nmer_cut,LibFile);

max_gap=min(max_gap,l-2*min_half_nmer);
% GET COUNTS FOR SAMPLE AND LIBRARY
if max_half_nmer==-1
    [ctsfr,ctsfl,bgsf,inx,all_len]=...
        get_counts_bgs_ng_v1(seq_colp5n,l,lib_nmer,count,File2x,...
        min_nmer,max_nmer,max_cut,Llib,Rlib,flip3to5,lib_model_dat1);
    all_len_nox=all_len;
else
    [ctsfr,ctsfl,bgsf,inx,all_len,all_len_nox]=...
        get_counts_bgs_gp_v1(seq_colp5n,l,lib_nmer,count,File2x,...
        min_half_nmer,max_half_nmer,max_gap,max_cut,Llib,Rlib,gap_mod,flip3to5,lib_model_dat1);
end

l_max_cut_arr         =length(max_cut_arr);
seeds_to_cont_arr     =cell(l_max_cut_arr,1);
seeds_to_cont2_arr    =cell(l_max_cut_arr,1);
seed_mot_acgt_arr     =cell(l_max_cut_arr,1);
pwm_all_arr           =cell(l_max_cut_arr,1);
int_pwm_seed_arr      =cell(l_max_cut_arr,1);
inx_pwm_seed_arr      =cell(l_max_cut_arr,1);
int_nu2x_c_arr2       =cell(l_max_cut_arr,1);
seqsx_inx_c_arr       =cell(l_max_cut_arr,1);
scores_sum20          =cell(l_max_cut_arr,1);
int_seed_pwm_score    =cell(l_max_cut_arr,1);
int_seed_calculated_arr    =cell(l_max_cut_arr,1);
count_seed_arr    =cell(l_max_cut_arr,1);
for i=1:l_max_cut_arr
    ix_cont=find(ctsfr>=max_cut_arr(i)); % INX OF CTS ABOVE max_cut_arr(i) IN RANDOM REGION
    if length(ix_cont)>=10 % ATLEAST 10 SEQUENCES ABOVE max_cut_arr(i)
        fprintf('\nMAX CUT ARR i=%d\n',i);
        ctsf=ctsfr(ix_cont)+ctsfl(ix_cont);
        inx_rev=get_rev_inx_v2(inx(ix_cont),max_half_nmer,l);
        inx_min=min([inx(ix_cont) inx_rev],[],2);
        [inx_u,inxx]=unique(inx_min);
        int_n2u=ctsf(inxx)./bgsf(ix_cont(inxx));
        int_n2u=int_n2u*l_seqs_lib2/l_seqs2;
        
        max_enrich=max(int_n2u);
        enrich_cutx=max(enrich_cut,max_enrich*enrich_perc_cut/100);
        inxx2=find(int_n2u>=enrich_cutx);
        int_n2u2=int_n2u(inxx2);
        int_nu2x=int_n2u2./(pow_factor.^(all_len_nox(ix_cont(inxx(inxx2)))));
        
        [~,int_nu2xs_ix]=sort(int_nu2x,'descend');
        max_weighted_enrich=int_nu2x(int_nu2xs_ix(1));
        inxx3=find(int_n2u2>=(max_weighted_enrich*weighted_enrich_perc_cut/100));
        
        seeds_to_cont_arr{i}=min(seeds_to_cont,length(inxx3));
        int_nu2x_c_arr2{i}=int_nu2x(inxx3(int_nu2xs_ix(1:seeds_to_cont_arr{i})));
        all_len_c=all_len(ix_cont(inxx(inxx2(inxx3(int_nu2xs_ix(1:seeds_to_cont_arr{i}))))));
        all_len_nox_c=all_len_nox(ix_cont(inxx(inxx2(inxx3(int_nu2xs_ix(1:seeds_to_cont_arr{i}))))));
        seqsx_inx_c_arr{i}=inx_u(inxx2(inxx3(int_nu2xs_ix(1:seeds_to_cont_arr{i}))));
        fprintf('\nGETTING PWMs\n');
        if max_cut_arr(i)>=max_cut_pwm
            seeds_to_cont2_arr{i}=min(seeds_to_cont2,seeds_to_cont_arr{i});
            if seeds_to_cont2_arr{i}>num_motifs_max
                [seed_mot_acgt_arr{i},pwm_all_arr{i},int_pwm_seed_arr{i},inx_pwm_seed_arr{i},~,...
                    int_seed_pwm_score{i},int_seed_calculated_arr{i},count_seed_arr{i}]=get_pwm_from_minseq_v1(...
                    seeds_to_cont2_arr{i},int_nu2x_c_arr2{i},all_len_c,all_len_nox_c,num_Ns,count,...
                    seqsx_inx_c_arr{i},seqs_uq,l,Llib,Rlib,num_motifs_max,max_half_nmer,File2x,...
                    lib_nmer,min_ct,pow_factor,flip3to5,l_seqs_lib2,l_seqs2,lib_model_dat1);
            end
        end
    end
end

clear bgsf ctsfl ctsfr inx all_len_c all_len_nox_c inx_min inx_rev inx_u inxx ix_cont...
    seq_colp5n seqs seqs_inx seq_m int_nu2xs_ix int_nu2x...
    int_n2u all_len all_len_nox ctsf seqs_inx_c Library_files Llibs...
    Rlibs TFs_file TFs_name TFs_round int_n2u2 seqs_uq count seq_col_tot bg...
    all_len_c all_len_nox_c seq_colp5n_tot;

try
  save([File1(1:end-4) '-' LibFile(1:end-4) '-res1x.mat'],'-mat7-binary');
catch
  save([File1(1:end-4) '-' LibFile(1:end-4) '-res1x.mat'],'-binary');
end
movefile([File1(1:end-4) '-' LibFile(1:end-4) '-res1x.mat'],[File1 '-' LibFile '-res1f.mat'],'f');

end
