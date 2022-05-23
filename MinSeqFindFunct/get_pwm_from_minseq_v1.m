% FUNCTION TO GET PWMs FROM MinSeqs 

function [seed_mot_acgt,pwm_all,int_pwm_seed,inx_pwm_seed,scores_pwm,...
    int_seed_pwm_score,int_seed_calculated,count_seed]=get_pwm_from_minseq_v1(...
    seeds_to_cont2,int_nu2x_c,all_len_c,all_len_nox_c,num_Ns,count,...
    seqsx_inx_c,seqs_uq,l,Llib,Rlib,num_motifs_max,max_half_nmer,...
    File2x,lib_nmer,min_ct,pow_factor,flip3to5,l_seqs_lib2,l_seqs2,lib_model_dat1,strnums_plus)
lLlib=length(Llib);
lRlib=length(Rlib);

int_nu2x_c2=int_nu2x_c(1:seeds_to_cont2);
all_len_c2=all_len_c(1:seeds_to_cont2);
all_len_nox_c2=all_len_nox_c(1:seeds_to_cont2);
seqsx_inx_c2=seqsx_inx_c(1:seeds_to_cont2);

[seq_col,~,~]=get_seq_col_v2(seqs_uq,Llib,Rlib);


aa=[1 0 2 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0];

max_iter_seed=1;
seed_mot_acgt=cell(num_motifs_max,max_iter_seed);
seed_col_all=cell(num_motifs_max,max_iter_seed);
pwm_all=cell(num_motifs_max,max_iter_seed);

int_pwm_seed=zeros(num_motifs_max,1);
int_seed_calculated=zeros(num_motifs_max,1);
count_seed=zeros(num_motifs_max,1);
inx_pwm_seed=zeros(num_motifs_max,1);
scores_pwm=zeros(seeds_to_cont2,num_motifs_max);
int_seed_pwm_score=zeros(num_motifs_max,1);
nm=0;
while (nm<num_motifs_max)
    nm=nm+1; itr=1;
    
    if nm>1
        int_nu2x_ct=int_nu2x_c2-max(scores_pwm(:,1:nm-1),[],2);
    else
        int_nu2x_ct=int_nu2x_c2-max(scores_pwm(:,1:nm),[],2); % IF FIRST PWM COME BACK WITHOUT ENOUGH COUNTS  
    end
    fprintf('\nPWM number:%d\n',nm);
    fprintf('Int Min=%1.3f Max=%1.3f\n',min(int_nu2x_ct),max(int_nu2x_ct));
    if max(int_nu2x_ct)<=0
        break;
    end
    [~,int_sx]=sort(int_nu2x_ct,'descend');
    inx_pwm_seed(nm)=int_sx(1);
    int_pwm_seed(nm)=int_nu2x_c2(int_sx(1));
    
    sub_len=0;
    cont_sub_len=1;
    while cont_sub_len
        seed_mot_acgt(nm,itr)=get_index_from_seqs_withN_v1(seqsx_inx_c2(inx_pwm_seed(nm)),...
            all_len_c2(inx_pwm_seed(nm)),max_half_nmer);
        num_Nsx=min(num_Ns,floor((l-length(seed_mot_acgt{nm,itr}))/2))-sub_len;
        num_Nsx=max(num_Nsx,0);
        seed_mot_acgt{nm,itr}=[repmat('N',1,num_Nsx) seed_mot_acgt{nm,itr} repmat('N',1,num_Nsx)];
        fprintf('%s\n',seed_mot_acgt{nm,itr});
        
        seed_col_all{nm,itr}=get_mot_to_seed_col(seed_mot_acgt{nm,itr});
        
        pwm_length=length(seed_mot_acgt{nm,itr});
        
        ctf=zeros(4*pwm_length,1);
        ctr=zeros(4*pwm_length,1);
        ctf_seed=0;
        ctr_seed=0;
        seed_mot_colp5nf_t=zeros(4*pwm_length,pwm_length);
        seed_mot_colp5nr_t=zeros(4*pwm_length,pwm_length);
        
        % PWM calculation 
        for i=lLlib+1:lLlib+l-pwm_length+1
            xx=double(seq_col(:,i*4-3:(i+pwm_length-1)*4));
            ixf=find(xx*seed_col_all{nm,itr}'>=pwm_length-1);
            ixr=find(xx*seed_col_all{nm,itr}(end:-1:1)'>=pwm_length-1);
            xxf=xx(ixf,:);
            xxr=xx(ixr,:);
            countf=count(ixf);
            countr=count(ixr);
            
            seed_col_t=seed_col_all{nm,itr};
            seed_col_tr=seed_col_t(end:-1:1);
            if i==lLlib+1
                seed_mot_acgt_t{1}=get_seed_col_to_mot(seed_col_t);
                seed_mot_colp5nf_seed=aa(cell2mat(seed_mot_acgt_t)-64);
                seed_mot_acgt_t{1}=get_seed_col_to_mot(seed_col_tr);
                seed_mot_colp5nr_seed=aa(cell2mat(seed_mot_acgt_t)-64);
            end
            if ~isempty(ixf)
                ctf_seed=ctf_seed+...
                    sum(countf(xxf*seed_col_t'==pwm_length));
            end
            if ~isempty(ixr)
                ctr_seed=ctr_seed+...
                    sum(countr(xxr*seed_col_tr'==pwm_length));
            end
            
            for pos=1:pwm_length
                for nucl=1:4
                    seed_col_t=seed_col_all{nm,itr};
                    seed_col_t(4*pos-3:4*pos)=0;
                    seed_col_t(4*pos-4+nucl)=1;
                    seed_col_tr=seed_col_t(end:-1:1);
                    if i==lLlib+1
                        seed_mot_acgt_t{1}=get_seed_col_to_mot(seed_col_t);
                        seed_mot_colp5nf_t(4*pos-4+nucl,:)=aa(cell2mat(seed_mot_acgt_t)-64);
                        seed_mot_acgt_t{1}=get_seed_col_to_mot(seed_col_tr);
                        seed_mot_colp5nr_t(4*pos-4+nucl,:)=aa(cell2mat(seed_mot_acgt_t)-64);
                    end
                    if ~isempty(ixf)
                        ctf(4*pos-4+nucl)=ctf(4*pos-4+nucl)+...
                            sum(countf(xxf*seed_col_t'==pwm_length));
                    end
                    if ~isempty(ixr)
                        ctr(4*pos-4+nucl)=ctr(4*pos-4+nucl)+...
                            sum(countr(xxr*seed_col_tr'==pwm_length));
                    end
                end
            end
        end
        
        % GET BG FOR PWM
        if flip3to5
            bgft=get_lib_prediction_v1(fliplr(seed_mot_colp5nf_t),pwm_length,...
                File2x,lib_nmer,l,fliplr(Rlib),fliplr(Llib),0,0,0,lib_model_dat1,strnums_plus);
            bgrt=get_lib_prediction_v1(fliplr(seed_mot_colp5nr_t),pwm_length,...
                File2x,lib_nmer,l,fliplr(Rlib),fliplr(Llib),0,0,0,lib_model_dat1,strnums_plus);
        else
            bgft=get_lib_prediction_v1(seed_mot_colp5nf_t,pwm_length,File2x,lib_nmer,l,'','',0,0,0,lib_model_dat1,strnums_plus);
            bgrt=get_lib_prediction_v1(seed_mot_colp5nr_t,pwm_length,File2x,lib_nmer,l,'','',0,0,0,lib_model_dat1,strnums_plus);
            bgft_seed=get_lib_prediction_v1(seed_mot_colp5nf_seed,pwm_length,File2x,lib_nmer,l,'','',0,0,0,lib_model_dat1,strnums_plus);
            bgrt_seed=get_lib_prediction_v1(seed_mot_colp5nr_seed,pwm_length,File2x,lib_nmer,l,'','',0,0,0,lib_model_dat1,strnums_plus);
        end
        
        ctt=ctf+ctr;
        bgt=bgft+bgrt;
        
        ctt_seed=ctf_seed+ctr_seed;
        bgt_seed=bgft_seed+bgrt_seed;
        
        pwm_all{nm,itr}=ones(4*pwm_length,1);
        cont_position=zeros(pwm_length,1);
        for pos=1:pwm_length
            if sum(ctt(4*pos-3:4*pos))>=min_ct &&  isempty(find(bgt(4*pos-3:4*pos)==0,1))
                cont_position(pos)=1;
                pwm_all{nm,itr}(4*pos-3:4*pos)=ctt(4*pos-3:4*pos)./bgt(4*pos-3:4*pos);
            end
        end
        int_seed_calculated(nm)=ctt_seed/bgt_seed*l_seqs_lib2/l_seqs2/(pow_factor^(all_len_nox_c2(inx_pwm_seed(nm))));
        count_seed(nm)=ctt_seed;
        ix_cont=find(cont_position);
        
        if isempty(ix_cont) && num_Nsx>0
            sub_len=sub_len+1;
        elseif isempty(ix_cont) && num_Nsx==0
            cont_sub_len=0;
        elseif (((ix_cont(end)-ix_cont(1))>(length(ix_cont)-1))||(ix_cont(end)<pwm_length-num_Nsx)||(ix_cont(1)>num_Nsx+1)) && sub_len<num_Ns && num_Nsx>0
            sub_len=sub_len+1;
        else
            cont_sub_len=0;
        end
    end
    if isempty(ix_cont)
        scores_pwm(inx_pwm_seed(nm),1)=int_nu2x_c2(inx_pwm_seed(nm));
        nm=nm-1;
    elseif (((ix_cont(end)-ix_cont(1))>(length(ix_cont)-1))||(ix_cont(end)<pwm_length-num_Nsx)||(ix_cont(1)>num_Nsx+1))
        scores_pwm(inx_pwm_seed(nm),1)=int_nu2x_c2(inx_pwm_seed(nm));
        nm=nm-1;
    else
        ix_start=min(ix_cont(1),num_Nsx+1);
        ix_stop=max(ix_cont(end),pwm_length-num_Nsx);
        pwm_all{nm,itr}=pwm_all{nm,itr}(4*ix_start-3:4*ix_stop);
        seed_mot_acgt{nm,itr}=seed_mot_acgt{nm,itr}(ix_start:ix_stop);
        
        % GET SCORES FROM PWM 
        scores_pwm(:,nm)=get_score_from_pwm_sequences_v3(seqsx_inx_c2,l,pwm_all(nm,itr),...
            int_seed_calculated(nm),inx_pwm_seed(nm),all_len_nox_c2,max_half_nmer,pow_factor);
        
        [seq_colx,~]=get_seq_col_v2(seed_mot_acgt(nm),[],[]);
        seed_pwm_score=get_score_from_pwm_sequences_v1(seq_colx,pwm_all(nm),1);
        int_seed_pwm_score(nm)=int_seed_calculated(nm)/seed_pwm_score* ...
            (pow_factor^(all_len_nox_c2(inx_pwm_seed(nm))));
        scores_pwm(inx_pwm_seed(nm),1)=int_nu2x_c2(inx_pwm_seed(nm));
    end
end
