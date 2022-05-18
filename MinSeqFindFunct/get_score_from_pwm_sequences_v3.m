% FUNCTION - SCORE SEQUENCES FROM PWM
function scores=get_score_from_pwm_sequences_v3(seqsx_inx_c,seq_l,pwms,int_seed,...
    inx_pwm_seed,all_len_nox,max_half_nmer,pow_factor)
seqsxu_c=get_index_from_seqs_withN_v1(seqsx_inx_c,seq_l,max_half_nmer);
[seq_colp5n,~,num_seqs]=get_seq_col_v2(seqsxu_c,repmat('N',1,seq_l),repmat('N',1,seq_l));

num_pwms=length(pwms);
if num_pwms==0
    scores=zeros(num_seqs,1);
else
    scores_p=cell(1,num_pwms);
    for nm=1:num_pwms
        pwm=pwms{nm};
        pwm_k=kron(pwm',ones(num_seqs,1));
        pwm_l=length(pwm)/4;
        loop_len=3*seq_l-pwm_l+1;
        scoresf=zeros(num_seqs,loop_len);
        scoresr=zeros(num_seqs,loop_len);
        % TIME CONSUMING LOOP
        for i=1:loop_len
            cols=4*i-3:4:4*(i+pwm_l-1);
            cts_pos=prod(seq_colp5n(:,cols)+seq_colp5n(:,cols+1)+seq_colp5n(:,cols+2)+seq_colp5n(:,cols+3),2);
            seq_colp5nk=seq_colp5n(:,4*i-3:4*(i+pwm_l-1)).*pwm_k;
            cols2=1:4:4*pwm_l;
            scoresf(:,i)=prod(seq_colp5nk(:,cols2)+seq_colp5nk(:,cols2+1)+...
                seq_colp5nk(:,cols2+2)+seq_colp5nk(:,cols2+3),2)./cts_pos;
        end
        pwm_k=kron(pwm(end:-1:1)',ones(num_seqs,1));
        % TIME CONSUMING LOOP
        for i=1:loop_len
            cols=4*i-3:4:4*(i+pwm_l-1);
            cts_pos=prod(seq_colp5n(:,cols)+seq_colp5n(:,cols+1)+seq_colp5n(:,cols+2)+seq_colp5n(:,cols+3),2);
            seq_colp5nk=seq_colp5n(:,4*i-3:4*(i+pwm_l-1)).*pwm_k;
            cols2=1:4:4*pwm_l;
            scoresr(:,i)=prod(seq_colp5nk(:,cols2)+seq_colp5nk(:,cols2+1)+...
                seq_colp5nk(:,cols2+2)+seq_colp5nk(:,cols2+3),2)./cts_pos;
        end
        scores_p{nm}=max([scoresf scoresr],[],2);
        scores_p{nm}=scores_p{nm}./(pow_factor.^all_len_nox);
        scores_p{nm}=scores_p{nm}.*int_seed(nm)/scores_p{nm}(inx_pwm_seed(nm));
    end
    scores=max(cell2mat(scores_p),[],2);
end
end
