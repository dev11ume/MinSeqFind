% FUNCTION TO SCORE SEQUENCES FROM PWMs
function scores=get_score_from_pwm_sequences_v1(seq_colp5n,pwms,scores_fac)
[num_seqs,seq_l]=size(seq_colp5n);
num_pwms=length(pwms);
if num_pwms==0
    scores=zeros(num_seqs,1);
else
    scores_p=cell(1,num_pwms);
    for nm=1:num_pwms
        pwm=pwms{nm};
        pwm_k=kron(pwm',ones(num_seqs,1));
        pwm_l=length(pwm)/4;
        loop_len=seq_l/4-pwm_l+1;
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
        scores_p{nm}=scores_p{nm}.*scores_fac(nm);
    end
    scores=max(cell2mat(scores_p),[],2);
end
end
