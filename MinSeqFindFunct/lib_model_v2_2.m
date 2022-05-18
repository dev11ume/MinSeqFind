% FUNCTION TO CREATE MARKOV MODEL OF LIBRARY POSITION SPECIFIC 
function [l,l_seqs1,l_seqs2,cont_seqs1,lib_model_dat1]=...
    lib_model_v2_2(File1,model_order_min,model_order_max,...
    lib_cut,interpol,gap_mod,flip3to5)
fprintf('Library Processing...\n');
warning off;
C=textscan_mod_v1(File1,'%s','\t');
seqs=C{1};

seqs=upper(seqs);
l=length(seqs{1});
l_seqs1=length(seqs);
cont_seqs1=zeros(l_seqs1,1);
% REMOVE SEQUENCE THAT CONTAINS 'Ns'
for i=1:l_seqs1
    if isempty(find(seqs{i}=='N', 1)); cont_seqs1(i)=1; end
    if l==20 % REMOVE SOME KNOWN CONTAMINATION SEQUENCES
        if sum(seqs{i}=='CGAATGATGGATTGCAACCG')>=18 ||...
                sum(seqs{i}=='TACCGATTACGTAATTTCGA')>=18 ||...
                sum(seqs{i}=='CGCGAATGACGTCAATCGGA')>=18 ||...
                sum(seqs{i}=='GCTAGATTGCGCAATCCAGT')>=18
            cont_seqs1(i)=0;
        end
    end
end
seqs=seqs(cont_seqs1==1);
l_seqs2=length(seqs);
if flip3to5
    seqs=fliplr(cell2mat(seqs));
end
[seq_colp5n,~]=get_seq_col_mat_only_v1(seqs,'','');

for model_order=model_order_min:model_order_max
    for i=1:l-model_order
        fprintf('Position %d\n',i);
        [M,fname]=get_M_matrix_lib_v2(seq_colp5n(:,i:i+model_order),interpol,model_order,lib_cut,File1,i,0,0);
        lib_model_dat1.(fname)=M;
    end
end

if gap_mod==1
    % GAPPED MODEL
    for model_order=model_order_min:model_order_max
        for i=1:l-model_order
            for i_gap_pos=1:model_order
                for i_gap_len=1:l-model_order-i
                    fprintf('Position %d\t%d\t%d\t\n',i,i_gap_pos,i_gap_len);
                    [M,fname]=get_M_matrix_lib_v2(seq_colp5n(:,[i:i+i_gap_pos-1 i+i_gap_pos+i_gap_len:i+i_gap_len+model_order]),...
                        interpol,model_order,lib_cut,File1,i,i_gap_pos,i_gap_len);
                    lib_model_dat1.(fname)=M;
                end
            end
        end
    end
    
    % ORDER LESS THAN MODEL_ORDER
    for model_orderi=0:model_order_max-1
        fprintf('model_orderi %d\n',model_orderi);
        for i=1:l-model_orderi
            [M,fname]=get_M_matrix_lib_v2(seq_colp5n(:,i:i+model_orderi),interpol,model_orderi,lib_cut,File1,i,0,0);
            lib_model_dat1.(fname)=M;
        end
    end
    
    % GAPPED MODEL
    for model_orderi=0:model_order_max-1
        for i_gap_pos=1:model_orderi
            i=1;%LEFT SIDE
            for i_gap_len=1:l-model_orderi-i
                [M,fname]=get_M_matrix_lib_v2(seq_colp5n(:,[i:i+i_gap_pos-1 i+i_gap_pos+i_gap_len:i+i_gap_len+model_orderi]),...
                    interpol,model_orderi,lib_cut,File1,i,i_gap_pos,i_gap_len);
                lib_model_dat1.(fname)=M;
            end
            ix=1;%RIGHT SIDE
            for i_gap_len=1:l-model_orderi-ix
                i=l-model_orderi-i_gap_len;
                [M,fname]=get_M_matrix_lib_v2(seq_colp5n(:,[i:i+i_gap_pos-1 i+i_gap_pos+i_gap_len:i+i_gap_len+model_orderi]),...
                    interpol,model_orderi,lib_cut,File1,i,i_gap_pos,i_gap_len);
                lib_model_dat1.(fname)=M;
            end
        end
    end
    
else
    
    for model_orderi=0:model_order_max-1
        fprintf('Shorter model: %d\n',model_orderi);
        i=1;%LEFT SIDE
        [M,fname]=get_M_matrix_lib_v2(seq_colp5n(:,i:i+model_orderi),interpol,model_orderi,lib_cut,File1,i,0,0);
        lib_model_dat1.(fname)=M;
        i=l-model_orderi;%RIGHT SIDE
        [M,fname]=get_M_matrix_lib_v2(seq_colp5n(:,i:i+model_orderi),interpol,model_orderi,lib_cut,File1,i,0,0);
        lib_model_dat1.(fname)=M;
    end
    
end
warning on;
fprintf('Done Processing Library\n');
end
