% FUNCTION TO LOAD MARKOV MODEL OF LIBRARY POSITION SPECIFIC 

function [l,l_seqs1,l_seqs2,cont_seqs1,lib_model_dat1]=...
    lib_model_load_var_v2(File1,model_order_min,model_order_max,...
    lib_cut,interpol,gap_mod,flip3to5)
fprintf('Library Loading ...\n');

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

File1=File1(1:end-4);

if interpol==1
    File1=[File1 '_i1_' num2str(lib_cut) ];
end

seqs=seqs(cont_seqs1==1);
l_seqs2=length(seqs);
if flip3to5
    seqs=fliplr(cell2mat(seqs));
end
for model_order=model_order_min:model_order_max
    for i=1:l-model_order
        fprintf('Position %d\n',i);
        fname=[File1 '_4_' num2str(model_order) '_order_model_' num2str(i)];
        M=dlmread([File1 '_4.' num2str(model_order) '.order.model.' num2str(i) '.txt'],...
            '\t');
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
                    fname=[File1 '_4_' num2str(model_order) ...
                        '_order_model_' num2str(i) '_' num2str(i_gap_pos) '_' ...
                        num2str(i_gap_len)];
                    M=dlmread([File1 '_4.' num2str(model_order) ...
                        '.order.model.' num2str(i) '_' num2str(i_gap_pos) '_' ...
                        num2str(i_gap_len) '.txt'],'\t');
                    lib_model_dat1.(fname)=M;
                end
            end
        end
    end
    
    % ORDER LESS THAN MODEL_ORDER
    for model_orderi=0:model_order_max-1
        fprintf('model_orderi %d\n',model_orderi);
        for i=1:l-model_orderi
            fname=[File1 '_4_' num2str(model_orderi) '_order_model_' num2str(i)];
            M=dlmread([File1 '_4.' num2str(model_orderi) '.order.model.' num2str(i) '.txt'],...
                '\t');
            lib_model_dat1.(fname)=M;
        end
    end
    
    % GAPPED MODEL
    for model_orderi=0:model_order_max-1
        for i_gap_pos=1:model_orderi
            i=1;%LEFT SIDE
            for i_gap_len=1:l-model_orderi-i
                fname=[File1 '_4_' num2str(model_orderi) ...
                    '_order_model_' num2str(i) '_' num2str(i_gap_pos) '_' ...
                    num2str(i_gap_len)];
                M=dlmread([File1 '_4.' num2str(model_orderi) ...
                    '.order.model.' num2str(i) '_' num2str(i_gap_pos) '_' ...
                    num2str(i_gap_len) '.txt'],'\t');
                lib_model_dat1.(fname)=M;
            end
            ix=1;%RIGHT SIDE
            for i_gap_len=1:l-model_orderi-ix
                i=l-model_orderi-i_gap_len;
                fname=[File1 '_4_' num2str(model_orderi) ...
                    '_order_model_' num2str(i) '_' num2str(i_gap_pos) '_' ...
                    num2str(i_gap_len)];
                M=dlmread([File1 '_4.' num2str(model_orderi) ...
                    '.order.model.' num2str(i) '_' num2str(i_gap_pos) '_' ...
                    num2str(i_gap_len) '.txt'],'\t');
                lib_model_dat1.(fname)=M;
            end
        end
    end
else
    for model_orderi=0:model_order_max-1
        fprintf('Shorter model: %d\n',model_orderi);
        i=1;%LEFT SIDE
        fname=[File1 '_4_' num2str(model_orderi) '_order_model_' num2str(i)];
        M=dlmread([File1 '_4.' num2str(model_orderi) '.order.model.' num2str(i) '.txt'],'\t');
        lib_model_dat1.(fname)=M;
        i=l-model_orderi;%RIGHT SIDE
        fname=[File1 '_4_' num2str(model_orderi) '_order_model_' num2str(i)];
        M=dlmread([File1 '_4.' num2str(model_orderi) '.order.model.' num2str(i) '.txt'],'\t');
        lib_model_dat1.(fname)=M;
    end
end
warning on;
fprintf('Done Processing Library\n');
end
