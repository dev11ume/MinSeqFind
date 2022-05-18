% CREATER MATRIX M FOR LIBRARY USED FOR ESTIMATION OF SEQUENCE ABUNDANCE IN LIBRARY 
function [M,fname]=get_M_matrix_lib_v2(seq_colp5n,interpol,model_order,lib_cut,File1,i,i_gap_pos,i_gap_len)

[l_seq,~]=size(seq_colp5n);
File1=File1(1:end-4);
fives = power(5,0:model_order);
mask_bits=dec2bin(0:2^(model_order+1)-1);
cts1=zeros(5^(model_order+1),1);
for j=1:2^(model_order+1)
    xx=bsxfun(@times,seq_colp5n,str2num(mask_bits(j,:)')');
    inx=sum(double(xx).*fives(ones(l_seq,1),:),2);
    cts1=cts1+histc(inx+1,1:5^(model_order+1));
end

if interpol==1
    for model_orderx=1:model_order
        ixt=5^(model_orderx):5^(model_orderx+1)-1;
        for j=1:model_order-model_orderx+1 % STARTING POSITION OF SEQ.
            ixt2=ixt*5^(j-1); 
            ctsi1=cts1(ixt2+1);
            ix_cut=find(ctsi1<lib_cut);
            if ~isempty(ix_cut)
                ixt2l=rem(ixt2(ix_cut),5^(model_orderx+j-1));
                ixt2r=floor(ixt(ix_cut)/5)*5^j;
                ixt2c=5^j*floor(ixt2l/5^j);
                ctsi1(ix_cut)=(cts1(ixt2r+1)./cts1(ixt2c+1)).*cts1(ixt2l+1);
            end
            cts1(ixt2+1)=ctsi1;
        end
    end
    File1=[File1 '_i1_' num2str(lib_cut) ];
end

seq_short=rem((0:5^(model_order+1)-1),5^(model_order));
M=[full(cts1) full(cts1(seq_short+1))];

if i_gap_len==0
    fname=[File1 '_4_' num2str(model_order) '_order_model_' num2str(i)];
    dlmwrite([File1 '_4.' num2str(model_order) '.order.model.' num2str(i) '.txt'],M,...
        'delimiter', '\t', 'precision', 8);
else
    fname=[File1 '_4_' num2str(model_order) ...
        '_order_model_' num2str(i) '_' num2str(i_gap_pos) '_' ...
        num2str(i_gap_len)];
    dlmwrite([File1 '_4.' num2str(model_order) ...
        '.order.model.' num2str(i) '_' num2str(i_gap_pos) '_' ...
        num2str(i_gap_len) '.txt'],M,'delimiter', '\t', 'precision', 8);
end

end
