% FUNCTION TO GET REVERSE COMPLEMENT INDEXES FOR K-MERS WITHOUT GAPS
function inx_rev=get_rev_inx_ng_v1(inx,tot_l)
l_inx=length(inx);
seq_5nrev=5-rem(floor((inx-1)*power(5,1-tot_l:0)),5);
seq_5nrev(seq_5nrev==5)=0;
fives = power(5,0:tot_l-1);
inx_rev=sum(seq_5nrev.*fives(ones(l_inx,1),:),2)+1;
end
