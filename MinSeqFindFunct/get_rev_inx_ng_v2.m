% FUNCTION TO GET REVERSE COMPLEMENT INDEXES FOR K-MERS WITHOUT GAPS - V2
function inx_rev=get_rev_inx_ng_v2(inx,tot_l)
l_inx=length(inx);
seq_5nrev=5-rem(floor((inx-1)*power(5,1-tot_l:0)),5);
seq_5nrev(seq_5nrev==5)=0;
fives = power(5,0:tot_l-1);
inx_rev=sum(seq_5nrev.*fives(ones(l_inx,1),:),2);
for i=1:tot_l
    ix=find(rem(inx_rev,5)==0);
    inx_rev(ix)=inx_rev(ix)/5;
end
inx_rev=inx_rev+1;
end
