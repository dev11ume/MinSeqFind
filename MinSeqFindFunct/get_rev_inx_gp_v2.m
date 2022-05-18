% FUNCTION TO GET REVERSE COMPLEMENT INDEXES FOR MINSEQS WITH GAPS 
function inx_rev=get_rev_inx_gp_v2(inx,max_half_nmer)
gap=floor((inx-1)/(5^(2*max_half_nmer)));
inx_rev=zeros(length(gap),1);
if ~isempty(find(gap~=0, 1))
    inxx=rem((inx(gap~=0)-1),(5^(2*max_half_nmer)));
    inx1=rem(inxx,5^max_half_nmer)+1;
    inx2=floor(inxx/(5^max_half_nmer))+1;
    inx1_rev=get_rev_inx_ng_v2(inx1,max_half_nmer)-1;
    inx2_rev=get_rev_inx_ng_v2(inx2,max_half_nmer)-1;
    inx_rev(gap~=0,1)=gap(gap~=0).*(5^(2*max_half_nmer))+inx2_rev+(5^(max_half_nmer))*inx1_rev+1;
end
if ~isempty(find(gap==0, 1))
    inx_rev(gap==0,1)=get_rev_inx_ng_v2(inx(gap==0),2*max_half_nmer);
end
end
