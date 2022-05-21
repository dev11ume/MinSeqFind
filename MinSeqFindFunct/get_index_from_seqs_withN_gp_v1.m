% FUNCTION- GET SEQUENES FROM INDEXES FOR SEQUENCES WITH N WITH GAPS
function seqs2x=get_index_from_seqs_withN_gp_v1(inx,tot_l,max_half_nmer)
mot='NACGT';

% NOTE- half_nmerj_pos CAN BE DERIVED FROM GAP AND INX
gap=floor((inx-1)/(5^(2*max_half_nmer)));
seqs2x=cell(length(inx),1);
if ~isempty(find(gap~=0, 1))
    inxx=rem((inx(gap~=0)-1),(5^(2*max_half_nmer)));
    inx1=rem(inxx,5^max_half_nmer);
    inx2=floor(inxx/(5^max_half_nmer));
    
    r1=rem(floor(inx1*power(5,1-max_half_nmer:0)),5)+1;
    seqs1=mot(r1(:,end:-1:1));
    
    r2=rem(floor(inx2*power(5,1-tot_l:0)),5)+1;
    seqs2=mot(r2(:,end:-1:1));
    
    max_gap=max(gap);
    for i=1:max_gap
        inxx=find(gap(gap~=0)>=i);
        seqs2(inxx,2:end)=seqs2(inxx,1:end-1);
        if i==1
            seqs2(inxx,1)='N';
        end
    end
    for i=1:max_half_nmer
        inxx=find(seqs1(:,i)~='N');
        if ~isempty(inxx)
          seqs2(inxx,i+1:end)=seqs2(inxx,i:end-1);
          seqs2(inxx,i)=seqs1(inxx,i);
        end
    end
    seqs2x(gap~=0)=cellstr(seqs2);
end

if ~isempty(find(gap==0, 1))
seqs2x(gap==0)=get_index_from_seqs_withN_ng_v1(inx(gap==0),tot_l);
end
end
