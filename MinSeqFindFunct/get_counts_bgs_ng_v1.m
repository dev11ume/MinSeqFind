% FUNCTION TO COUNT AND FIND BG FOR NON-GAPPED K-MER SEQUENCES 

function [ctsfr,ctsfl,bgsf,inx,all_len]...
    =get_counts_bgs_ng_v1(seq_colp5n,l,lib_nmer,count,File2x,...
    min_nmer,max_nmer,max_cut,Llib,Rlib,flip3to5,lib_model_dat1,strnums_plus)
fprintf('\nCounting Nmers and getting BG...\n');

% Assuming no Ns IN seq_colp5n
seq_colp5nr=5-seq_colp5n(:,end:-1:1);
lLlib=length(Llib);
lRlib=length(Rlib);

aa=[1 0 2 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0];

% COUNTING half_nmer GAP half_nmer
[l_seq,~]=size(seq_colp5n);

cts=sparse(5^max_nmer,1);
cts1=sparse(5^max_nmer,1);
cts2=sparse(5^max_nmer,1);
bgs=sparse(5^max_nmer,1);

all_len=sparse(5^max_nmer,1);

for nmer=min_nmer:max_nmer
    fprintf('nmer:%d\n',nmer);
    fivesi = power(5,0:nmer-1);
    tot_l=nmer;
    ctsi=sparse(5^max_nmer,1);
    ctsi1=sparse(5^max_nmer,1);
    ctsi2=sparse(5^max_nmer,1);
    
    % COUNTING WITH LEFT CONSTANT REGION
    for i=1:lLlib
        inx=sum(seq_colp5n(:,i:i+nmer-1).*fivesi(ones(l_seq,1),:),2)+1;
        ctsi1=ctsi1+sparse(inx,1,count,5^max_nmer,1);
    end
    % COUNTING RANDOM REGION
    for i=1+lLlib:l-tot_l+1+lLlib
        inx=sum(seq_colp5n(:,i:i+nmer-1).*fivesi(ones(l_seq,1),:),2)+1;
        ctsi=ctsi+sparse(inx,1,count,5^max_nmer,1);
    end
    % COUNTING WITH RIGHT CONSTANT REGION
    for i=l+lLlib-tot_l+1+1:l+lLlib+lRlib-tot_l+1
        inx=sum(seq_colp5n(:,i:i+nmer-1).*fivesi(ones(l_seq,1),:),2)+1;
        ctsi2=ctsi2+sparse(inx,1,count,5^max_nmer,1);
    end
    
    % COUNTING REVERSE COMPLEMENT
    % COUNTING WITH LEFT CONSTANT REGION
    for i=1:lRlib
        inx=sum(seq_colp5nr(:,i:i+nmer-1).*fivesi(ones(l_seq,1),:),2)+1;
        ctsi1=ctsi1+sparse(inx,1,count,5^max_nmer,1);
    end
    % COUNTING RANDOM REGION
    for i=1+lRlib:l-tot_l+1+lRlib
        inx=sum(seq_colp5nr(:,i:i+nmer-1).*fivesi(ones(l_seq,1),:),2)+1;
        ctsi=ctsi+sparse(inx,1,count,5^max_nmer,1);
    end
    % COUNTING WITH RIGHT CONSTANT REGION
    for i=l+lRlib-tot_l+1+1:l+lLlib+lRlib-tot_l+1
        inx=sum(seq_colp5nr(:,i:i+nmer-1).*fivesi(ones(l_seq,1),:),2)+1;
        ctsi2=ctsi2+sparse(inx,1,count,5^max_nmer,1);
    end
    
    %Find palidromes and divide counts
    inx=find(ctsi);
    inx_pal=get_palindrome_inx_ng_v1(inx,tot_l);
    ctsi(inx_pal)=ctsi(inx_pal)/2;
    inx1=find(ctsi1);
    inx_pal=get_palindrome_inx_ng_v1(inx1,tot_l);
    ctsi1(inx_pal)=ctsi1(inx_pal)/2;
    inx2=find(ctsi2);
    inx_pal=get_palindrome_inx_ng_v1(inx2,tot_l);
    ctsi2(inx_pal)=ctsi2(inx_pal)/2;
    
    %Removing low count sequences
    inx2=find(ctsi(inx)<max_cut);
    ctsi(inx(inx2))=0;
    
    % GET BG FOR SEQUENCES WITH >0 COUNTS
    bg_inx=find(ctsi);
    if ~isempty(bg_inx)
        seqs_bg=get_index_from_seqs_withN_ng_v1(bg_inx,tot_l);
        seq_colp5n_bg=aa(cell2mat(seqs_bg)-64);
        if flip3to5 % FLIPPED AND EXCHANGED RLIB LLIB
            bg=get_lib_prediction_v1(fliplr(seq_colp5n_bg),tot_l,...
                File2x,lib_nmer,l,fliplr(Rlib),fliplr(Llib),0,0,0,lib_model_dat1,strnums_plus);
        else
            bg=get_lib_prediction_v1(seq_colp5n_bg,tot_l,File2x,lib_nmer,l,Llib,Rlib,0,0,0,lib_model_dat1,strnums_plus);
        end
        
        bg_inx_rev=get_rev_inx_ng_v1(bg_inx,tot_l);
        ix_to_add=find(bg_inx_rev~=bg_inx); % BECAUSE WE REMOVED REV-COMPL
        seqs_bg=get_index_from_seqs_withN_ng_v1(bg_inx_rev,tot_l);
        seq_colp5n_bg=aa(cell2mat(seqs_bg)-64);
        if flip3to5 % FLIPPED AND EXCHANGED RLIB LLIB
            bg_rev=get_lib_prediction_v1(fliplr(seq_colp5n_bg),tot_l,...
                File2x,lib_nmer,l,fliplr(Rlib),fliplr(Llib),0,0,0,lib_model_dat1,strnums_plus);
        else
            bg_rev=get_lib_prediction_v1(seq_colp5n_bg,tot_l,File2x,lib_nmer,l,Llib,Rlib,0,0,0,lib_model_dat1,strnums_plus);
        end
        
        cts=cts+ctsi;
        bgs=bgs+sparse(bg_inx,1,bg,5^max_nmer,1);
        bgs=bgs+sparse(bg_inx(ix_to_add),1,bg_rev(ix_to_add),5^max_nmer,1);
        all_len=all_len+sparse(bg_inx,1,tot_l,5^max_nmer,1);
        
        % CONSIDERING ONLY INX THOSE WHICH HAVE ENOUGH COUNTS FOR RANDOM REGION
        cts1=cts1+sparse(bg_inx,1,ctsi1(bg_inx),5^max_nmer,1);
        cts2=cts2+sparse(bg_inx,1,ctsi2(bg_inx),5^max_nmer,1);
    end
end

inx=find(cts);

all_len=full(all_len(inx));

ctsfr=full(cts(inx));
ctsf1=full(cts1(inx));
ctsf2=full(cts2(inx));
ctsfl=ctsf1+ctsf2;
bgsf=full(bgs(inx));

fprintf('Done\n');
end
