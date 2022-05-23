% FUNCTION TO COUNT AND FIND BG FOR GAPPED MinSeqs

function [ctsfr,ctsfl,bgsf,inx,all_len,all_len_nox]...
    =get_counts_bgs_gp_v1(seq_colp5n,l,lib_nmer,count,File2x,...
    min_half_nmer,max_half_nmer,max_gap,max_cut,Llib,Rlib,gap_mod,flip3to5,lib_model_dat1,strnums_plus)
fprintf('\nCounting Nmers and getting BG...\n');

% NO Ns IN seq_colp5n
seq_colp5nr=5-seq_colp5n(:,end:-1:1);
lLlib=length(Llib);
lRlib=length(Rlib);

aa=[1 0 2 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0];

% COUNTING half_nmer GAP half_nmer
[l_seq,~]=size(seq_colp5n);
max_gap=min(max_gap,l-2*min_half_nmer);

m_size=(5^(2*max_half_nmer))*(max_gap+1);

ct_times=sparse(m_size,1);
cts=sparse(m_size,1);
cts1=sparse(m_size,1);
cts2=sparse(m_size,1);
bgs=sparse(m_size,1);

all_len=sparse(m_size,1);
all_len_nox=sparse(m_size,1);

for half_nmeri=min_half_nmer:max_half_nmer
    fivesi = power(5,0:half_nmeri-1);
    for half_nmerj=half_nmeri:max_half_nmer
        fivesj = power(5,0:half_nmerj-1);
        for j=0:min(l-half_nmeri-half_nmerj,max_gap)
            fprintf('Nmer-Gap-Nmer=%d\t%d\t%d\n',half_nmeri,j,half_nmerj);
            tot_l=half_nmeri+half_nmerj+j;
            
            ctsi=sparse(m_size,1);
            ctsi1=sparse(m_size,1);
            ctsi2=sparse(m_size,1);
            
            % COUNTING WITH LEFT CONSTANT REGION
            for i=1:lLlib
                inx=j*(5^(2*max_half_nmer))+sum(seq_colp5n(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                    +(5^(min(max_half_nmer,j*1000+half_nmeri)))*...
                    sum(seq_colp5n(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq,1),:),2)+...
                    +1;
                ctsi1=ctsi1+sparse(inx,1,count,m_size,1);
            end
            % COUNTING RANDOM REGION
            for i=1+lLlib:l-tot_l+1+lLlib
                inx=j*(5^(2*max_half_nmer))+sum(seq_colp5n(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                    +(5^(min(max_half_nmer,j*1000+half_nmeri)))*...
                    sum(seq_colp5n(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq,1),:),2)+...
                    +1;
                ctsi=ctsi+sparse(inx,1,count,m_size,1);
            end
            % COUNTING WITH RIGHT CONSTANT REGION
            for i=lLlib+l-tot_l+2:l+lLlib+lRlib-tot_l+1
                inx=j*(5^(2*max_half_nmer))+sum(seq_colp5n(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                    +(5^(min(max_half_nmer,j*1000+half_nmeri)))*...
                    sum(seq_colp5n(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq,1),:),2)+...
                    +1;
                ctsi2=ctsi2+sparse(inx,1,count,m_size,1);
            end
            
            % COUNTING REVERSE COMPLEMENT
            % COUNTING WITH LEFT CONSTANT REGION
            for i=1:lRlib
                inx=j*(5^(2*max_half_nmer))+sum(seq_colp5nr(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                    +(5^(min(max_half_nmer,j*1000+half_nmeri)))*...
                    sum(seq_colp5nr(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq,1),:),2)+...
                    +1;
                ctsi1=ctsi1+sparse(inx,1,count,m_size,1);
            end
            % COUNTING RANDOM REGION
            for i=1+lRlib:l-tot_l+1+lRlib
                inx=j*(5^(2*max_half_nmer))+sum(seq_colp5nr(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                    +(5^(min(max_half_nmer,j*1000+half_nmeri)))*...
                    sum(seq_colp5nr(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq,1),:),2)+...
                    +1;
                ctsi=ctsi+sparse(inx,1,count,m_size,1);
            end
            % COUNTING WITH RIGHT CONSTANT REGION
            for i=lRlib+l-tot_l+2:l+lRlib+lLlib-tot_l+1
                inx=j*(5^(2*max_half_nmer))+sum(seq_colp5nr(:,i:i+half_nmeri-1).*fivesi(ones(l_seq,1),:),2)+...
                    +(5^(min(max_half_nmer,j*1000+half_nmeri)))*...
                    sum(seq_colp5nr(:,i+half_nmeri+j:i+half_nmeri+j+half_nmerj-1).*fivesj(ones(l_seq,1),:),2)+...
                    +1;
                ctsi2=ctsi2+sparse(inx,1,count,m_size,1);
            end
            
            %Find palidromes and divide counts
            inx=find(ctsi);
            if (j==0 && rem(half_nmeri+half_nmerj,2)==0) || half_nmeri==half_nmerj
                inx_pal=get_palindrome_inx_gp_v1(inx,max_half_nmer);
                ctsi(inx_pal)=ctsi(inx_pal)/2;
                inx1=find(ctsi1);
                inx_pal=get_palindrome_inx_gp_v1(inx1,max_half_nmer);
                ctsi1(inx_pal)=ctsi1(inx_pal)/2;
                inx2=find(ctsi2);
                inx_pal=get_palindrome_inx_gp_v1(inx2,max_half_nmer);
                ctsi2(inx_pal)=ctsi2(inx_pal)/2;
            end
            
            %Removing low count sequences
            inx2=find(ctsi(inx)<max_cut);
            ctsi(inx(inx2))=0;
            
            % GET BG FOR SEQUENCES WITH >0 COUNTS
            bg_inx=find(ctsi);
            if ~isempty(bg_inx)
                  
                seqs_bg=get_index_from_seqs_withN_gp_v1(bg_inx,tot_l,max_half_nmer);
                seq_colp5n_bg=aa(cell2mat(seqs_bg)-64);
                if flip3to5 % FLIPPED AND EXCHANGED RLIB LLIB, CHANGED half_nmeri
                    bg=get_lib_prediction_v1(fliplr(seq_colp5n_bg),tot_l,...
                        File2x,lib_nmer,l,fliplr(Rlib),fliplr(Llib),j,half_nmerj,gap_mod,lib_model_dat1,strnums_plus);
                else
                    bg=get_lib_prediction_v1(seq_colp5n_bg,tot_l,File2x,...
                        lib_nmer,l,Llib,Rlib,j,half_nmeri,gap_mod,lib_model_dat1,strnums_plus);
                end
                
                % FOR REVERSE
                bg_inx_rev=get_rev_inx_gp_v2(bg_inx,max_half_nmer);
                ix_to_add=find(bg_inx_rev~=bg_inx); % BECAUSE WE REMOVED REV-COMPL
                seqs_bg=get_index_from_seqs_withN_gp_v1(bg_inx_rev,tot_l,max_half_nmer);
                seq_colp5n_bg=aa(cell2mat(seqs_bg)-64);
                if flip3to5 % FLIPPED AND EXCHANGED RLIB LLIB, CHANGED half_nmerj
                    bg_rev=get_lib_prediction_v1(fliplr(seq_colp5n_bg),tot_l,...
                        File2x,lib_nmer,l,fliplr(Rlib),fliplr(Llib),j,half_nmeri,gap_mod,lib_model_dat1,strnums_plus);
                else
                    bg_rev=get_lib_prediction_v1(seq_colp5n_bg,tot_l,File2x,...
                        lib_nmer,l,Llib,Rlib,j,half_nmerj,gap_mod,lib_model_dat1,strnums_plus);
                end
                
                cts=cts+ctsi;
                
                bgs=bgs+sparse(bg_inx,1,bg,m_size,1);
                bgs=bgs+sparse(bg_inx(ix_to_add),1,bg_rev(ix_to_add),m_size,1);
                
                all_len=all_len+sparse(bg_inx,1,tot_l,m_size,1);
                all_len_nox=all_len_nox+sparse(bg_inx,1,half_nmeri+half_nmerj,m_size,1);
                % CONSIDERING ONLY INX THOSE WHICH HAVE ENOUGH COUNTS FOR RANDOM REGION
                cts1=cts1+sparse(bg_inx,1,ctsi1(bg_inx),m_size,1);
                cts2=cts2+sparse(bg_inx,1,ctsi2(bg_inx),m_size,1);
                % NOTE 4-0-3 AND 3-0-4 WILL BE COUNTED TWICE AND WILL BE NORMALIZED BY ct_times
                ct_times=ct_times+sparse(bg_inx,1,1,m_size,1);
            end
        end
    end
end

inx=find(cts);

ctsfr=full(cts(inx))./full(ct_times(inx));

all_len=full(all_len(inx))./full(ct_times(inx));
all_len_nox=full(all_len_nox(inx))./full(ct_times(inx));

ctsf1=full(cts1(inx))./full(ct_times(inx));
ctsf2=full(cts2(inx))./full(ct_times(inx));
bgsf=full(bgs(inx))./full(ct_times(inx));
ctsfl=ctsf1+ctsf2;

fprintf('Done\n');
end
