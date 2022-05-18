% FUNCTION FOR INITIAL DATA PROCESSING 
function [l_seqs_raw,l_seqs2,cont_seqs1,cont_seqs2,seq_colp5n,count,seqs_uq,to_contf]=...
    init_data_process_v1(File1,Llib,Rlib,nmer_cut,LibFile)
fprintf('Initial data processing...');

warning off;
C=textscan_mod_v1(File1,'%s','\t');
warning on;
seqs=C{1};
seqs=upper(seqs);
if length(seqs)>3000000
  seqs=seqs(1:3000000); % SO THAT SYSTEM DOESN'T GO OUT OF RAM --  
end

l=length(seqs{1});
l_seqs_raw=length(seqs);
cont_seqs1=zeros(l_seqs_raw,1);
for i=1:l_seqs_raw
    % THROWING SEQUENCES CONTAINING Ns
    if isempty(find(seqs{i}=='N', 1)); cont_seqs1(i)=1; end
    % REMOVE SOME KNOWN CONTAMINATION SEQUENCES
    if sum(seqs{i}=='CGAATGATGGATTGCAACCG')>=18 ||...
            sum(seqs{i}=='TACCGATTACGTAATTTCGA')>=18 ||...
            sum(seqs{i}=='CGCGAATGACGTCAATCGGA')>=18 ||...
            sum(seqs{i}=='GCTAGATTGCGCAATCCAGT')>=18
        cont_seqs1(i)=0;
    end
end
if strcmpi(Llib,'x')
    Llib='';
end
if strcmpi(Rlib,'x')
    Rlib='';
end
seqs=seqs(cont_seqs1==1);

[seqs_uq,~,seqs_inx_c]=unique(seqs);
l_seqs_uq=length(seqs_uq);
count=histc(seqs_inx_c,1:l_seqs_uq);

% NOTE CAN'T CONVERT >20MERS INTO NUMBERS, PRECISION WILL KILL IT
[seq_colp5n,~]=get_seq_mat_v1(seqs_uq,[],[]);
seq_colp5n_rev=5-seq_colp5n(:,end:-1:1);
map_it='ACGT';
seqs_uq_rev=map_it(seq_colp5n_rev);

ctsc=containers.Map(seqs_uq(1:l_seqs_uq),1:l_seqs_uq);
cbig=find(count>=nmer_cut);

throw_seq=zeros(length(cbig),1);
for i=1:length(cbig)
    seq_revi=seqs_uq_rev(cbig(i),:);
    if isKey(ctsc,seq_revi)
        if count(ctsc(seq_revi))<=(count(cbig(i))/nmer_cut)
            throw_seq(i)=1;
        end
    else
        throw_seq(i)=1;
    end
end

ix_throw=find(throw_seq);
cont_seqs2=ones(l_seqs_uq,1);
ix2=cbig(ix_throw);
for i=1:length(ix_throw)
    seqi=seqs_uq{ix2(i)};
    for j=1:l
        seq_ix1=[seqi(1:j-1) 'A' seqi(:,[j+1:end])];
        seq_ix2=[seqi(1:j-1) 'C' seqi(:,[j+1:end])];
        seq_ix3=[seqi(1:j-1) 'G' seqi(:,[j+1:end])];
        seq_ix4=[seqi(1:j-1) 'T' seqi(:,[j+1:end])];
        if isKey(ctsc,seq_ix1)
            cont_seqs2(ctsc(seq_ix1))=0;
        end
        if isKey(ctsc,seq_ix2)
            cont_seqs2(ctsc(seq_ix2))=0;
        end
        if isKey(ctsc,seq_ix3)
            cont_seqs2(ctsc(seq_ix3))=0;
        end
        if isKey(ctsc,seq_ix4)
            cont_seqs2(ctsc(seq_ix4))=0;
        end
    end
end

count=count(cont_seqs2==1);
seqs_uq=seqs_uq(cont_seqs2==1);
l_seqs2=sum(count);
[seq_colp5n,~]=get_seq_mat_v1(seqs_uq,Llib,Rlib);

cont_seqs=zeros(l_seqs_raw,1);
ix1=find(cont_seqs1==1);
ix2=find(cont_seqs2==1);
cont_seqs(ix1(ix2))=1;
to_contf=find(cont_seqs);

fprintf('Done\n');
end
