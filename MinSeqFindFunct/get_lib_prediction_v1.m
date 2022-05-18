% FUNCTION TO GET PREDICTED ABUNDANCE OF SEQUENCES IN THE LIBRARY 
function bg0=get_lib_prediction_v1(seq_colp5n,pwm_l,File2x,lib_nmer,l,Llib,Rlib,gapm,half_nmeri,gap_mod,lib_model_dat1)

if gap_mod==0
    gapm=0; half_nmeri=0;
end
aa=[1 0 2 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0];
l_Llib=length(Llib);
l_Rlib=length(Rlib);

Llib_seq_col=aa(Llib-64);
Rlib_seq_col=aa(Rlib-64);

lib_nmer=lib_nmer+1; % IF 5TH ORDER MARKOV MODEL IS USED, THE LENGTH OF FINAL SEQ IS 6MER

C=textscan_mod_v1([File2x num2str(lib_nmer-1) '.order.model.' num2str(1) '.txt'],'%f%f','\t');
count_6mer2=C{1};

fname=[strrep(File2x,'.','_') num2str(lib_nmer-1) '_order_model_' num2str(1)];
C=lib_model_dat1.(fname);
count_6mer=C(:,1);

tot_count=max(count_6mer);

[l_seq,~]=size(seq_colp5n);

bg0=zeros(l_seq,1);
for i=1:l+l_Llib+l_Rlib-pwm_l+1
    % CHECK IF SEQ_COL MATCHES CONSTANT REGION IF IN CONSTANT REGION
    match=ones(l_seq,1);
    for k=1:min(l_Llib-i+1,pwm_l)
        match=match & ( (seq_colp5n(:,k)==0) |  (seq_colp5n(:,k)==Llib_seq_col(i+k-1)) );
    end
    for k=max(l+l_Llib-i+1+1,1):pwm_l
        match=match & ( (seq_colp5n(:,k)==0) |  (seq_colp5n(:,k)==Rlib_seq_col(k+i-l-l_Llib-1)) );
    end
    ixx=find(match);
    l_seqx=length(ixx);
    
    if l_seqx>0
        ix3=[i:i+half_nmeri-1 i+half_nmeri+gapm:i+pwm_l-1]; % NOT INCLUDING Ns, BUT CONSTANT
        ix4=ix3((ix3>=l_Llib+1) & (ix3<=l_Llib+l)); % NOT INCLUDING Ns & CONSTANT
        
        % IF ONLY IN CONSTANT REGION
        if isempty(ix4)
            bg1x=tot_count*ones(l_seqx,1);
            bg0(ixx)=bg0(ixx)+bg1x;
        else
            lib_nmeri=lib_nmer;
            if length(ix4)<lib_nmer
                lib_nmeri=length(ix4);
            end
            bg1x=zeros(l_seqx,1);
            bg2x=zeros(l_seqx,1);
            for ind=1:length(ix4)-lib_nmeri+1
                strt=ix4(ind)-l_Llib;
                File_name=[File2x num2str(lib_nmeri-1) '.order.model.' num2str(strt) '.txt'];
                ix4_diff=ix4(ind+1:ind+lib_nmeri-1)-ix4(ind:ind+lib_nmeri-2);
                pos_gapi=find(ix4_diff>1);
                if ~isempty(pos_gapi) % IF NOT IN CONTINUATION
                    if length(pos_gapi)>1
                        disp('error multi gaps');
                    end
                    gapx=ix4_diff(pos_gapi)-1;
                    File_name=[File2x num2str(lib_nmeri-1) '.order.model.' ...
                        num2str(strt) '_' num2str(pos_gapi) '_' num2str(gapx) '.txt'];
                end
                
                fname=[strrep(File2x,'.','_') num2str(lib_nmeri-1) '_order_model_' num2str(strt)];
                
                C=lib_model_dat1.(fname);
                count_6mer=C(:,1);
                Lib_count_5meri=C(:,2);
                
                Lib_prob_5meri=count_6mer./Lib_count_5meri;
                fivesi = power(5,0:lib_nmeri-1);
                
                indn=sum(seq_colp5n(ixx,ix4(ind:ind+lib_nmeri-1)-i+1).*fivesi(ones(l_seqx,1),:),2)+1;
                if ind==1
                    bg1x(:)=Lib_count_5meri(indn);
                    bg1x(bg1x==0)=1;
                end
                bg1x=bg1x.*Lib_prob_5meri(indn);
            end
            bg0(ixx)=bg0(ixx)+bg1x;
        end
    end
end

end
