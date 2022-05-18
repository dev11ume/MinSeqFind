% AUTHOR - DEVESH BHIMSARIA
% PROGRAM TO GET MinSeqs BY ANALYSIS OF SEQUENCING DATA 

% FILE CONTAINING LIST OF ALL SAMPLES (SEE README FILE FOR MORE INFORMATION)
list_file='list_NRs_bhimsaria_1_demo.txt';

% SETTING UP PARAMETERS TO RUN THE MINSEQ FUNCTION

% MINSEQ RELATED PARAMETERS
param.min_half_nmer=4; % MINIMUM MONOMER SIZE IN MINSEQ 
param.max_half_nmer=8; % MAXIMUM MONOMER SIZE IN MINSEQ; SET max_half_nmer=-1 TO RUN 6 TO 15 MER WITHOUT GAPS OR SPACERS
param.min_nmer=6; % USED ONLY IF max_half_nmer=-1; MINIMUM K-MER SIZE TO RUN WITHOUT GAPS
param.max_nmer=15;  % USED ONLY IF max_half_nmer=-1; MAXIMUM K-MER SIZE TO RUN WITHOUT GAPS
param.max_gap=100; % MAXIMUM GAP BETWEEN MONOMERS
param.pow_factor=2; % POWER FACTOR (SEE SUPPLEMENTARY FILE 1)
param.max_cut_arr=50; % CUTOFF FOR THRESHOLDING - MINIMUM NUMBER OF SEQUENCES TO CONSIDER THE MINSEQ
param.flip3to5=0; % 1= FLIP SEQUENCES 3' TO 5'
param.seeds_to_cont=100000; % NUMBER OF MINSEQS TO SAVE 

% LIBRARY RELATED PARAMETERS
param.lib_nmer=5; % ORDER OF MARKOV MODEL
param.gap_mod=0; % 1=GAPPED MARKOV MODEL , 0=NON-GAPPED MARKOV MODEL
param.interpol=1; % 1=INTERPOLATED MARKOV MODEL
param.lib_cut=50; % CUTOFF FOR THRESHOLDING - MINIMUM NUMBER OF SEQUENCES TO CONSIDER THE LIBRARY INTERPOLATED MARKOV MODEL

% PWM RELATED PARAMETERS
param.max_cut_pwm=50; % CUTOFF FOR THRESHOLDING - MINIMUM NUMBER OF SEQUENCES TO CONSIDER MINSEQ FOR PWM
param.min_ct=20; % CUTOFF FOR THRESHOLDING - MINIMUM NUMBER OF SEQUENCES TO CONSIDER FLANKING POSITIONS FOR PWM
param.seeds_to_cont2=1000; % NUMBER OF MINSEQS TO CONTINUE FOR PWM CALCULATION
param.num_motifs_max=20; % NUMBER OF PWM MOTIFS TO EXTRACT
param.num_Ns=3; % FLANKING LENGTH FOR PWMs

% THRESHOLDING TO REMOVE SOME MINSEQS
param.enrich_cut=0; % ABSOLUTE ENRICHMENT CUT OF MINSEQS 
param.enrich_perc_cut=0; % PERCENTAGE OF MAX ENRICHMENT TO CUT FOR MINSEQS 
param.weighted_enrich_perc_cut=0;% PERCENTAGE OF MAX WEIGHTED ENRICHMENT TO CUT FOR MINSEQS 

seq_to_print=100000;% NUMBER OF MINSEQS TO PRINT TO FILE

path_files='MinSeqFindFunct';
addpath(path_files);

C=textscan_db1(list_file,'%s%s%*d%s%s%s','\t');
TFs_file=C{1};Library_files=C{2};TF_sample=C{3};Llibs=C{4};Rlibs=C{5};

dir_nm=[path_files '-op'];

% CREATE DIRECTORIES FOR OUTPUT FILES
mkdir(dir_nm);
mkdir([dir_nm '/TXT']);
mkdir([dir_nm '/CHEN']);
mkdir([dir_nm '/MinSeqs']);

for ii=1:length(TFs_file)
    % CREATE DIRECTORY FOR OUTPUT MAT FILE
    mkdir([dir_nm '/' TFs_file{ii} '_' Library_files{ii}]);
    
    % PARAMETERS & FILE NAMES TO PASS TO THE FUNCTION
    File1=[TFs_file{ii} '.txt'];
    LibFile=[Library_files{ii} '.txt'];
    Llib=Llibs{ii};
    Rlib=Rlibs{ii};
    if Llib=='x'
        Llib='';
    end
    if Rlib=='x'
        Rlib='';
    end
    
    % RUN MINSEQ PROGRAM TO GET MINSEQS
    % minseqfind_1(File1,LibFile,Llib,Rlib,param);
    
    % movefile([TFs_file{ii} '-' Library_files{ii} '-res1f.mat'],[dir_nm '/' TFs_file{ii} '_' Library_files{ii} '/resf.mat']);   
    write_pwms_file_v1(dir_nm, [TFs_file{ii} '_' Library_files{ii}]);
    write_pwms_file_chen_v1(dir_nm, [TFs_file{ii} '_' Library_files{ii}]);
    write_minseq_to_file_v1(dir_nm, [TFs_file{ii} '_' Library_files{ii}],1,seq_to_print);
end
rmpath(path_files);
