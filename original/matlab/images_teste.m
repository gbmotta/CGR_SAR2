%% 
tic
%Define the directory containing the FASTA files
fasta_dir = 'E:\FASTA\Gabriel\unique_fasta\';

% Load the file containing the names of the FASTA files you want to read
file_names_unique = importdata('E:\FASTA\Gabriel\unique_fasta\fasta_names.txt');

mat_fasta = struct('Header', {}, 'Sequence', {});

% Loop through the file names and read each corresponding FASTA file
for m = 1:length(file_names_unique)
    file_path_unique = [fasta_dir, file_names_unique{m}];
    fasta_data_unique = fastaread(file_path_unique);
    mat_fasta= [mat_fasta;fasta_data_unique];

end

clearvars fasta_data_unique 

toc
%%
tic

[tf, loc] = ismember(test_data.VirusName, {mat_fasta.Header});
mat_fasta_clean = mat_fasta(loc,:);
cell_fasta = (struct2cell(mat_fasta_clean))';
seqsCGR_teste = zeros(64,64,1,height(test_data));

clearvars mat_fasta loc mat_fasta_clean

toc
%%
tic

for y=1:length(cell_fasta)
    %disp(y)
    p  = CGR(cell_fasta{y,2});
 
    k=64;
    count = histcounts2(p(1,:),p(2,:),k);
    c = rot90(count);

    b=8; %Número de bits
    CGRI = floor(((2^b-1)/max(max(c)))*c);
    seqsCGR_teste(:,:,1,y) = CGRI;
    
end

%save(nome_arquivo_mat, 'seqsCGR_img', '-v7.3');
toc