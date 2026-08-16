tic

% Diretorio de trabalho
fasta_dir = 'E:\FASTA\Gabriel\unique_fasta\';

% Arquivo txt com os FASTA desejados
file_names_unique = importdata('E:\FASTA\Gabriel\unique_fasta\fasta_names.txt');

mat_fasta = struct('Header', {}, 'Sequence', {});

% Loop para selecionar somentes os FASTA desejados
for i = 1:length(file_names_unique)
    file_path_unique = [fasta_dir, file_names_unique{i}];
    fasta_data_unique = fastaread(file_path_unique);
    mat_fasta= [mat_fasta;fasta_data_unique];
   
end

%Salvar um arquivo mat que sera utilizado na seleção das variantes
table_mat_fasta = struct2table(mat_fasta);
save('table_mat_fasta', 'selected_fasta', '-v7.3');

toc