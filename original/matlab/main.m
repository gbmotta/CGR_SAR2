%%
%Iniciar a contagem do tempo de execução do script, ao final, terá um "toc"
tic

cd E:\FASTA\Gabriel\unique_fasta\

nome_arquivo = 'merged';
nome_arquivo_csv = strcat(nome_arquivo,'_labels.csv');
nome_arquivo_fasta = strcat(nome_arquivo,'.fasta');

cd E:\FASTA\Gabriel\
dataset1 = readtable(nome_arquivo_csv, 'Delimiter', ',');

dataset1.VirusName = string(dataset1.VirusName);
dataset1.Variant = string(dataset1.Variant);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

table_variant_count= [];
table_variant_coun_median= [];

table_variant_count = table(categories(categorical(dataset1.Variant)),countcats(categorical(dataset1.Variant)),'VariableNames',{'Variant', 'Count'});

table_variant_count_median = table_variant_count(table_variant_count.Count >= median(table_variant_count.Count),:); 

table_variant_count_median(strcmp(table_variant_count_median.Variant, 'UNKNOWN'),:) = [];

clean_dataset1 = ismember(dataset1.Variant, cellstr(table_variant_count_median.Variant));

clean_dataset2 = dataset1(clean_dataset1,:);

categories_variant = unique(table_variant_count_median.Variant);

positions = cell(length(categories_variant), 1);

table_rand_seqs = [];
minimum_variant_count = min(table_variant_count_median.Count);

for i = 1:length(categories_variant)
    category = categories_variant(i);
    positions = find(clean_dataset2.Variant == category);
    random_position = randperm(minimum_variant_count)';
    random_sequence = positions(random_position);
    table_rand_seqs = [table_rand_seqs;random_sequence];
end

table_final = clean_dataset2(table_rand_seqs,:);

% Loop through each row of the data
for i = 1:size(table_final, 1)
    % Split the string into words
    variant_type = strsplit(table_final.Variant{i});
    
    % Extract the first two words and join them into a string
    variant_name = join(variant_type(1:2), ' ');
    
    % Update the table with the extracted information
    table_final.Variant_New{i} = variant_name;
end

table_final.Variant = string(table_final.Variant_New);
table_final = removevars(table_final, 'Variant_New');

toc

%%

tic

LabelExperiment_CGR = [];

LabelExperiment_CGR(table_final.Variant == 'VOC Alpha') = 0;
LabelExperiment_CGR(table_final.Variant == 'VOC Beta') = 1;
LabelExperiment_CGR(table_final.Variant == 'VOC Delta') = 2;
LabelExperiment_CGR(table_final.Variant == 'VOC Gamma') = 3;
LabelExperiment_CGR(table_final.Variant == 'VOI Epsilon') = 4;
LabelExperiment_CGR(table_final.Variant == 'VOI Eta') = 5;
LabelExperiment_CGR(table_final.Variant == 'VOI Iota') = 6;
LabelExperiment_CGR(table_final.Variant == 'VOI Kappa') = 7;
LabelExperiment_CGR(table_final.Variant == 'VOI Lambda') = 8;
LabelExperiment_CGR(table_final.Variant == 'VOI Mu') = 9;
LabelExperiment_CGR(table_final.Variant == 'VOI Zeta') = 10;
LabelExperiment_CGR(table_final.Variant == 'VOC Omicron') = 11;

LabelExperiment_CGR = categorical(LabelExperiment_CGR);

toc

%%
tic
cd E:\FASTA\Gabriel\unique_fasta

unique_fasta_directory = 'E:\FASTA\Gabriel\unique_fasta\';
file_selected_fasta = 'selected_fasta.mat';
selected_fasta_path = fullfile(unique_fasta_directory, file_selected_fasta);

if exist(selected_fasta_path, 'file') == 2

    load selected_fasta.mat;
    
else
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
    if exist(selected_fasta_path, 'file') == 0
        table_mat_fasta = struct2table(mat_fasta);
        save('table_mat_fasta', 'selected_fasta', '-v7.3');
    end
end
toc

%%
tic

[tf, loc] = ismember(table_final.VirusName, table_mat_fasta.Header);
mat_fasta_clean = table_mat_fasta(loc,:);
cell_fasta = (table2cell(mat_fasta_clean));
seqsCGR_img = {}; 
seqsCGR_train_1 = zeros(64,64,1,height(mat_fasta_clean));

toc

%%

tic

for i=1:length(cell_fasta)
    p  = CGR(cell_fasta{i,2});
 
    k=64;
    count = histcounts2(p(1,:),p(2,:),k);
    c = rot90(count);

    b=8; %Número de bits
    CGRI = floor(((2^b-1)/max(max(c)))*c);
    seqsCGR_train_1(:,:,1,i) = CGRI;
    
end

cd E:\FASTA\Gabriel\Results_Experimento1

%Caso queira salvar os arquivos mat 
%save(nome_arquivo_mat, 'seqsCGR_img', '-v7.3');
toc

%%
%Para kfold
kfoldsize = 10; %numero de kfolds - por exemplo, k=10 gera uma divisão de 90%/10%
cvobj = cvpartition (LabelExperiment_CGR,'k',kfoldsize);

layer0 = imageInputLayer([64 64 1],'Normalization','rescale-zero-one','Name','GenomeImage');

layer1 = convolution2dLayer(7,32,'Padding','same');
layer2 = reluLayer;
layer3 = maxPooling2dLayer(2,'Stride',2);

layer4 = convolution2dLayer(5,64,'Padding','same');
layer5 = reluLayer;
layer6 = maxPooling2dLayer(2,'Stride',2);

layer7 = convolution2dLayer(5,64,'Padding','same');
layer8 = reluLayer;
layer9 = maxPooling2dLayer(2,'Stride',2);

layer10 = fullyConnectedLayer(64);
layer10a = dropoutLayer(0.4);

layer11 = fullyConnectedLayer(32);
layer11a = dropoutLayer(0.4);

layer12 = fullyConnectedLayer(length(categories_variant));
layer13 = softmaxLayer;
layer14 = classificationLayer;

layers = [layer0 ...
              layer1  layer2  layer3 ...
              layer4  layer5  layer6 ...
              layer7  layer8  layer9 ...
              layer10  layer10a ...
              layer11 layer11a layer12 layer13 layer14];

%Opções e Laço
vnet = cell(1:kfoldsize);
vinfo = cell(1:kfoldsize);
teste_dados = cell(1:kfoldsize);
teste_cat = cell(1:kfoldsize);
val_med = 0;
loss_med = 0;

for i=1:kfoldsize

VectorTrainingData = seqsCGR_train_1(:,:,:,cvobj.training(i));
VectorTrainingLabel = LabelExperiment_CGR(cvobj.training(i));

VectorValidationData = seqsCGR_train_1(:,:,:,cvobj.test(i));
VectorValidationLabel = LabelExperiment_CGR(cvobj.test(i));

teste_dados{i,1} = VectorValidationData;
teste_cat{i,1} = VectorValidationLabel;

options = trainingOptions('adam', ...
        'MiniBatchSize',128,...
        'MaxEpochs',10, ...
        'ValidationData',{VectorValidationData VectorValidationLabel}, ...
        'ValidationFrequency',8, ...
        'Verbose',true, ...
        'InitialLearnRate',0.001, ...
        'Shuffle','every-epoch');
        
        
        %'Plots','training-progress',...
        %'ExecutionEnvironment', 'cpu');
        %'ResetInputNormalization',false);
        %'GradientDecayFactor',0.75);
        
[net,info] = trainNetwork(VectorTrainingData,VectorTrainingLabel,layers,options);  
vnet{i}=net;
vinfo{i} = info;

val_med(i) = vinfo{i}.FinalValidationAccuracy;
loss_med(i) = vinfo{i}.FinalValidationLoss;
end

%%
Val_final = (sum(val_med)/kfoldsize)
stdev = std(val_med)

erro_final = (sum(loss_med)/kfoldsize)
std_erro = std(loss_med)

Resultados_erro_acu = table([Val_final], [stdev], [erro_final], [std_erro], 'VariableNames',{'Acurácia Validação' 'Desvio Padrão Acurácia da Validação' 'Erro Validação' 'Desvio Padrão do Erro da Validação'})

toc