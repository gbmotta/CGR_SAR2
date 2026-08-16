% % 
% There are some commented lines as %#ok<NASGU>,  %#ok<AGROW>,
% %#ok<SAGROW>, %#ok<PREALL> 
% that refers to removing the warning for
% pre allocate variables, as it's usuless as we can see here: 
% https://www.mathworks.com/matlabcentral/answers/394671-no-performance-improvement-with-preallocating-for-object-arrays

% Here we'll load the principal csv, read it and prepare for further 
% modifications.

tic

%Change the working directory.
cd E:\FASTA\Gabriel

nome_arquivo = 'merged'; %The name of the principal csv file

nome_arquivo_csv = strcat(nome_arquivo,'_labels.csv'); % The csv file

% Read the principal csv and put in dataset1.
dataset_csv = readtable(nome_arquivo_csv, 'Delimiter', ',');

% Transform the columns in strings to the next analysis.
dataset_csv.VirusName = string(dataset_csv.VirusName);
dataset_csv.Variant = string(dataset_csv.Variant);   

toc

%%
tic
% Here we'll create the final dataset that will be compared with the fasta
% files to create the images and thier annotations. We'll clean the csv
% dataset and get the samples of variants of our interest for the CNN
% analysis.

% Create a table with the variants and how many samples of each exists in
% the dataset_csv
table_variant_count = table(categories(categorical(dataset_csv.Variant)), ...
    countcats(categorical(dataset_csv.Variant)), ...
    'VariableNames',{'Variant', 'Count'});

% Create a table with median of the samples of variants and remove the
% samples with variants classified as UNKNOWN
table_variant_median = table_variant_count ...
(table_variant_count.Count >= median(table_variant_count.Count),:); 
table_variant_median(strcmp(table_variant_median.Variant, ...
'UNKNOWN'),:) = [];

% Here we create a clean csv logical data set with only variants with a number of
% samples above the median of all variants and without samples marked with
% variant UNKNOWN
clean_dataset_csv = ismember(dataset_csv.Variant, ...
cellstr(table_variant_median.Variant));

% Here we create a clean csv dataset named clean_csv_dataset2 with the
% logical values that we get in the line above,
clean_csv_dataset2 = dataset_csv(clean_dataset_csv,:);

% Here we'll see which variants left from the median analysis.
categories_samples = unique(table_variant_median.Variant);

% Create empty variables for increase speed in the for loop below. 
positions = cell(length(categories_samples), 1);
table_rand_seqs = [];

% Creation of a dataset that will be used to generate the images with
% random samples of all the selected variants.
for i = 1:length(categories_samples)
    positions = find(clean_csv_dataset2.Variant == categories_samples(i));
    random_position = randperm(length(positions))';
    random_sequence = positions(random_position(1:min...
    (table_variant_median.Count)));
    table_rand_seqs = [table_rand_seqs;random_sequence]; %#ok<AGROW>
end

% The dataset with random samples.
table_final = clean_csv_dataset2(table_rand_seqs,:);

% Remove unused variables.
clearvars clean_dataset_csv dataset_csv table_rand_seqs clean_csv_dataset2 ...
positions random_sequence random_position
toc
toc


%%

% Here we'll divide the dataset in two, one with 10% of the samples,
% choosen randomly tha will be used as a test dataset in future analysis
% and the validation dataset, with 90% of the samples, that will be usued
% to validate the CNN.

tic

test_percentage = 0.1;
test_data = []; %#ok<NASGU>
validation_data = []; %#ok<NASGU>
indices = randperm(height(table_final), round(height(table_final)...
*test_percentage));
test_data = table_final(indices,:); 
validation_data = setdiff(table_final,test_data);
% Remove unused variables.
clearvars table_final indices
toc

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
    mat_fasta= [mat_fasta;fasta_data_unique]; %#ok<AGROW> 

end

clearvars fasta_data_unique 
%%

tic

[tf, loc] = ismember(validation_data.VirusName, {mat_fasta.Header});
mat_fasta_clean = mat_fasta(loc,:);
cell_fasta = (struct2cell(mat_fasta_clean))';
seqsCGR_train = zeros(64,64,1,height(validation_data));

clearvars mat_fasta loc

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
    seqsCGR_train(:,:,1,y) = CGRI;
    
end

%save(nome_arquivo_mat, 'seqsCGR_img', '-v7.3');
toc

%% 
tic
% Get unique variant labels and their indices
[variant_labels, ~, variant_indices] = unique(validation_data.Variant);

% Initialize LabelExperiment_CGR with zeros
LabelExperiment_CGR = zeros(size(validation_data.Variant)); %#ok<PREALL> 

% Set the corresponding indices
LabelExperiment_CGR = variant_indices - 1;

% Convert to categorical
LabelExperiment_CGR = categorical(LabelExperiment_CGR);

toc

%% 
tic
% k-fold
kfoldsize = 10; 
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

layer12 = fullyConnectedLayer(length(categories_samples));
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
validation_dados = cell(1:kfoldsize);
validation_cat = cell(1:kfoldsize);
val_med = 0;
loss_med = 0;

for i=1:kfoldsize

VectorTrainingData = seqsCGR_train(:,:,:,cvobj.training(i));
VectorTrainingLabel = LabelExperiment_CGR(cvobj.training(i));

VectorValidationData = seqsCGR_train(:,:,:,cvobj.test(i));
VectorValidationLabel = LabelExperiment_CGR(cvobj.test(i));

validation_dados{i,1} = VectorValidationData;
validation_cat{i,1} = VectorValidationLabel;

options = trainingOptions('adam', ...
        'MiniBatchSize',128,...
        'MaxEpochs',40, ...
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

val_med(i) = vinfo{i}.FinalValidationAccuracy; %#ok<SAGROW> 
loss_med(i) = vinfo{i}.FinalValidationLoss; %#ok<SAGROW> 
end

Val_final = (sum(val_med)/kfoldsize);
stdev = std(val_med);

erro_final = (sum(loss_med)/kfoldsize);
std_erro = std(loss_med);

Resultados_erro_acu_realm = table(Val_final, stdev, erro_final, std_erro, 'VariableNames',{'Acurácia Validação' 'Desvio Padrão Acurácia da Validação' 'Erro Validação' 'Desvio Padrão do Erro da Validação'});

toc