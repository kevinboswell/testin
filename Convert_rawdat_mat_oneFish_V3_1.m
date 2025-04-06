% Reads all raw BEM .dat data belonging to one fish and consolidates them 
% into one .mat file.
% Each .dat file corresponds to the TS of one frequency at one orientation
% angle.
%This assumes that .dat files are in the format "XXX_freq_angle.dat"

%Camilo Roa, FIU, 2023.
% 2025, Fixed problems reading frequency and angle from filename

%% 
clearvars
close all

%% Select folder for one fish

fishFolder=uigetdir(pwd,'Fish folder'); %Select where all .dat files are for a specific fish

d = dir([fishFolder, '/*.dat']); %List .dat files
n=length(d);                     %Count number of .dat files

% Max and minimum angles based on models generated from BEMM. Currently set as: 
% Max angle (44 degrees) and Min angle (-46 degrees) in 2 degree increments
% Max freq (200 kHz) and min freq (30 kHz) in 2 kHz increments

angle_min = -46;
angle_max = 44;
Angle=angle_min:2:angle_max;

freq_min = 30000;
freq_max = 200000;
Frequency=freq_min:2000:freq_max;

formatSpec = '%s';  %Specify reading format as string
TS=zeros(length(Frequency),length(Angle)); %Initialize TS matrix

%Loop all files in folder
for ii=1:n
    fileID=fopen(fullfile(fishFolder,d(ii).name),'r');    %Open .dat file
    A=fscanf(fileID,formatSpec); %Read file
    fclose(fileID); %Close .dat file
    
    [~,fileName]=fileparts(d(ii).name); %Extract name of .dat file
    strParts=strsplit(fileName,'_');
    idxAngle=Angle==str2double(strParts(end)); %Find index for current angle
    idxFreq=Frequency==str2double(strParts(end-1)); %Find index for current frequency
    dataNum=str2double(split(A(2:end-1),')(')); %Convert data from string to numerical
    
    idxAnglesData=imag(dataNum)==0; %Find all real values
    anglesData=round(dataNum(idxAnglesData)*180/pi);
    idxSigmaData=find(anglesData==180+str2double(strParts(end)));
    TS(idxFreq,idxAngle)=20*log10(abs(dataNum(idxSigmaData+361))); %TS for angle and freq
    
    clear fileID A fishFreqAngle
end


%% === Prepare the first dataset: mean TS over frequency range, vs angle ===
TS_linear_angles = 20.^(TS(:, idxAngle) / 10); % Linear TS for selected angles
TS_linear_selected_freq = TS_linear_angles(idxFreq, :); % Subset frequencies

nFreq = sum(idxFreq);

mean_TS_angle_linear = mean(TS_linear_selected_freq, 1);
std_TS_angle_linear = std(TS_linear_selected_freq, 0, 1);
sem_TS_angle_linear = std_TS_angle_linear ./ sqrt(nFreq);

mean_TS_angle_dB = 20 * log10(mean_TS_angle_linear);
std_TS_angle_dB = 20 * log10(mean_TS_angle_linear + std_TS_angle_linear) - mean_TS_angle_dB;
sem_TS_angle_dB = 20 * log10(mean_TS_angle_linear + sem_TS_angle_linear) - mean_TS_angle_dB;

angle_min = -10;
angle_max = 20;
idxAngle = (Angle >= angle_min) & (Angle <= angle_max);
selectedAngles = Angle(idxAngle);
nAngles = sum(idxAngle);


%% === Prepare the second dataset: mean TS over angle range, vs frequency ===
TS_linear_freq = 20.^(TS(:, idxAngle) / 10);

mean_TS_freq_linear = mean(TS_linear_freq, 2);
std_TS_freq_linear = std(TS_linear_freq, 0, 2);
sem_TS_freq_linear = std_TS_freq_linear ./ sqrt(nAngles);

mean_TS_freq_dB = 20 * log10(mean_TS_freq_linear);
std_TS_freq_dB = 20 * log10(mean_TS_freq_linear + std_TS_freq_linear) - mean_TS_freq_dB;
sem_TS_freq_dB = 20 * log10(mean_TS_freq_linear + sem_TS_freq_linear) - mean_TS_freq_dB;

selectedFrequencies = Frequency / 1000;


%% === Ensure row vectors to avoid fill() errors ===
selectedAngles = selectedAngles(:)';
mean_TS_angle_dB = mean_TS_angle_dB(:)';
sem_TS_angle_dB = sem_TS_angle_dB(:)';
std_TS_angle_dB = std_TS_angle_dB(:)';

selectedFrequencies = selectedFrequencies(:)';
mean_TS_freq_dB = mean_TS_freq_dB(:)';
sem_TS_freq_dB = sem_TS_freq_dB(:)';
std_TS_freq_dB = std_TS_freq_dB(:)';

idxFreq = (Frequency / 1000 >= freq_min) & (Frequency / 1000 <= freq_max);
idxAngle = (Angle >= angle_min) & (Angle <= angle_max);

% if ~any(idxFreq)
%     uialert(fig, 'No frequencies found in the selected range.', 'Frequency Error');
%     return;
% end
% 
% if ~any(idxAngle)
%     uialert(fig, 'No angles found in the selected range.', 'Angle Error');
%     return;
% end




%% Save TS to .mat
newFileName=extractBefore(fileName,'_');

newMatPath=fullfile(fishFolder,'MatFiles');
if ~exist(newMatPath,'dir')
        mkdir(newMatPath);
end

newFullFileName=fullfile(newMatPath,newFileName);
save(newFullFileName,"TS","Angle","Frequency", "fileName")

%% Save image
newImagPath=fullfile(fishFolder,'ImageFiles');
if ~exist(newImagPath,'dir')
        mkdir(newImagPath);
end

newImagFileName=fullfile(newImagPath,[newFileName,'.jpg']);

%% === Create Figure 1  ---TS by Frequency and Angle ---
figure;
imagesc(Angle,Frequency/1000,TS)
set(gca,'YDir','normal')
xlabel('Incidence angle [deg]');ylabel('Frequency [kHz]')
clim([-70 -20])
h=colorbar;
h.Label.Position(1)=4;
ylabel(h,'TS [dB]','Fontsize',14,'Rotation',270)
colormap("jet")
%set(gca, 'Visible', 'off')        
%set(gca,'position',[0 0 1 1],'units','normalized')
saveas(gcf,newImagFileName);


%% === Create Figure 2 -- Combined statisticial subplots ===
figure;

% --- Subplot 1: Mean TS vs Angle ---
subplot(2, 1, 1);
hold on;

% Shaded SEM region
%  fill([selectedAngles, fliplr(selectedAngles)], ...
%            [mean_TS_angle_dB - sem_TS_angle_dB, fliplr(mean_TS_angle_dB + sem_TS_angle_dB)], ...
%            [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

%  STD shaded region 
  fill([selectedAngles, fliplr(selectedAngles)], ...
            [mean_TS_angle_dB - std_TS_angle_dB, fliplr(mean_TS_angle_dB + std_TS_angle_dB)], ...
            [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);



plot(selectedAngles, mean_TS_angle_dB, 'b-', 'LineWidth', 2);

xlabel('Incidence Angle [deg]');
ylabel('Mean TS [dB]');
title(sprintf('Mean TS vs Angle\n(%.1f–%.1f kHz)', freq_min, freq_max));
ylim([-60 -35]);
grid on;
hold off;


% --- Subplot 2: Mean TS vs Frequency ---
subplot(2, 1, 2);
hold on;

% Shaded SEM region
%fill([selectedFrequencies, fliplr(selectedFrequencies)], ...
%     [mean_TS_freq_dB - sem_TS_freq_dB, fliplr(mean_TS_freq_dB + sem_TS_freq_dB)], ...
%     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Light blue

% STD shaded region
 fill([selectedFrequencies, fliplr(selectedFrequencies)], ...
      [mean_TS_freq_dB - std_TS_freq_dB, fliplr(mean_TS_freq_dB + std_TS_freq_dB)], ...
      [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3); % Light red

plot(selectedFrequencies, mean_TS_freq_dB, 'b-', 'LineWidth', 2);

xlabel('Frequency [kHz]');
ylabel('Mean TS [dB]');
title(sprintf('Mean TS vs Frequency\n(%.0f° to %.0f°)', min(selectedAngles), max(selectedAngles)));
ylim([-60 -35]);
grid on;
hold off;


% === Save the combined figure ===
combinedFigName = fullfile(newImagPath, [newFileName, '_TS_Angle_and_Frequency_Subplots.jpg']);
saveas(gcf, combinedFigName);
