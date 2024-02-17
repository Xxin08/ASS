% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Auther: Xin Xiong
%  Information Engineering University, Institute of Surveying and Mapping
%  No.62, Kexue avenue, Zhengzhou City, Henan Provice, China
%  versionASS for multimodal image matching V1.0
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
warning('off');
addpath('.\Images');

%% data sets
disp('Visible-SAR');
str = 'VS_7';
% str = 'VS_9';
% str = 'VS_4';

disp(str);
str_ref = [str '_1.tif'];
str_sen= [str '_2.tif'];
str_CP = [str '_CPp.txt'];

image_1=imread(str_ref);
image_2=imread(str_sen);


%% Convert input image format
[~,~,L1]=size(image_1);
[~,~,L2]=size(image_2);
if(L1==3)
    image_11=rgb2gray(image_1);
else
    image_11=image_1;
end
if(L2==3)
    image_22=rgb2gray(image_2);
else
    image_22=image_2;
end

%Converted to floating point data
image_11=im2single(image_11);
image_22=im2single(image_22);

%% set parameters
t1=clock;
o = 8;
Nr = 3;
PS = 84;
K = 0.008; % Number of extracted feature points = K * image_width * image_height
cform='affine';% 'similarity','affine','perspective'

%% ASS Matching method
tic;
[kps1, U1, M1] = ASS_detector(image_11, o, PS, K);
[kps2, U2, M2] = ASS_detector(image_22, o, PS, K);
disp(['ASS detector',num2str(toc),'s']);
tic;
ssf1 = ASS_descriptor(M1, kps1, o, PS, Nr, 1);
ssf2 = ASS_descriptor(M2, kps2, o, PS, Nr, 2);
ssf2.des=ASS_transdesR(ssf1.des,ssf2.des,2*o);
disp(['ASS descriptor',num2str(toc),'s']);

%% Matching
tic;
[solution,rmse,cor1,cor2]=ASS_match(ssf1.des,double(ssf1.kps),ssf2.des,double(ssf2.kps),cform);
disp(['NNDR matching',num2str(toc),'s']);

t2=clock;
disp(['Total time',num2str(etime(t2,t1)),'s']);

%% Matching result
appendimages(image_1,image_2,cor1,cor2);

%% Registration result
image_fusion(image_1,image_2,inv(solution));















