%Change directory to the Spatial modelling package directory
cd("/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My Drive/Marmoset_shared_folder/Shiny_plots_Erin/SpatialModelling/")

addpath(genpath('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/code/gpml-matlab-v3.6-2015-07-07'))

%Load CS7
[OBJ1,section] = LoadCS5('3D');
[D,Locations,XYZ,CellType,ShotsCS5] = LoadShots('CS5');
[OutputCS5] = loadCS5Scaffold(D,Locations,ShotsCS5);

%Load CS6
[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');
[OutputCS6] = loadCS6Scaffold(D,Locations,ShotsCS6);

%Edit file paths to the outputs from genFinalMappingsCorrelations.R script
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/FinalDylanMapping/Correlation_int_locations1_wall_rr2.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/FinalDylanMapping/Correlation_int_locations2_wall_rr2.csv'
File3 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/FinalDylanMapping/Correlation_int_set1_wall_rr2.csv'

%Loop through each cell and use a GP to estimate the embryo wide
%correlation
[targInd,C] = ProjectDataCorrelation(File1,File2,File3,OutputCS6);
[M,V] = loadProbsCorrelation(OutputCS6,targInd,C,OBJ2);

[targIndCS5,CCS5] = ProjectDataCorrelationCS5(File1,File2,File3,OutputCS5);
[MCS5,VCS5] = loadProbsCorrelationCS5(OutputCS5,targIndCS5,CCS5,OBJ1);

S2 = importdata(File2);

%Subset specific cells types for visulastion
indPLAXA = find(strcmp(S2.textdata(2:end,3),"PLAXA")==1)
indPLAXA1 = find(strcmp(S2.textdata(2:end,3),"PLAXA1")==1)

indConv = find(strcmp(S2.textdata(2:end,3),"Conventional")==1)
indConv1 = find(strcmp(S2.textdata(2:end,3),"Conventional1")==1)
indConv2 = find(strcmp(S2.textdata(2:end,3),"Conventional2")==1)

indHyp = find(strcmp(S2.textdata(2:end,3),"HYPO")==1)
indTr = find(strcmp(S2.textdata(2:end,3),"cmTSCPAVS")==1)
indTb1 = find(strcmp(S2.textdata(2:end,3),"OKAEP5esc")==1)
indTb2 = find(strcmp(S2.textdata(2:end,3),"newTSP3")==1)

[OBJ1b,a1,b1,OutputCS5b] = transformCS5(OBJ1,OutputCS5,'all');

OBJ1d = transformOBJ(OBJ1,pi/2,[0,1,0]);
OBJ1e = transformOBJ(OBJ1d,'flipX');
OBJ1f = transformOBJ(OBJ1e,'flipZ');
OBJ1test = transformOBJ(OBJ1f,6.0879,[0.6144,0.7715,0.1653]);

S1 = importdata(File1);
S2 = importdata(File2);
C1 = importdata(File3);

%Also load the blastocyst
OBJ0b=importdata('./Data/SpatialData/Blastocyst_slice_2.mat'); %read_wobj('Medium_resolution_CS5_version10-fordylan-310.obj')
meh = load('./Data/SpatialData/lBlHR.mat')
OBJ0b = meh.O4;
OBJ0b.vertices = Quaternion3(pi*1.28,[0,1,0],OBJ0b.vertices);

R5 = C1.data(find(strcmp(S1.textdata(2:end,3),'Tb_CS3')==1),:);
R3 = C1.data(find(strcmp(S1.textdata(2:end,3),'Hyp_CS3')==1),:);
R1 = C1.data(find(strcmp(S1.textdata(2:end,3),'Epi_CS3')==1),:);
R0 = C1.data(find(strcmp(S1.textdata(2:end,3),'ICM_CS3')==1),:);
Rm = C1.data(find(strcmp(S1.textdata(2:end,3),'cMor_CS2')==1),:); 

%%%%% 
%First look at Niave
subplot(3,5,1);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R5(:,indPLAXA1)))./(size(R5(:,indPLAXA1),1)*size(R5(:,indPLAXA1),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R3(:,indPLAXA1)))./(size(R3(:,indPLAXA1),1)*size(R3(:,indPLAXA1),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R1(:,indPLAXA1)))./(size(R1(:,indPLAXA1),1)*size(R1(:,indPLAXA1),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.45,.66] )
material shiny
title('PLAXA')

subplot(3,5,2);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R5(:,indConv1)))./(size(R5(:,indConv1),1)*size(R5(:,indConv1),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R3(:,indConv1)))./(size(R3(:,indConv1),1)*size(R3(:,indConv1),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R1(:,indConv1)))./(size(R1(:,indConv1),1)*size(R1(:,indConv1),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny
title('Conv')
 
subplot(3,5,6);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,3)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,8)),1)','FaceColor','interp','LineStyle','none');
% f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,2)),1)','FaceColor','interp','LineStyle','none');
% f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,1)),1)','FaceColor','interp','LineStyle','none');
%  f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,3)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1 b1])
%view([[355.6210,-0.0545]])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.45,.66] )  
material shiny
 
subplot(3,5,7);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indConv1,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indConv1,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indConv1,3)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indConv1,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indConv1,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indConv1,8)),1)','FaceColor','interp','LineStyle','none');
% f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,2)),1)','FaceColor','interp','LineStyle','none');
% f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,1)),1)','FaceColor','interp','LineStyle','none');
%  f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,3)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1 b1])
%view([[355.6210,-0.0545]])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.4,.65] )  
material shiny
 
 
subplot(3,5,11);
f5 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,2)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1test.objects(24).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,8)),1)','FaceColor','interp','LineStyle','none');

% f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,2)),1)','FaceColor','interp','LineStyle','none');
% f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,1)),1)','FaceColor','interp','LineStyle','none');
%  f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,3)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view([375.7415,-26.7022])

% view([a1 b1])
%view([[355.6210,-0.0545]])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.45,.65] )  
material shiny
 
subplot(3,5,12);
f5 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indConv1,2)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1test.objects(24).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indConv1,8)),1)','FaceColor','interp','LineStyle','none');

% f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,2)),1)','FaceColor','interp','LineStyle','none');
% f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,1)),1)','FaceColor','interp','LineStyle','none');
%  f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indPLAXA1,3)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view([375.7415,-26.7022])

% view([a1 b1])
%view([[355.6210,-0.0545]])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.4,.65] )  
material shiny
 
 
 
 

subplot(3,5,3);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R5(:,indTr)))./(size(R5(:,indTr),1)*size(R5(:,indTr),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R3(:,indTr)))./(size(R3(:,indTr),1)*size(R3(:,indTr),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R1(:,indTr)))./(size(R1(:,indTr),1)*size(R1(:,indTr),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.56] )
material shiny
title('PAVS')

subplot(3,5,8);
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,3)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,8)),1)','FaceColor','interp','LineStyle','none');

% f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,2)),1)','FaceColor','interp','LineStyle','none');
% f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,1)),1)','FaceColor','interp','LineStyle','none');
%  f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,3)),1)','FaceColor','interp','LineStyle','none');

axis equal
axis off
view([a1 b1])
%view([[355.6210,-0.0545]])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.56] )
material shiny 

subplot(3,5,13);
f5 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,2)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1test.objects(24).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,8)),1)','FaceColor','interp','LineStyle','none');

axis equal
axis off
view([375.7415,-26.7022])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.56] )
material shiny
 
 
 
subplot(3,5,4);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R5(:,indTb2)))./(size(R5(:,indTb2),1)*size(R5(:,indTb2),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R3(:,indTb2)))./(size(R3(:,indTb2),1)*size(R3(:,indTb2),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R1(:,indTb2)))./(size(R1(:,indTb2),1)*size(R1(:,indTb2),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.4] )
material shiny
title('TPS3')

subplot(3,5,5);
f1 = patch('Faces',OBJ0b.objects(4).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R5(:,indTb1)))./(size(R5(:,indTb1),1)*size(R5(:,indTb1),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',OBJ0b.objects(12).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R3(:,indTb1)))./(size(R3(:,indTb1),1)*size(R3(:,indTb1),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ0b.objects(8).data.vertices,'Vertices',OBJ0b.vertices,'FaceVertexCData',(sum(sum(R1(:,indTb1)))./(size(R1(:,indTb1),1)*size(R1(:,indTb1),2))).*ones(length(OBJ0b.vertices),1),'FaceColor','interp','LineStyle','none');
%view([164.2550,-40.2050])
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.15,.45] )
material shiny
title('TBLC')
  
 
subplot(3,5,9);
 
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb2,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb2,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb2,3)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb2,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb2,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb2,8)),1)','FaceColor','interp','LineStyle','none');

% f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,2)),1)','FaceColor','interp','LineStyle','none');
% f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,1)),1)','FaceColor','interp','LineStyle','none');
%  f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,3)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1 b1])
%view([[355.6210,-0.0545]])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.4] )
material shiny
 
subplot(3,5,10);
 
f1 = patch('Faces',OBJ1b.objects(4).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb1,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb1,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb1,3)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(16).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb1,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(12).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb1,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb1,8)),1)','FaceColor','interp','LineStyle','none');
% f1 = patch('Faces',OBJ1b.objects(20).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,2)),1)','FaceColor','interp','LineStyle','none');
% f5 = patch('Faces',OBJ1b.objects(24).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,1)),1)','FaceColor','interp','LineStyle','none');
%  f2 = patch('Faces',OBJ1b.objects(8).data.vertices,'Vertices',  OBJ1b.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTr,3)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view([a1 b1])
%view([[355.6210,-0.0545]])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.15,.45] )
material shiny
 
 
subplot(3,5,14);
f5 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb2,2)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1test.objects(24).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb2,8)),1)','FaceColor','interp','LineStyle','none');
 
axis equal
axis off
%view([375.7415,-26.7022])
view([360.6926,-20.7542])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.4] )
material shiny
 
subplot(3,5,15);
f5 = patch('Faces',OBJ1test.objects(20).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb1,2)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1test.objects(24).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb1,8)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ1test.objects(16).data.vertices,'Vertices',  OBJ1test.vertices,'FaceVertexCData',mean(cell2mat(MCS5(indTb1,5)),1)','FaceColor','interp','LineStyle','none');

axis equal
axis off
%view([375.7415,-26.7022])
view([360.6926,-20.7542])
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.15,.45] )
material shiny
print('-dpng',['./Mapping.png'],'-r800')

%Now preimplantation
load('MissingSections/eBlHR2.mat')
O5.vertices = Quaternion3(pi*0.95,[0,1,0],O5.vertices); 
load('MissingSections/cMHR.mat')
load('/Users/christopherpenfold/Library/CloudStorage/GoogleDrive-cap76@cam.ac.uk/My Drive/Marmoset_shared_folder/3D_Modellinng/LowRes/OM.mat')

subplot(3,5,1);
f1 = patch('Faces',OM.objects(4).data.vertices,'Vertices',OM.vertices,'FaceVertexCData',(sum(sum(Rm(:,indPLAXA1)))./(size(Rm(:,indPLAXA1),1)*size(Rm(:,indPLAXA1),2))).*ones(length(OM.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny
 
subplot(3,5,2);
f1 = patch('Faces',OM.objects(4).data.vertices,'Vertices',OM.vertices,'FaceVertexCData',(sum(sum(Rm(:,indConv1)))./(size(Rm(:,indConv1),1)*size(Rm(:,indConv1),2))).*ones(length(OM.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.45,.65] )
material shiny

subplot(3,5,3);
f1 = patch('Faces',OM.objects(4).data.vertices,'Vertices',OM.vertices,'FaceVertexCData',(sum(sum(Rm(:,indTr)))./(size(Rm(:,indTr),1)*size(Rm(:,indTr),2))).*ones(length(OM.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.56] )
material shiny

subplot(3,5,4);
f1 = patch('Faces',OM.objects(4).data.vertices,'Vertices',OM.vertices,'FaceVertexCData',(sum(sum(Rm(:,indTb2)))./(size(Rm(:,indTb2),1)*size(Rm(:,indTb2),2))).*ones(length(OM.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.4] )
material shiny
  
subplot(3,5,5);
f1 = patch('Faces',OM.objects(4).data.vertices,'Vertices',OM.vertices,'FaceVertexCData',(sum(sum(Rm(:,indTb1)))./(size(Rm(:,indTb1),1)*size(Rm(:,indTb1),2))).*ones(length(OM.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.15,.45] )
material shiny
  
%print('-dpng',['~/Desktop/Mapping2.png'],'-r800')

subplot(3,5,6);
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R0(:,indPLAXA1)))./(size(R0(:,indPLAXA1),1)*size(R0(:,indPLAXA1),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indPLAXA1)))./(size(R5(:,indPLAXA1),1)*size(R5(:,indPLAXA1),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indPLAXA1)))./(size(R5(:,indPLAXA1),1)*size(R5(:,indPLAXA1),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny
 
subplot(3,5,7);
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R0(:,indConv1)))./(size(R0(:,indConv1),1)*size(R0(:,indConv1),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indConv1)))./(size(R5(:,indConv1),1)*size(R5(:,indConv1),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indConv1)))./(size(R5(:,indConv1),1)*size(R5(:,indConv1),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny
 
subplot(3,5,8);
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R0(:,indTr)))./(size(R0(:,indTr),1)*size(R0(:,indTr),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indTr)))./(size(R5(:,indTr),1)*size(R5(:,indTr),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indTr)))./(size(R5(:,indTr),1)*size(R5(:,indTr),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.56]  )
material shiny
 
subplot(3,5,9);
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R0(:,indTb2)))./(size(R0(:,indTb2),1)*size(R0(:,indTb2),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indTb2)))./(size(R5(:,indTb2),1)*size(R5(:,indTb2),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indTb2)))./(size(R5(:,indTb2),1)*size(R5(:,indTb2),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.3,.4] )
material shiny

subplot(3,5,10);
f1 = patch('Faces',O5.objects(4).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R0(:,indTb1)))./(size(R0(:,indTb1),1)*size(R0(:,indTb1),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f0 = patch('Faces',O5.objects(10).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indTb1)))./(size(R5(:,indTb1),1)*size(R5(:,indTb1),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',O5.objects(8).data.vertices,'Vertices',O5.vertices,'FaceVertexCData',(sum(sum(R5(:,indTb1)))./(size(R5(:,indTb1),1)*size(R5(:,indTb1),2))).*ones(length(O5.vertices),1),'FaceColor','interp','LineStyle','none');
axis equal
axis off
camlight('left')
material dull 
colormap(parula);
set(gca,'clim',[.15,.45] )
material shiny

print('-dpng',['./PreMapping.png'],'-r800')

%Check out the CS6
[OBJ2b,a1,b1] = transformCS6(OBJ2,'all');

subplot(3,5,1);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,3)),1)','FaceColor','interp','LineStyle','none'); 
f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,6)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny
 
subplot(3,5,2);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,3)),1)','FaceColor','interp','LineStyle','none'); 
f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,6)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny 
 
subplot(3,5,3);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,3)),1)','FaceColor','interp','LineStyle','none'); 
f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,6)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny
 
subplot(3,5,4);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,3)),1)','FaceColor','interp','LineStyle','none'); 
f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,6)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.3,.4] )
material shiny
 
subplot(3,5,5);
f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,2)),1)','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,3)),1)','FaceColor','interp','LineStyle','none'); 
f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,6)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.15,.45] )
material shiny
 
subplot(3,5,6);
% f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,2)),1)','FaceColor','interp','LineStyle','none');
%f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,3)),1)','FaceColor','interp','LineStyle','none'); 
%f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,6)),1)','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indPLAXA1,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny
 
subplot(3,5,7);
%f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,2)),1)','FaceColor','interp','LineStyle','none');
%f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,3)),1)','FaceColor','interp','LineStyle','none'); 
% f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,6)),1)','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indConv1,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny

subplot(3,5,8);%
%f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,2)),1)','FaceColor','interp','LineStyle','none');
%f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,3)),1)','FaceColor','interp','LineStyle','none'); 
%f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,6)),1)','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTr,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.4,.65] )
material shiny
 
subplot(3,5,9);
%f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,2)),1)','FaceColor','interp','LineStyle','none');
%f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,3)),1)','FaceColor','interp','LineStyle','none'); 
%f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,6)),1)','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb2,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.3,.4] )
material shiny
 
subplot(3,5,10);
%f1 = patch('Faces',OBJ2b.objects(20).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,1)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(4).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,2)),1)','FaceColor','interp','LineStyle','none');
%f2 = patch('Faces',OBJ2b.objects(16).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,3)),1)','FaceColor','interp','LineStyle','none'); 
%f5 = patch('Faces',OBJ2b.objects(8).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,4)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(32).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,5)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(24).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,6)),1)','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2b.objects(12).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,7)),1)','FaceColor','interp','LineStyle','none');
f5 = patch('Faces',OBJ2b.objects(28).data.vertices,'Vertices',  OBJ2b.vertices,'FaceVertexCData',mean(cell2mat(M(indTb1,8)),1)','FaceColor','interp','LineStyle','none');
axis equal
axis off
view(a1,b1)
camlight('left')
colormap(parula);
set(gca,'clim',[.15,.45] )
material shiny
 
print('-dpng',['./MappingCS6.png'],'-r800')