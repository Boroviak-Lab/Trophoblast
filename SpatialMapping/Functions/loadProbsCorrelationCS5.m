function [ConvCS5,ConvfCS5] = loadProbsCorrelationCS5(OutputCS5,targInd,C,OBJ1)

inds = find(targInd~=0);
XYZ = OutputCS5.cleanX(inds,:);
Tissue = OutputCS5.cleanAnotaton(inds);
C = C.data(targInd(inds),:);
%
%keyboard
XtrainCS5 = XYZ;
XtrainCS5(:,1:3) = XtrainCS5(:,1:3) / 400;
%XtestCS6 = D2(:,1:3);
%XtestCS6(:,1:3) =XtestCS6(:,1:3)/400;

ind1 = find(strcmp(Tissue,'Am_CS5')==1);
ind2 = find(strcmp(Tissue,'EmDisc_CS5')==1 | strcmp(Tissue,'EmDisc_Gast_CS5')==1| strcmp(Tissue,'EmDisc_gast_CS5')==1);
ind4 = find(strcmp(Tissue,'PGC_CS5')==1);
ind3 = find(strcmp(Tissue,'VE_CS5')==1);
ind5 = find(strcmp(Tissue,'ExMes_CS5')==1);
%ind6 = find(strcmp(Tissue,'VE_CS5')==1);
ind7 = find(strcmp(Tissue,'SYS_CS5')==1);
ind8 = find(strcmp(Tissue,'Tb_CS5')==1);

%GP priors
pg1 = {@priorGauss,0,1};  pg2 = {@priorGauss,0,1};  pg3 = {@priorGauss,log(0.5),1};  pc = {@priorClamped};
prior1.mean = {[]}; prior1.cov  = {pg1;pg2}; prior1.lik  = {pg3};
im1 = {@infPrior,@infExact,prior1};                % inference method


%keyboard
k = 1;
%k = 1;
for i = 1:size(C,2)
    
    Ytrain1 = C(ind1,i);
    Ytrain2 = C(ind2,i);
    Ytrain3 = C(ind3,i);
    Ytrain4 = C(ind4,i);
    Ytrain5 = C(ind5,i);
  %  Ytrain6 = C(ind6,i);
    Ytrain7 = C(ind7,i);
    Ytrain8 = C(ind8,i);
    

    
    Ymean = ([Ytrain1;Ytrain2]);
    Ymean = mean(Ymean(find(Ymean~=0)));
    
    Ymean1 = ([Ytrain1]);
    Ymean1 = mean(Ymean1(find(Ymean1~=0)));
    
    Ymean2 = ([Ytrain2]);
    Ymean2 = mean(Ymean2(find(Ymean2~=0)));

    Ymean3 = ([Ytrain3]);
    Ymean3 = mean(Ymean3(find(Ymean3~=0)));
 
    Ymean4 = ([Ytrain4]);
    Ymean4 = mean(Ymean4(find(Ymean4~=0)));
    Ymean5 = ([Ytrain5]);
    Ymean5 = mean(Ymean5(find(Ymean5~=0)));
    Ymean7 = ([Ytrain7]);
    Ymean7 = mean(Ymean7(find(Ymean7~=0)));
    Ymean8 = ([Ytrain8]);
    Ymean8 = mean(Ymean8(find(Ymean8~=0)));
    
    
    
    Ytrainin1 = Ytrain1-Ymean1; %1 ./ ( 1 + exp(- k*(Ytrain1-Ymean) ));
    Ytrainin2 = Ytrain2-Ymean2; %1 ./ ( 1 + exp(- k*(Ytrain2-Ymean) ));
    Ytrainin3 = Ytrain3-Ymean3; %1 ./ ( 1 + exp(- k*(Ytrain2-Ymean) ));
    Ytrainin4 = Ytrain4-Ymean4; %1 ./ ( 1 + exp(- k*(Ytrain2-Ymean) ));
    Ytrainin5 = Ytrain5-Ymean5; %1 ./ ( 1 + exp(- k*(Ytrain2-Ymean) ));
    Ytrainin7 = Ytrain7-Ymean7; %1 ./ ( 1 + exp(- k*(Ytrain2-Ymean) ));
    Ytrainin8 = Ytrain8-Ymean8; %1 ./ ( 1 + exp(- k*(Ytrain2-Ymean) ));
  
    
        %First do Tb
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin1(:,1));       
        par  = {@meanConst,@covSEiso,'likGauss',XtrainCS5(ind1,:), Ytrainin1(:,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  XtrainCS5(ind1,:), Ytrainin1(:,1), [OBJ1.vertices]/400);        
        ConvCS5{i,1} = m_1'+Ymean1;
        ConvfCS5{i,1} = s_1';  
                
        %substrain = find(LabBinCS6==uLab(j));      
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin2(:,1));       
        par  = {@meanConst,@covSEiso,'likGauss',XtrainCS5(ind2,:), Ytrainin2(:,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  XtrainCS5(ind2,:), Ytrainin2(:,1), [OBJ1.vertices]/400);        
        ConvCS5{i,2} = m_1'+Ymean2;
        ConvfCS5{i,2} = s_1';          
    
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin3(:,1));       
        par  = {@meanConst,@covSEiso,'likGauss',XtrainCS5(ind3,:), Ytrainin3(:,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  XtrainCS5(ind3,:), Ytrainin3(:,1), [OBJ1.vertices]/400);        
        ConvCS5{i,3} = m_1'+Ymean3;
        ConvfCS5{i,3} = s_1';             
        
        
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin3(:,1));       
        par  = {@meanConst,@covSEiso,'likGauss',XtrainCS5(ind4,:), Ytrainin4(:,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  XtrainCS5(ind4,:), Ytrainin4(:,1), [OBJ1.vertices]/400);        
        ConvCS5{i,4} = m_1'+Ymean4;
        ConvfCS5{i,4} = s_1';             
        
        
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin5(:,1));       
        par  = {@meanConst,@covSEiso,'likGauss',XtrainCS5(ind5,:), Ytrainin5(:,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  XtrainCS5(ind5,:), Ytrainin5(:,1), [OBJ1.vertices]/400);        
        ConvCS5{i,5} = m_1'+Ymean5;
        ConvfCS5{i,5} = s_1';             
        
        
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin3(:,1));       
        par  = {@meanConst,@covSEiso,'likGauss',XtrainCS5(ind7,:), Ytrainin7(:,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  XtrainCS5(ind7,:), Ytrainin7(:,1), [OBJ1.vertices]/400);        
        ConvCS5{i,7} = m_1'+Ymean7;
        ConvfCS5{i,7} = s_1';             
        
        
        hyp1.cov  = [log(2); log(0.1) ]; 
        hyp1.lik  = [log(0.1/2)];  
        hyp1.mean = mean(Ytrainin8(:,1));       
        par  = {@meanConst,@covSEiso,'likGauss',XtrainCS5(ind8,:), Ytrainin8(:,1)};
        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  XtrainCS5(ind8,:), Ytrainin8(:,1), [OBJ1.vertices]/400);        
        ConvCS5{i,8} = m_1'+Ymean8;
        ConvfCS5{i,8} = s_1';             
                
        
        
end

return

%First get the most varied cells and 
for i = 1:315
    m1(i,1) = max([ConvCS6{i,1},ConvCS6{i,2}]) - min(max([ConvCS6{i,1};ConvCS6{i,2}]));
end


listType = find(strcmp(Idents2.textdata(2:end,2),"Amnioid_#2")==1);
m2 = m1(listType);
[aa bb] =  sort(m2,'descend')
h = figure(1)
for i = 1:min(length(listType),25)
subplot(5,5,i);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),1}','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),2}','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2C.objects(24).data.vertices,'Vertices',OBJ2C.vertices,'FaceVertexCData',Conv{i,3}','FaceColor','interp','LineStyle','none');
%f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',ConvCS6{listType(i),3}','FaceColor','interp','LineStyle','none');
view([-239.7545,27.0685])
%view([-92.8171,-38.7904])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
end
%savefig('/Volumes/GoogleDrive/Shared drives/Munger et al./Chris'' Sequencing Dungeon/Updated Figures/3DFigures/Am_reps.fig')
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
%print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmD_All.png'])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 30])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/AmB_reps_k=4.png'])


clf
listType = find(strcmp(Idents2.textdata(2:end,2),"Amnioid_#1")==1);
m2 = m1(listType);
[aa bb] =  sort(m2,'descend')
h = figure(1)
for i = 1:min(length(listType),25)
subplot(5,5,i);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),1}','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),2}','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2C.objects(24).data.vertices,'Vertices',OBJ2C.vertices,'FaceVertexCData',Conv{i,3}','FaceColor','interp','LineStyle','none');
%f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',ConvCS6{listType(i),3}','FaceColor','interp','LineStyle','none');
view([-239.7545,27.0685])
%view([-92.8171,-38.7904])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 30])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/Am_reps_k=4.png'])


clf
listType = find(strcmp(Idents2.textdata(2:end,2),"EmDisc_#1")==1);
m2 = m1(listType);
[aa bb] =  sort(m2,'descend')
h = figure(1)
for i = 1:min(length(listType),25)
subplot(5,5,i);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),1}','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),2}','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2C.objects(24).data.vertices,'Vertices',OBJ2C.vertices,'FaceVertexCData',Conv{i,3}','FaceColor','interp','LineStyle','none');
%f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',ConvCS6{listType(i),3}','FaceColor','interp','LineStyle','none');
%view([-239.7545,27.0685])
view([-92.8171,-38.7904])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('right')
material shiny 
material dull 

end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 30])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmD_reps_k=4.png'])




material([1 1 1])


clf
listType = find(strcmp(Idents2.textdata(2:end,2),"EmDisc_#2")==1);
m2 = m1(listType);
[aa bb] =  sort(m2,'descend')
h = figure(1)
for i = 1:min(length(listType),25)
subplot(5,5,i);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),1}','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',ConvCS6{listType(bb(i)),2}','FaceColor','interp','LineStyle','none');
%f5 = patch('Faces',OBJ2C.objects(24).data.vertices,'Vertices',OBJ2C.vertices,'FaceVertexCData',Conv{i,3}','FaceColor','interp','LineStyle','none');
%f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',ConvCS6{listType(i),3}','FaceColor','interp','LineStyle','none');
%view([-239.7545,27.0685])
view([-92.8171,-38.7904])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('right')
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 30])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmDisc_reps_k=4.png'])



%for i = 1:315
%    Xtrainin = Xtrain2;
%    Ytrainin = MD{1}(:,i);%sum(MD{Dataset}(:,CL{i}),2) / sum(DL.data(CL{i}));    
%    for j = 1:length(uLab)
%        substrain = find(LabBin==uLab(j));      
%        hyp1.cov  = [log(2); log(0.1) ]; 
%        hyp1.lik  = [log(0.1/2)];  
%        hyp1.mean = mean(Ytrainin(substrain,1));       
%        par  = {@meanConst,@covSEiso,'likGauss',Xtrainin(substrain,:), Ytrainin(substrain,1)};
%        hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});
%        [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrainin(substrain,:), Ytrainin(substrain,1), [OBJ2.vertices]/400);        
%        Conv{i,j} = m_1';
%        Convf{i,j} = s_1';  
%end
%end



%Averages.......
listType = find(strcmp(Idents2.textdata(2:end,2),"EmDisc_#1")==1);
X_EmD_CS5_1 = ConvCS5{listType(1),1};X_EmD_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_EmD_CS5_1 = X_EmD_CS5_1 + ConvCS5{listType(i),1};
X_EmD_CS5_2 = X_EmD_CS5_2 + ConvCS5{listType(i),2};
end
X_EmD_CS5_1 = X_EmD_CS5_1/i;
X_EmD_CS5_2 = X_EmD_CS5_2/i;
X_EmD_CS6_1 = ConvCS6{listType(1),1};X_EmD_CS6_2 = ConvCS6{listType(1),2};X_EmD_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_EmD_CS6_1 = X_EmD_CS6_1 + ConvCS6{listType(i),1};
X_EmD_CS6_2 = X_EmD_CS6_2 + ConvCS6{listType(i),2};
X_EmD_CS6_3 = X_EmD_CS6_3 + ConvCS6{listType(i),3};
end
X_EmD_CS6_1 = X_EmD_CS6_1/i;
X_EmD_CS6_2 = X_EmD_CS6_2/i;
X_EmD_CS6_3 = X_EmD_CS6_3/i;
X_EmD_CS7_1 = ConvCS7{listType(1),1};X_EmD_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_EmD_CS7_1 = X_EmD_CS7_1 + ConvCS7{listType(i),1};
X_EmD_CS7_2 = X_EmD_CS7_2 + ConvCS7{listType(i),2};
end
X_EmD_CS7_1 = X_EmD_CS7_1/i;
X_EmD_CS7_2 = X_EmD_CS7_2/i;

save('X_EmD_CS6_1_k=4.mat','X_EmD_CS6_1')
save('X_EmD_CS6_2_k=4.mat','X_EmD_CS6_2')
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_EmD_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_EmD_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_EmD_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_EmD_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_EmD_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_EmD_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[.2 .7] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmD_All.png'])

listType = find(strcmp(Idents2.textdata(2:end,2),"EmDisc_#2")==1);
X_EmDisc_CS5_1 = ConvCS5{listType(1),1};X_EmDisc_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_EmDisc_CS5_1 = X_EmDisc_CS5_1 + ConvCS5{listType(i),1};
X_EmDisc_CS5_2 = X_EmDisc_CS5_2 + ConvCS5{listType(i),2};
end
X_EmDisc_CS5_1 = X_EmDisc_CS5_1/i;
X_EmDisc_CS5_2 = X_EmDisc_CS5_2/i;
X_EmDisc_CS6_1 = ConvCS6{listType(1),1};X_EmDisc_CS6_2 = ConvCS6{listType(1),2};X_EmDisc_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_EmDisc_CS6_1 = X_EmDisc_CS6_1 + ConvCS6{listType(i),1};
X_EmDisc_CS6_2 = X_EmDisc_CS6_2 + ConvCS6{listType(i),2};
X_EmDisc_CS6_3 = X_EmDisc_CS6_3 + ConvCS6{listType(i),3};
end
X_EmDisc_CS6_1 = X_EmDisc_CS6_1/i;
X_EmDisc_CS6_2 = X_EmDisc_CS6_2/i;
X_EmDisc_CS6_3 = X_EmDisc_CS6_3/i;
X_EmDisc_CS7_1 = ConvCS7{listType(1),1};X_EmDisc_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_EmDisc_CS7_1 = X_EmDisc_CS7_1 + ConvCS7{listType(i),1};
X_EmDisc_CS7_2 = X_EmDisc_CS7_2 + ConvCS7{listType(i),2};
end
X_EmDisc_CS7_1 = X_EmDisc_CS7_1/i;
X_EmDisc_CS7_2 = X_EmDisc_CS7_2/i;

save('X_EmDisc_CS6_1_k=4.mat','X_EmDisc_CS6_1')
save('X_EmDisc_CS6_2_k=4.mat','X_EmDisc_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_EmDisc_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_EmDisc_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_EmDisc_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_EmDisc_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_EmDisc_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_EmDisc_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[.2 .7] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/EmDisc_All.png'])

listType = find(strcmp(Idents2.textdata(2:end,2),"Amnioid_#1")==1);
X_Am_CS5_1 = ConvCS5{listType(1),1};X_Am_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_Am_CS5_1 = X_Am_CS5_1 + ConvCS5{listType(i),1};
X_Am_CS5_2 = X_Am_CS5_2 + ConvCS5{listType(i),2};
end
X_Am_CS5_1 = X_Am_CS5_1/i;
X_Am_CS5_2 = X_Am_CS5_2/i;
X_Am_CS6_1 = ConvCS6{listType(1),1};X_Am_CS6_2 = ConvCS6{listType(1),2};X_Am_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_Am_CS6_1 = X_Am_CS6_1 + ConvCS6{listType(i),1};
X_Am_CS6_2 = X_Am_CS6_2 + ConvCS6{listType(i),2};
X_Am_CS6_3 = X_Am_CS6_3 + ConvCS6{listType(i),3};
end
X_Am_CS6_1 = X_Am_CS6_1/i;
X_Am_CS6_2 = X_Am_CS6_2/i;
X_Am_CS6_3 = X_Am_CS6_3/i;
X_Am_CS7_1 = ConvCS7{listType(1),1};X_Am_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_Am_CS7_1 = X_Am_CS7_1 + ConvCS7{listType(i),1};
X_Am_CS7_2 = X_Am_CS7_2 + ConvCS7{listType(i),2};
end
X_Am_CS7_1 = X_Am_CS7_1/i;
X_Am_CS7_2 = X_Am_CS7_2/i;
save('X_Am_CS6_1_k=4.mat','X_Am_CS6_1')
save('X_Am_CS6_2_k=4.mat','X_Am_CS6_2')
clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_Am_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_Am_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_Am_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_Am_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_Am_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_Am_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[.2 .7] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/Am_All.png'])

listType = find(strcmp(Idents2.textdata(2:end,2),"Amnioid_#2")==1);
X_AmB_CS5_1 = ConvCS5{listType(1),1};X_AmB_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_AmB_CS5_1 = X_AmB_CS5_1 + ConvCS5{listType(i),1};
X_AmB_CS5_2 = X_AmB_CS5_2 + ConvCS5{listType(i),2};
end
X_AmB_CS5_1 = X_AmB_CS5_1/i;
X_AmB_CS5_2 = X_AmB_CS5_2/i;
X_AmB_CS6_1 = ConvCS6{listType(1),1};X_AmB_CS6_2 = ConvCS6{listType(1),2};X_AmB_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_AmB_CS6_1 = X_AmB_CS6_1 + ConvCS6{listType(i),1};
X_AmB_CS6_2 = X_AmB_CS6_2 + ConvCS6{listType(i),2};
X_AmB_CS6_3 = X_AmB_CS6_3 + ConvCS6{listType(i),3};
end
X_AmB_CS6_1 = X_AmB_CS6_1/i;
X_AmB_CS6_2 = X_AmB_CS6_2/i;
X_AmB_CS6_3 = X_AmB_CS6_3/i;
X_AmB_CS7_1 = ConvCS7{listType(1),1};X_AmB_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_AmB_CS7_1 = X_AmB_CS7_1 + ConvCS7{listType(i),1};
X_AmB_CS7_2 = X_AmB_CS7_2 + ConvCS7{listType(i),2};
end
X_AmB_CS7_1 = X_AmB_CS7_1/i;
X_AmB_CS7_2 = X_AmB_CS7_2/i;
save('X_AmB_CS6_1_k=4.mat','X_AmB_CS6_1')
save('X_AmB_CS6_2_k=4.mat','X_AmB_CS6_2')
clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_AmB_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_AmB_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_AmB_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_AmB_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[.2 .7] )
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_AmB_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_AmB_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[.2 .7] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/Amnioid_bead_All.png'])




%Am lineagesFGF_noMEF
listType = find(strcmp(Idents2.textdata(2:end,2),"FGF_noMEF")==1);
X_FGFNo_CS5_1 = ConvCS5{listType(1),1};X_FGFNo_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_FGFNo_CS5_1 = X_FGFNo_CS5_1 + ConvCS5{listType(i),1};
X_FGFNo_CS5_2 = X_FGFNo_CS5_2 + ConvCS5{listType(i),2};
end
X_FGFNo_CS5_1 = X_FGFNo_CS5_1/i;
X_FGFNo_CS5_2 = X_FGFNo_CS5_2/i;
X_FGFNo_CS6_1 = ConvCS6{listType(1),1};X_FGFNo_CS6_2 = ConvCS6{listType(1),2};X_FGFNo_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_FGFNo_CS6_1 = X_FGFNo_CS6_1 + ConvCS6{listType(i),1};
X_FGFNo_CS6_2 = X_FGFNo_CS6_2 + ConvCS6{listType(i),2};
X_FGFNo_CS6_3 = X_FGFNo_CS6_3 + ConvCS6{listType(i),3};
end
X_FGFNo_CS6_1 = X_FGFNo_CS6_1/i;
X_FGFNo_CS6_2 = X_FGFNo_CS6_2/i;
X_FGFNo_CS6_3 = X_FGFNo_CS6_3/i;
X_FGFNo_CS7_1 = ConvCS7{listType(1),1};X_FGFNo_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_FGFNo_CS7_1 = X_FGFNo_CS7_1 + ConvCS7{listType(i),1};
X_FGFNo_CS7_2 = X_FGFNo_CS7_2 + ConvCS7{listType(i),2};
end
X_FGFNo_CS7_1 = X_FGFNo_CS7_1/i;
X_FGFNo_CS7_2 = X_FGFNo_CS7_2/i;
save('X_FGFNo_CS6_1_k=4.mat','X_FGFNo_CS6_1')
save('X_FGFNo_CS6_2_k=4.mat','X_FGFNo_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_FGFNo_CS5_2' - X_AmB_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_FGFNo_CS6_1' - X_AmB_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_FGFNo_CS6_2' - X_AmB_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_FGFNo_CS6_3' - X_AmB_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_FGFNo_CS7_2' - X_AmB_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_FGFNo_CS7_1' - X_AmB_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/FGF_vs_AmB_All.png'])



listType = find(strcmp(Idents2.textdata(2:end,2),"ActA_noMEF")==1);
X_ActANo_CS5_1 = ConvCS5{listType(1),1};X_ActANo_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_ActANo_CS5_1 = X_ActANo_CS5_1 + ConvCS5{listType(i),1};
X_ActANo_CS5_2 = X_ActANo_CS5_2 + ConvCS5{listType(i),2};
end
X_ActANo_CS5_1 = X_ActANo_CS5_1/i;
X_ActANo_CS5_2 = X_ActANo_CS5_2/i;
X_ActANo_CS6_1 = ConvCS6{listType(1),1};X_ActANo_CS6_2 = ConvCS6{listType(1),2};X_ActANo_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_ActANo_CS6_1 = X_ActANo_CS6_1 + ConvCS6{listType(i),1};
X_ActANo_CS6_2 = X_ActANo_CS6_2 + ConvCS6{listType(i),2};
X_ActANo_CS6_3 = X_ActANo_CS6_3 + ConvCS6{listType(i),3};
end
X_ActANo_CS6_1 = X_ActANo_CS6_1/i;
X_ActANo_CS6_2 = X_ActANo_CS6_2/i;
X_ActANo_CS6_3 = X_ActANo_CS6_3/i;
X_ActANo_CS7_1 = ConvCS7{listType(1),1};X_ActANo_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_ActANo_CS7_1 = X_ActANo_CS7_1 + ConvCS7{listType(i),1};
X_ActANo_CS7_2 = X_ActANo_CS7_2 + ConvCS7{listType(i),2};
end
X_ActANo_CS7_1 = X_ActANo_CS7_1/i;
X_ActANo_CS7_2 = X_ActANo_CS7_2/i;
save('X_ActANo_CS6_1_k=4.mat','X_ActANo_CS6_1')
save('X_ActANo_CS6_2_k=4.mat','X_ActANo_CS6_2')



clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_ActANo_CS5_2' - X_AmB_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_ActANo_CS6_1' - X_AmB_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_ActANo_CS6_2' - X_AmB_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_ActANo_CS6_3' - X_AmB_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_ActANo_CS7_2' - X_AmB_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_ActANo_CS7_1' - X_AmB_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/ActANo_vs_AmB_All.png'])



listType = find(strcmp(Idents2.textdata(2:end,2),"BMP_noMEF")==1);
X_BMPNo_CS5_1 = ConvCS5{listType(1),1};X_BMPNo_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_BMPNo_CS5_1 = X_BMPNo_CS5_1 + ConvCS5{listType(i),1};
X_BMPNo_CS5_2 = X_BMPNo_CS5_2 + ConvCS5{listType(i),2};
end
X_BMPNo_CS5_1 = X_BMPNo_CS5_1/i;
X_BMPNo_CS5_2 = X_BMPNo_CS5_2/i;
X_BMPNo_CS6_1 = ConvCS6{listType(1),1};X_BMPNo_CS6_2 = ConvCS6{listType(1),2};X_BMPNo_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_BMPNo_CS6_1 = X_BMPNo_CS6_1 + ConvCS6{listType(i),1};
X_BMPNo_CS6_2 = X_BMPNo_CS6_2 + ConvCS6{listType(i),2};
X_BMPNo_CS6_3 = X_BMPNo_CS6_3 + ConvCS6{listType(i),3};
end
X_BMPNo_CS6_1 = X_BMPNo_CS6_1/i;
X_BMPNo_CS6_2 = X_BMPNo_CS6_2/i;
X_BMPNo_CS6_3 = X_BMPNo_CS6_3/i;
X_BMPNo_CS7_1 = ConvCS7{listType(1),1};X_BMPNo_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_BMPNo_CS7_1 = X_BMPNo_CS7_1 + ConvCS7{listType(i),1};
X_BMPNo_CS7_2 = X_BMPNo_CS7_2 + ConvCS7{listType(i),2};
end
X_BMPNo_CS7_1 = X_BMPNo_CS7_1/i;
X_BMPNo_CS7_2 = X_BMPNo_CS7_2/i;

save('X_BMPNo_CS6_1_k=4.mat','X_BMPNo_CS6_1')
save('X_BMPNo_CS6_2_k=4.mat','X_BMPNo_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_BMPNo_CS5_2' - X_AmB_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_BMPNo_CS6_1' - X_AmB_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_BMPNo_CS6_2' - X_AmB_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_BMPNo_CS6_3' - X_AmB_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_BMPNo_CS7_2' - X_AmB_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_BMPNo_CS7_1' - X_AmB_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/BMP_vs_AmB_All.png'])



listType = find(strcmp(Idents2.textdata(2:end,2),"ActA_MEF")==1);
X_ActA_CS5_1 = ConvCS5{listType(1),1};X_ActA_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_ActA_CS5_1 = X_ActA_CS5_1 + ConvCS5{listType(i),1};
X_ActA_CS5_2 = X_ActA_CS5_2 + ConvCS5{listType(i),2};
end
X_ActA_CS5_1 = X_ActA_CS5_1/i;
X_ActA_CS5_2 = X_ActA_CS5_2/i;
X_ActA_CS6_1 = ConvCS6{listType(1),1};X_ActA_CS6_2 = ConvCS6{listType(1),2};X_ActA_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_ActA_CS6_1 = X_ActA_CS6_1 + ConvCS6{listType(i),1};
X_ActA_CS6_2 = X_ActA_CS6_2 + ConvCS6{listType(i),2};
X_ActA_CS6_3 = X_ActA_CS6_3 + ConvCS6{listType(i),3};
end
X_ActA_CS6_1 = X_ActA_CS6_1/i;
X_ActA_CS6_2 = X_ActA_CS6_2/i;
X_ActA_CS6_3 = X_ActA_CS6_3/i;
X_ActA_CS7_1 = ConvCS7{listType(1),1};X_ActA_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_ActA_CS7_1 = X_ActA_CS7_1 + ConvCS7{listType(i),1};
X_ActA_CS7_2 = X_ActA_CS7_2 + ConvCS7{listType(i),2};
end
X_ActA_CS7_1 = X_ActA_CS7_1/i;
X_ActA_CS7_2 = X_ActA_CS7_2/i;
save('X_ActA_CS6_1_k=4.mat','X_ActA_CS6_1')
save('X_ActA_CS6_2_k=4.mat','X_ActA_CS6_2')



%Emb lineages

%Emb lineages
listType = find(strcmp(Idents2.textdata(2:end,2),"BMP_MEF")==1);
X_BMP_CS5_1 = ConvCS5{listType(1),1};X_BMP_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_BMP_CS5_1 = X_BMP_CS5_1 + ConvCS5{listType(i),1};
X_BMP_CS5_2 = X_BMP_CS5_2 + ConvCS5{listType(i),2};
end
X_BMP_CS5_1 = X_BMP_CS5_1/i;
X_BMP_CS5_2 = X_BMP_CS5_2/i;
X_BMP_CS6_1 = ConvCS6{listType(1),1};X_BMP_CS6_2 = ConvCS6{listType(1),2};X_BMP_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_BMP_CS6_1 = X_BMP_CS6_1 + ConvCS6{listType(i),1};
X_BMP_CS6_2 = X_BMP_CS6_2 + ConvCS6{listType(i),2};
X_BMP_CS6_3 = X_BMP_CS6_3 + ConvCS6{listType(i),3};
end
X_BMP_CS6_1 = X_BMP_CS6_1/i;
X_BMP_CS6_2 = X_BMP_CS6_2/i;
X_BMP_CS6_3 = X_BMP_CS6_3/i;
X_BMP_CS7_1 = ConvCS7{listType(1),1};X_BMP_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_BMP_CS7_1 = X_BMP_CS7_1 + ConvCS7{listType(i),1};
X_BMP_CS7_2 = X_BMP_CS7_2 + ConvCS7{listType(i),2};
end
X_BMP_CS7_1 = X_BMP_CS7_1/i;
X_BMP_CS7_2 = X_BMP_CS7_2/i;
clf
save('X_BMP_CS6_1_k=4.mat','X_BMP_CS6_1')
save('X_BMP_CS6_2_k=4.mat','X_BMP_CS6_2')
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_BMP_CS5_2' - X_EmDisc_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_BMP_CS6_1' - X_EmDisc_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_BMP_CS6_2' - X_EmDisc_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_BMP_CS6_3' - X_EmDisc_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_BMP_CS7_2' - X_EmDisc_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_BMP_CS7_1' - X_EmDisc_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/BMP_vs_EmDisc_All.png'])




listType = find(strcmp(Idents2.textdata(2:end,2),"SB43_MEF")==1);
X_SB43_CS5_1 = ConvCS5{listType(1),1};X_SB43_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_SB43_CS5_1 = X_SB43_CS5_1 + ConvCS5{listType(i),1};
X_SB43_CS5_2 = X_SB43_CS5_2 + ConvCS5{listType(i),2};
end
X_SB43_CS5_1 = X_SB43_CS5_1/i;
X_SB43_CS5_2 = X_SB43_CS5_2/i;
X_SB43_CS6_1 = ConvCS6{listType(1),1};X_SB43_CS6_2 = ConvCS6{listType(1),2};X_SB43_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_SB43_CS6_1 = X_SB43_CS6_1 + ConvCS6{listType(i),1};
X_SB43_CS6_2 = X_SB43_CS6_2 + ConvCS6{listType(i),2};
X_SB43_CS6_3 = X_SB43_CS6_3 + ConvCS6{listType(i),3};
end
X_SB43_CS6_1 = X_SB43_CS6_1/i;
X_SB43_CS6_2 = X_SB43_CS6_2/i;
X_SB43_CS6_3 = X_SB43_CS6_3/i;
X_SB43_CS7_1 = ConvCS7{listType(1),1};X_SB43_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_SB43_CS7_1 = X_SB43_CS7_1 + ConvCS7{listType(i),1};
X_SB43_CS7_2 = X_SB43_CS7_2 + ConvCS7{listType(i),2};
end
X_SB43_CS7_1 = X_SB43_CS7_1/i;
X_SB43_CS7_2 = X_SB43_CS7_2/i;
save('X_SB43_CS6_1_k=4.mat','X_SB43_CS6_1')
save('X_SB43_CS6_2_k=4.mat','X_SB43_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_SB43_CS5_2' - X_EmDisc_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_SB43_CS6_1' - X_EmDisc_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_SB43_CS6_2' - X_EmDisc_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_SB43_CS6_3' - X_EmDisc_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_SB43_CS7_2' - X_EmDisc_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_SB43_CS7_1' - X_EmDisc_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/SB43_vs_EmDisc_All.png'])



listType = find(strcmp(Idents2.textdata(2:end,2),"CHIR_MEF")==1);
X_CHIR_CS5_1 = ConvCS5{listType(1),1};X_CHIR_CS5_2 = ConvCS5{listType(1),2};
for i = 2:length(listType)
X_CHIR_CS5_1 = X_CHIR_CS5_1 + ConvCS5{listType(i),1};
X_CHIR_CS5_2 = X_CHIR_CS5_2 + ConvCS5{listType(i),2};
end
X_CHIR_CS5_1 = X_CHIR_CS5_1/i;
X_CHIR_CS5_2 = X_CHIR_CS5_2/i;
X_CHIR_CS6_1 = ConvCS6{listType(1),1};X_CHIR_CS6_2 = ConvCS6{listType(1),2};X_CHIR_CS6_3 = ConvCS6{listType(1),3};
for i = 2:length(listType)
X_CHIR_CS6_1 = X_CHIR_CS6_1 + ConvCS6{listType(i),1};
X_CHIR_CS6_2 = X_CHIR_CS6_2 + ConvCS6{listType(i),2};
X_CHIR_CS6_3 = X_CHIR_CS6_3 + ConvCS6{listType(i),3};
end
X_CHIR_CS6_1 = X_CHIR_CS6_1/i;
X_CHIR_CS6_2 = X_CHIR_CS6_2/i;
X_CHIR_CS6_3 = X_CHIR_CS6_3/i;
X_CHIR_CS7_1 = ConvCS7{listType(1),1};X_CHIR_CS7_2 = ConvCS7{listType(1),2};
for i = 2:length(listType)
X_CHIR_CS7_1 = X_CHIR_CS7_1 + ConvCS7{listType(i),1};
X_CHIR_CS7_2 = X_CHIR_CS7_2 + ConvCS7{listType(i),2};
end
X_CHIR_CS7_1 = X_CHIR_CS7_1/i;
X_CHIR_CS7_2 = X_CHIR_CS7_2/i;
save('X_CHIR_CS6_1_k=4.mat','X_CHIR_CS6_1')
save('X_CHIR_CS6_2_k=4.mat','X_CHIR_CS6_2')

clf
subplot(1,3,1);
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',X_CHIR_CS5_2' - X_EmDisc_CS5_2','FaceColor','interp','LineStyle','none');
view([-59.3039,41.1722])
%colorbar
set(gca,'clim',[-0.1 .1] )
axis equal
axis off
camlight('headlight')
subplot(1,3,2);
f1 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',X_CHIR_CS6_1' - X_EmDisc_CS6_1','FaceColor','interp','LineStyle','none');
hold on
f2 = patch('Faces',OBJ2B.objects(8).data.vertices,'Vertices',OBJ2B.vertices,'FaceVertexCData',X_CHIR_CS6_2' - X_EmDisc_CS6_2','FaceColor','interp','LineStyle','none');
f7 = patch('Faces',OBJ2E.objects(12).data.vertices,'Vertices',OBJ2E.vertices,'FaceVertexCData',X_CHIR_CS6_3' - X_EmDisc_CS6_3','FaceColor','interp','LineStyle','none');
view([ 165.7876,21.7629])
set(gca,'clim',[-0.1 .1] )
%colorbar
axis equal
axis off
camlight('headlight')
subplot(1,3,3);
f1 = patch('Faces',OBJ3.objects(20).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_CHIR_CS7_2' - X_EmDisc_CS7_2','FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ3.objects(16).data.vertices,'Vertices',OBJ3.vertices,'FaceVertexCData',X_CHIR_CS7_1' - X_EmDisc_CS7_1','FaceColor','interp','LineStyle','none');
%colorbar
set(gca,'clim',[-0.1 .1] )
view([24.5784,1.2879])
axis equal
axis off
camlight('headlight')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 30 10])
print('-dpng',['/Users/christopherpenfold/Desktop/Thorsten/FINAL/Clara_Amnioids_Figuers/3D/CHIR_vs_EmDisc_All.png'])


return

D1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/ESCs_3DPCAcoordinates.csv');
 Distancematrix = dist(D1.data');
 IDs = D1.textdata(2:end,1:3);
 esc = find(strcmp(IDs(:,2),'ESC_primed')==1);
 naive = find(strcmp(IDs(:,2),'ESC_naive')==1);

  n1 = find(strcmp(IDs(:,2),'EmDisc_CS5')==1);
  n2 = find(strcmp(IDs(:,2),'EmDisc_CS6')==1);


  
 for i = 1:length(CellIDs)
     IndsMatch(i,1) = find(strcmp( ['X' CellIDs{i}],D1.textdata(2:end,1))==1);
 end
 Dhat1 = ( max(max(Distancematrix)) - mean(Distancematrix(IndsMatch,esc),2) ) ./max(max(Distancematrix));
 Dhat2 = ( max(max(Distancematrix)) - mean(Distancematrix(IndsMatch,naive),2) ) ./max(max(Distancematrix));
 Dhat3 = ( max(max(Distancematrix)) - min(Distancematrix(IndsMatch,esc)')' ) ./max(max(Distancematrix));
 Dhat4 = ( max(max(Distancematrix)) - min(Distancematrix(IndsMatch,naive)')' ) ./max(max(Distancematrix)); 
 Dhat5 = ( max(max(Distancematrix)) - Distancematrix(IndsMatch,esc) ) ./max(max(Distancematrix));
 Dhat6 = ( max(max(Distancematrix)) - Distancematrix(IndsMatch,naive) ) ./max(max(Distancematrix)); 
 

% D1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/ESCs_3DPCAcoordinates.csv');
% Distancematrix = dist(D1.data');
% IDs = D1.textdata(2:end,1:3);
% 
% esc = find(strcmp(IDs(:,2),'ESC_primed')==1);
% naive = find(strcmp(IDs(:,2),'ESC_naive')==1);
% 
% INVIVO = setdiff(1:1:length(IDs),[esc;naive]);
% 
% for i = 1:length(IDs)    
%     type = IDs(i,2);
%     %indexclude = find(strcmp(IDs(:,2),type)==1);    
%     %type = IDs(find(Distancematrix(i,:)<15),2);
%     Closest = find(Distancematrix(i,:)<10);
%     N1(i,1) = length(intersect(Closest,esc));    
%     N2(i,1) = length(intersect(Closest,naive));
%     N3(i,1) = length(Closest);    
% 
%     Closest = find(Distancematrix(i,:)<12);
%     N1(i,2) = length(intersect(Closest,esc));    
%     N2(i,2) = length(intersect(Closest,naive));
%     N3(i,2) = length(Closest);    
%   
%     Closest = find(Distancematrix(i,:)<15);
%     N1(i,3) = length(intersect(Closest,esc));   
%     N2(i,3) = length(intersect(Closest,naive));
%     N3(i,3) = length(Closest);    
%     
%     Closest = find(Distancematrix(i,:)<18);
%     N1(i,4) = length(intersect(Closest,esc));   
%     N2(i,4) = length(intersect(Closest,naive));
%     N3(i,4) = length(Closest);    
% 
%     Closest = find(Distancematrix(i,:)<20);
%     N1(i,5) = length(intersect(Closest,esc));   
%     N2(i,5) = length(intersect(Closest,naive));
%     N3(i,5) = length(Closest);        
% end
% 
% %imagesc([N1./N3,N2./N3])
% %set(gca,'YTick',1:1:size(IDs,1),'YTickLabel',IDs(:,2))
% 
% 
% %DMAP1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/CCAFigForESCs/DimRed/ESC_all.csv')
% %targ = DMAP1.textdata(1,2:end);
% %targ = strrep(targ,'"','');
% %targ = strrep(targ,'RNA.','');
% 
% %estarg = DMAP1.textdata(2:end,1);
% %estarg = strrep(estarg,'RNA.','');
% 
% %Targ = targ;
% %Targ = strrep(Targ,'EmDisc_CS5_','');
% %Targ = strrep(Targ,'ExMes_CS5_','');
% %Targ = strrep(Targ,'VE_CS5_','');
% %Targ = strrep(Targ,'Tb_CS5_','');
% %Targ = strrep(Targ,'Am_CS5_','');
% for i = 1:length(CellIDs)
%     try
%     P1(i,:) = N1(find(strcmp(IDs(:,1),[CellIDs{i}])==1),:)./N3(find(strcmp(IDs(:,1),[CellIDs{i}])==1),:);
%     P2(i,:) = N2(find(strcmp(IDs(:,1),[CellIDs{i}])==1),:)./N3(find(strcmp(IDs(:,1),[CellIDs{i}])==1),:);
%     catch
%             P1(i,:) = N1(find(strcmp(IDs(:,1),['X',CellIDs{i}])==1),:)./N3(find(strcmp(IDs(:,1),['X',CellIDs{i}])==1),:);
%     P2(i,:) = N2(find(strcmp(IDs(:,1),['X',CellIDs{i}])==1),:)./N3(find(strcmp(IDs(:,1),['X',CellIDs{i}])==1),:);
%     end
% end    
%     

set(0, 'DefaultTextInterpreter', 'none')
%for i = 1:size(P1,1)

    
%     Ytrain2 = P2(find(LabBin==1),3);
%     Ytrain = P1(find(LabBin==1),3);
%     %Ytrain = (Pt(find(isnan(Pt(:,i))~=1),i) );
%     Xtrain = Output.Xtrain(find(LabBin==1),:);
 % 
 
 
  pg1 = {@priorGauss,0,2};  
  pg2 = {@priorGauss,0,2};  
  pg3 = {@priorGauss,log(0.5),2};  
  pc = {@priorClamped};
  prior1.mean = {[]};  
  prior1.cov  = {pg1;pg2;pg2;pg2};
  prior1.lik  = {pg3};
  im1 = {@infPrior,@infExact,prior1};
%  pg1 = {@priorGauss,0,1};  pg2 = {@priorGauss,0,1};  pg3 = {@priorGauss,log(0.5),1};  pc = {@priorClamped};
%  prior1.mean = {[]}; prior1.cov  = {pg1;pg2}; prior1.lik  = {pg3};
%  im1 = {@infPrior,@infExact,prior1};                % inference method
%  
%  hyp1.cov  = [var(Ytrain); log(2)];
%  hyp1.lik  = [var(Ytrain)/2];    
%  hyp1.mean = mean(Ytrain);       
% 
% par  = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain};
% hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
% 
% %[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain, Xtrain(:,1:3));
% [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain, Output.Xtest(Output.ind1,:));
%  hyp1.cov  = [var(Ytrain2); log(2)];
%  hyp1.lik  = [var(Ytrain2)/2];    
%  hyp1.mean = mean(Ytrain2);       
% 
% [m_2 s_2] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain2, Output.Xtest(Output.ind1,:));


%scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,[0,0,0],'fill')
%hold on
%scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill')


%h = figure(1)
% subplot(1,2,1);
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% hold on
% scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill')
% %t%itle([estarg{i}])
% view([ 48, 18])
% alpha = 0.05;
% %set(h1,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
% 
% yl = ylim;
% xl = xlim;
% zl = zlim;
% set(gca,'clim',[0 0.2] )
% 
% set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% set(gca,'fontsize', 12)
% set(gca,'TickLength',[0 0])
% cb=colorbar;
% %cb.Position = cb.Position + [0.1 -0.02 0 0]
% %cy=get(cb,'YTick');
% %set(cb,'YTick',[0 0.1)])
% set(gca,'fontsize', 12)
% %axis equal
% 
% subplot(1,2,2);
% %h = figure(2)
% scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_2,'fill')
% hold on
% scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain2,'fill')
% %title([estarg{i}])
% view([ 48, 18])
% alpha = 0.05;
% %set(h1,'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
% 
% yl = ylim;
% xl = xlim;
% zl = zlim;
% set(gca,'clim',[0 0.2] )
% 
% set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% set(gca,'fontsize', 12)
% set(gca,'TickLength',[0 0])
% cb=colorbar;
% %cb.Position = cb.Position + [0.1 -0.02 0 0]
% %cy=get(cb,'YTick');
% %set(cb,'YTick',[0 0.1)])
% set(gca,'fontsize', 12)
% axis equal
% 
% clf

% 
% 
% Dcorr = importdata("~/Desktop/Thorsten/FINAL/CCAFigForPaper_wCS6_v2//DimRed/ESC_corr_wCS6.csv")
% 
% 
% targs = Dcorr.textdata(1,2:end);
% targs = strrep(targs,'"','');
% targs = strrep(targs,'RNA.','');
% escs = Dcorr.textdata(2:end,1);
% 
% CoC = zeros(length(escs),length(CellIDs4));
% for i = 1:length(CellIDs4)
%     idnss = find(contains(targs,CellIDs4{i}));
%     try
%     CoC(:,i) = Dcorr.data(:,idnss);
%     catch
%     CoC(:,1) =0;
%     end
% end
% 
% 
% 
% %   
% % 
% % for i = 1:length(escs)
% %     
% %     train = CoC(i,find(LabBin==1))';
% %     Xtrain = Output.Xtrain(find(LabBin==1),:);
% %     
% %     hyp1.cov  = [var(Ytrain); log(2)]; hyp1.lik  = [var(Ytrain)/2];  hyp1.mean = mean(Ytrain);       
% % 
% %     par  = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain};
% %     hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
% % 
% %     [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain, Output.Xtest(Output.ind1,:));
% %     
% %     
% %     h = figure(1)
% %     scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% %     hold on
% %     scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill','MarkerEdgeColor','k')
% %     title([escs{i}])
% %     view([ -78.8,14.9])
% %     alpha = 0.05;
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.2 0.9] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
% % 
% %    %pause 
% %     print(h,['~/Desktop/Sections/CS6_' escs{i} '.pdf'],'-dpdf','-painters');
% %     clf
% % end
% 
% %print(h, '-dpdf', ['~/Desktop/' estarg{i} '.pdf'],'-r500','-painters')
% %clf
% %end
% 
% %plot(m_1,Ytrain,'o')
% 
% %[L1 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
% %[L2 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
% 
% %for i = 1:25
% %subplot(5,5,i);scatter3(Output.Xtrain(:,1),Output.Xtrain(:,2),Output.Xtrain(:,3),50,Pt(:,i),'fill')
% %
% %end
% % h = figure(5)
% % scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% % hold on
% % scatter3(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3),95,'fill')
% % text(Output.Xtrain(find(Output.LabBin==5),1),Output.Xtrain(find(Output.LabBin==5),2),Output.Xtrain(find(Output.LabBin==5),3), CellIDs4(find(Output.LabBin==5)));
% % savefig(h,['~/Desktop/FIG/CS5_Tb'])
% % clf
% % 
% 
% %escs(1:45)
% %escs(46:61)
% 
% 
%   Ytrain1 = reshape(CoC([1:45],find(LabBin==1))',45*46,1);
%   Ytrain2 = reshape(CoC([46:61],find(LabBin==1))',16*46,1);  
%   
%   Xtrain1 = repmat(Output.Xtrain(find(LabBin==1),:),45,1);
%   Xtrain2 = repmat(Output.Xtrain(find(LabBin==1),:),16,1);   
%       
%   hyp1.cov  = [var(Ytrain1); log(2)]; hyp1.lik  = [var(Ytrain1)/2];  hyp1.mean = mean(Ytrain1);       
%   par  = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain1};
%   hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
%   [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain1(:,1:3), Ytrain1, [OBJ2.vertices]/400);
%         
%   Output.Conv = m_1;
% 
% %    h = figure(2)
% %    subplot(2,3,1);
% %    scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% %    hold on
% %    %scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill','MarkerEdgeColor','k')
% %     title(['Conv'])
% %     view([-121.7479,14.1676])
% %     alpha = 0.05;
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.4 0.7] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
% %     subplot(2,3,4);
% %     scatter3(Xtrain1(:,1)+randn(2070,1)*0.05,Xtrain1(:,2)+randn(2070,1)*0.05,Xtrain1(:,3)+randn(2070,1)*0.05,250,Ytrain1,'fill','MarkerEdgeColor','k')
% %     view([-121.7479,14.1676])
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.4 0.7] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
%    %pause 
%     %print(h,['~/Desktop/Sections/CS6_' escs{i} '.pdf'],'-dpdf','-painters');
%     %clf
%     
%     hyp1.cov  = [var(Ytrain2); log(2)]; hyp1.lik  = [var(Ytrain2)/2];  hyp1.mean = mean(Ytrain2);       
%     par  = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain2};
%     hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
%     [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain2(:,1:3), Ytrain2, [OBJ2.vertices]/400);
%     
%     
%     %h = figure(1)
% %     subplot(2,3,2);
% %     scatter3(Output.Xtest(Output.ind1,1),Output.Xtest(Output.ind1,2),Output.Xtest(Output.ind1,3),5,m_1,'fill')
% %     hold on
% %     %scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),250,Ytrain,'fill','MarkerEdgeColor','k')
% %     title(['Naive'])
% %     view([-121.7479,14.1676])
% %     alpha = 0.05;
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.4 0.7] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
% % 
% %     
% %     subplot(2,3,5); scatter3(Xtrain2(:,1)+randn(736,1)*0.05,Xtrain2(:,2)+randn(736,1)*0.05,Xtrain2(:,3)+randn(736,1)*0.05,250,Ytrain2,'fill','MarkerEdgeColor','k')
% %     view([-121.7479,14.1676])
% %     yl = ylim;
% %     xl = xlim;
% %     zl = zlim;
% %     set(gca,'clim',[0.4 0.7] )
% %     set(gca,'ZTick',linspace(zl(1),zl(2),5),'YTick',linspace(yl(1),yl(2),5),'XTick',linspace(xl(1),xl(2),5),'xticklabel',num2str(''),'yticklabel',num2str(''),'zticklabel',num2str(''))
% %     set(gca,'fontsize', 12)
% %     set(gca,'TickLength',[0 0])
% %     cb=colorbar;
% %     set(gca,'fontsize', 12)
%     
% Output.Naive = m_1;
% 
%     %print(h,['~/Desktop/Sections/CS6_' escs{i} '.pdf'],'-dpdf','-painters');
%     %clf    




set(0, 'DefaultTextInterpreter', 'none')

%pg1 = {@priorGauss,0,2};  
%pg2 = {@priorGauss,0,2};  
%pg3 = {@priorGauss,log(0.5),2};  
%pc = {@priorClamped};
%prior1.mean = {[]};  
%prior1.cov  = {pg1;pg2;pg2;pg2};
%prior1.lik  = {pg3};
%im1 = {@infPrior,@infExact,prior1};
pg1 = {@priorGauss,0,1};  pg2 = {@priorGauss,0,1};  pg3 = {@priorGauss,log(0.5),1};  pc = {@priorClamped};
prior1.mean = {[]}; prior1.cov  = {pg1;pg2}; prior1.lik  = {pg3};
im1 = {@infPrior,@infExact,prior1};                % inference method

Dcorr = importdata("~/Desktop/Thorsten/FINAL/CCAFigForPaper_wCS6_v2//DimRed/ESC_corr_wCS6.csv");


targs = Dcorr.textdata(1,2:end);
targs = strrep(targs,'"','');
targs = strrep(targs,'RNA.','');
escs = Dcorr.textdata(2:end,1);

CoC = zeros(length(escs),length(CellIDs4));
for i = 1:length(CellIDs4)
    idnss = find(contains(targs,CellIDs4{i}));
    try
    CoC(:,i) = Dcorr.data(:,idnss);
    catch
    CoC(:,1) =0;
    end
end

Dhat5 = Dhat5';
Dhat6 = Dhat6';

for i = 1:size(CoC)
    weightedDist(i,1:3) = sum(repmat(CoC(i,:)',1,3).*Xtrain,1);
end


%keyboard
uLab = unique(LabBin);
for i = 1:length(uLab)
    
    
  Ytrain1 = reshape(CoC([1:45],find(LabBin==uLab(i)))',45*length(find(LabBin==uLab(i))),1);
  Ytrain2 = reshape(CoC([46:61],find(LabBin==uLab(i)))',16*length(find(LabBin==uLab(i))),1);  
  
  Xtrain1 = repmat(Output.Xtrain(find(LabBin==uLab(i)),:),45,1);
  Xtrain1 = Xtrain1 + randn(size(Xtrain1,1),3)*0.05;
  Xtrain2 = repmat(Output.Xtrain(find(LabBin==uLab(i)),:),16,1);
  Xtrain2 = Xtrain2 + randn(size(Xtrain2,1),3)*0.05;   
      
  hyp1.cov  = [log(2); var(Ytrain1)]; hyp1.lik  = [var(Ytrain2)/2];  hyp1.mean = mean(Ytrain1);       
  par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain1};
  hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
  [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain1(:,1:3), Ytrain1, [OBJ2.vertices]/400);
        
       
  H{i,1} = hyp2;
  Output.Conv{i,1} = m_1;

    hyp1.cov  = [log(2); var(Ytrain2)]; hyp1.lik  = [var(Ytrain2)/2];  hyp1.mean = mean(Ytrain2);       
    par  = {@meanConst,@covSEiso,'likGauss',Xtrain2(:,1:3), Ytrain2};
    hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
    [m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain2(:,1:3), Ytrain2, [OBJ2.vertices]/400);
    
    
    Output.Naive{i,1} = m_1;
    
      H{i,2} = hyp2;
      
      
  
 
Ytrain3 = Dhat1(find(LabBin==uLab(i)));
Ytrain4 = Dhat2(find(LabBin==uLab(i)));
Ytrain5 = Dhat3(find(LabBin==uLab(i)));
Ytrain6 = Dhat4(find(LabBin==uLab(i)));
Xtrain3 = Output.Xtrain(find(LabBin==uLab(i)),:);
Xtrain4 = Xtrain3;
Xtrain5 = Xtrain3;
Xtrain6 = Xtrain3;

hyp1.cov  = [log(2); var(Ytrain3)]; hyp1.lik  = [var(Ytrain3)/2];  hyp1.mean = mean(Ytrain3);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain3};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain3(:,1:3), Ytrain3, [OBJ2.vertices]/400);

 Output.Conv{i,2} = m_1;
      

hyp1.cov  = [log(2); var(Ytrain4)]; hyp1.lik  = [var(Ytrain4)/2];  hyp1.mean = mean(Ytrain4);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain4};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain4(:,1:3), Ytrain4, [OBJ2.vertices]/400);

Output.Naive{i,2} = m_1;

hyp1.cov  = [log(2); var(Ytrain5)]; hyp1.lik  = [var(Ytrain5)/2];  hyp1.mean = mean(Ytrain5);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain5};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain5(:,1:3), Ytrain5, [OBJ2.vertices]/400);

Output.Conv{i,3} = m_1;


hyp1.cov  = [log(2); var(Ytrain6)]; hyp1.lik  = [var(Ytrain6)/2];  hyp1.mean = mean(Ytrain6);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain6};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain6(:,1:3), Ytrain6, [OBJ2.vertices]/400);

Output.Naive{i,3} = m_1;
    

  Ytrain7 = reshape(Dhat5(:,find(LabBin==uLab(i)))',45*length(find(LabBin==uLab(i))),1);
  Ytrain8 = reshape(Dhat6(:,find(LabBin==uLab(i)))',16*length(find(LabBin==uLab(i))),1);  

  
  
hyp1.cov  = [log(2); var(Ytrain7)]; hyp1.lik  = [var(Ytrain7)/2];  hyp1.mean = mean(Ytrain5);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain7};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain1(:,1:3), Ytrain7, [OBJ2.vertices]/400);

Output.Conv{i,4} = m_1;


hyp1.cov  = [log(2); var(Ytrain8)]; hyp1.lik  = [var(Ytrain8)/2];  hyp1.mean = mean(Ytrain8);       
par  = {@meanConst,@covSEiso,'likGauss',Xtrain1(:,1:3), Ytrain8};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain2(:,1:3), Ytrain8, [OBJ2.vertices]/400);
 
Output.Naive{i,4} = m_1;

      
      %f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',OBJ1.vertices,'FaceVertexCData',m_1,'FaceColor','interp','LineStyle','none');
%hold on

end

Output.weightedDist = weightedDist;
Output.CoC = CoC;
Output.Xtrain = Xtrain;