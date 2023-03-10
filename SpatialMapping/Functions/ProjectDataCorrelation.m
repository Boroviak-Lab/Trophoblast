function [targInd,C1] = ProjectDataCorrelation(File1,File2,File3,OutputCS6)

S1 = importdata(File1);
S2 = importdata(File2);

C1 = importdata(File3);

targs = S1.textdata(2:end,2); %C1.textdata(2:end,2);
regs = S2.textdata(2:end,3);

targs = strrep(targs,'P1_','')
targs = strrep(targs,'P2_','')

%regs = strrep(regs,'P1_','')
%regs = strrep(regs,'P2_','')


for i = 1:length( OutputCS6.cleanShot )
     OutputCS6.cleanShot{i} = ['E15C2_' OutputCS6.cleanShot{i}];
end

targInd = zeros(length(OutputCS6.cleanShot),1);


%keyboard
%X = NaN*ones(length(targs),1);
%Y = NaN*ones(length(targs),1);
%Z = NaN*ones(length(targs),1);

%CellType = cell(length(targs),1);

for i = 1:length(OutputCS6.cleanShot)
    try
        targInd(i,1) = find(strcmp(OutputCS6.cleanShot{i} , targs)==1);        
        %X(i,1) = OutputCS6.cleanX(targInd(i,1),1);
        %Y(i,1) = OutputCS6.cleanX(targInd(i,1),2);
        %Z(i,1) = OutputCS6.cleanX(targInd(i,1),3);              
        %CellType{i,1} = OutputCS6.cleanAnotaton{ targInd(i,1) };
    end
end