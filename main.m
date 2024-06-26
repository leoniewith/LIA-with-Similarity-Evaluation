clc
close all
clear
rng('shuffle')
R = [550 300 3];
allhv = [];
allreflag1 = [];
allreflag2 = [];
allreflag3 = [];
allcrflag1 = [];
allcrflag2 = [];
allcrflag3 = [];
allsiflag1 = [];
allsiflag2 = [];
allsiflag3 = [];
PortPosition = [0 50];
MaxRobotNo = 3;
global XWind YWind RWind
XWind = 59.27;
YWind = 21.39;
RWind = 15;
load('x.mat')
load('y.mat')
load('r.mat')
XVessel = x(11:20);
YVessel = y(11:20);
RVessel = r(11:20);
reflag1 = [];
reflag2 = [];
reflag3 = [];
crflag1 = [];
crflag2 = [];
crflag3 = [];
siflag1 = [];
siflag2 = [];
siflag3 = [];
DistanceWeight = zeros(1,length(XVessel));
for i = 1:length(XWind)
    distance = sqrt((XVessel-XWind(i)).^2+(YVessel-YWind(i)).^2);
    DistanceWeight(find(distance<RWind(i))) = DistanceWeight(find(distance<RWind(i))) + distance(find(distance<RWind(i)))./RWind(i);
end
GenotypeLength = [4 4 2 2];
PopulationSize = 100;
MaxGeneration = 100;
MemorySize = 5;
ArtiSize = 10;
Pc = 0.9;
Pm = 0.1;
N = PopulationSize/2+10;
VarietyNumber = length(XVessel);
FinalPop = [];
AllPop = [];
UCB_std = 0;
cu_re_1 = cell(1, sum(GenotypeLength)*VarietyNumber);
cu_re_0 = cell(1, sum(GenotypeLength)*VarietyNumber);
dis_factor = 0.9;
InitPop = randi([0 1],PopulationSize,VarietyNumber*sum(GenotypeLength)+VarietyNumber*4);
DecodePop = DecodeFunction(InitPop,XVessel,YVessel,RVessel,GenotypeLength,PopulationSize,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
pa = DecodePop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
IsDomi = IDAf(pa);
DomiIndex = find(IsDomi==1);
DomiPop = DecodePop(DomiIndex,:);
hillcoff = hillco(DomiPop(:,end-2:end)');
useindex = find(hillcoff>0);
DomiPop = DomiPop(useindex,:);
iter = 1;
while iter <= MaxGeneration
    newpop = DecodePop;
    % Actor-critic-inspired crossover
    crossPop = [];
    for i = 1:PopulationSize
        index1 = randi(PopulationSize);
        index2 = randi(PopulationSize);
        parent1 = newpop(index1,:);
        parent2 = newpop(index2,:);
        tempPa1 = parent1((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
        tempGene1 = parent1(1:sum(GenotypeLength)*VarietyNumber);
        tempPa2 = parent2((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
        tempGene2 = parent2(1:sum(GenotypeLength)*VarietyNumber);
        newIndi = [];
        for j = 1:sum(GenotypeLength)*VarietyNumber
            newGene1 = [tempGene1(1:j-1) 1-tempGene1(j) tempGene1(j+1:end)];
            decodeGene1 = DecodeFunction(newGene1,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
            decodePa1 = decodeGene1((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
            rela1 =  IDAf([tempPa1; decodePa1]);
            stp1 = find(newpop(:, j) == tempGene1(j));
            stp2 = find(newpop(:, j) == 1 - tempGene1(j));
            stp1Pa = newpop(stp1, (sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
            stp2Pa = newpop(stp2, (sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
            rela =  IDAf([stp1Pa; stp2Pa]);
            if length(stp1) ~= 0
                cureward1 = sum(rela(1:length(stp1))) / length(stp1);
            else
                cureward1 = 0;
            end
            if length(stp2) ~= 0
                cureward2 = sum(rela((length(stp1)+1):PopulationSize)) / length(stp2);
            else
                cureward2 = 0;
            end
            if i == 1
                if tempGene1(j) == 1
                    cu_re_1{j} = [cu_re_1{j}, cureward1];
                    cu_re_0{j} = [cu_re_0{j}, cureward2];
                else
                    cu_re_0{j} = [cu_re_0{j}, cureward1];
                    cu_re_1{j} = [cu_re_1{j}, cureward2];
                end
            end
            newGene2 = [tempGene2(1:j-1) 1-tempGene2(j) tempGene2(j+1:end)];
            decodeGene2 = DecodeFunction(newGene2,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
            decodePa2 = decodeGene2((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
            rela2 =  IDAf([tempPa2; decodePa2]);
            if (rela1(1) == rela2(1)) && (rela1(1) == 1) && rand < Pc
                if tempGene1(j) == tempGene2(j)
                    newIndi = [newIndi tempGene1(j)];
                else
                    cuallreward1 = 0;
                    cuallreward0 = 0;
                    for ir = 1 : (length(cu_re_1{j}) - 1)
                        cuallreward1 = cuallreward1 + cu_re_1{j}(ir) * dis_factor ^ (length(cu_re_1{j}) - ir);
                        cuallreward0 = cuallreward0 + cu_re_0{j}(ir) * dis_factor ^ (length(cu_re_0{j}) - ir);
                    end
                    if tempGene1(j) == 1
                        cureward1 = cu_re_1{j}(end);
                        cu_re1 = cuallreward1;
                        cureward2 = cu_re_0{j}(end);
                        cu_re2 = cuallreward0;
                    else
                        cureward1 = cu_re_0{j}(end);
                        cu_re1 = cuallreward0;
                        cureward2 = cu_re_1{j}(end);
                        cu_re2 = cuallreward1;
                    end
                    if cureward1 + cu_re1 / iter / MaxGeneration > cureward2 + cu_re2 / iter / MaxGeneration
                        newIndi = [newIndi tempGene1(j)];
                    else
                        newIndi = [newIndi tempGene2(j)];
                    end
                end
            elseif rela1(1) == 1 && rand < Pc
                newIndi = [newIndi tempGene1(j)];
            elseif rela2(1) == 1 && rand < Pc
                newIndi = [newIndi tempGene2(j)];
            else
                rela = IDAf([tempPa1; tempPa2]);
                if rela(2) == 0
                    newIndi = [newIndi tempGene1(j)];
                else
                    newIndi = [newIndi tempGene2(j)];
                end
            end
        end
        crossPop = [crossPop; newIndi];
    end
    DecodePop = DecodeFunction(crossPop,XVessel,YVessel,RVessel,GenotypeLength,PopulationSize,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
    % Mutation process regulation
    renumber1 = 0;
    renumber2 = 0;
    renumber3 = 0;
    repath = [];
    nonrepath = [];
    nonreno = [];
    crnumber1 = 0;
    crnumber2 = 0;
    crnumber3 = 0;
    crpath = [];
    noncrpath = [];
    noncrno = [];
    sinumber1 = 0;
    sinumber2 = 0;
    sinumber3 = 0;
    sipath = [];
    nonsipath = [];
    nonsino = [];
    for itry = 1:PopulationSize
        [flag,positionindex,inter_index] = trytest2(DecodePop(itry,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
        renumber1 = renumber1+flag(1);
        if flag(1) ~= 0
            temp1 = MultiMutationFunction1(positionindex,VarietyNumber,GenotypeLength,1,DecodePop(itry,1:sum(GenotypeLength)*length(XVessel)));
            temp1(:,sum(GenotypeLength)*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3) = 0;
            decodetemp = DecodeFunction(temp1,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
            [flag,positionindex,~] = trytest2(decodetemp(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
            renumber2 = renumber2+flag(1);
            repath = [repath;DecodePop(itry,:)];
        else
            nonrepath = [nonrepath;DecodePop(itry,:)];
            nonreno = [nonreno itry];
        end
        if flag(2) ~= 0
            crnumber1 = crnumber1+1;
            temp1 = MultiMutationFunction2(VarietyNumber,GenotypeLength,1,DecodePop(itry,1:sum(GenotypeLength)*length(XVessel)));
            temp1(:,sum(GenotypeLength)*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3) = 0;
            decodetemp = DecodeFunction(temp1,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
            [flag,positionindex,~] = trytest2(decodetemp(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
            if flag(2) ~= 0
                crnumber2 = crnumber2+1;
            end
            crpath = [crpath;DecodePop(itry,:)];
        else
            noncrpath = [noncrpath;DecodePop(itry,:)];
            noncrno = [noncrno itry];
        end
        if inter_index == 1
            sinumber1 = sinumber1 + inter_index;
            temp1 = MultiMutationFunction1(positionindex,VarietyNumber,GenotypeLength,1,DecodePop(itry,1:sum(GenotypeLength)*length(XVessel)));
            temp1(:,sum(GenotypeLength)*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3) = 0;
            decodetemp = DecodeFunction(temp1,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
            [flag,positionindex,inter_index] = trytest2(decodetemp(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
            sinumber2 = sinumber2 + inter_index;
            sipath = [sipath;DecodePop(itry,:)];
        else
            nonsipath = [nonsipath;DecodePop(itry,:)];
            nonsino = [nonsino itry];
        end
    end
    newrepath = [];
    for ttry = 1:size(repath,1)
        redun = repath(ttry,:);
        for t = 1:sum(GenotypeLength)*VarietyNumber
            gene = redun(t);
            mugene = 1-redun(t);
            pb1 = size(nonrepath,1)/PopulationSize;
            if isempty(nonrepath)
                pab1 = 0;
            else
                pab1 = length(find(nonrepath(:,t)==mugene))/(size(nonrepath,1)+1);
            end
            pb2 = size(repath,1)/PopulationSize;
            pab2 = length(find(repath(:,t)==mugene))/(size(repath,1)+1);
            pb1a =  (pb1*pab1)/(pb1*pab1+pb2*pab2);
            if rand <= pb1a
                redun(t) = mugene;
            end
        end
        redunnew = DecodeFunction(redun,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
        [flagnew,positionindex,~] = trytest2(redunnew(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
        renumber3 = renumber3+flagnew(1);
        newrepath = [newrepath;redunnew];
    end
    newcrpath = [];
    for ttry = 1:size(crpath,1)
        crdun = crpath(ttry,:);
        for t = 1:sum(GenotypeLength)*VarietyNumber
            gene = crdun(t);
            mugene = 1-crdun(t);
            pb1 = size(noncrpath,1)/PopulationSize;
            if isempty(noncrpath)
                pab1 = 0;
            else
                pab1 = length(find(noncrpath(:,t)==mugene))/(size(noncrpath,1)+1);
            end
            pb2 = size(crpath,1)/PopulationSize;
            pab2 = length(find(crpath(:,t)==mugene))/(size(crpath,1)+1);
            pb1a =  (pb1*pab1)/(pb1*pab1+pb2*pab2);
            if rand <= pb1a
                crdun(t) = mugene;
            end
        end
        crdunnew = DecodeFunction(crdun,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
        [flagnew,positionindex,~] = trytest2(crdunnew(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
        if flagnew(2) ~= 0
            crnumber3 = crnumber3+1;
        end
        newcrpath = [newcrpath;crdunnew];
    end
    newsipath = [];
    for ttry = 1:size(sipath,1)
        sidun = sipath(ttry,:);
        for t = 1:sum(GenotypeLength)*VarietyNumber
            gene = sidun(t);
            mugene = 1-sidun(t);
            pb1 = size(nonsipath,1)/PopulationSize;
            if isempty(nonsipath)
                pab1 = 0;
            else
                pab1 = length(find(nonsipath(:,t)==mugene))/(size(nonsipath,1)+1);
            end
            pb2 = size(sipath,1)/PopulationSize;
            pab2 = length(find(sipath(:,t)==mugene))/(size(sipath,1)+1);
            pb1a =  (pb1*pab1)/(pb1*pab1+pb2*pab2);
            if rand <= pb1a
                sidun(t) = mugene;
            end
        end
        sidunnew = DecodeFunction(sidun,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
        [flagnew,positionindex,inter_index] = trytest2(sidunnew(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
        sinumber3 = sinumber3 + inter_index;
        newsipath = [newsipath;sidunnew];
    end
    reno = setdiff([1:PopulationSize],nonreno);
    crno = setdiff([1:PopulationSize],noncrno);
    sino = setdiff([1:PopulationSize],nonsino);
    for ttry = 1:PopulationSize
        cadi_replace = [];
        if ismember(ttry, reno) == 1
            cadi_replace = [cadi_replace; newrepath(find(reno == ttry), :)];
        end
        if ismember(ttry, crno) == 1
            cadi_replace = [cadi_replace; newcrpath(find(crno == ttry), :)];
        end
        if ismember(ttry, sino) == 1
            cadi_replace = [cadi_replace; newsipath(find(sino == ttry), :)];
        end
        size_cadi = size(cadi_replace);
        if size_cadi(1) ~= 0
            DecodePop(ttry,:) = cadi_replace(randi([1,size_cadi(1)]), :);
        end
    end    
    reflag1 = [reflag1 renumber1];
    reflag2 = [reflag2 renumber2];
    reflag3 = [reflag3 renumber3];
    crflag1 = [crflag1 crnumber1];
    crflag2 = [crflag2 crnumber2];
    crflag3 = [crflag3 crnumber3];
    siflag1 = [siflag1 sinumber1];
    siflag2 = [siflag2 sinumber2];
    siflag3 = [siflag3 sinumber3];
    DecodePop = DecodeFunction(DecodePop,XVessel,YVessel,RVessel,GenotypeLength,PopulationSize,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
    % Evaluation
    pa = DecodePop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
    % Selection process regulation of memory cell
    IsDomi = IDAf(pa);
    DomiIndex = find(IsDomi==1);
    DomiPop1=[DomiPop;DecodePop(DomiIndex,:)];
    pa = DomiPop1(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
    IsDomi = IDAf(pa);
    DomiIndex = find(IsDomi==1);
    DomiPop2=DomiPop1(DomiIndex,:);
    hillcoff = hillco(DomiPop2(:,end-2:end)');
    useindex = find(hillcoff>0);
    DomiPop = DomiPop2(useindex,:);
    tempsize = size(DomiPop);
    if tempsize(1) == 0
        stop
    end
    FinalPop = [FinalPop;DomiPop];
    temphv = hypervolume(DomiPop(:,161:163), R, 100000);
    allhv = [allhv temphv];
    AllPop = [AllPop;DecodePop];
    randorder = randperm(PopulationSize);
    DecodePop = [DomiPop(1:min(size(DomiPop,1),MemorySize),:);DecodePop(randorder(1:PopulationSize-min(size(DomiPop,1),MemorySize)),:)];
    % Similarity evaluation
    std_reward_0 = zeros(length(cu_re_0), 3);
    better_0_index = [];
    for i = 1:length(cu_re_0)
        [H,p_value,Z] = Mann_Kendall(cu_re_0{i}(end-min(10, length(cu_re_0{i})-1):end));
        std_reward_0(i, 1) = H;
        std_reward_0(i, 2) = p_value;
        std_reward_0(i, 3) = Z;
        if H == 1 && p_value <= 0.05 && Z > 0
            better_0_index = [better_0_index i];
        end
    end
    std_reward_1 = zeros(length(cu_re_1), 3);
    better_1_index = [];
    for i = 1:length(cu_re_1)
        [H,p_value,Z] = Mann_Kendall(cu_re_1{i}(end-min(10, length(cu_re_1{i})-1):end));
        std_reward_1(i, 1) = H;
        std_reward_1(i, 2) = p_value;
        std_reward_1(i, 3) = Z;
        if H == 1 && p_value <= 0.05 && Z > 0
            better_1_index = [better_1_index i];
        end
    end
    mask_similarity = [];
    if isempty(better_0_index) && isempty(better_1_index)
        for ttry = 1:PopulationSize
            temp_mask = [];
            for im = 1:sum(GenotypeLength)
                mask_genes = [0:9] * sum(GenotypeLength) + im;
                remain_genes = setdiff([1:length(XVessel) * sum(GenotypeLength)],mask_genes);
                remain_similarity = [];
                for jtry = 1:PopulationSize
                    if ttry ~= jtry
                        remain_similarity = [remain_similarity (1 ./ (sum(DecodePop(ttry, remain_genes) == DecodePop(jtry, remain_genes)) / length(remain_genes))) ^ 2];
                    end
                end
                temp_mask = [temp_mask mean(remain_similarity)];
            end
            mask_similarity = [mask_similarity mean(temp_mask)];
        end
    else
        for ttry = 1:PopulationSize
            temp_mask = [];
            mask_genes = better_0_index;
            remain_genes = setdiff([1:length(XVessel) * sum(GenotypeLength)],mask_genes);
            remain_similarity = [];
            for jtry = 1:PopulationSize
                if ttry ~= jtry
                    remain_similarity = [remain_similarity (1 ./ (sum(DecodePop(ttry, remain_genes) == DecodePop(jtry, remain_genes)) / length(remain_genes))) ^ 2];
                end
            end
            temp_mask = [temp_mask mean(remain_similarity)];
            mask_genes = better_1_index;
            remain_genes = setdiff([1:length(XVessel) * sum(GenotypeLength)],mask_genes);
            remain_similarity = [];
            for jtry = 1:PopulationSize
                if ttry ~= jtry
                    remain_similarity = [remain_similarity (1 ./ (sum(DecodePop(ttry, remain_genes) == DecodePop(jtry, remain_genes)) / length(remain_genes))) ^ 2];
                end
            end
            temp_mask = [temp_mask mean(remain_similarity)];
            mask_similarity = [mask_similarity mean(temp_mask)];
        end
    end
    genecoorpop = DecodePop(:, 1:120);
    [geneR,geneP] = corrcoef(genecoorpop);
    genecoor = geneP < 0.05;
    genecoor = double(genecoor);
    coorwithothers = sum(genecoor) / 120;
    coorthre = quantile(coorwithothers,0.75,2);  
    secondMask =  find(coorwithothers > coorthre);
    secondRemain = setdiff([1:length(XVessel) * sum(GenotypeLength)],secondMask);
    second_mask_similarity = [];
    for ttry = 1:PopulationSize
        seconde_remain_similarity = [];
        for jtry = 1:PopulationSize
            if ttry ~= jtry
                seconde_remain_similarity = [seconde_remain_similarity (1 ./ (sum(DecodePop(ttry, secondRemain) == DecodePop(jtry, secondRemain)) / length(secondRemain))) ^ 2];
            end
        end
        second_mask_similarity = [second_mask_similarity mean(seconde_remain_similarity)];
    end
    [Y1,~] = mapminmax(mask_similarity,0,1);
    [Y2,~] = mapminmax(second_mask_similarity,0,1);
    third_matrix = [];
    [Ytemp,~] = mapminmax(DecodePop(:, 161)',0,1);
    third_matrix = [third_matrix Ytemp'];
    [Ytemp,~] = mapminmax(DecodePop(:, 162)',0,1);
    third_matrix = [third_matrix Ytemp'];
    [Ytemp,~] = mapminmax(DecodePop(:, 163)',0,1);
    third_matrix = [third_matrix Ytemp'];
    dist_matrix = [];
    for ttry = 1:PopulationSize
        temp = [];
        for jtry = 1:PopulationSize
            if ttry ~= jtry
                tempdist = sqrt((third_matrix(ttry, 1) -  third_matrix(jtry, 1)) ^ 2 + (third_matrix(ttry, 2) -  third_matrix(jtry, 2)) ^ 2 + (third_matrix(ttry, 3) -  third_matrix(jtry, 3)) ^ 2);
                temp = [temp; tempdist];
            end
        end
        dist_matrix = [dist_matrix; mean(temp) max(temp) min(temp)];
    end
    dist_matrix(:,4) = max(dist_matrix(:,1)) - dist_matrix(:,1);
    dist_matrix(:,5) = max(dist_matrix(:,2)) - dist_matrix(:,2);
    dist_matrix(:,6) = max(dist_matrix(:,3)) - dist_matrix(:,3);
    dist_matrix(:,7) = IDAf(dist_matrix(:, 4:6));
    [Ytemp,~] = mapminmax(dist_matrix(:,4)',0,1);
    third_matrix = [third_matrix Ytemp'];
    third_matrix = [third_matrix mean(third_matrix(:,1:4)')'];
    tempY3 = max(third_matrix(:,5)) - third_matrix(:,5);
    [Y3,~] = mapminmax(tempY3',0,1);
    wy1 = 1/mean(Y1);
    wy2 = 1/mean(Y2);
    wy3 = 1/mean(Y3);
    final_Y = wy1 * Y1 + wy2 * Y2 + wy3 * Y3;
    select_prb = final_Y ./ sum(final_Y);
    selectResult = [];
    chooseCount = PopulationSize;
    for it = 1:chooseCount
        sumP = 0;
        randNum = rand();
        for i = 1:length(select_prb)
            sumP = sumP + select_prb(i);
            if sumP > randNum
                selectResult(it) = i;
                break;
            end
        end
    end
    length(unique(selectResult))
    DecodePop = DecodePop(selectResult,:);
    randorder = randperm(PopulationSize);
    DecodePop = [DomiPop(1:min(size(DomiPop,1),MemorySize),:);DecodePop(randorder(1:PopulationSize-min(size(DomiPop,1),MemorySize)),:)];
    iter = iter+1
end
DomiPop = DecodeFunction(DomiPop,XVessel,YVessel,RVessel,GenotypeLength,size(DomiPop,1),VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
Epa = DomiPop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
load train
sound(y,Fs)
% Output solutions in memory cell
f1 = Epa(:,1);
f2 = Epa(:,2);
f3 = Epa(:,3);
figure(2)
axis([min(f1)-50 max(f1)+50 min(f2)-50 max(f2)+50 min(f3)-0.5 max(f3)+0.5])
hold on
h = plot3(f1,f2,f3,'o','MarkerFaceColor',[214/255, 39/255, 40/255],'MarkerEdgeColor',[214/255, 39/255, 40/255]); grid on;
hXLabel = xlabel('Objective Function 1');
hYLabel = ylabel('Objective Function 2');
hZLabel = zlabel('Objective Function 3');
title('Pareto-optimal Solution Space')
hold on
set(h,'LineWidth',1.5)
set(gca, 'Box', 'on', ...
    'LineWidth', 1,...
    'TickDir', 'in', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1])
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set(gca, 'FontSize', 13)
set([hXLabel, hYLabel, hZLabel], 'FontSize', 14)
VesselPro = 1;
VesselPose(1,:) = [-2*VesselPro 0.8*VesselPro];
VesselPose(2,:) = [-VesselPro -0.8*VesselPro];
VesselPose(3,:) = [VesselPro -0.8*VesselPro];
VesselPose(4,:) = [2*VesselPro 0.8*VesselPro];
VesselPose(5,:) = [-2*VesselPro 0.8*VesselPro];
VesselPose1(1,:) = [-VesselPro 5*VesselPro];
VesselPose1(2,:) = [-VesselPro 0.8*VesselPro];
VesselPose1(3,:) = [VesselPro 2.5 * VesselPro];
VesselPose1(4,:) = [-VesselPro 5*VesselPro];
for iplot = 1:size(Epa,1)
    figure
    grid on
    axis([0 100 0 100])
    hold on
    axis square;
    hold on
    for i = 1:length(XWind)
        theta = 0:pi/50:2*pi;
        x0 = XWind(i)+RWind(i)*cos(theta);
        y0 = YWind(i)+RWind(i)*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+12*cos(theta);
        y0 = YWind(i)+12*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+9*cos(theta);
        y0 = YWind(i)+9*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+6*cos(theta);
        y0 = YWind(i)+6*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+3*cos(theta);
        y0 = YWind(i)+3*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+0.5*cos(theta);
        y0 = YWind(i)+0.5*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
    end
    for i = 1:length(XVessel)
        fill([VesselPose(1,1)+XVessel(i),VesselPose(2,1)+XVessel(i),VesselPose(3,1)+XVessel(i),VesselPose(4,1)+XVessel(i),VesselPose(5,1)+XVessel(i)], [VesselPose(1,2)+YVessel(i),VesselPose(2,2)+YVessel(i),VesselPose(3,2)+YVessel(i),VesselPose(4,2)+YVessel(i),VesselPose(5,2)+YVessel(i)],[31/255,119/255,180/255])
        theta = 0:pi/50:2*pi;
        x0 = XVessel(i)+RVessel(i)*cos(theta);
        y0 = YVessel(i)+RVessel(i)*sin(theta);
        plot(x0,y0,'-.','Color',[31/255,119/255,180/255],'LineWidth',1.2);
        for ixp = 1:4
            plot([VesselPose(ixp,1)+XVessel(i) VesselPose(ixp+1,1)+XVessel(i)],[VesselPose(ixp,2)+YVessel(i) VesselPose(ixp+1,2)+YVessel(i)],'Color',[31/255,119/255,180/255])
            hold on
        end
        for ixp = 1:3
            plot([VesselPose1(ixp,1)+XVessel(i) VesselPose1(ixp+1,1)+XVessel(i)],[VesselPose1(ixp,2)+YVessel(i) VesselPose1(ixp+1,2)+YVessel(i)],'Color',[31/255,119/255,180/255],'LineWidth',1.2)
            hold on
        end
    end
    h=plot(PortPosition(1),PortPosition(2),'s','Color','k','MarkerFaceColor','k','MarkerSize',5);
    hold on
    hXLabel=xlabel('x');
    hYLabel=ylabel('y');
    result = DomiPop(iplot,VarietyNumber*sum(GenotypeLength)+1:end-3);
    catetype=[];
    rou=[];
    cita=[];
    dijige=[];
    for i=1:length(result)
        if mod(i,4)==1
            rou=[rou result(i)];
        end
        if mod(i,4)==2
            cita=[cita result(i)];
        end
        if mod(i,4)==3
            catetype=[catetype result(i)];
        end
        if mod(i,4)==0
            dijige=[dijige result(i)];
        end
    end
    for i=1:length(XVessel)
        x111(i)=XVessel(i)+rou(i)*cos(cita(i));
        y111(i)=YVessel(i)+rou(i)*sin(cita(i));
    end
    colorarray = 'kkkkkkkkkkkkkkkkkkkkkkk';
    for i=1:MaxRobotNo
        cate=find(catetype==i);
        if isempty(cate)
            continue
        end
        [vv,v]=sort(dijige(cate));
        po=[PortPosition;[x111(cate(v))' y111(cate(v))'];PortPosition];
        for j = 1:size(po',2)-1
            arrow([po(j,1),po(j,2)],[po(j+1,1),po(j+1,2)],'Color','k','Type','Line','Width',2);
            hold on
        end
        plot(x111(cate(v))',y111(cate(v))','.','Color',[255/255,127/255,14/255],'MarkerSize',10);
    end
    set(h,'LineWidth',1.5)
    set(gca, 'Box', 'on', ...
        'LineWidth', 1,...
        'XGrid', 'off', 'YGrid', 'off', ...
        'TickDir', 'in', 'TickLength', [.015 .015], ...
        'XMinorTick', 'off', 'YMinorTick', 'off', ...
        'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1])
    set(gca, 'FontName', 'Helvetica')
    set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
    set(gca, 'FontSize', 15)
    set([hXLabel, hYLabel], 'FontSize', 16)
end
Epa = DomiPop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
f1 = Epa(:,1);
f2 = Epa(:,2);
f3 = Epa(:,3);
plot3(f1,f2,f3,'ro','MarkerFaceColor','r'); grid on;
legend(['All solution count ',num2str(size(AllPop,1))],['Optimal solution count ',num2str(size(DomiPop,1))])
figure(21)
grid on
subplot(3,1,1)
hold on
plot([1:MaxGeneration],reflag1,'d-','Color',[0.1,0.1,0.1],'MarkerFaceColor',[0.1,0.1,0.1],'MarkerSize',3)
hold on
plot([1:MaxGeneration],reflag2,'d-','Color',[0.1,0.7,0.2],'MarkerFaceColor',[0.1,0.7,0.2],'MarkerSize',3)
hold on
plot([1:MaxGeneration],reflag3,'d-','Color',[0.8,0.1,0.6],'MarkerFaceColor',[0.8,0.1,0.6],'MarkerSize',3)
hold on
legend(['No mutation ',num2str(sum(reflag1))],['Random multi-point mutation ',num2str(sum(reflag2))],['Mutation process regulation ',num2str(sum(reflag3))])
xlabel('Number of Iterations')
ylabel('Number of Redundant Patrol')
subplot(3,1,2)
grid on
plot([1:MaxGeneration],crflag1,'d-','Color',[0.1,0.1,0.1],'MarkerFaceColor',[0.1,0.1,0.1],'MarkerSize',3)
hold on
plot([1:MaxGeneration],crflag2,'d-','Color',[0.1,0.7,0.2],'MarkerFaceColor',[0.1,0.7,0.2],'MarkerSize',3)
hold on
plot([1:MaxGeneration],crflag3,'d-','Color',[0.8,0.1,0.6],'MarkerFaceColor',[0.8,0.1,0.6],'MarkerSize',3)
legend(['No mutation ',num2str(sum(crflag1))],['Random multi-point mutation ',num2str(sum(crflag2))],['Mutation process regulation ',num2str(sum(crflag3))])
xlabel('Number of Iteration')
ylabel('Number of Overlappint Patrol')
subplot(3,1,3)
hold on
plot([1:MaxGeneration],siflag1,'d-','Color',[0.1,0.1,0.1],'MarkerFaceColor',[0.1,0.1,0.1],'MarkerSize',3)
hold on
plot([1:MaxGeneration],siflag2,'d-','Color',[0.1,0.7,0.2],'MarkerFaceColor',[0.1,0.7,0.2],'MarkerSize',3)
hold on
plot([1:MaxGeneration],siflag3,'d-','Color',[0.8,0.1,0.6],'MarkerFaceColor',[0.8,0.1,0.6],'MarkerSize',3)
hold on
legend(['No mutation ',num2str(sum(siflag1))],['Random multi-point mutation ',num2str(sum(siflag2))],['Mutation process regulation ',num2str(sum(siflag3))])
xlabel('Number of Iterations')
ylabel('Number of Self-intersected Patrol')

figure(20)
subplot(3,1,1)
grid on
hold on
plot([1:MaxGeneration],reflag1,'d-','MarkerFaceColor','default','MarkerSize',3)
hold on
plot([1:MaxGeneration],reflag2,'rd-','MarkerFaceColor','default','MarkerSize',3)
hold on
plot([1:MaxGeneration],reflag3,'kd-','MarkerFaceColor','default','MarkerSize',3)
hold on
legend(['No mutation'],['Random multi-point mutation'],['Mutation process regulation'])
xlabel('Number of Iterations')
ylabel('Number of Redundant Patrol')
title('Redundant Patrol')
subplot(3,1,2)
grid on
hold on
plot([1:MaxGeneration],crflag1,'d-','MarkerSize',3)
hold on
plot([1:MaxGeneration],crflag2,'rd-','MarkerSize',3)
hold on
plot([1:MaxGeneration],crflag3,'kd-','MarkerSize',3)
legend(['No mutation'],['Random multi-point mutation'],['Mutation process regulation'])
xlabel('Number of Iteration')
ylabel('Number of Overlappint Patrol')
title('Overlapping Patrol')
subplot(3,1,3)
grid on
hold on
plot([1:MaxGeneration],siflag1,'d-','MarkerFaceColor','default','MarkerSize',3)
hold on
plot([1:MaxGeneration],siflag2,'rd-','MarkerFaceColor','default','MarkerSize',3)
hold on
plot([1:MaxGeneration],siflag3,'kd-','MarkerFaceColor','default','MarkerSize',3)
hold on
legend(['No mutation'],['Random multi-point mutation'],['Mutation process regulation'])
xlabel('Number of Iterations')
ylabel('Number of Self-intersected Patrol')
title('Self-intersected Patrol')
allreflag1 = [allreflag1; reflag1];
allreflag2 = [allreflag2; reflag2];
allreflag3 = [allreflag3; reflag3];
allcrflag1 = [allcrflag1; crflag1];
allcrflag2 = [allcrflag2; crflag2];
allcrflag3 = [allcrflag3; crflag3];
allsiflag1 = [allsiflag1; siflag1];
allsiflag2 = [allsiflag2; siflag2];
allsiflag3 = [allsiflag3; siflag3];
save(['res',sprintf('%d.mat', allhv(end))])