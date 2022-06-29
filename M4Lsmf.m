function [] = M4Lsmf( nTypes,instanseIdx,Gcell,Rcell,thetaCell,max_iter )
LM=Rcell{2};
alpha = 10^7;
beta = 10^6;
[r_thetacell, c_thetacell] = size(thetaCell);
G_enum = cell(size(Gcell));
G_denom = cell(size(Gcell));
Wr = zeros(length(Rcell),1);
Wh = ones(r_thetacell,c_thetacell);
threshold = 0.00001;


Y=LM;
Rcell{2} = Y;

for ii = 1:length(Gcell)
    G_enum{ii}=0;
    G_denom{ii} =0;
end
Srcell = cell(size(Rcell));
Shcell = cell(r_thetacell,c_thetacell);

for iter = 1:max_iter
    
    mus=zeros(length(Rcell),1);
    
    for rr = 1:length(instanseIdx)
        i = fix(instanseIdx{rr}/nTypes)+1;
        j = mod(instanseIdx{rr},nTypes);
        if j ==0
            i = i-1;
            j = 4;
        end
        Gmatii = Gcell{i};
        Gmatjj = Gcell{j};
        
        Rmat = Rcell{rr};
        
        Smat = pinv(Gmatii'*Gmatii)*Gmatii'*Rmat*Gmatjj*pinv(Gmatjj'*Gmatjj);
        Smat(isnan(Smat))=0;
        Srcell{rr}=Smat;
        
        result = Rcell{rr}-Gmatii*Srcell{rr}*Gmatjj';
        R = sum(sum(result.^2));
        mus(rr)=R;
        
    end
    
    kmus=zeros(r_thetacell, c_thetacell);
    
    for ii = 1:length(Gcell)
        Gmatii = Gcell{ii};
        temp_sh = 0;
        temp_wh = 0;
        
        for jj = 1:c_thetacell
            if ~isempty(thetaCell{ii,jj})
                temp_sh = temp_sh+Wh(ii,jj).*(Gmatii'*thetaCell{ii,jj}*Gmatii);
                temp_wh = temp_wh+Wh(ii,jj);
            end
        end
        Shmat = pinv(Gmatii'*Gmatii)*(temp_sh/temp_wh)*pinv(Gmatii'*Gmatii);
        Shmat(isnan(Shmat)) = 0;
        Shcell{ii}  = Shmat;
    end
    
    for ii = 1:r_thetacell
        for jj = 1:c_thetacell
            if ~isempty(thetaCell{ii,jj})
                Gmatii = Gcell{ii};
                result2= thetaCell{ii,jj}-Gmatii*Shcell{ii}*Gmatii';
                R2 = sum(sum(result2.^2));
                kmus(ii,jj) = R2;
            end
        end
    end
    
    [~,Wr]=getOptimalWrWeights(mus,alpha);
    fprintf('the Wr is finished!\n');
    
    [~,Wh]=getOptimalWhWeights(kmus,beta);
    fprintf('the Wh is finished!\n');
    
    for rr = 1:length(instanseIdx)
        i = fix(instanseIdx{rr}/nTypes)+1;
        j = mod(instanseIdx{rr},nTypes);
        if j ==0
            i = i-1;
            j = 4;
        end
        temp1 = Rcell{rr}*Gcell{j}*Srcell{rr}';
        temp1(isnan(temp1))=0;
        t = abs(temp1);
        temp1p = (t+temp1)/2;
        temp1n = (t-temp1)/2;
        
        temp2 = Srcell{rr}*Gcell{j}'*Gcell{j}*Srcell{rr}';
        temp2(isnan(temp2))=0;
        t = abs(temp2);
        temp2p = (t+temp2)/2;
        temp2n = (t-temp2)/2;
        
        temp3 = Rcell{rr}'*Gcell{i}*Srcell{rr};
        temp3(isnan(temp3))=0;
        t = abs(temp3);
        t(t<0)=0;
        temp3p = (t+temp3)/2;
        temp3n = (t-temp3)/2;
        
        temp4 = Srcell{rr}'*Gcell{i}'*Gcell{i}*Srcell{rr};
        temp4(isnan(temp4))=0;
        t = abs(temp4);
        t(t<0)=0;
        temp4p = (t+temp4)/2;
        temp4n = (t-temp4)/2;
        
        
        G_enum{i} = G_enum{i}+Wr(rr).* temp1p+Wr(rr).*Gcell{i}*temp2n;
        G_denom{i}= G_denom{i}+Wr(rr).*temp1n+Wr(rr).*Gcell{i}*temp2p;
        
        G_enum{j} = G_enum{j}+ Wr(rr).*temp3p+Wr(rr).*Gcell{j}*temp4n;
        G_denom{j}= G_denom{j}+Wr(rr).*temp3n+Wr(rr).*Gcell{j}*temp4p;
        
    end
    
    for ii = 1:r_thetacell
        for jj= 1:c_thetacell
            if ~isempty(thetaCell{ii,jj})
                temp5 = thetaCell{ii,jj}*Gcell{ii}*Shcell{ii}';
                temp5(isnan(temp5))=0;
                t = abs(temp5);
                t(t<0)=0;
                temp5p = (t+temp5)/2;
                temp5n = (t-temp5)/2;
                
                temp6 = Gcell{ii}*Shcell{ii}*Gcell{ii}'*Gcell{ii}*Shcell{ii}';
                temp6(isnan(temp6)) = 0;
                t = abs(temp6);
                t(t<0)=0;
                temp6p = (t+temp6)/2;
                temp6n = (t-temp6)/2;
                
                G_enum{ii} = G_enum{ii}+2.*Wh(ii,jj).*temp5p+2.*Wh(ii,jj).*temp6n;
                G_denom{ii} = G_denom{ii}+2.*Wh(ii,jj).*temp5n+2.*Wh(ii,jj).*temp6p;
            end
        end
    end
    
    for ii = 1:length(Gcell)
        G_denom{ii}=G_denom{ii}+eps;
        factor = sqrt(G_enum{ii}./G_denom{ii});
        Gcell{ii}=Gcell{ii}.*factor;
        Gcell{ii}(isnan(Gcell{ii}))=0;
        Gcell{ii}(isinf(Gcell{ii}))=0;
    end
     
    result = Rcell{2}-Gcell{1}*Srcell{2}*Gcell{3}';
    R = sum(sum(result.^2));
    if R<threshold
        break;
    end
   fprintf('the iteration is: %d£¬RSS£º%f\n',iter,R);
end
newF=Gcell{1}*Srcell{2}*Gcell{3}';

end

