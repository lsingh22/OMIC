function [] = read_multifil(filename, varargin)

%Read in coils from the coils file.
%Plot each coil and keep track of the current

%Assumptions: rectangular cross section = four vertices 

if nargin > 2 
    disp('ERROR: too many inputs!');
        return;
        
elseif nargin == 0
    disp('ERROR: coils filename required!');
    return;
    
else
    try
        data = textread(filename, '%s', 'delimiter', '\n');
    catch
        error('Coils File Not Found In Directory');
    end
        
    lines = size(data);     %data variable is an array of each line of filename
    numCoil = 0;
    
    i = 4;      %i,j,k,l,n,o are dummy variables
    j = 1;
    while i < lines(1)
        str = strjoin(data(i));
        C = strsplit(str,' ');   %splits elements in each line
        S = size(C);
        if S(2) == 6 || S(2) == 7
            numCoil = numCoil + 1;
            coilLine(1,j) = i;
            j = j + 1;
        end
        i = i + 1;
    end
    
    firstLine = zeros(1,numCoil);
    lastLine = zeros(1,numCoil);
    dataNum = zeros(numCoil,1);
    firstLine(1,1) = 4;
    k = 1;
    while k <= numCoil
        lastLine(1,k) = coilLine(1,k);
        if k ~= 1
            firstLine(1,k) = coilLine(1,k-1) + 1;

        end
        dataNum(k,1) = lastLine(1,k) - firstLine(1,k) + 1;
        k = k + 1;
    end
        
    X = zeros(max(dataNum),numCoil);   
    Y = zeros(max(dataNum),numCoil);
    Z = zeros(max(dataNum),numCoil);
    I = zeros(1,numCoil);   

    l = 1;
    m = 1;
    while m <= numCoil
        while l <= dataNum(m,1)
            str = strjoin(data(firstLine(1,m) + l - 1));
            C = strsplit(str,' ');
            X(l,m) = str2double(cell2mat(C(1,1)));
            Y(l,m) = str2double(cell2mat(C(1,2)));
            Z(l,m) = str2double(cell2mat(C(1,3))); 
            if(l == 1)
                %I(1,m) = str2double(cell2mat(C(1,4)));
            end
            l = l + 1;
        end
        m = m + 1;
        l = 1;
    end
    
    assignin('base','nice',X);
    n = 1;
    o = 1;
    
    Nseg = 128;

    
    for ii = 1:numCoil/8
       xx = zeros(Nseg + 1,5); 
       yy = zeros(Nseg + 1,5); 
       zz = zeros(Nseg + 1,5); 
       %xx = zeros(Nseg, 5); yy = xx; zz = xx; 
        for jj = 1: Nseg +1 
            xx( jj, 1 ) = X( (jj-1)*5 + 1 , ii);
            xx( jj, 2 ) = X( (jj-1)*5 + 2 , ii);
            xx( jj, 3 ) = X( (jj-1)*5 + 3 , ii);
            xx( jj, 4 ) = X( (jj-1)*5 + 4 , ii); 
            xx( jj, 5 ) = X( (jj-1)*5 + 5 , ii);
            
            yy( jj, 1 ) = Y( (jj-1)*5 + 1 , ii);
            yy( jj, 2 ) = Y( (jj-1)*5 + 2 , ii);
            yy( jj, 3 ) = Y( (jj-1)*5 + 3 , ii);
            yy( jj, 4 ) = Y( (jj-1)*5 + 4 , ii); 
            yy( jj, 5 ) = Y( (jj-1)*5 + 5 , ii);
            
            zz( jj, 1 ) = Z( (jj-1)*5 + 1 , ii);
            zz( jj, 2 ) = Z( (jj-1)*5 + 2 , ii);
            zz( jj, 3 ) = Z( (jj-1)*5 + 3 , ii);
            zz( jj, 4 ) = Z( (jj-1)*5 + 4 , ii); 
            zz( jj, 5 ) = Z( (jj-1)*5 + 5 , ii); 
        end
            % xx(Nseg +1, :) = xx(1,:);
            % yy(Nseg +1, :) = yy(1,:);
            % zz(Nseg +1, :) = zz(1,:);
            
%         coil = 6;
%         
%         
%         if (ii==coil)
%             surf(xx,yy,zz,'FaceColor','r'); hold on;
%         end  
        
        %if (ii==1 || ii==2 || ii==3 || ii==4 || ii==5 || ii==6)
           surf(xx,yy,zz,'FaceColor','r'); hold on;
        %end  
       % surf(xx,yy,zz,'FaceColor','blue'); hold on;
        
%         if (mod(ii,2) == 0)
%             surf(xx,yy,zz,'FaceColor','red'); hold on;
%         else
%             surf(xx,yy,zz,'FaceColor','blue'); hold on;
%         end
    end
    axis equal;
%assignin('base','nice',xx);
end 