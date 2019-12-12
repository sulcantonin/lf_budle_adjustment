%% Light Field Feature detection
% Inputs :      LF minimum num. of feature rays per light field feature
%            disparity_threshold: threshold for deviation the disparity
%            scale_difference:


% (c) O. Johannsen, A. Sulc, B. Goldluecke, University of Konstanz
% License: Creative Commons BY-SA 4.0,
% Please cite our paper if you use our code in your work:
% O. Johannsen, A. Sulc, B. Goldluecke: On Linear Structure from Motion for Light Field Cameras
function [ M,Dc ] = lf_extract_features( LF, varargin)

% parsing input parameters
p = inputParser;
addRequired (p,'LF',@isnumeric);
addParameter(p,'support',2,@isnumeric);
addParameter(p,'disparity_threshold',0.01,@isnumeric);
addParameter(p,'scale_difference',0.5,@isnumeric);
addParameter(p,'cross',false,@islogical);

parse(p,LF,varargin{:});
LF = p.Results.LF;
support = p.Results.support;
disparity_threshold = p.Results.disparity_threshold;
scale_difference = p.Results.scale_difference;
cross = p.Results.cross;

T = size( LF,1);
S = size( LF,2);

Sc = (S+1)/2;
Tc = (T+1)/2;

F = cell(T,S);
D = cell(T,S);

% detecting SIFT features, vl-feat is necessary (http://www.vlfeat.org/)
%{
for t=1:T
    I = im2single(rgb2gray(squeeze(LF(t,Sc,:,:,:))));
    [F{t,Sc},D{t,Sc}] = vl_sift(I);
end

for s=1:S
    I = im2single(rgb2gray(squeeze(LF(Tc,s,:,:,:))));
    [F{Tc,s},D{Tc,s}] = vl_sift(I);
end
%}
for t=1:T
    for s =1:S
        if cross && (s ~= Sc && t ~= Tc)
            continue
        end
        
        I = im2single(rgb2gray(squeeze(LF(t,s,:,:,:))));
        [F{t,s},D{t,s}] = vl_sift(I);
    end
end



%% matching and putting together features in center view together with features
% in other subaperture views into one structure M

% center view descriptor
Dc = D{Tc,Sc};
% center view features
Fc = F{Tc,Sc};
% number of features
Nc = size(F{Tc,Sc},2);
% each cell is one set of feature rays which should represent the same
% light field feature (the same 3D point)
M = cell(Nc,1);
for t = 1:T
    for s = 1:S
        % skip center view
        if s == Sc && t == Tc || isempty(D{t,s})
            continue
        end
        % matches of current view with Center view
        m = vl_ubcmatch(Dc,D{t,s});
        for i=1:size(m,2)
            m1 = m(1,i);
            m2 = m(2,i);
            
            xc = Fc(1,m1); % x of cv
            yc = Fc(2,m1); % y of cv
            scalec = Fc(3,m1); % scale of cv
            x = F{t,s}(1,m2);
            y = F{t,s}(2,m2);
            scale = F{t,s}(3,m2); % scale
            
            % scale of the same 3D point should be same
            if abs(scalec-scale) > scale_difference
                continue
            end
            
            % adding center view if it wasn't done yet
            if isempty(M{m1})
                M{m1} = [Sc,Tc,xc,yc];
            end
            % adding current (t,s,) view in pixel coordinates
            M{m1} = cat(1,M{m1},[s t x y]);
        end
    end
end


%% Calulating (rough) disparity for future remove,
% filling up disparity tables for respective directions
dtable = cell(Nc,1);
for i=1:Nc
    M_ = M{i};
    if isempty(M_)
        continue
    end
    t = M_(:,1); s = M_(:,2);
    y = M_(:,3); x = M_(:,4);
    
    Mc = length(t);
    dtable{i} = zeros(Mc,1);
    
    % center view
    cv = t == Tc & s == Sc;
    yc = y(cv);
    xc = x(cv);
    
    for j=1:Mc
        dy = (y(j)-yc) ./ (t(j)-Tc);
        dx = (x(j)-xc) ./ (s(j)-Sc);
        
        if ~ isfinite( dy )
            dy = 0;
        end
        
        if ~ isfinite( dx )
            dx = 0;
        end
        
        dtable{i}(j) = (dy + dx) / ((dy == 0) + (dx == 0));
    end
end

%% Eliminating features in other subaperture viwes whose dmap differs
% by scale disparity_threshold from median or > 5
for i=1:Nc
    d = dtable{i};
    d = d(d ~= 0);
    dcv = 0;
    
    if ~isempty(d)
        % center view disparity
        dcv = nanmedian( d );
    end
    % center view is always the first element
    dtable{i}(1) = dcv;
    dtable{i}(:) = dtable{i} .* ( ...
        (abs(dtable{i} - dcv) < disparity_threshold) & ...
        (abs(dtable{i}) <= 5) ) ;
end

for i=1:Nc
    d = dtable{i};
    M{i} = M{i}(d~=0,:);
end


%% Eliminating those features which doesn't have enough subaperture views
% (2 rays are needed to triangulate a 3D point)
cnt = arrayfun(@(i) size(M{i},1),1:length(M));
keep = find(cnt >= max(2,support));
M = M(keep);
Dc = Dc(:,keep);