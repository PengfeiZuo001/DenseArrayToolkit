iEvent = 20;
evid = eventIDs{iEvent};
[gather, matchIndex] = getCommonEventGather(DataStruct, evid);
[gatherReconstructed, d1_otg, reconGrid] = rankReduction(gather, gridStruct, RankReductionParam);
% Convert offsets to h = distance - min(distance)
px = linspace(-0.01,0.01,20);
py = linspace(-0.01,0.01,20);
t=(0:size(d1_otg,1)-1)*dt;
hx = reconGrid.x;
hy = reconGrid.y;
Param.hx  = hx;       % required by radon_op
Param.hy  = hy;
Param.px  = px;    % just for illustration if radon_op uses velocity
Param.py  = py;
Param.nt = length(t);
Param.dt = dt;
Param.type = 1;     % might indicate forward/backward radon in radon_op
npx = length(px);
npy = length(py);
% Preallocate transform model "ma"
ma = zeros(Param.nt, npx, npy);
N1 = 10;
N2 = 1;
%% 4.3 Perform Radon Transform on Z
try
    % Initialize the model with zeros
    mi_z = yc_pcg(@radon3d_op, Param, d1_otg, ma, N1, N2, 1);
    d1_otg_radon = radon3d_op(mi_z, Param, 1);  % forward modeling from the found model
    figure;
    imagesc([reshape(d1_otg,550,11*14) reshape(d1_otg_radon,550,11*14)])
    caxis([-0.02 0.02])
catch ME
    warnMsg = sprintf('[RadonTransform] Event %s, Z-component failed: %s', ...
        eventID, ME.message);
    warning(warnMsg);
    commonEventGather = appendHistory(commonEventGather, warnMsg);
    DataStruct(matchIndex) = commonEventGather;
end