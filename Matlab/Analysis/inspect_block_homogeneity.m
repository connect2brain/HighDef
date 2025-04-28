ID = 'sub-003';
exp_id = 'R2L';
responses = readtable(sprintf('B:/Projects/2023-01 HighDef/Results/%s_%s_raw.csv', ID, exp_id));



allBlocks = unique(responses.Block)';
allowedBlocks = [];
for i = allBlocks
    if max(responses.Intensity_percentMSO(responses.Block == i)) == max(responses.Intensity_percentMSO)
        allowedBlocks = [allowedBlocks, i];
    end
end

nAllowedBlocks = length(allowedBlocks);
blockColors = [linspace(0, 1, nAllowedBlocks); fliplr(linspace(0, 1, nAllowedBlocks)); [linspace(0,1,floor(nAllowedBlocks/2)) linspace(1,0,ceil(nAllowedBlocks/2))]]';


labels = {};
figure;
hold on
for i = 1:length(allowedBlocks)
    bl = allowedBlocks(i);
    blockMask = responses.Block == bl;
    thresholdInThisBlock = quantile(responses.CsE_FDI_in_uV(blockMask), 0.75);
    mask = blockMask & responses.CsE_FDI_in_uV > thresholdInThisBlock;
    scatter3(responses.p1(mask), responses.p2(mask), responses.p3(mask), 10, blockColors(i,:), 'o', "filled")
    scatter3(median(responses.p1(mask)), median(responses.p2(mask)), median(responses.p3(mask)), 50, blockColors(i,:), 'o', "filled")
    labels{2*i-1} = sprintf('Block %d', bl);
    labels{2*i} = '';
end

legend(labels)
axis square
box on
