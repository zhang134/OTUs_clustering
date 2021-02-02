% calculate the MCC for OTU clustering methods
% input: similarity matrix: matrix, label file: label, similarity
% threshold: similarity
% output: mcc
% example: mccIs = mcc(matrix, label_cdhit97,0.97)
function mccIs = mcc(matrix, label, similarity)
label_cell = label_to_cell(label);
TP = 0; % same OTU, distance <= threshold
FP = 0; % same OTU, distance >  threshold
TN = 0; % different OTUs, distance > threshold
FN = 0; % different OTUs, distance <= threshold
label_num = length(label_cell);
% TP and FP
for i = 1:label_num
    otu = label_cell{i};
    otu_size = length(otu);
    fprintf("i=%d, total = %d\n", i, label_num);
    if otu_size > 1
        for jj = 1 : otu_size - 1
            seqID1 = otu(jj);
            for kk = jj + 1 : otu_size
                seqID2 = otu(kk);
                if matrix(seqID1,seqID2) >= similarity
                    TP = TP + 1;
                else
                    FP = FP + 1;
                end
            end
        end
    end
end

% FN and TN
for i = 1:label_num-1
    otu1 = label_cell{i};
    otu1_size = length(otu1);
    fprintf("i=%d, total = %d\n", i, label_num);
    for j = i + 1 : label_num
        otu2 = label_cell{j};
        otu2_size = length(otu2);
        for k = 1:otu1_size
            seqqID1 = otu1(k);
            for v = 1:otu2_size
                seqqID2 = otu2(v);
                if matrix(seqqID1,seqqID2) < similarity
                    TN = TN + 1;
                else
                    FN = FN + 1;
                end
            end
        end
    end
end       
fenzi = TP*TN-FP*FN;
fenmu = sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN));
mccIs = fenzi/fenmu;
end

% convert two columns label to cell
function out=label_to_cell(label)
all_label = label(:,2);
unique_label=unique(all_label);
label_num = length(unique_label);
out = cell(label_num,1);
for i = 1:label_num
    label_is = unique_label(i);
    indexxx = find(all_label==label_is);
    out{i} = indexxx;
end
end
    