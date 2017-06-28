function [indexOnset, indexOffset] = detectPM(tEqui, voltPM, ...
    sogliaDerPos, sogliaDerNeg, sogliaGrowth, sogliaWin)
%detect onset and offset of calcium events based on growth and derivatives
%Input parameters:
%tEqui: time array
%voltPM: denoised fluorescence trace
%sogliaDerPos: threshold for the first right derivative at onset
%sogliaDerNeg: threshold for the first right derivative at offset
%sogliaGrowth: threshold for the growth
%sogliaWin: maximum duration of en event onset
%Output parameters:
%indexOnset: time index of all detected events onset
%indexOffset: time index of all detected events offset

%compute first derivative
DvoltPM=(voltPM(2:end)-voltPM(1:end-1))/(tEqui(2)-tEqui(1));

%putative onset time
indexUp = DvoltPM>sogliaDerPos;
indexUp = diff(indexUp);
Ind = find(indexUp == 1);
%putative offset time
indexDown = DvoltPM<sogliaDerNeg;
indexDown = diff(indexDown);
iI = 1;

indexOnset = [];
indexOffset = [];

for j=1:length(Ind)
    i1 = Ind(j);
    if iI>1 && Ind(j)<= indexOffset(iI-1)
        j=j+1;
        continue
    end
    if isempty(find(indexDown(i1:end)==1,1,'first'))
        lastUp = length(voltPM);
    else
        lastUp = i1 + find(indexUp(i1:end)==-1,1,'first');
        firstDown = i1 + find(indexDown(i1:end)==1,1,'first');
        if firstDown-lastUp <= sogliaWin
            onset = i1;
            [~,k] = max(voltPM(i1:firstDown));
            offset =  min(i1+k,length(voltPM));
            if voltPM(offset)-voltPM(onset)>sogliaGrowth
                indexOnset(iI) = min(onset+1,length(voltPM));
                indexOffset(iI) = max(1,offset-1);
                iI = iI+1;
            end
        end
    end
    
end