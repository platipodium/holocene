% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz  <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later
%
% calculate area based average of RGR
%
tii=1;
n=length(wei);
rgr_ar=zeros(1,length(tim));

% loop over time window
for ti=breakPoints
  % load data for each time window
  load([outputDirectory 'mat/AllPop' tag '_' num2str(tii) '.mat']);
        % variables: poptime=tm, ymv=ymv,trgr=tirgr,rgr=rgrv,nreg=nregions
  %%ii=find(~isnan(rgr))' %'
  rgr(isnan(rgr))=0;

  if length(trgr)~= n,
      fprintf('dim mismatch %d: %d\n',length(trgr),n);
      ti=9E9;
  else
    tarea=sum(areav{tii}(1:nregv(tii))); % total area

    rgrm = areav{tii}(1:nregv(tii))*rgr/tarea; % weighed average

    % down-values near the edges of the time window
    i0= (ti-tim0+ddt(1))/dt+1;
    ii=i0:(i0+n-1);
    rgr_ar(ii)=rgr_ar(ii)+fliplr(wei.*rgrm); % stores into vector
    %%fprintf('%d: %d %1.0f\t nreg %d %d %d\t%1.2f %1.2f\n',tii,ti,tarea*1E-4,nreg,size(rgr,1),size(areav{tii},2),1E3*mean(rgrm),mean(areav{tii}(1:nregv(tii))));

    tii=tii+1;
    nregv2(tii)=nreg;

  end
end %ti
