% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz  <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later
%
% moving average with distance decaying filter
%
function [arg]=movweighavg(times,values,window,offset)

values=reshape(values,size(times));
ji=1;
nut=length(times);
for it = 1:nut
  adt=abs(times-times(it));
  ind=find(adt<window*0.5);
  weigh=1./(adt(ind)+offset); % offset in years
  mavg(ji)=sum(values(ind).*weigh)/sum(weigh);
  ji=ji+1;
end
arg=mavg;

return
