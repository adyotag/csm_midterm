classdef Read_output < handle
   properties
       f = NaN;
   end
   methods
       function self = Read_output(f)
          self.f = f;
       end
       function r = nsad(self)
            filename = 'NodalStressAndDisp_' + extractBefore(self.f, ".") + '.mat';
            load(filename);
            r=farr;
       end
       function r = ips(self)
            filename = 'IPStresses_' + extractBefore(self.f, ".") + '.mat';
            load(filename);
            r=farr;
       end
       
   end    
end