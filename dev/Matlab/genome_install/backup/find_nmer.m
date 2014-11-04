% DNA sequence tool.
%   Returns location in nmers array of recieved seq.
%   Faster than matlab(find) because nmer vector is highly structured.
function [nmer,err] = find_nmer(seq)
   steps = length(seq);
   nt    = zeros(1,steps);
   err   = false;

   for i = 1:steps
       switch(upper(seq(i)))
           case 'A'
               nt(i) = 0;
           case 'T'
               nt(i) = 1;
           case 'C'
               nt(i) = 2;
           case 'G'
               nt(i) = 3;
           otherwise
               err = true;
       end;
   end;
   nmer = 1;
   for i = 1:steps
       nmer = nmer + nt(i)*4^(i-1);
   end;
end
