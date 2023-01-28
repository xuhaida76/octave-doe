# Copyright (C) 2012 - Maria Christopoulou
# Copyright (C) 2001 - Per A. Brodtkorb
#
# This file is released under the terms of the CeCILL_V2 license : http:#www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function [y, id] = ryates(ef)
    # Reverse Yates algorithm to give estimated responses
    #
    # Calling Sequence
    #    y = ryates(ef)
    #    [y, id] = ryates(ef)
    #
    # Parameters
    #    ef : a n-by-1 vector of doubles, containing average response, main effects and interaction effects. ef(1) is the average response and ef(2:$) contain the main effects and interaction effects.
    #    y : a n-by-1 vector of doubles, the estimated response given the effects.
    #    id : a (n-1)-by-k matrix of doubles, identification matrix of main and interaction effects.
    #
    # Description
    #  Reverse Yates algorithm to give estimated responses
    #
	# This is the reverse function of scidoe_yates.
	#
    # This function is adapted from the Wafo toolbox:
	#
    # http:#www.maths.lth.se/matstat/wafo/documentation/wafodoc/wafo/wstats/ryates.html
    #
    # Examples
    # # "Statistics for experimenters Design,Innovation and Discovery",
    # # Second Edition, George E.P Box, J. Stuart Hunter, William Hunter, Wiley Series, 2005
    # # Tables 5.6, 5.7 - Pilot Plant Investigation Experiment
    #    D = scidoe_ff2n(3);
    #    y = [60 72 54 68 52 83 45 80];
    #    [ef,id] = scidoe_yates(y);
    #    y1 = scidoe_ryates(ef)
    #
    # Bibliography 
    # "Statistics for experimenters Design,Innovation and Discovery",
    # Second Edition, George E.P Box, J. Stuart Hunter, William Hunter, Wiley Series, 2005
    #
    # Authors
    # Copyright (C) 2012 - Maria Christopoulou
    # Copyright (C) 2001 - Per A. Brodtkorb

    # Check number of input and output arguments
    #[lhs,rhs] = argn();
    #apifun_checkrhs("scidoe_ryates",rhs,1)
    #apifun_checklhs("scidoe_ryates",lhs,1:2)
    #
    # Check type
    #apifun_checktype("scidoe_ryates",ef,"ef",1,"constant")
    #
    # Check whether ef is a column vector
    #apifun_checkveccol("scidoe_ryates",ef,"ef",1)
    #
    # Check if length of ef is power of two
    n = size(ef,1);   # Number of runs
    k = log2(n);      # Number of variables.
    #assert_checkalmostequal(k,round(k))

    # Reverse yates algorithm
    y = ef*(n/2);
    y(1,:) = y(1,:)*2;
    if (nargout>1) 
        [y,id] = yates(y(end:-1:1,:));
    else
        y = yates(y(end:-1:1,:));
    end
    y = y(end:-1:1,:)/2;
    y(end,:) = y(end,:)*2;

endfunction

%!test
%# Box, G.E.P, Hunter, W.G. and Hunter, J.S. (1978)
%# Statistics for experimenters, John Wiley & Sons, pp 342
%# Table 5.2 - Pilot Plant Investigation, page 177

%! y = [60 72 54 68 52 83 45 80]';
%! [ef, id] = yates(y);
%! y1 = ryates(ef);
%! assert(y,y1);