## Copyright (C) 2013 - Michael Baudin
## Copyright (C) 2012 - Maria Christopoulou
##
## This file must be used under the terms of the CeCILL.
## This source file is licensed as described in the file COPYING, which
## you should have received as part of this distribution.  The terms
## are also available at
## http:##www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function H = fracfact(gen)
    ## Fractional Factorial Design
    ##
    ## Calling Sequence
    ##    H = scidoe_fracfact(gen)
    ##
    ## Parameters
    ##    gen : a 1-by-n matrix of strings, consisting of lowercase ,uppercase letters or operators "-" and "+", indicating the factors of the experiment
    ##    H : a m-by-n matrix of doubles, the fractional factorial design. m is 2^k, where k is the number of letters in gen, and n is the total number of entries in "gen" 
    ##
    ## Description
    ## The function computes a fractional factorial design.
    ## Fractional factorial designs are used instead of full factorial 
    ## designs, because they need fewer runs.
    ##
    ## In "gen" we define the main factors of the experiment and the 
    ## factors whose levels are the products of the main factors.
    ## For example, if 
    ## <screen>
    ## gen = "a b ab"
    ## </screen>
    ## then "a" and "b" are the main factors, while 
    ## the 3rd factor is the product of the first two.
    ## If we input uppecase letters in gen, we get the same result.
    ## We can also use the operators "+" and "-" in gen.
    ##
    ## For example, if 
    ## <screen>
    ## gen = "a b -ab"
    ## </screen>
    ## then the 3rd factor is the opposite of the 
    ## product of "a" and "b".
    ## 
    ## The output matrix includes a two level full factorial design, 
    ## built by the main factors of gen, and the products of the main factors.
    ## The columns of H follow the sequence of "gen". 
    ## For example, if 
    ## <screen>
    ## gen = "a b ab c"
    ## </screen>
    ## then columns H(:,1), H(:,2) and H(4,:) 
    ## include the two level full factorial design and H(:,3) includes 
    ## the products of the main factors.
    ##
    ## Examples
    ## ## Fractional factorial design, where the 3rd factors 
    ## ## is the product of the first two
    ## H = scidoe_fracfact("a b ab")
    ## ## The same with uppercase letters
    ## H = scidoe_fracfact("A B AB")
    ## ## Use of operators
    ## H = scidoe_fracfact("a b -ab c +abc")
    ##
    ## Bibliography
    ## http:##www.itl.nist.gov/div898/handbook/pri/section3/pri3342.htm
    ## http:##www.mathworks.com/help/toolbox/stats/fracfact.html
    ## 
    ## Authors
    ## Copyright (C) 2012 - Maria Christopoulou
    ## Copyright (C) 2013 - Michael Baudin

    ## [lhs,rhs] = argn()
    ## apifun_checkrhs("scidoe_fracfact",rhs,1)
    ## apifun_checklhs("scidoe_fracfact",lhs,1)
    ##
    ## apifun_checktype("scidoe_fracfact",gen,"gen",1,"string")
    ##
    ## Recognize letters and combinations
    A = strsplit(gen,{"-" " " "+"});
    C = cellfun("length",A);
    ## Indices of single letters (main factors)
    I = find(C==1)
    ## Indices of letter combinations. 
    ## We need them to fill out H2 properly
    J = find(C~=1) 
    ##
    ## Check if there are "-" or "+" operators in gen
    U = strsplit(gen,{" "})
    ## If R1 is either null or not, the result is not 
    ## changed, since it is a multiplication by 1.
    R1 = strmatch("+",U)
    R2 = strmatch("-",U)
    ## Fill in design with two level factorial design
    H1 = 2*ff2n(length(I))-1
    H2 = zeros(size(H1,1),length(C))
    H2(:,I) = H1
    ##
    ## Recognize combinations and fill in the rest of matrix H2
    ## with the proper products
    ##
    for k=J(1,:)
        ## For Lowercase letters
        xx=int16(A{k})-96
        ## For Uppercase letters
        if(xx<0) 
            xx=int16(A{k})-64
        end
        H2(:,k) = prod(H1(:,xx),2)
    end
    ## Update design if gen includes "-" operator
    if (~isempty(R2)) 
        H2(:,R2) = (-1).*H2(:,R2);
    end
    ## Fractional Factorial Design
    H = H2;

endfunction

%!test

%! B_exp = [
%!   -1   -1    1;
%!    1   -1   -1;
%!   -1    1   -1;
%!    1    1    1;
%! ];
%! B_exp=sortdesign(B_exp)
%! B_comp = fracfact("a b ab")
%! B_comp = sortdesign(B_comp)
%! assert(B_comp,B_exp);

%! B_exp = [
%!	   -1   -1    1;
%!	    1   -1   -1;
%!	   -1    1   -1;
%!      1    1    1;
%!	];
%! B_exp=sortdesign(B_exp);
%! B_comp=fracfact("A B AB");
%! B_comp = sortdesign(B_comp)
%! assert(B_comp,B_exp)

%!test
%! B_exp=[
%!	   -1   -1   -1   -1    1;
%!	    1   -1   -1   -1   -1;
%!	   -1    1   -1   -1   -1;
%!	    1    1   -1   -1    1;
%!	   -1   -1    1   -1   -1;
%!	    1   -1    1   -1    1;
%!	   -1    1    1   -1    1;
%!	    1    1    1   -1   -1;
%!	   -1   -1   -1    1   -1;
%!	    1   -1   -1    1    1;
%!	   -1    1   -1    1    1;
%!	    1    1   -1    1   -1;
%!	   -1   -1    1    1    1;
%!	    1   -1    1    1   -1;
%!	   -1    1    1    1   -1;
%!	    1    1    1    1    1;
%!	];
%! B_exp=sortdesign(B_exp);
%! B_comp=fracfact("a b c d abcd");
%! B_comp = sortdesign(B_comp)
%! assert(B_comp,B_exp)

%!test
%! B_exp=[
%!	   -1   -1   -1   -1    1;
%!	    1   -1   -1   -1   -1;
%!	   -1    1   -1   -1   -1;
%!	    1    1   -1   -1    1;
%!	   -1   -1    1   -1   -1;
%!	    1   -1    1   -1    1;
%!	   -1    1    1   -1    1;
%!	    1    1    1   -1   -1;
%!	   -1   -1   -1    1   -1;
%!	    1   -1   -1    1    1;
%!	   -1    1   -1    1    1;
%!	    1    1   -1    1   -1;
%!	   -1   -1    1    1    1;
%!	    1   -1    1    1   -1;
%!	   -1    1    1    1   -1;
%!	    1    1    1    1    1;
%!	];
%! B_exp=sortdesign(B_exp);
%! B_comp=fracfact("a b c d -abcd");
%! B_comp = sortdesign(B_comp);
%! assert(B_comp,B_exp)
