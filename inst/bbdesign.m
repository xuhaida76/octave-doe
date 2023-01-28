## Copyright (C) 2012 - Michael Baudin
## Copyright (C) 2012 - Maria Christopoulou
## Copyright (C) 2009 - Yann Collette
##
## This file must be used under the terms of the CeCILL.
## This source file is licensed as described in the file COPYING, which
## you should have received as part of this distribution.  The terms
## are also available at
## http:##www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function bbdesign(varargin)

    ## Box-Benhken design of experiments
    ##
    ## Calling Sequence
    ##    H = bbdesign(nbvar)
    ##    H = bbdesign(nbvar,"center",nbcenter)
    ##
    ## Parameters
    ##    nbvar : a 1-by-1 matrix of doubles, integer value, nbvar >= 3. The number of variables of the experiment
    ##    nbcenter : a 1-by-1 matrix of doubles, integer value, positive. The number of repetitions of the central point in the design (default nbcenter depends on nbvar).
    ##    H : a m-by-nbvar matrix of doubles, the design of experiments in the range [-1,1], where m=nbvar*4+nbcenter
    ##
    ## Description
    ## Creates a Box-Benkhen Design of Experiments. 
    ## This type of design is sufficient for quadratic models. 
    ##
    ## This is a three-level design, where each variable 
    ## can only have one of the following values: -1, 0, +1. 
    ##
    ## By default, scidoe_bbdesign produces a number of center 
    ## points which depends on nbvar. 
    ## If nbvar>16, then nbcenter is equal to nbvar. 
    ## If nbvar<=16, the nbcenter = T(nbvar), where T is from the following table:
    ##
    ## <literal>[0 0 3 3 6 6 6 8 9 10 12 12 13 14 15 16]</literal>
    ##
    ## To have only one center point, use the 
    ## <literal>H=scidoe_bbdesign(nbvar,"center",1)</literal> calling sequence.
    ##
    ## In this design, there is no point at the corners of the 
    ## cube [-1,+1]. 
    ## This property can be advantageous when the corners of the cube 
    ## are too expensive, or even impossible, to compute. 
    ##
    ## On output, the rows in H are sorted in increasing order.
    ##
    ## This function does not entirely reproduce Matlab's bbdesign. 
    ## For nbvar = 9,16, the function gives different results than 
    ## the known bibliography.
    ## 
    ##
    ## Examples
    ## ## Create a Box-Benkhen Design with three 
    ## ## variables one repetition of center:
    ## H = scidoe_bbdesign(3)
    ## ## Plot the design
    ## h = scf();
    ## param3d(H(:,1),H(:,2),H(:,3))
    ## h.children.children.mark_mode="on";
    ## h.children.children.line_mode="off";
    ## h.children.children.mark_size=1;
    ## scidoe_plotcube(3)
    ## xtitle("Box-Benkhen Design","X1","X2","X3")
    ##
    ## ## With three variables and the 
    ## ## center repeated two times:
    ## H = scidoe_bbdesign(3,2)
    ##
    ## ## Print the number of experiments
    ## for nbvar = 3 : 10
    ##   H = bbdesign(nbvar);
    ##   m = size(H,"r");
    ##   mprintf("nbvar=%d, Num. Experiments=%d\n",..
    ##      nbvar,m)
    ## end
    ##
    ## Bibliography
    ## http:##en.wikipedia.org/wiki/Box%E2%80%93Behnken_design
    ## Design and analysis of experiments, Montgomery, 5th Edition, Wiley, 2001
    ## Plans d'expériences Applications à l'entreprise, Gilbert Saporta, Jean-Jacques Droesbeke, Jeanne Fine, Editions Technip, 1997
    ##
    ## Authors
    ## Copyright (C) 2012 - Michael Baudin
    ## Copyright (C) 2012 - Maria Christopoulou
    ## Copyright (C) 2009 - Yann Collette

    ## Check number of input and output arguments:
    ## [lhs, rhs] = argn()
    ## apifun_checkrhs("scidoe_bbdesign",rhs,[1 3])
    ## apifun_checklhs("scidoe_bbdesign",lhs,1)
    ##
  ##

    p=inputParser();
  
    nbvar = varargin{1}
    ## Check type, size, content of nbvar
    ## apifun_checktype("scidoe_bbdesign",nbvar,"nbvar",1,"constant")
    ## apifun_checkscalar("scidoe_bbdesign",nbvar,"nbvar",1)
    ## apifun_checkgreq("scidoe_bbdesign",nbvar,"nbvar",1,3)
    ## apifun_checkflint("scidoe_bbdesign",nbvar,"nbvar",1)
    ##
    ## Input and check of nbcenter
    
    if (nargin==1) 
        tableOfExtraPoints = [0 0 3 3 6 6 6 8 9 10 12 12 13 14 15 16];
        if (nbvar<=16) 
            nbcenter = tableOfExtraPoints(nbvar)
        else
            nbcenter = nbvar
        end
    end

    if (nargin==3) 
    ##    nbcenter = apifun_argindefault (varargin,3,[]) 
    ##    default.center = nbcenter;
        ##Set key value pairs
      ##    options = apifun_keyvaluepairs(default)
      p.addParameter("center",[]);
      p.parse(varargin{2:end});
      nbcenter = p.Results.center
    end
    ##
    ## Check type, content of nbcenter
    ## apifun_checktype("scidoe_bbdesign",nbcenter,"nbcenter",3,"constant")
    ## apifun_checkscalar("scidoe_bbdesign",nbcenter,"nbcenter",3)
    ## apifun_checkgreq("scidoe_bbdesign",nbcenter,"nbcenter",3,1)
    ## apifun_checkflint("scidoe_bbdesign",nbcenter,"nbcenter",3);
    ##
    ## Box Behnken Designs for Nbvar = 3:7
    ## Plans d'expériences Applications à l'entreprise, Gilbert Saporta, Jean-Jacques Droesbeke, Jeanne Fine, Editions Technip, 1997
    ## pages 254 - 255
    ##
    switch(nbvar)
        
    case 3 

        H =[
            -1 -1  0;
            1 -1  0;
            -1  1  0;
            1  1  0;
            -1  0 -1;
            1  0 -1;
            -1  0  1;
            1  0  1;
            0 -1 -1;
            0  1 -1;
            0 -1  1;
            0  1  1;
        ];

       
    case 4 

        H =[
            -1 -1  0  0;
            1 -1  0  0;
            -1  1  0  0;
            1  1  0  0;
            0  0 -1 -1;
            0  0  1 -1;
            0  0 -1  1;
            0  0  1  1;
            -1  0  0 -1;
            1  0  0 -1;
            -1  0  0  1;
            1  0  0  1;
            0 -1 -1  0;
            0  1 -1  0;
            0 -1  1  0;
            0  1  1  0;
            -1  0 -1  0;
            1  0 -1  0;
            -1  0  1  0;
            1  0  1  0;
            0 -1  0 -1;
            0  1  0 -1;
            0 -1  0  1;
            0  1  0  1;
        ];

        
    case 5 

        H = [
             -1 -1  0  0  0;
             1 -1  0  0  0;
             -1  1  0  0  0;
             1  1  0  0  0;
             0  0 -1 -1  0;
             0  0  1 -1  0;
             0  0 -1  1  0;
             0  0  1  1  0;
             0 -1  0  0 -1;
             0  1  0  0 -1;
             0 -1  0  0  1;
             0  1  0  0  1;
             -1  0 -1  0  0;
             1  0 -1  0  0;
             -1  0  1  0  0;
             1  0  1  0  0;
             0  0  0 -1 -1;
             0  0  0  1 -1;
             0  0  0 -1  1;
             0  0  0  1  1;
             0 -1 -1  0  0;
             0  1 -1  0  0;
             0 -1  1  0  0;
             0  1  1  0  0;
             -1  0  0 -1  0;
             1  0  0 -1  0;
             -1  0  0  1  0;
             1  0  0  1  0;
             0 0  -1  0 -1;
             0  0  1  0 -1;
             0  0 -1  0  1;
             0  0  1  0  1;
             -1  0  0  0 -1;
             1  0  0  0 -1;
             -1  0  0  0  1;
             1  0  0  0  1;
             0 -1  0 -1  0;
             0  1  0 -1  0;
             0 -1  0  1  0;
             0  1  0  1  0;
        ];

        ## nbvar = 6
    case 6 

        H = [
             -1 -1  0 -1  0  0;
             1 -1  0 -1  0  0;
             -1  1  0 -1  0  0;
             1  1  0 -1  0  0;
             -1 -1  0  1  0  0;
             1 -1  0  1  0  0;
             -1  1  0  1  0  0;
             1  1  0  1  0  0;
             0 -1 -1  0 -1  0;
             0  1 -1  0 -1  0;
             0 -1  1  0 -1  0;
             0  1  1  0 -1  0;
             0 -1 -1  0  1  0;
             0  1 -1  0  1  0;
             0 -1  1  0  1  0;
             0  1  1  0  1  0;
             0  0 -1 -1  0 -1;
             0  0  1 -1  0 -1;
             0  0 -1  1  0 -1;
             0  0  1  1  0 -1;
             0  0 -1 -1  0  1;
             0  0  1 -1  0  1;
             0  0 -1  1  0  1;
             0  0  1  1  0  1;
             -1  0  0 -1 -1  0;
             1  0  0 -1 -1  0;
             -1  0  0  1 -1  0;
             1  0  0  1 -1  0;
             -1  0  0 -1  1  0;
             1  0  0 -1  1  0;
             -1  0  0  1  1  0;
             1  0  0  1  1  0;
             0 -1  0  0 -1 -1;
             0  1  0  0 -1 -1;
             0 -1  0  0  1 -1;
             0  1  0  0  1 -1;
             0 -1  0  0 -1  1;
             0  1  0  0 -1  1;
             0 -1  0  0  1  1;
             0  1  0  0  1  1;
             -1 -1  0  0  0 -1;
             1 -1  0  0  0 -1;
             -1  1  0  0  0 -1;
             1  1  0  0  0 -1;
             -1 -1  0  0  0  1;
             1 -1  0  0  0  1;
             -1  1  0  0  0  1;
             1  1  0  0  0  1;
        ];

        ## nbvar = 7
    case 7 

        H = [
             0  0  0 -1 -1 -1  0;
             0  0  0  1 -1 -1  0;
             0  0  0 -1  1 -1  0;
             0  0  0  1  1 -1  0;
             0  0  0 -1 -1  1  0;
             0  0  0  1 -1  1  0;
             0  0  0 -1  1  1  0;
             0  0  0  1  1  1  0;
             -1  0  0  0  0 -1 -1;
             1  0  0  0  0 -1 -1;
             -1  0  0  0  0  1 -1;
             1  0  0  0  0  1 -1;
             -1  0  0  0  0 -1  1;
             1  0  0  0  0 -1  1;
             -1  0  0  0  0  1  1;
             1  0  0  0  0  1  1;
        0 -1  0  0 -1  0 -1;
        0  1  0  0 -1  0 -1;
        0 -1  0  0  1  0 -1;
        0  1  0  0  1  0 -1;
        0 -1  0  0 -1  0  1;
        0  1  0  0 -1  0  1;
        0 -1  0  0  1  0  1;
        0  1  0  0  1  0  1;
        -1 -1  0 -1  0  0  0;
        1 -1  0 -1  0  0  0;
        -1  1  0 -1  0  0  0;
        1  1  0 -1  0  0  0;
        -1 -1  0  1  0  0  0;
        1 -1  0  1  0  0  0;
        -1  1  0  1  0  0  0;
        1  1  0  1  0  0  0;
        0  0 -1 -1  0  0 -1;
        0  0  1 -1  0  0 -1;
        0  0 -1  1  0  0 -1;
        0  0  1  1  0  0 -1;
        0  0 -1 -1  0  0  1;
        0  0  1 -1  0  0  1;
        0  0 -1  1  0  0  1;
        0  0  1  1  0  0  1;
        -1  0 -1  0 -1  0  0;
        1  0 -1  0 -1  0  0;
        -1  0  1  0 -1  0  0;
        1  0  1  0 -1  0  0;
        -1  0 -1  0  1  0  0;
        1  0 -1  0  1  0  0;
        -1  0  1  0  1  0  0;
        1  0  1  0  1  0  0;
        0 -1 -1  0  0 -1  0;
        0  1 -1  0  0 -1  0;
        0 -1  1  0  0 -1  0;
        0  1  1  0  0 -1  0;
        0 -1 -1  0  0  1  0;
        0  1 -1  0  0  1  0;
        0 -1  1  0  0  1  0;
        0  1  1  0  0  1  0;
        ];

    end
    ##
    ## If nbvar > 7, then the design is computed by the following:
    ## Depending on the number of variables nbvar, we compute two level design
    ## which is populated with the real doe.
    if (nbvar==[8 13 14 15] | nbvar>16 ) 
        H = bbdesign2(nbvar)
    elseif (nbvar==[9]) 
        H = bbdesign3(nbvar)
    elseif(nbvar==[10 11 12 16])
        H = bbdesign4(nbvar)
    end
    C = zeros(nbcenter,nbvar);
    H = [H;C]

#    H = sortdesign(H)
endfunction

function H = bbdesign2(nbvar)
    H_fact = 2*ff2n(2)-1;
    ## Number of factorial points per block:
    ## number of rows in H_fact
    fppb = 2^2
    nb_lines = 2*(nbvar-1)*nbvar;
    H = zeros(nb_lines,nbvar);
    Index = 0;
    for i=1:nbvar-1
        for j=i+1:nbvar
            Index = Index + 1;
            startindex = max([1 ((Index-1)*fppb+1)])
            stopindex = Index*fppb
            rows = startindex:stopindex;
            H(rows,[i j]) = H_fact;
        end
    end
endfunction

function H = bbdesign3(nbvar)
    H_fact=2*ff2n(3)-1;
    fppb = 2^3
    nb_lines = nbvar*fppb
    H = zeros(nb_lines,nbvar);
    Index = 0;
    for i=1:nbvar-2
        for j=i+1:nbvar-1
            for k=j+1:nbvar
                Index = Index + 1;
                if (Index<=nbvar)
                    startindex = max([1 ((Index-1)*fppb+1)])
                    stopindex = Index*fppb
                    rows = startindex:stopindex;
                    H(rows,[i j k]) = H_fact;
                end
            end
        end
    end
endfunction
function H = bbdesign4(nbvar)
    H_fact=2*ff2n(4)-1;
    fppb = 2^4
    nb_lines = nbvar*fppb
    H = zeros(nb_lines,nbvar);
    Index = 0;
    for i=1:nbvar-3
        for j=i+1:nbvar-2
            for k=j+1:nbvar-1
                for l=k+1:nbvar
                    Index = Index + 1;
                    if (Index<=nbvar)
                        startindex = max([1 ((Index-1)*fppb+1)])
                        stopindex = Index*fppb
                        rows = startindex:stopindex;
                        H(rows,[i j k l]) = H_fact;
                    end
                end
            end
        end
    end
endfunction



%!test
%!
%!
%! nb_var = 3;
%! C = bbdesign(3,"center",1);
%! E = [
%!    -1 -1  0;
%!     1 -1  0;
%!    -1  1  0;
%!     1  1  0;
%!    -1  0 -1;
%!     1  0 -1;
%!    -1  0  1;
%!     1  0  1;
%!     0 -1 -1;
%!     0  1 -1;
%!     0 -1  1;
%!     0  1  1;
%!     0  0  0;
%! ];
%! E = sortdesign(E);
%! assert(C,E);
%!
%!test
%!# Test with both input parameters
%!
%! C = bbdesign(5,"center",3);
%! E =   [
%!     -1.  -1.    0.    0.    0.;
%!      1.  -1.    0.    0.    0.;
%!      -1.  1.    0.    0.    0.;
%!      1.   1.    0.    0.    0.;
%!      -1.  0.   -1.    0.    0.;
%!      1.   0.   -1.    0.    0.
%!      -1.  0.    1.    0.    0.
%!      1.   0.    1.    0.    0.
%!      -1.  0.    0.   -1.    0.
%!      1.   0.    0.   -1.    0.
%!      -1.  0.    0.    1.    0.
%!      1.   0.    0.    1.    0.
%!      -1.  0.    0.    0.   -1.
%!      1.   0.    0.    0.   -1.
%!      -1.  0.    0.    0.    1.
%!      1.   0.    0.    0.    1.
%!      0.  -1.   -1.    0.    0.
%!      0.   1.   -1.    0.    0.
%!      0.  -1.    1.    0.    0.
%!      0.   1.    1.    0.    0.
%!      0.  -1.    0.   -1.    0.
%!      0.   1.    0.   -1.    0.
%!      0.  -1.    0.    1.    0.
%!      0.   1.    0.    1.    0.
%!      0.  -1.    0.    0.   -1.
%!      0.   1.    0.    0.   -1.
%!      0.  -1.    0.    0.    1.
%!      0.   1.    0.    0.    1.
%!      0.   0.   -1.   -1.    0.
%!      0.   0.    1.   -1.    0.
%!      0.   0.   -1.    1.    0.
%!      0.   0.    1.    1.    0.
%!      0.   0.   -1.    0.   -1.
%!      0.   0.    1.    0.   -1.
%!      0.   0.   -1.    0.    1.
%!      0.   0.    1.    0.    1.
%!      0.   0.    0.   -1.   -1.
%!      0.   0.    0.    1.   -1.
%!      0.   0.    0.   -1.    1.
%!      0.   0.    0.    1.    1.
%!      0.   0.    0.    0.    0.
%!      0.   0.    0.    0.    0.
%!      0.   0.    0.    0.    0.
%! ];
%! E = sortdesign(E);
%!
%! assert(C,E);
%!
%!# Print the number of experiments With variable point at center.
%!# Source: http://en.wikipedia.org/wiki/Box%E2%80%93Behnken_design
%! expected = [
%!     0
%!     0
%!     15
%!     27
%!     46	
%!     54	
%!     62	
%!     120	
%!     130	
%!     170	
%!     188	
%!     204	
%!     0
%!     0
%!     0
%!     396	
%! ];
%! for nbvar = 3 : 16
%!    H = scidoe_bbdesign(nbvar);
%!    m = size(H,"r");
%!    printf("nbvar=%d, Num. Experiments=%d (expected=%d)\n",..
%!    nbvar,m,expected(nbvar));
%!    if (nbvar==9 | nbvar==16) then
%!        // Skip this case:
%!        // TODO : fix this bug !
%!        continue
%!    end
%!    if (expected(nbvar)<>0) then
%!        assert_checkequal(m,expected(nbvar));
%!    end
%! end
%!
%!
%!# Print the number of experiments.
%!# With 1 point at center.
%!# Source: http://en.wikipedia.org/wiki/Box%E2%80%93Behnken_design
%! expected = [
%!     0
%!     0
%!     13	
%!     25	
%!     41	
%!     49	
%!     57	
%!    113
%!    121
%!    161
%!    177
%!    193
%!    0
%!    0
%!    0
%!    385
%! ];
%! for nbvar = 3 : 16
%!    H = scidoe_bbdesign(nbvar,"center",1);
%!    m = size(H,"r");
%!    printf("nbvar=%d, Num. Experiments=%d (expected=%d)\n",
%!              nbvar,m,expected(nbvar))
%!    if (nbvar==9 | nbvar==16) then
%!        // Skip this case:
%!        // TODO : fix this bug !
%!        continue
%!    end
%!    if (expected(nbvar)<>0) then
%!        assert_checkequal(m,expected(nbvar));
%!    end
%! end
%!
%!// Table from the R-rsm package.
%!// bbd(3), with the first row deleted, 
%!// so that the center point appears three times.
%! E = [
%!     0 -1 -1
%!     1 -1  0
%!     0  0  0
%!    -1  1  0
%!    -1  0  1
%!     0  0  0
%!     1  1  0
%!     1  0  1
%!     0  1  1
%!     0 -1  1
%!     0  1 -1
%!     1  0 -1
%!     0  0  0
%!    -1  0 -1
%!    -1 -1  0
%! ];
%! C = bbdesign(3);
%! E = sortdesign(E);
%!
%! assert(C,E);
%!
%# Table from the R-rsm package.
%#  bbd(4), with 9 zero center points  deleted, changing the number of rows from%# 
%# 36 down to 27.
%! E=[
%!    0  0  1 -1
%!    0  0 -1 -1
%!    1 -1  0  0
%!    0  0 -1  1
%!    0  0  0  0
%!   -1 -1  0  0
%!   -1  1  0  0
%!    1  1  0  0
%!    0  0  1  1
%!   -1  0  0  1
%!    1  0  0 -1
%!    1  0  0  1
%!    0  1 -1  0
%!   -1  0  0 -1
%!    0  1  1  0
%!    0 -1  1  0
%!    0 -1 -1  0
%!    0  1  0 -1
%!    0 -1  0  1
%!    0  1  0  1
%!    0  0  0  0
%!   -1  0 -1  0
%!    0  0  0  0
%!    1  0 -1  0
%!    1  0  1  0
%!   -1  0  1  0
%!    0 -1  0 -1
%! ];
%! C = bbdesign(4);
%! E = sortdesign(E);
%! assert(C,E);
%!# Table from the R-rsm package.
%!# bbd(5), with 2 zero center points 
%!# deleted, changing the number of rows from 
%!# 48 down to 46.
%! E = [
%!    -1  0 -1  0  0
%!     0  0  0 -1 -1
%!     0 -1  0  0 -1
%!     0  0  0  1 -1
%!     0  0  1 -1  0
%!     0  0  0  1  1
%!     0  1  0  0  1
%!     1 -1  0  0  0
%!    -1 -1  0  0  0
%!     1  0 -1  0  0
%!     1  1  0  0  0
%!     0  0 -1 -1  0
%!    -1  1  0  0  0
%!     0  0  0  0  0
%!     1  0  1  0  0
%!     0  0  0  0  0
%!    -1  0  1  0  0
%!     0  0  0 -1  1
%!     0 -1  0  0  1
%!     0  0 -1  1  0
%!     0  1  0  0 -1
%!     0  0  1  1  0
%!     0 -1  1  0  0
%!     0  0  0  0  0
%!     0  1  1  0  0
%!     1  0  0 -1  0
%!     0  0  1  0  1
%!     0  0  0  0  0
%!     0 -1  0  1  0
%!     0  0  1  0 -1
%!     0  0 -1  0  1
%!     1  0  0  1  0
%!     0 -1  0 -1  0
%!     1  0  0  0 -1
%!     0  1  0 -1  0
%!     1  0  0  0  1
%!     0  1  0  1  0
%!     0 -1 -1  0  0
%!     0  0  0  0  0
%!     0  0 -1  0 -1
%!    -1  0  0 -1  0
%!     0  0  0  0  0
%!     0  1 -1  0  0
%!    -1  0  0  0  1
%!    -1  0  0  0 -1
%!    -1  0  0  1  0
%! ];
%! C = bbdesign(5);
%! E = sortdesign(E);
%! assert(C,E);
%!
%!# Table from the R-rsm package.
%!# bbd(6), with two center points added, 
%!# changing the number of rows from 52 to 54
%! if (%f) 
%!    E = [
%!    0  0  0  0  0  0
%!    0  0  0  0  0  0
%!    1  0  0  0 -1 -1
%!    0  0  0  0  0  0
%!    0 -1 -1  0  1  0
%!    -1  0  1  0  0  1
%!    -1  0  0  1  1  0
%!    0  0  0  0  0  0
%!    1 -1  0  1  0  0
%!    -1  0  0  0 -1  1
%!    0 -1 -1  0 -1  0
%!    0  0  1 -1  0  1
%!    1  0 -1  0  0 -1
%!    -1  0  0  1 -1  0
%!    0  1  1  0 -1  0
%!    -1  0  0  0  1  1
%!    0  0 -1  1  0  1
%!    0  0 -1  1  0 -1
%!    -1  0 -1  0  0  1
%!    -1  1  0 -1  0  0
%!    1  0  0  0  1  1
%!    1  1  0 -1  0  0
%!    -1  0 -1  0  0 -1
%!    0  0  1 -1  0 -1
%!    1 -1  0 -1  0  0
%!    1  0  0  1 -1  0
%!    0 -1  1  0 -1  0
%!    -1  0  0  0 -1 -1
%!    1  0  1  0  0 -1
%!    0  0 -1 -1  0 -1
%!    0 -1  1  0  1  0
%!    0  0  1  1  0  1
%!    0  0  1  1  0 -1
%!    0  0  0  0  0  0
%!    0  0 -1 -1  0  1
%!    0  1 -1  0 -1  0
%!    1  0  0  0 -1  1
%!    -1 -1  0 -1  0  0
%!    -1  0  0  0  1 -1
%!    -1  0  1  0  0 -1
%!    1  0  0  0  1 -1
%!    1  0  0 -1  1  0
%!    1  1  0  1  0  0
%!    1  0 -1  0  0  1
%!    1  0  0  1  1  0
%!    0  1  1  0  1  0
%!    -1 -1  0  1  0  0
%!    0  1 -1  0  1  0
%!    -1  0  0 -1 -1  0
%!    1  0  0 -1 -1  0
%!    0  0  0  0  0  0
%!    1  0  1  0  0  1
%!    -1  0  0 -1  1  0
%!    -1  1  0  1  0  0
%!    ];
%!    C = bbdesign(6);
%!    E = sortdesign(E);
%!
%!    assert(C,E);
%!end
%!
%!# Table from the R-rsm package.
%!# bbd(7), with two center points added, 
%!# changing the number of rows from 60 to 62
%!E = [
%!    0  0  0  0  0  0  0
%!    0  0  0  0  0  0  0
%!    0  0  0  1 -1  1  0
%!    0 -1  0  0 -1  0  1
%!    1  0 -1  0  1  0  0
%!   -1  1  0  1  0  0  0
%!    1  0  0  0  0  1 -1
%!    0  0  0  1 -1 -1  0
%!    0  0  0  0  0  0  0
%!    1  0  0  0  0  1  1
%!   -1  1  0 -1  0  0  0
%!    1 -1  0 -1  0  0  0
%!    0  1 -1  0  0  1  0
%!    1  1  0 -1  0  0  0
%!    0  0 -1  1  0  0 -1
%!    0  0 -1 -1  0  0  1
%!    0  0  0 -1 -1  1  0
%!    1  0  1  0  1  0  0
%!    0  0  1  1  0  0  1
%!    0  0  0  0  0  0  0
%!   -1  0  1  0  1  0  0
%!    0  1  1  0  0  1  0
%!    0  0  0  0  0  0  0
%!    0 -1  0  0  1  0  1
%!    0 -1 -1  0  0 -1  0
%!    0  0 -1 -1  0  0 -1
%!   -1 -1  0 -1  0  0  0
%!    0  0  1  1  0  0 -1
%!   -1  0  0  0  0 -1 -1
%!    0 -1  0  0 -1  0 -1
%!    0 -1 -1  0  0  1  0
%!   -1  0 -1  0 -1  0  0
%!    0 -1  1  0  0 -1  0
%!    0  1  1  0  0 -1  0
%!    0  0  0 -1 -1 -1  0
%!    0  1  0  0  1  0  1
%!    0  0  0  0  0  0  0
%!   -1 -1  0  1  0  0  0
%!    0  0 -1  1  0  0  1
%!    1 -1  0  1  0  0  0
%!    1  0  1  0 -1  0  0
%!    0  0  0  1  1  1  0
%!    1  0  0  0  0 -1  1
%!    0  1  0  0 -1  0 -1
%!    0  0  0 -1  1  1  0
%!    0 -1  1  0  0  1  0
%!    0  1  0  0 -1  0  1
%!    0  0  1 -1  0  0  1
%!    0  0  0 -1  1 -1  0
%!   -1  0  0  0  0  1  1
%!   -1  0  1  0 -1  0  0
%!    0 -1  0  0  1  0 -1
%!    0  1 -1  0  0 -1  0
%!   -1  0  0  0  0 -1  1
%!    1  0 -1  0 -1  0  0
%!    0  0  0  1  1 -1  0
%!    1  0  0  0  0 -1 -1
%!    0  0  1 -1  0  0 -1
%!    0  1  0  0  1  0 -1
%!    1  1  0  1  0  0  0
%!   -1  0 -1  0  1  0  0
%!   -1  0  0  0  0  1 -1
%! ];
%! C = bbdesign(7);
%! E = sortdesign(E);
%! assert(C,E);
