# Copyright (C) 2012-2013 - Michael Baudin
# Copyright (C) 2012 - Maria Christopoulou
# Copyright (C) 2009 - Yann Collette
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution.  The terms
# are also available at
# http:#www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function H = scidoe_ccdesign(varargin) # In progress
    # A Central Composite Design of Experiments
    #
    # Calling Sequence
    #     H = scidoe_ccdesign(nbvar)
    #     H = scidoe_ccdesign(nbvar,"Name",value,...)
    #
    # Parameters
    #     nbvar : a 1-by-1 matrix of doubles, positive integer, the number of variables of the experiment
    #     Name : a 1-by-1 matrix of strings. The options are "type", "alpha" and "center".
    #     value : the value
    # Description
    # This function produces a Central Composite Design (CCD) of experiments.
    # A CCD is composed by :
    # <itemizedlist>
    #   <listitem>
    #     <para>
    #       a factorial block design, 
    #     </para>
    #   </listitem>
    #   <listitem>
    #     <para>
    #       a block of axial (star) points and 
    #     </para>
    #   </listitem>
    #   <listitem>
    #     <para>
    #       a matrix of center points (zeros).
    #     </para>
    #   </listitem>
    # </itemizedlist>
    #
    # The available options are the following.
    # <itemizedlist>
    #   <listitem>
    #     <para>
    #       "type" : "circumscribed" (default), "inscribed" or "faced". 
    #     </para>
    #   </listitem>
    #   <listitem>
    #     <para>
    #       "alpha" : "orthogonal" (default) or "rotatable".
    #     </para>
    #   </listitem>
    #   <listitem>
    #     <para>
    #       "center" : a 1-by-2 row vector of doubles, the number of center points in each block of the design.
    #     </para>
    #   </listitem>
    # </itemizedlist>
    #
    #
    # The "circumscribed" and "inscribed" can be rotatable designs, 
    # but "faced" cannot. For "faced" CCD, alpha =1.
    # If "type" is specified, while "alpha" is not, then default value 
    # is "orthogonal".
    #
    # If the design is "orthogonal", 
    # <screen>
    # alpha = sqrt(nbvar*(1+(nao/na))/(1+(nco/nc))), 
    # </screen>
    # where 
    # <itemizedlist>
    #   <listitem>
    #     <para>
    #       nc : the number of factorial points, 
    #     </para>
    #   </listitem>
    #   <listitem>
    #     <para>
    #       nco : the number of center points added to the factorial design, 
    #     </para>
    #   </listitem>
    #   <listitem>
    #     <para>
    #       na : the number of axial points, 
    #     </para>
    #   </listitem>
    #   <listitem>
    #     <para>
    #       nao : the number of center points added to the axial points.
    #     </para>
    #   </listitem>
    # </itemizedlist>
    #
    # If the design is "rotatable", 
    # <screen>
    # alpha = nc^1/4, 
    # </screen>
    # where nc=2^nbvar are the factorial points.
    #
    # If the design is "faced", then alpha=1.
    #
    # The user can specify the number of center points in each block of the 
    # design (factorial and axial).
    # If "Parameter" is [2 3], then there are 2 center points in the 
    # factorial block and 3 points in the axial points block.
    # Default value is [4 4] (a total of 8 center points in the CCD).
    #
    # There is no strict order in the input of the above options.
    # For example, 
    # <screen>
    # scidoe_ccdesign(2,"type","inscribed","center",[2 3]) 
    # </screen>
    # will output 
    # the same as 
    # <screen>
    # scidoe_ccdesign(2,"center",[2 3],"type","inscribed")
    # </screen>
    #
    # Examples
    # # Circumscribed orthogonal design
    # H = scidoe_ccdesign(2)
    # # Inscribed orthogonal design
    # H = scidoe_ccdesign(3,"type","inscribed")
    # # Faced CCD
    # H = scidoe_ccdesign(2,"type","faced")
    # # Inscribed rotatable CC Design
    # H = scidoe_ccdesign(3,"type","inscribed","alpha","rotatable")
    # # Orthogonal CCD with 2 2 center points in the factorial block 
    # # and 3 center points in the star points block
    # H = scidoe_ccdesign(3,"center",[2 3])
    # 
    # Bibliography
    #    George E. P. Box, J. Stuart Hunter, William G. Hunter, "Statistics for experimenters Design,Innovation and Discovery", Second Edition, 2005
    #    http:#en.wikipedia.org/wiki/Central_composite_design
    #    http:#rgm2.lab.nig.ac.jp/RGM2/func.php?rd_id=rsm:ccd
    #    http:#www.mathworks.com/help/toolbox/stats/f56635.html
    #    http:#www.itl.nist.gov/div898/handbook/pri/section3/pri3361.htm
    #
    # Authors
    # Copyright (C) 2012-2013 - Michael Baudin
    # Copyright (C) 2012 - Maria Christopoulou
    # Copyright (C) 2009 - Yann Collette

    [lhs, rhs] = argn()
    apifun_checkrhs("scidoe_ccdesign",rhs,1:7)
    apifun_checklhs("scidoe_ccdesign",lhs,1)

    nbvar = varargin(1)
    # Check type, size, content of nbvar
    apifun_checktype("scidoe_ccdesign",nbvar,"nbvar",1,"constant")
    apifun_checkscalar("scidoe_ccdesign",nbvar,"nbvar",1)
    apifun_checkgreq("scidoe_ccdesign",nbvar,"nbvar",1,2)
    apifun_checkflint("scidoe_ccdesign",nbvar,"nbvar",1)
    #
    # Set default CCD, alpha value and number of center points
    default.type = "circumscribed"
    default.alpha = "orthogonal"
    default.center = [4 4];
    #
    options = apifun_keyvaluepairs(default,varargin(2:$))
    #
    typevalue = options.type
    alphavalue = options.alpha
    centervalue = options.center
    #
    apifun_checktype("scidoe_ccdesign",typevalue,"typekey",3,"string")
    apifun_checkoption("scidoe_ccdesign",typevalue,"typekey",3,["circumscribed" "inscribed" "faced"])
    apifun_checktype("scidoe_ccdesign",alphavalue,"alphakey",5,"string")
    apifun_checkoption("scidoe_ccdesign",alphavalue,"alphakey",5,["orthogonal" "rotatable"])
    apifun_checktype("scidoe_ccdesign",centervalue,"centerkey",7,"constant")

    # Orthogonal Design
    if (alphavalue == "orthogonal") 
        [H2,a] = scidoe_star(nbvar,"alpha","orthogonal","center",centervalue);
    end        
    #
    # Rotatable Design
    if (alphavalue == "rotatable") 
        [H2,a] = scidoe_star(nbvar,"alpha","rotatable")
    end 
    #
    #
    # Inscribed CCD
    if (typevalue == "inscribed") then
        H1 = 2*scidoe_ff2n(nbvar)-1;
        H1 = H1./a; # Scale down the factorial points
        H2 = scidoe_star(nbvar)
    end
    #
    # Faced CCD
    if (typevalue == "faced") then
        H2 = scidoe_star(nbvar); # Value of alpha is always 1 in Faced ccd
        H1 = 2*scidoe_ff2n(nbvar)-1;
    end
    #
    # Circumscribed Design
    if (typevalue == "circumscribed") then
        H1 = 2*scidoe_ff2n(nbvar)-1;
    end
    #
    # Center points
    C1 = zeros(centervalue(1),nbvar);
    C2 = zeros(centervalue(2),nbvar);
    #
    # Central Composite Design
    H=[H1;C1;H2;C2];

endfunction

# <-- JVM NOT MANDATORY -->
#
# The expected designs are generated with 
# the ccd function of the "rsm" package of R
#
# 
# Test with one parameter
B_comp = scidoe_ccdesign(2);
#
B_exp = [-1.000000 -1.000000
          1.000000 -1.000000
         -1.000000  1.000000
          1.000000  1.000000
          0.000000  0.000000
          0.000000  0.000000
          0.000000  0.000000
          0.000000  0.000000
         -1.414214  0.000000
          1.414214  0.000000
          0.000000 -1.414214
          0.000000  1.414214
          0.000000  0.000000
          0.000000  0.000000
          0.000000  0.000000
          0.000000  0.000000
          ];
B_comp=scidoe_sortdesign(B_comp);
B_exp=scidoe_sortdesign(B_exp);
assert_checkalmostequal(B_exp,B_comp,1.e-4)
#
# Test with one parameter
B_comp = scidoe_ccdesign(5);

B_exp = [
     -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
      1.000000 -1.000000 -1.000000 -1.000000 -1.000000
     -1.000000  1.000000 -1.000000 -1.000000 -1.000000
      1.000000  1.000000 -1.000000 -1.000000 -1.000000
     -1.000000 -1.000000  1.000000 -1.000000 -1.000000
      1.000000 -1.000000  1.000000 -1.000000 -1.000000
     -1.000000  1.000000  1.000000 -1.000000 -1.000000
      1.000000  1.000000  1.000000 -1.000000 -1.000000
     -1.000000 -1.000000 -1.000000  1.000000 -1.000000
      1.000000 -1.000000 -1.000000  1.000000 -1.000000
     -1.000000  1.000000 -1.000000  1.000000 -1.000000
      1.000000  1.000000 -1.000000  1.000000 -1.000000
     -1.000000 -1.000000  1.000000  1.000000 -1.000000
      1.000000 -1.000000  1.000000  1.000000 -1.000000
     -1.000000  1.000000  1.000000  1.000000 -1.000000
      1.000000  1.000000  1.000000  1.000000 -1.000000
     -1.000000 -1.000000 -1.000000 -1.000000  1.000000
      1.000000 -1.000000 -1.000000 -1.000000  1.000000
     -1.000000  1.000000 -1.000000 -1.000000  1.000000
      1.000000  1.000000 -1.000000 -1.000000  1.000000
     -1.000000 -1.000000  1.000000 -1.000000  1.000000
      1.000000 -1.000000  1.000000 -1.000000  1.000000
     -1.000000  1.000000  1.000000 -1.000000  1.000000
      1.000000  1.000000  1.000000 -1.000000  1.000000
     -1.000000 -1.000000 -1.000000  1.000000  1.000000
      1.000000 -1.000000 -1.000000  1.000000  1.000000
     -1.000000  1.000000 -1.000000  1.000000  1.000000
      1.000000  1.000000 -1.000000  1.000000  1.000000
     -1.000000 -1.000000  1.000000  1.000000  1.000000
      1.000000 -1.000000  1.000000  1.000000  1.000000
     -1.000000  1.000000  1.000000  1.000000  1.000000
      1.000000  1.000000  1.000000  1.000000  1.000000
      0.000000  0.000000  0.000000  0.000000  0.000000
      0.000000  0.000000  0.000000  0.000000  0.000000
      0.000000  0.000000  0.000000  0.000000  0.000000
      0.000000  0.000000  0.000000  0.000000  0.000000
     -2.494438  0.000000  0.000000  0.000000  0.000000
      2.494438  0.000000  0.000000  0.000000  0.000000
      0.000000 -2.494438  0.000000  0.000000  0.000000
      0.000000  2.494438  0.000000  0.000000  0.000000
      0.000000  0.000000 -2.494438  0.000000  0.000000
      0.000000  0.000000  2.494438  0.000000  0.000000
      0.000000  0.000000  0.000000 -2.494438  0.000000
      0.000000  0.000000  0.000000  2.494438  0.000000
      0.000000  0.000000  0.000000  0.000000 -2.494438
      0.000000  0.000000  0.000000  0.000000  2.494438
      0.000000  0.000000  0.000000  0.000000  0.000000
      0.000000  0.000000  0.000000  0.000000  0.000000
      0.000000  0.000000  0.000000  0.000000  0.000000
      0.000000  0.000000  0.000000  0.000000  0.000000
      ];
B_comp=scidoe_sortdesign(B_comp);
B_exp=scidoe_sortdesign(B_exp);
assert_checkalmostequal(B_comp,B_exp,1.e-4);
#
# Test with parameters
B_comp = scidoe_ccdesign(3,"type","inscribed");

B_exp = [-0.5477226 -0.5477226 -0.5477226
          0.5477226 -0.5477226 -0.5477226
         -0.5477226  0.5477226 -0.5477226
          0.5477226  0.5477226 -0.5477226
         -0.5477226 -0.5477226  0.5477226
          0.5477226 -0.5477226  0.5477226
         -0.5477226  0.5477226  0.5477226
          0.5477226  0.5477226  0.5477226
          0.0000000  0.0000000  0.0000000
          0.0000000  0.0000000  0.0000000
          0.0000000  0.0000000  0.0000000
          0.0000000  0.0000000  0.0000000
         -1.0000000  0.0000000  0.0000000
          1.0000000  0.0000000  0.0000000
          0.0000000 -1.0000000  0.0000000
          0.0000000  1.0000000  0.0000000
          0.0000000  0.0000000 -1.0000000
          0.0000000  0.0000000  1.0000000
          0.0000000  0.0000000  0.0000000
          0.0000000  0.0000000  0.0000000
          0.0000000  0.0000000  0.0000000
          0.0000000  0.0000000  0.0000000
          ];
B_comp=scidoe_sortdesign(B_comp);
B_exp=scidoe_sortdesign(B_exp);
assert_checkalmostequal(B_comp,B_exp,1.e-4);
#
# 
B_comp = scidoe_ccdesign(3,"type","faced");

B_exp = [-1 -1 -1
          1 -1 -1
         -1  1 -1
          1  1 -1
         -1 -1  1
          1 -1  1
         -1  1  1
          1  1  1
          0  0  0
          0  0  0
          0  0  0
          0  0  0
         -1  0  0
          1  0  0
          0 -1  0
          0  1  0
          0  0 -1
          0  0  1
          0  0  0
          0  0  0
          0  0  0
          0  0  0
          ];
B_comp=scidoe_sortdesign(B_comp);
B_exp=scidoe_sortdesign(B_exp);
assert_checkalmostequal(B_comp,B_exp,1.e-4);
#
# Test "alpha" parameter
B_comp = scidoe_ccdesign(3,"alpha","rotatable");

B_exp = [-1.000000 -1.000000 -1.000000
          1.000000 -1.000000 -1.000000
         -1.000000  1.000000 -1.000000
          1.000000  1.000000 -1.000000
         -1.000000 -1.000000  1.000000
          1.000000 -1.000000  1.000000
         -1.000000  1.000000  1.000000
          1.000000  1.000000  1.000000
          0.000000  0.000000  0.000000
          0.000000  0.000000  0.000000
          0.000000  0.000000  0.000000
          0.000000  0.000000  0.000000
         -1.681793  0.000000  0.000000
          1.681793  0.000000  0.000000
          0.000000 -1.681793  0.000000
          0.000000  1.681793  0.000000
          0.000000  0.000000 -1.681793
          0.000000  0.000000  1.681793
          0.000000  0.000000  0.000000
          0.000000  0.000000  0.000000
          0.000000  0.000000  0.000000
          0.000000  0.000000  0.000000];
B_comp=scidoe_sortdesign(B_comp);
B_exp=scidoe_sortdesign(B_exp);
assert_checkalmostequal(B_exp,B_comp,1.e-4);
#
#
# Orthogonal CCD with 2 center points in the factorial block and 3 center points in the star points block
# ccd(2,alpha='orthogonal',n0=c(2,3),randomize=FALSE)
B_exp=[
-1.000000 -1.000000
 1.000000 -1.000000
-1.000000  1.000000
 1.000000  1.000000
 0.000000  0.000000
 0.000000  0.000000
-1.527525  0.000000
 1.527525  0.000000
 0.000000 -1.527525
 0.000000  1.527525
 0.000000  0.000000
 0.000000  0.000000
 0.000000  0.000000
 ];
B_comp = scidoe_ccdesign(2,"alpha","orthogonal","center",[2 3]);
B_comp=scidoe_sortdesign(B_comp);
B_exp=scidoe_sortdesign(B_exp);
assert_checkalmostequal(B_exp,B_comp,1.0e-4)
#
#
# Rotatable and inscribed design
#ccd(2,alpha='rotatable',randomize=FALSE,inscribed=TRUE)
#
B_exp=[
-0.7071068 -0.7071068
 0.7071068 -0.7071068
-0.7071068  0.7071068
 0.7071068  0.7071068
 0.0000000  0.0000000
 0.0000000  0.0000000
 0.0000000  0.0000000
 0.0000000  0.0000000
-1.0000000  0.0000000
 1.0000000  0.0000000
 0.0000000 -1.0000000
 0.0000000  1.0000000
 0.0000000  0.0000000
 0.0000000  0.0000000
 0.0000000  0.0000000
 0.0000000  0.0000000];
B_comp=scidoe_ccdesign(2,"alpha","rotatable","type","inscribed")
B_comp=scidoe_sortdesign(B_comp);
B_exp=scidoe_sortdesign(B_exp);
assert_checkalmostequal(B_exp,B_comp,1.0e-4)