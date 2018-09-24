=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut

#!/usr/bin/perl -w
eval unpack u=>q{_=7-E('-T<FEC=#L*"B`@("!M>2`H)&=T9BD@/2!`05)'5CL*("`@(`H@("`@;7D@)6%N;F\["B`@("!O<&5N_*$9)3$4L(B1G=&8B*2!\?"!D:64@(D%B;W)T:6YG+BX@0V%N)W0@;W!E;B`D9W1F(#H@)"%<;B(["B`@("!W_:&EL92AM>2`D;&EN93T\1DE,13XI>PH@("`@("`@(&-H;VUP("1L:6YE.PH@("`@("`@(&YE>'0@:68H)&QI_;F4]?B]>7",O*3L*("`@("`@("!M>2!`87)R87D@/2!S<&QI="@O7'0O+"1L:6YE*3L*("`@("`@("!N97AT_(&EF*"1A<G)A>5LR72!N92`B=')A;G-C<FEP="(I.PH@("`@("`@(&UY("@D8VAR+"1S=&%R="PD96YD+"1S_=')A;F0L)&YA;64I(#T@*"1A<G)A>5LP72PD87)R87E;,UTL)&%R<F%Y6S1=+"1A<G)A>5LV72PD87)R87E;_.%TI.PH@("`@("`@("1C:'(@/2`B8VAR(B`N("1C:'(@:68H)&-H<B%^+V-H<B\I.PH@("`@("`@("1C:'(]_?G,O8VAR350O8VAR32\["B`@("`@("`@;7D@)&ED(#T@)#$@:68H)&YA;64]?B]T<F%N<V-R:7!T7%]I9"!<_(BA<=RLI7")<.R\I.PH@("`@("`@(&UY("1G;B`]("1N86UE.PH@("`@("`@("1G;CU^<R\H+BHI9V5N95Q?_;F%M92!<(B\O.PH@("`@("`@("1G;CU^<R]<(EP[*"XJ*2\O.PH@("`@("`@("1G;CU^<R]G96YE7%]I9"!<_(B\O.PH@("`@("`@(&EF*"$D9VXI>PH@("`@("`@(`EP<FEN="`B0V%N)W0@9FEN9"!G96YE('-Y;6)O;"!F_;W(@:60@/2`D:61<;B(["B`@("`@("`@"65X:70["B`@("`@("`@?0H@("`@("`@("1A;FYO>R(D:60B?2`]_("(D9VY<="1S=')A;F0B.PH@("`@?0H@("`@8VQO<V4H1DE,12D["B`@("`*("`@(&UY("1O=71F;B`]("1G_=&8["B`@("`D;W5T9FX]?G,O7"YG=&8O7"Y04TDM<VEG;6%<+FUA<'!I;F=<+G1X="\["B`@("!O<&5N*$]5_5"PB/B1O=71F;B(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N("1O=71F;EQN(CL*"69O<F5A8V@@_;7D@)&ED*'-O<G0@:V5Y<R`E86YN;RE["B`@("`@("`@;F5X="!I9B@D:60@97$@(B(I.PH@("`@("`@(&EF_*"$D86YN;WLD:61]*7L*("`@("`@("`)<')I;G0@(B1I9"!H87,@;F\@;6%P<&EN9R!R97-U;'1S+EQN(CL*_("`@("`@("!]"B`@("`@("`@<')I;G0@3U54("(D:61<="(@+B`D86YN;WLD:61]("X@(EQN(CL*("`@('T*/("`@(&-L;W-E*$]55"D[}


