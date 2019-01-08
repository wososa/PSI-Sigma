=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w
eval unpack u=>q{_=7-E('-T<FEC=#L*(`EU<V4@4W1A=&ES=&EC<SHZ0F%S:6,@<7<H.F%L;"D["@D*("`@(&UY("@D9&(L)&)A_;2PD='EP92D@/2!`05)'5CL*("`@(`H@("`@;7D@)&=E;F]M92`](")(=6UA;B(["B`@("!I9B@D8F%M/7XO_4V5Q=6EN<R\I>PH@("`@"21G96YO;64@/2`B25,B.PH@("`@?0H@("`@;7D@)71M<#L*("`@(&UY("5U;FEQ_=64["B`@("!M>2`H)&=C:'(L)'-T87)T+"1E;F0I(#T@*"(M(BPP+#`I.PH@("`@;7D@)&5X;VYS.PH@("`@_;7D@)&UA>"`](#`["B`@("!M>2`E8F5D.PH@("`@"B`@("!M>2`D<W1A<G1T:6UE(#T@=&EM93L*"B`@("!M_>2`E8FEN.PH@("`@;7D@)&)I;G-I>F4@/2`Q,#`["B`@("!M>2!`8VAR;VUO<V]M97,["B`@("!M>2`D8VAR_9F]R;6%T(#T@,#L*:68H)&=E;F]M92!E<2`B25,B*7L*"7!U<V@H0&-H<F]M;W-O;65S+"))4R(I.PH))&-H_<F9O<FUA="`](#$["GUE;'-E>PH)9F]R*&UY("1I(#T@,3L@)&D@/#T@,C([("1I*RLI>PH)"7!U<V@H0&-H_<F]M;W-O;65S+"1I*3L*"7T*"7!U<V@H0&-H<F]M;W-O;65S+")8(BD["@EP=7-H*$!C:')O;6]S;VUE<RPB_62(I.PH)<'5S:"A`8VAR;VUO<V]M97,L(DTB*3L*"7!U<V@H0&-H<F]M;W-O;65S+")-5"(I.PH);W!E;BA)_3E!55"P@)RU\)RPB<V%M=&]O;',@=FEE=R`B("X@)&)A;2`N("(@?"!H96%D("UN(#$B*2!O<B!D:64@)"$[_"@EW:&EL92`H;7D@)&EN<'5T(#T@/$E.4%54/BD@>PH)"6-H;VUP("1I;G!U=#L*"0EM>2!`87)R87D@/2!S_<&QI="@O7'0O+"1I;G!U="D["@D);7D@*"1N86UE+"1C:'(L)&9L86<L)'-S+"1C:6=A<BD@/2`H)&%R<F%Y_6S!=+"1A<G)A>5LR72PD87)R87E;,5TL)&%R<F%Y6S-=+"1A<G)A>5LU72D["@D):68H)&-H<CU^+UYC:'(O_*7L*"0D))&-H<F9O<FUA="`](#$["@D)?0H)?0H)8VQO<V4H24Y0550I.PH)<')I;G0@(F-H<F9O<FUA="`]_("(@+B`D8VAR9F]R;6%T("X@(EQN(CL*?0H)<WES=&5M*")M:V1I<B!<7W=O<V]S871M<"(I.PH*9F]R96%C_:"!M>2`D:6-H<BA`8VAR;VUO<V]M97,I>PH);7D@)6%N;F\["B`@("!M>2`E:6YT<F]N<SL*("`@('!R:6YT_(")$;VEN9RXN+B!C:')O;6]S;VUE("1I8VAR7&XB.PH@("`@;W!E;BA&24Q%+"(D9&(B*2!\?"!D:64@(D%B_;W)T:6YG+BX@0V%N)W0@;W!E;B`D9&(@.B`D(5QN(CL*("`@('=H:6QE*&UY("1L:6YE/3Q&24Q%/BE["B`@_("`@("`@8VAO;7`@)&QI;F4["B`@("`@("`@;F5X="!I9B@D;&EN93U^+UY<(R\I.PH@("`@("`@(&YE>'0@_:68H)&QI;F4@97$@(B(I.PH)"6UY("@D8VAR+"1I,7,L)&DQ92PD:3)S+"1I,F4L)'1E<RPD=&5E+"1A;FYO_+"1A<RPD864L)&YA;64L)&=N*2`]('-P;&ET*"]<="\L)&QI;F4I.PH@("`@("`@("1C:'(]?G,O8VAR+R\[_"@D);F5X="!I9B@D8VAR(&YE("1I8VAR*3L*("`@("`@("!N97AT(&EF*"1N86UE(7XO7%]27%\O*3L*"0EM_>2`D8S$@/2`D87,@*R`Q.PH)"6UY("1C,B`]("1A92`M(#$["@D);7D@)&EN=&5R=F%L(#T@<W!R:6YT9B@B_)60B+"@D864M)&%S*2\U*3L*"0EM>2`D8S,@/2`H)&%S*R1I;G1E<G9A;"TQ*3L*"0EM>2`D8S0@/2`H)&%S_*R@D:6YT97)V86PJ,BDM,2D["@D);7D@)&,U(#T@*"1A<RLH)&EN=&5R=F%L*C,I+3$I.PH)"6UY("1L;V,@_/2`B)&DR<UQT)&DR92(["B`@("`@("`@)&%N;F][)&-H<GU[)&,Q?7LD;&]C?2LK.PH@("`@("`@("1A;FYO_>R1C:')]>R1C,GU[)&QO8WTK*SL*("`@("`@("`D86YN;WLD8VAR?7LD8S-]>R1L;V-]*RL["B`@("`@("`@_)&%N;F][)&-H<GU[)&,T?7LD;&]C?2LK.PH@("`@("`@("1A;FYO>R1C:')]>R1C-7U[)&QO8WTK*SL@("`*_("`@("`@("`D:6YT<F]N<WLD8VAR?7LD;&]C?7LD8S%](#T@,#L*("`@("`@("`D:6YT<F]N<WLD8VAR?7LD_;&]C?7LD8S)](#T@,#L*("`@("`@("`D:6YT<F]N<WLD8VAR?7LD;&]C?7LD8S-](#T@,#L*("`@("`@("`D_:6YT<F]N<WLD8VAR?7LD;&]C?7LD8S1](#T@,#L*("`@("`@("`D:6YT<F]N<WLD8VAR?7LD;&]C?7LD8S5]_(#T@,#L*("`@('T*("`@(&-L;W-E*$9)3$4I.PH*"6UY("1C;W5N="`](#`["@EM>2`D861D(#T@(B(["@DC_(T9O<B!);&QU;6EN82!S:&]R="!R96%D<PH):68H)'1Y<&4@/3T@,2E["@D))&%D9"`]("(M1B`R-38@+68@_,B`M<2`R-34B.PH)?0H):68H)&-H<F9O<FUA="`]/2`Q("8F("1I8VAR(7XO7F-H<B\I>PH)"21I8VAR(#T@_(F-H<B(@+B`D:6-H<CL*"7T*"6]P96XH24Y0550L("<M?"<L(G-A;71O;VQS('9I97<@(B`N("1A9&0@+B`B_("(@+B`D8F%M("X@(B`B("X@)&EC:'(@+B`B(BD@;W(@9&EE("0A("X@(EQR(B`N("(H97)R;W(I)&%D9"`D_8F%M("1I8VAR7'(B.PH)=VAI;&4@*&UY("1I;G!U="`](#Q)3E!55#XI('L*"0EC:&]M<"`D:6YP=70["@D)_;F5X="!I9B@D:6YP=70@97$@(B(I.PH)"21C;W5N="LK.PH)"6UY($!A<G)A>2`]('-P;&ET*"]<="\L)&EN_<'5T*3L*"0EM>2`H)&YA;64L)&-H<BPD9FQA9RPD<W,L)&-I9V%R*2`]("@D87)R87E;,%TL)&%R<F%Y6S)=_+"1A<G)A>5LQ72PD87)R87E;,UTL)&%R<F%Y6S5=*3L*"0EI9B@A)&-I9V%R*7L*"0D)<')I;G0@(G!R;V)L_96UA=&EC(&EN<'5T(#T@)&EN<'5T7&XB.PH)"0EN97AT.PH)"7T*"0DC(T9O<B!);&QU;6EN82!S:&]R="!R_96%D<PH)"6EF*"1T>7!E(#T](#$I>PH)"0EN97AT(&EF*"1C:6=A<CU^+UM.4T@J72\I.PH)"7T*"0DC(T9O_<B!-:6Y)3TX@;&]N9R!R96%D<PH)"6EF*"1T>7!E(#T](#(I>PH)"0EN97AT(&EF*"1C:6=A<B!E<2`B*B(I_.PH)"7T*"0DD8VAR/7YS+V-H<B\O.PH)"6UY($!#(#T@<W!L:70H+UM!+5I=+RPD8VEG87(I.PH)"6UY($!,_(#T@<W!L:70H+ULP+3E=*B\L)&-I9V%R*3L*"0EM>2`D;F,@/2!S8V%L87(@0$,["@D);7D@)&YL(#T@<V-A_;&%R($!,.PH)"0H)"6EF*"1T>7!E(#T](#$I>PH)"0EF;W(H;7D@)&D],#LD:3PD;F,[)&DK*RE["@D)"0EI_9B@A)$Q;)&DK,5TI>PH)"0D)"7!R:6YT(")I(#T@)&E<;B(["@D)"0D)<')I;G0@(F,@/2`D;F-<;B(["@D)_"0D)<')I;G0@(FP@/2`D;FQ<;B(["@D)"0D)<')I;G0@(D-;)&E=(#T@(B`N("1#6R1I72`N(")<;B(["@D)_"0D)<')I;G0@(DQ;)&E=(#T@(B`N("1,6R1I72`N(")<;B(["@D)"0D)<')I;G0@(F-I9V%R(#T@)&-I9V%R_7&XB.PH)"0D)"65X:70["@D)"0E]"@D)"0EI9B@D3%LD:2LQ72!N92`B32(I>PH)"0D)"21S<RL])$-;)&E=_(&EF*"1,6R1I*S%=(&5Q(").(B!\?"`D3%LD:2LQ72!E<2`B1"(I.PH)"0D)?65L<V5["@D)"0D);7D@)&YB_(#T@)$-;)&E=.PH)"0D)"69O<BAM>2`D:CTP.R1J/"1N8CLD:BLK*7L*"0D)"0D);7D@)'!O<R`]("1S<RLD_:CL*"0D)"0D):68H<V-A;&%R(&ME>7,@)7LD86YN;WLD8VAR?7LD<&]S?7T@/3T@,"E["@D)"0D)"0EN97AT_.PH)"0D)"0E]96QS97L*"0D)"0D)"69O<F5A8V@@;7D@)&QO8RAK97ES("5[("1A;FYO>R1C:')]>R1P;W-]_('TI>PH)"0D)"0D)"6YE>'0@:68H(21A;FYO>R1C:')]>R1P;W-]>R1L;V-]*3L*"0D)"0D)"0DD:6YT<F]N_<WLD8VAR?7LD;&]C?7LD<&]S?2LK.PH)"0D)"0D)?0H)"0D)"0E]"@D)"0D)?0H)"0D)"21S<RL])$-;)&E=_.PH)"0D)?0H)"0E]"@D)?0H)"6EF*"1T>7!E(#T](#(I>PH)"0EM>2`D;&]N9TX@/2`P.PH)"0EF;W(H;7D@_)&D],#LD:3PD;F,[)&DK*RE["@D)"0EI9B@D3%LD:2LQ72!E<2`B3B(I>PH)"0D)"6EF*"1#6R1I72`^(#4P_*7L*"0D)"0D))&QO;F=.(#T@,3L*"0D)"0D);&%S=#L*"0D)"0E]"@D)"0E]"@D)"7T*"0D);F5X="!I9B@D_;&]N9TX@/3T@,2D["@D)"69O<BAM>2`D:3TP.R1I/"1N8SLD:2LK*7L*"0D)"6EF*"$D3%LD:2LQ72E["@D)_"0D)<')I;G0@(FD@/2`D:5QN(CL*"0D)"0EP<FEN="`B8R`]("1N8UQN(CL*"0D)"0EP<FEN="`B;"`]("1N_;%QN(CL*"0D)"0EP<FEN="`B0ULD:5T@/2`B("X@)$-;)&E=("X@(EQN(CL*"0D)"0EP<FEN="`B3%LD:5T@_/2`B("X@)$Q;)&E=("X@(EQN(CL*"0D)"0EP<FEN="`B8VEG87(@/2`D8VEG87)<;B(["@D)"0D)97AI=#L*_"0D)"7T*"0D)"6EF*"1,6R1I*S%=(&5Q(").(B!\?"`D3%LD:2LQ72!E<2`B22(@?'P@)$Q;)&DK,5T@97$@_(E,B*7L*"0D)"0DD<W,K/21#6R1I72!I9B@D3%LD:2LQ72!E<2`B3B(I.PH)"0D)?65L<V5["@D)"0D);7D@_)&YB(#T@)$-;)&E=.PH)"0D)"69O<BAM>2`D:CTP.R1J/"1N8CLD:BLK*7L*"0D)"0D);7D@)'!O<R`]("1S_<RLD:CL*"0D)"0D):68H<V-A;&%R(&ME>7,@)7LD86YN;WLD8VAR?7LD<&]S?7T@/3T@,"E["@D)"0D)"0EN_97AT.PH)"0D)"0E]96QS97L*"0D)"0D)"69O<F5A8V@@;7D@)&QO8RAK97ES("5[("1A;FYO>R1C:')]>R1P_;W-]('TI>PH)"0D)"0D)"6YE>'0@:68H(21A;FYO>R1C:')]>R1P;W-]>R1L;V-]*3L*"0D)"0D)"0DD:6YT_<F]N<WLD8VAR?7LD;&]C?7LD<&]S?2LK.PH)"0D)"0D)?0H)"0D)"0E]"@D)"0D)?0H)"0D)"21S<RL])$-;_)&E=.PH)"0D)?0H)"0E]"@D)?0H)?0H)8VQO<V4H24Y0550I.PH*("`@(&EF*"1I8VAR(7XO7F-H<B\I>PH@_("`@"21I8VAR(#T@(F-H<B(@+B`D:6-H<CL*("`@('T*"6]P96XH3U54+"(^7W=O<V]S871M<"\B("X@)&)A_;2`N("(N(B`N("1I8VAR("X@(BYT>'0B*2!\?"!D:64@(D%B;W)T:6YG+BX@0V%N)W0@;W!E;B!O=71P=71T_;7!F:6QE+G1X="`Z("0A7&XB.PH)9F]R96%C:"!M>2`D8VAR*'-O<G0@:V5Y<R`E:6YT<F]N<RE["@D)9F]R_96%C:"!M>2`D;&]C*'-O<G0@:V5Y<R`E>R`D:6YT<F]N<WLD8VAR?2!]*7L*"0D);7D@*"1S<RPD964I(#T@_<W!L:70H+UQT+RPD;&]C*3L*"0D);7D@0'9A;'5E<SL*"0D);7D@)'!A<W,@/2`P.PH)"0EF;W)E86-H(&UY_("1P;W,H<V]R="!K97ES("5[("1I;G1R;VYS>R1C:')]>R1L;V-]('TI>PH)"0D);7D@)'8@/2`P.PH)"0D)_:68H(21I;G1R;VYS>R1C:')]>R1L;V-]>R1P;W-]*7L*"0D)"0DD=B`](#`["@D)"0E]96QS97L*"0D)"0DD_=B`]("1I;G1R;VYS>R1C:')]>R1L;V-]>R1P;W-].PH)"0D)"21P87-S*RL["@D)"0E]"0H)"0D)<'5S:"A`_=F%L=65S+"1V*3L*"0D)?0H)"0EN97AT(&EF*"1P87-S(#P@,R!\?"!S8V%L87(@0'9A;'5E<R`]/2`P*3L*_"0D);7D@)&UE9&EA;B`]('-P<FEN=&8H(B5D(BQM961I86XH0'9A;'5E<RDI.PH)"0EP<FEN="!/550@(B1C_:')<="1L;V-<="(@+B`D;65D:6%N("XB7&XB.PH@("`@"7T*("`@('T*("`@(&-L;W-E*$]55"D["GT*"6UY_("1O=71F;B`]("1B86T["B`@("`D;W5T9FX]?G,O7"Y!;&EG;F5D*"XJ*2\O.PH@("`@)&]U=&9N/7YS+UPN_8F%M*"XJ*2\O.PH@("`@)&]U=&9N("X]("(N25(N;W5T+G1A8B(["@ES>7-T96TH(F-A="!<7W=O<V]S871M_<%PO(B`N("1B86T@+B`B+F-H<EPJ+G1X="`^("(@+B`D;W5T9FXI.PH)<WES=&5M*")R;2!<7W=O<V]S871M_<%PO(B`N("1B86T@+B`B+F-H<EPJ+G1X="(I.PH@("`@;7D@)'-T;W!T:6UE(#T@=&EM93L*("`@(&UY("1M_:6YS(#T@*"1S=&]P=&EM92TD<W1A<G1T:6UE*2`O(#8P.PH@("`@<')I;G0@(E-P96YT("1M:6YS(&UI;G-<#;B([}
