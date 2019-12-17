=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w

eval unpack u=>q{_=7-E('-T<FEC=#L*(`EU<V4@4W1A=&ES=&EC<SHZ0F%S:6,@<7<H.F%L;"D["@D*("`@(&UY("@D9&(L)&)A_;2PD='EP92D@/2!`05)'5CL*("`@(`H@("`@)&)A;3U^<R\H+BHI7"\O+SL*("`@(`H@("`@;7D@)&=E;F]M_92`](")(=6UA;B(["B`@("!I9B@D8F%M/7XO4V5Q=6EN<R\I>PH@("`@"21G96YO;64@/2`B25,B.PH@("`@_?0H@("`@;7D@)71M<#L*("`@(&UY("5U;FEQ=64["B`@("!M>2`H)&=C:'(L)'-T87)T+"1E;F0I(#T@*"(M_(BPP+#`I.PH@("`@;7D@)&5X;VYS.PH@("`@;7D@)&UA>"`](#`["B`@("!M>2`E8F5D.PH@("`@"B`@("!M_>2`D<W1A<G1T:6UE(#T@=&EM93L*"B`@("!M>2`E8FEN.PH@("`@;7D@)&)I;G-I>F4@/2`Q,#`["B`@("!M_>2!`8VAR;VUO<V]M97,["B`@("!M>2`D8VAR9F]R;6%T(#T@,#L*("`@(&UY("1R96%D9F]R;6%T(#T@(G-I_;F=L92UE;F0B.PH@("`@"@EI9B@D9V5N;VUE(&5Q("))4R(I>PH)"7!U<V@H0&-H<F]M;W-O;65S+"))4R(I_.PH)"21C:')F;W)M870@/2`Q.PH)?65L<V5["@D)9F]R*&UY("1I(#T@,3L@)&D@/#T@,C([("1I*RLI>PH)_"0EP=7-H*$!C:')O;6]S;VUE<RPD:2D["@D)?0H)"7!U<V@H0&-H<F]M;W-O;65S+")8(BD["@D)<'5S:"A`_8VAR;VUO<V]M97,L(EDB*3L*"0EP=7-H*$!C:')O;6]S;VUE<RPB32(I.PH)"7!U<V@H0&-H<F]M;W-O;65S_+")-5"(I.PH)"6]P96XH24Y0550L("<M?"<L(G-A;71O;VQS('9I97<@(B`N("1B86T@+B`B('P@:&5A9"`M_;B`Q(BD@;W(@9&EE("0A.PH)"7=H:6QE("AM>2`D:6YP=70@/2`\24Y0550^*2!["@D)"6-H;VUP("1I;G!U_=#L*"0D);7D@0&%R<F%Y(#T@<W!L:70H+UQT+RPD:6YP=70I.PH)"0EM>2`H)&YA;64L)&-H<BPD9FQA9RPD_<W,L)&-I9V%R*2`]("@D87)R87E;,%TL)&%R<F%Y6S)=+"1A<G)A>5LQ72PD87)R87E;,UTL)&%R<F%Y6S5=_*3L*"0D):68H)&-H<CU^+UYC:'(O*7L*"0D)"21C:')F;W)M870@/2`Q.PH)"0E]"@D)?0H)"6-L;W-E*$E._4%54*3L*"0EP<FEN="`B8VAR;VUO<V]M92!F;W)M870@/2`B("X@)&-H<F9O<FUA="`N(")<;B(["@E]"@D*_"6]P96XH24Y0550L("(M?"(L(G-A;71O;VQS('9I97<@+68@,B`B("X@)&)A;2`N("(@?"!H96%D("UN(#4B_*2!O<B!D:64@)"$["@EW:&EL92`H;7D@)&EN<'5T(#T@/$E.4%54/BD@>PH)"6-H;VUP("1I;G!U=#L*"0DD_<F5A9&9O<FUA="`](")P86ER960M96YD(CL*"7T*"6-L;W-E*$E.4%54*3L*"7!R:6YT(")R96%D(&9O<FUA_="`]("(@+B`D<F5A9&9O<FUA="`N(")<;B(["@EI9B@D<F5A9&9O<FUA="!E<2`B<VEN9VQE+65N9"(I>PH)_"7!R:6YT(");5V%R;FYI;F==.B!4:&4@=7-E(&]F('-I;F=L92UE;F0@<VAO<G0@<F5A9',@:7,@;F]T(')E_8V]M;65N9&5D(&9O<B!A;'1E<FYA=&EV92!S<&QI8VEN9R!A;F%L>7-I<RY<;B([(`H)?0H)"@ES>7-T96TH_(FUK9&ER(%Q?=V]S;W-A=&UP(BD["@IF;W)E86-H(&UY("1I8VAR*$!C:')O;6]S;VUE<RE["@EM>2`E86YN_;SL*("`@(&UY("5I;G1R;VYS.PH@("`@<')I;G0@(D1O:6YG+BXN(&-H<F]M;W-O;64@)&EC:')<;B(["B`@_("!O<&5N*$9)3$4L(B1D8B(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N("1D8B`Z("0A7&XB.PH@_("`@=VAI;&4H;7D@)&QI;F4]/$9)3$4^*7L*("`@("`@("!C:&]M<"`D;&EN93L*("`@("`@("!N97AT(&EF_*"1L:6YE/7XO7EPC+RD["B`@("`@("`@;F5X="!I9B@D;&EN92!E<2`B(BD["@D);7D@*"1C:'(L)&DQ<RPD_:3%E+"1I,G,L)&DR92PD=&5S+"1T964L)&%N;F\L)&%S+"1A92PD;F%M92PD9VXI(#T@<W!L:70H+UQT+RPD_;&EN92D["B`@("`@("`@)&-H<CU^<R]C:'(O+SL*"0EN97AT(&EF*"1C:'(@;F4@)&EC:'(I.PH@("`@("`@_(&YE>'0@:68H)&YA;64A?B]<7U)<7R\I.PH)"6UY("1C,2`]("1A<R`K(#$["@D);7D@)&,R(#T@)&%E("T@_,3L*"0EM>2`D:6YT97)V86P@/2!S<')I;G1F*"(E9"(L*"1A92TD87,I+S4I.PH)"6UY("1C,R`]("@D87,K_)&EN=&5R=F%L+3$I.PH)"6UY("1C-"`]("@D87,K*"1I;G1E<G9A;"HR*2TQ*3L*"0EM>2`D8S4@/2`H)&%S_*R@D:6YT97)V86PJ,RDM,2D["@D);7D@)&QO8R`]("(D:3)S7'0D:3)E(CL*("`@("`@("`D86YN;WLD8VAR_?7LD8S%]>R1L;V-]*RL["B`@("`@("`@)&%N;F][)&-H<GU[)&,R?7LD;&]C?2LK.PH@("`@("`@("1A;FYO_>R1C:')]>R1C,WU[)&QO8WTK*SL*("`@("`@("`D86YN;WLD8VAR?7LD8S1]>R1L;V-]*RL["B`@("`@("`@_)&%N;F][)&-H<GU[)&,U?7LD;&]C?2LK.R`@(`H@("`@("`@("1I;G1R;VYS>R1C:')]>R1L;V-]>R1C,7T@_/2`P.PH@("`@("`@("1I;G1R;VYS>R1C:')]>R1L;V-]>R1C,GT@/2`P.PH@("`@("`@("1I;G1R;VYS>R1C_:')]>R1L;V-]>R1C,WT@/2`P.PH@("`@("`@("1I;G1R;VYS>R1C:')]>R1L;V-]>R1C-'T@/2`P.PH@("`@_("`@("1I;G1R;VYS>R1C:')]>R1L;V-]>R1C-7T@/2`P.PH@("`@?0H@("`@8VQO<V4H1DE,12D["@H);7D@_)&-O=6YT(#T@,#L*"6UY("1A9&0@/2`B(CL*"2,C1F]R($EL;'5M:6YA('-H;W)T(')E861S"@EI9B@D='EP_92`]/2`Q*7L*"0DD861D(#T@(BU&(#(U-B`M9B`R("UQ(#(U-2(@:68H)')E861F;W)M870@97$@(G!A:7)E_9"UE;F0B*3L*"0DD861D(#T@(BU&(#(U-B`M<2`R-34B(&EF*"1R96%D9F]R;6%T(&5Q(")S:6YG;&4M96YD_(BD["@E]"@EI9B@D8VAR9F]R;6%T(#T](#$@)B8@)&EC:'(A?B]>8VAR+RE["@D))&EC:'(@/2`B8VAR(B`N_("1I8VAR.PH)?0H):68H)&-H<F9O<FUA="`]/2`P*7L*"0DD:6-H<CU^<R]C:'(O+SL*"7T*"6]P96XH24Y0_550L("<M?"<L(G-A;71O;VQS('9I97<@(B`N("1A9&0@+B`B("(@+B`D8F%M("X@(B`B("X@)&EC:'(@+B`B_(BD@;W(@9&EE("0A("X@(EQR(B`N("(H97)R;W(I)&%D9"`D8F%M("1I8VAR7'(B.PH)=VAI;&4@*&UY("1I_;G!U="`](#Q)3E!55#XI('L*"0EC:&]M<"`D:6YP=70["@D);F5X="!I9B@D:6YP=70@97$@(B(I.PH)"21C_;W5N="LK.PH)"6UY($!A<G)A>2`]('-P;&ET*"]<="\L)&EN<'5T*3L*"0EM>2`H)&YA;64L)&-H<BPD9FQA_9RPD<W,L)&-I9V%R*2`]("@D87)R87E;,%TL)&%R<F%Y6S)=+"1A<G)A>5LQ72PD87)R87E;,UTL)&%R<F%Y_6S5=*3L*"0EI9B@A)&-I9V%R*7L*"0D)<')I;G0@(G!R;V)L96UA=&EC(&EN<'5T(#T@)&EN<'5T7&XB.PH)_"0EN97AT.PH)"7T*"0DC(T9O<B!);&QU;6EN82!S:&]R="!R96%D<PH)"6EF*"1T>7!E(#T](#$I>PH)"0EN_97AT(&EF*"1C:6=A<CU^+UM.4T@J72\I.PH)"7T*"0DC(T9O<B!-:6Y)3TX@;&]N9R!R96%D<PH)"6EF*"1T_>7!E(#T](#(I>PH)"0EN97AT(&EF*"1C:6=A<B!E<2`B*B(I.PH)"7T*"0DD8VAR/7YS+V-H<B\O.PH)"6UY_($!#(#T@<W!L:70H+UM!+5I=+RPD8VEG87(I.PH)"6UY($!,(#T@<W!L:70H+ULP+3E=*B\L)&-I9V%R*3L*_"0EM>2`D;F,@/2!S8V%L87(@0$,["@D);7D@)&YL(#T@<V-A;&%R($!,.PH)"0H)"6EF*"1T>7!E(#T](#$I_>PH)"0EF;W(H;7D@)&D],#LD:3PD;F,[)&DK*RE["@D)"0EI9B@A)$Q;)&DK,5TI>PH)"0D)"7!R:6YT(")I_(#T@)&E<;B(["@D)"0D)<')I;G0@(F,@/2`D;F-<;B(["@D)"0D)<')I;G0@(FP@/2`D;FQ<;B(["@D)"0D)_<')I;G0@(D-;)&E=(#T@(B`N("1#6R1I72`N(")<;B(["@D)"0D)<')I;G0@(DQ;)&E=(#T@(B`N("1,6R1I_72`N(")<;B(["@D)"0D)<')I;G0@(F-I9V%R(#T@)&-I9V%R7&XB.PH)"0D)"65X:70["@D)"0E]"@D)"0EI_9B@D3%LD:2LQ72!N92`B32(I>PH)"0D)"21S<RL])$-;)&E=(&EF*"1,6R1I*S%=(&5Q(").(B!\?"`D3%LD_:2LQ72!E<2`B1"(I.PH)"0D)?65L<V5["@D)"0D);7D@)&YB(#T@)$-;)&E=.PH)"0D)"69O<BAM>2`D:CTP_.R1J/"1N8CLD:BLK*7L*"0D)"0D);7D@)'!O<R`]("1S<RLD:CL*"0D)"0D):68H<V-A;&%R(&ME>7,@)7LD_86YN;WLD8VAR?7LD<&]S?7T@/3T@,"E["@D)"0D)"0EN97AT.PH)"0D)"0E]96QS97L*"0D)"0D)"69O<F5A_8V@@;7D@)&QO8RAK97ES("5[("1A;FYO>R1C:')]>R1P;W-]('TI>PH)"0D)"0D)"6YE>'0@:68H(21A;FYO_>R1C:')]>R1P;W-]>R1L;V-]*3L*"0D)"0D)"0DD:6YT<F]N<WLD8VAR?7LD;&]C?7LD<&]S?2LK.PH)"0D)_"0D)?0H)"0D)"0E]"@D)"0D)?0H)"0D)"21S<RL])$-;)&E=.PH)"0D)?0H)"0E]"@D)?0H)"6EF*"1T>7!E_(#T](#(I>PH)"0EM>2`D;&]N9TX@/2`P.PH)"0EF;W(H;7D@)&D],#LD:3PD;F,[)&DK*RE["@D)"0EI9B@D_3%LD:2LQ72!E<2`B3B(I>PH)"0D)"6EF*"1#6R1I72`^(#4P*7L*"0D)"0D))&QO;F=.(#T@,3L*"0D)"0D)_;&%S=#L*"0D)"0E]"@D)"0E]"@D)"7T*"0D);F5X="!I9B@D;&]N9TX@/3T@,2D["@D)"69O<BAM>2`D:3TP_.R1I/"1N8SLD:2LK*7L*"0D)"6EF*"$D3%LD:2LQ72E["@D)"0D)<')I;G0@(FD@/2`D:5QN(CL*"0D)"0EP_<FEN="`B8R`]("1N8UQN(CL*"0D)"0EP<FEN="`B;"`]("1N;%QN(CL*"0D)"0EP<FEN="`B0ULD:5T@/2`B_("X@)$-;)&E=("X@(EQN(CL*"0D)"0EP<FEN="`B3%LD:5T@/2`B("X@)$Q;)&E=("X@(EQN(CL*"0D)"0EP_<FEN="`B8VEG87(@/2`D8VEG87)<;B(["@D)"0D)97AI=#L*"0D)"7T*"0D)"6EF*"1,6R1I*S%=(&5Q(")._(B!\?"`D3%LD:2LQ72!E<2`B22(@?'P@)$Q;)&DK,5T@97$@(E,B*7L*"0D)"0DD<W,K/21#6R1I72!I9B@D_3%LD:2LQ72!E<2`B3B(I.PH)"0D)?65L<V5["@D)"0D);7D@)&YB(#T@)$-;)&E=.PH)"0D)"69O<BAM>2`D_:CTP.R1J/"1N8CLD:BLK*7L*"0D)"0D);7D@)'!O<R`]("1S<RLD:CL*"0D)"0D):68H<V-A;&%R(&ME>7,@_)7LD86YN;WLD8VAR?7LD<&]S?7T@/3T@,"E["@D)"0D)"0EN97AT.PH)"0D)"0E]96QS97L*"0D)"0D)"69O_<F5A8V@@;7D@)&QO8RAK97ES("5[("1A;FYO>R1C:')]>R1P;W-]('TI>PH)"0D)"0D)"6YE>'0@:68H(21A_;FYO>R1C:')]>R1P;W-]>R1L;V-]*3L*"0D)"0D)"0DD:6YT<F]N<WLD8VAR?7LD;&]C?7LD<&]S?2LK.PH)_"0D)"0D)?0H)"0D)"0E]"@D)"0D)?0H)"0D)"21S<RL])$-;)&E=.PH)"0D)?0H)"0E]"@D)?0H)?0H)8VQO_<V4H24Y0550I.PH*("`@(&EF*"1I8VAR(7XO7F-H<B\I>PH@("`@"21I8VAR(#T@(F-H<B(@+B`D:6-H<CL*_("`@('T*("`@(`H);W!E;BA/550L(CY?=V]S;W-A=&UP+R(@+B`D8F%M("X@(BXB("X@)&EC:'(@+B`B+G1X_="(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N(&]U='!U='1M<&9I;&4N='AT(#H@)"%<;B(["@EF_;W)E86-H(&UY("1C:'(H<V]R="!K97ES("5I;G1R;VYS*7L*"0EF;W)E86-H(&UY("1L;V,H<V]R="!K97ES_("5[("1I;G1R;VYS>R1C:')]('TI>PH)"0EM>2`H)'-S+"1E92D@/2!S<&QI="@O7'0O+"1L;V,I.PH)"0EM_>2!`=F%L=65S.PH)"0EM>2`D<&%S<R`](#`["@D)"69O<F5A8V@@;7D@)'!O<RAS;W)T(&ME>7,@)7L@)&EN_=')O;G-[)&-H<GU[)&QO8WT@?2E["@D)"0EM>2`D=B`](#`["@D)"0EI9B@A)&EN=')O;G-[)&-H<GU[)&QO_8WU[)'!O<WTI>PH)"0D)"21V(#T@,#L*"0D)"7UE;'-E>PH)"0D)"21V(#T@)&EN=')O;G-[)&-H<GU[)&QO_8WU[)'!O<WT["@D)"0D))'!A<W,K*SL*"0D)"7T)"@D)"0EP=7-H*$!V86QU97,L)'8I.PH)"0E]"@D)"6YE_>'0@:68H)'!A<W,@/"`S('Q\('-C86QA<B!`=F%L=65S(#T](#`I.PH)"0EM>2`D;65D:6%N(#T@<W!R:6YT_9B@B)60B+&UE9&EA;BA`=F%L=65S*2D["@D)"7!R:6YT($]55"`B)&-H<EQT)&QO8UQT(B`N("1M961I86X@_+B)<;B(["B`@("`)?0H@("`@?0H@("`@8VQO<V4H3U54*3L*?0H);7D@)&]U=&9N(#T@)&)A;3L*("`@("1O_=71F;CU^<R]<+D%L:6=N960H+BHI+R\["B`@("`D;W5T9FX]?G,O7"YB86TH+BHI+R\["B`@("`D;W5T9FX@_+CT@(BY)4BYO=70N=&%B(CL*"7-Y<W1E;2@B8V%T(%Q?=V]S;W-A=&UP7"\B("X@)&)A;2`N("(N8VAR7"HN_='AT(#X@(B`N("1O=71F;BD["@ES>7-T96TH(G)M(%Q?=V]S;W-A=&UP7"\B("X@)&)A;2`N("(N8VAR7"HN_='AT(BD["B`@("!M>2`D<W1O<'1I;64@/2!T:6UE.PH@("`@;7D@)&UI;G,@/2`H)'-T;W!T:6UE+21S=&%RL='1I;64I("\@-C`["B`@("!P<FEN="`B4W!E;G0@)&UI;G,@;6EN<UQN(CL};

print $@;
