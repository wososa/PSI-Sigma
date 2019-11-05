=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w
eval unpack u=>q{_=7-E('-T<FEC=#L*"B`@("!M>2`H)&=T9BPD<6-H<BPD;F]V96QJ=6YC=&EO;F-R:71E<FEA*2`]($!!4D=6_.PH@("`@"B`@("!M>2`H)&=C:'(L)'-T87)T+"1E;F0I(#T@*"(M(BPP+#`I.PH@("`@;7D@)&5X;VYS.PH@_("`@;7D@)&-O=6YT(#T@,#L*("`@(&UY("1M87@@/2`P.PH@("`@"B`@("!M>2`D<W1A<G1T:6UE(#T@=&EM_93L*"B`@("!M>2`E86YN;SL*("`@(&UY("5E>&]N86YN;SL*("`@(&]P96XH1DE,12PB)&=T9B(I('Q\(&1I_92`B06)O<G1I;F<N+B!#86XG="!O<&5N("1G=&8@.B`D(5QN(CL*("`@('=H:6QE*&UY("1L:6YE/3Q&24Q%_/BE["B`@("`@("`@8VAO;7`@)&QI;F4["B`@("`@("`@;F5X="!I9B@D;&EN93U^+UY<(R\I.PH@("`@("`@_(&UY($!A<G)A>2`]('-P;&ET*"]<="\L)&QI;F4I.PH@("`@("`@(&YE>'0@:68H)&%R<F%Y6S)=(&YE(")G_96YE(B`F)B`D87)R87E;,ET@;F4@(F5X;VXB*3L*("`@("`@("!M>2`H)&-H<BPD8V%T+"1S=&%R="PD96YD_+"1S=')A;F0L)&YA;64I(#T@*"1A<G)A>5LP72PD87)R87E;,5TL)&%R<F%Y6S-=+"1A<G)A>5LT72PD87)R_87E;-ETL)&%R<F%Y6SA=*3L*("`@("`@("`D8VAR(#T@(F-H<B(@+B`D8VAR(&EF*"1C:'(A?B]C:'(O*3L*_("`@("`@("!N97AT(&EF*"1C:'(@;F4@)'%C:'(I.PH@("`@("`@(",C<F5M;W9E('1H:7,@9F]R('1O;6%T_;PH@("`@("`@(",D;F%M93U^<R]<+BA<9"LI7")<.R]<(EP[+V<["B`@("`@("`@:68H)&%R<F%Y6S)=(&5Q_(")G96YE(BE["B`@("`@("`@"21N86UE/7YS+R@N*BEG96YE7%]N86UE(%PB+R\["B`@("`@("`@"21N86UE_/7YS+UPB7#L@*"XJ*2\O.PH@("`@("`@(`DD;F%M93U^<R]G96YE7%]I9"!<(B\O.PH@("`@("`@(`DD86YN_;WLD8VAR?7LB)'-T87)T7'0D96YD(GTK*SL*("`@("`@("`))&%N;F][(B1C:'(B?7LB)'-T87)T7'0D96YD_(GT@/2`D;F%M93L*("`@("`@("!]"B`@("`@("`@:68H)&%R<F%Y6S)=(&5Q(")E>&]N(BE["@D)"6UY("1%_3E-4(#T@)&YA;64["B`)"0DD14Y35#U^<R\H+BHI=')A;G-C<FEP=%Q?:60@7"(O+SL*("`@(`D))$5.4U0]_?G,O7")<.R`H+BHI+R\["@D)"6EF*"1%3E-4/7XO7%\O*7L*"0D)"21%3E-4/7YS+UQ?+UPN+V<["@D)"7T*_("`@("`@("`):68H)&YA;64]?B]T<F%N<V-R:7!T7V)I;W1Y<&4O*7L*("`@("`@("`)"21N86UE/7YS+R@N_*BET<F%N<V-R:7!T7V)I;W1Y<&4@7"(O+SL*("`@("`@("`)"21N86UE/7YS+UPB7#L@*"XJ*2\O.PH@("`@_("`@(`D))&-A="`]("1N86UE.PH@("`@("`@(`E]"B`@("`@("`@"6EF*"1N86UE/7XO=')A;G-C<FEP=%]T_>7!E+RE["@D)"0DD;F%M93U^<R\H+BHI=')A;G-C<FEP=%]T>7!E(%PB+R\["B`@("`@("`@"0DD;F%M93U^_<R]<(EP[("@N*BDO+SL*("`@("`@("`)"21C870@/2`D;F%M93L*("`@("`@("`)?0H@("`@("`@(`EM>2`D_;&%B96P@/2`M,3L*("`@("`@("`))&QA8F5L(#T@,2!I9B@D8V%T(&5Q(")N;VYS96YS95]M961I871E9%]D_96-A>2(I.PH@("`@("`@(`DC)&5X;VYA;FYO>R1C:')]>R(D<W1A<G1<="1E;F0B?7LD14Y35'TK*SL*("`@_("`@("`):68H(21E>&]N86YN;WLB)&-H<B)]>R(D<W1A<G1<="1E;F0B?7LD14Y35'TI>PH@("`@("`@(`D)_)&5X;VYA;FYO>R(D8VAR(GU[(B1S=&%R=%QT)&5N9")]>R1%3E-4?2`]("1L86)E;#L*("`@("`@("`)?65L_<V5["B`@("`@("`@"0DD97AO;F%N;F][(B1C:'(B?7LB)'-T87)T7'0D96YD(GU[)$5.4U1](#T@)&QA8F5L_(&EF*"1L86)E;"`]/2`Q*3L*("`@("`@("`)?0H@("`@("`@('T)"B`@("!]"B`@("!C;&]S92A&24Q%*3L*_"B`@("!M>2!`9FEL97,["B`);W!E;BA&24Q%+")G<F]U<&$N='AT(BD@?'P@9&EE(")!8F]R=&EN9RXN($-A_;B=T(&]P96X@9W)O=7!A+G1X="`Z("0A7&XB.PH@("`@=VAI;&4H;7D@)&QI;F4]/$9)3$4^*7L*("`@("`@_("!C:&]M<"`D;&EN93L*("`@("`@("!N97AT(&EF*"1L:6YE(&5Q("(B*3L*("`@("`@("!M>2`D86-C97-S_:6]N(#T@)&QI;F4["B`@("`@("`@)&%C8V5S<VEO;CU^<R]<+D%L:6=N961<+G-O<G1E9$)Y0V]O<F1<+F]U_=%PN8F%M+R\["@D))&%C8V5S<VEO;CU^<R]<+G-O<G1E9%PN;W5T7"YB86TO+SL*"0DD86-C97-S:6]N/7YS_+UPN8F%M+R\["@D))&%C8V5S<VEO;CU^<R]<+B0O+SL*"0DD86-C97-S:6]N+CT@(BY32BYO=70N=&%B(CL*_("`@("`@("!P=7-H*$!F:6QE<RPD86-C97-S:6]N*3L*("`@("`@("`*("`@('T*("`@(&-L;W-E*$9)3$4I_.PH@"6]P96XH1DE,12PB9W)O=7!B+G1X="(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N(&=R;W5P_8BYT>'0@.B`D(5QN(CL*("`@('=H:6QE*&UY("1L:6YE/3Q&24Q%/BE["B`@("`@("`@8VAO;7`@)&QI;F4[_"B`@("`@("`@;F5X="!I9B@D;&EN92!E<2`B(BD["B`@("`@("`@;7D@)&%C8V5S<VEO;B`]("1L:6YE.PH@_("`@("`@("1A8V-E<W-I;VX]?G,O7"Y!;&EG;F5D7"YS;W)T961">4-O;W)D7"YO=71<+F)A;2\O.PH)"21A_8V-E<W-I;VX]?G,O7"YS;W)T961<+F]U=%PN8F%M+R\["@D))&%C8V5S<VEO;CU^<R]<+F)A;2\O.PH)"21A_8V-E<W-I;VX]?G,O7"XD+R\["@D))&%C8V5S<VEO;BX]("(N4THN;W5T+G1A8B(["B`@("`@("`@<'5S:"A`_9FEL97,L)&%C8V5S<VEO;BD["B`@("!]"B`@("!C;&]S92A&24Q%*3L*("`*("`@(&UY("532CL*("`@(&UY_("1C;W5N='-J(#T@,#L*("`@(&9O<F5A8V@@;7D@)&IF;BA`9FEL97,I>PH@("`@"6YE>'0@:68H+7H@)&IF_;BD["B`@("`)(W!R:6YT("(D:F9N+BXN7&XB.PH@("`@"6UY("1C;W5N="`](#`["@D);W!E;BA&24Q%+"`B_)&IF;B(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N("1J9FY<;B(["B`@("`@("`@=VAI;&4H;7D@_)&QI;F4]/$9)3$4^*7L*("`@("`@("`)8VAO;7`@)&QI;F4["B`@("`@("`@"6UY($!A<G)A>2`]('-P;&ET_*"]<="\L)&QI;F4I.PH@("`@("`@("`@("!M>2`H)&-H<BPD<W1A<G0L)&5N9"PD;G5M*2`]("@D87)R87E;_,%TL)&%R<F%Y6S%=+"1A<G)A>5LR72PD87)R87E;-ETI.PH@("`@("`@("`@("`D8VAR(#T@(F-H<B(@+B`D_8VAR(&EF*"1C:'(A?B]C:'(O*3L*("`@("`@("`@("`@;F5X="!I9B@D8VAR(&YE("1Q8VAR*3L*("`@("`@_("`@("`@;F5X="!I9B@D;G5M(#P@)&YO=F5L:G5N8W1I;VYC<FET97)I82D["B`@("`@("`@("`@("132GLD_8VAR?7LB)'-T87)T7'0D96YD(GTK*SL*("`@("`@("`@("`@)&-O=6YT*RL["B`@("`@("`@?0H@("`@("`@_(&-L;W-E*$9)3$4I.PH@("`@("`@("1C;W5N='-J*RL["B`@("`@("`@(W!R:6YT(")N=6UB97(@;V8@=F%L_:60@:G5N8W1I;VX@/2`D8V]U;G1<;B(["B`@("!]"B`@("!P<FEN="`B3G5M8F5R(&]F('9A;&ED("Y32BYO_=70N=&%B(&9I;&5S(#T@)&-O=6YT<VI<;B(["@H@("`@;7D@)71M<#L*("`@(&UY("5B960["B`@("!M>2`E_9&(["B`@("!M>2`E=6YI<75E.PH@("`@;7D@)'1O=&%L7V5X;VX@/2`P.PH@("`@;W!E;BA&24Q%+"(D9W1F_(BD@?'P@9&EE(")!8F]R=&EN9RXN($-A;B=T(&]P96X@)&=T9B`Z("0A7&XB.PH@("`@=VAI;&4H;7D@)&QI_;F4]/$9)3$4^*7L*("`@("`@("!C:&]M<"`D;&EN93L*("`@("`@("!M>2`D;G5M(#T@<V-A;&%R(&ME>7,@_)71M<#L*("`@("`@("`C:68H96]F*7L*("`@("`@("`C?65L<V5["B`@("`@("`@(PEN97AT(&EF*"1L:6YE_(7XO7B1Q8VAR+RD["B`@("`@("`@(WT*("`@("`@("!N97AT(&EF*"1L:6YE(7XO7'1E>&]N7'0O("8F("1L_:6YE(7XO7'1T<F%N<V-R:7!T7'0O*3L*("`@("`@("!M>2!`87)R87D@/2!S<&QI="@O7'0O+"1L:6YE*3L*_("`@("`@("!M>2`H)&-H<BPD8V%T+"1S<RPD964L)'-T<F%N9"PD86YN;RD@/2`H)&%R<F%Y6S!=+"1A<G)A_>5LR72PD87)R87E;,UTL)&%R<F%Y6S1=+"1A<G)A>5LV72PD87)R87E;.%TI.PH)"21C:'(@/2`B8VAR(B`N_("1C:'(@:68H)&-H<B%^+V-H<B\I.PH)"6YE>'0@:68H)&-H<B!N92`D<6-H<BD["B`@("`@("`@)&=C:'(@_/2`D8VAR(&EF*"1G8VAR(&5Q("(M(BD["B`@("`@("`@:68H)'-T87)T(#T](#`@)B8@)&5N9"`]/2`P("8F_("1C870@97$@(G1R86YS8W)I<'0B*7L*("`@("`@("`))&=C:'(@/2`D8VAR.PH@("`@("`@(`DD<W1A<G0@_/2`D<W,["B`@("`@("`@"21E;F0@/2`D964["B`@("`@("`@?0H*"0EI9B@D8V%T(&5Q(")E>&]N(B`F)B`D_9V-H<B!E<2`D8VAR*7L*"0D);7D@)'1I9"`]("(B.PH)"0EM>2`D=&ED(#T@)&%N;F\["B`)"0DD=&ED/7YS_+R@N*BET<F%N<V-R:7!T7%]I9"!<(B\O.PH@("`@"0DD=&ED/7YS+UPB7#L@*"XJ*2\O.PH)"0EI9B@D=&ED_/7XO7%\O*7L*"0D)"21T:60]?G,O7%\O7"XO9SL*"0D)?0H@("`@("`@"0EP=7-H*$![)'1M<'LD=&ED?7TL_(B1C:')<="1S<UQT)&5E(BD["B`@("`@("`)"6EF*"$D=6YI<75E>R(D8VAR7'0D<W-<="1E92)]*7L*("`@_("`@(`D)"21U;FEQ=65[(B1C:')<="1S<UQT)&5E(GTK*SL*("`@("`@(`D)"21T;W1A;%]E>&]N*RL["B`@_("`@("`)"7T*"0E]"@H@("`@("`@(&EF*"1C870@97$@(G1R86YS8W)I<'0B('Q\(&5O9BE["B`@("`@("`@_"6EF*"1S<R`^("1E;F0@?'P@96]F*7L*("`@("`@("`)"21C;W5N="LK.PH@("`@("`@(`D);7D@)&YU;2`]_('-C86QA<B!K97ES("5T;7`["@D)"0EI9B@D;G5M(#X](#$I>PH)"0D)"21M87@@/2`D;G5M(&EF*"1N=6T@_/B`D;6%X*3L*"0D)"0EM>2`D:60@/2`B)&=C:'(B("X@(E\B("X@(B1S=&%R="(@+B`B7R(@+B`B)&5N9"([_"@D)"0D)<G5N*"1I9"P@)'1O=&%L7V5X;VXL("5T;7`I.PH)"0D)?0H@("`@("`@(`D))71M<"`]("@I.PH@_("`@("`@(`D))75N:7%U92`]("@I.PH@("`@("`@(`D))&=C:'(@/2`D8VAR.PH@("`@("`@(`D))'-T87)T_(#T@)'-S.PH@("`@("`@(`D))&5N9"`]("1E93L*("`@("`@("`)"21T;W1A;%]E>&]N(#T@,#L*("`@("`@_("`)?65L<V5["B`@("`@("`@"0EI9B@D964@/B`D96YD("8F("1G8VAR(&5Q("1C:'(I>PH@("`@("`@(`D)_"21E;F0@/2`D964["@D)"0E]"@D)"7T*"0D);F5X=#L*"0E]"B`@("`@("`@"@H)?0H*("`@(&UY("1S=&]P_=&EM92`]('1I;64["B`@("!M>2`D;6EN<R`]('-P<FEN=&8H(B4N,F8B+"@H)'-T;W!T:6UE+21S=&%R='1I_;64I("\@-C`I*3L*("`@(`H);W!E;BA/550L(CXD<6-H<BYB960B*2!\?"!D:64@(D%B;W)T:6YG+BX@0V%N_)W0@;W!E;B`D<6-H<BYB960@.B`D(5QN(CL*"69O<F5A8V@@;7D@)&)E9&]U='!U="AS;W)T(&ME>7,@)6)E_9"E["@D)<')I;G0@3U54("1B961[)&)E9&]U='!U='T@+B`B7&XB.PH)?0H)8VQO<V4H3U54*3L*"0H);W!E_;BA/550L(CXD<6-H<BYD8B(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N)'%C:'(N9&(@.B`D(5QN_(CL*"69O<F5A8V@@;7D@)&1B;W5T<'5T*'-O<G0@:V5Y<R`E9&(I>PH)"7!R:6YT($]55"`D9&)[)&1B;W5T_<'5T?2`N(")<;B(["@E]"@EC;&]S92A/550I.PH@("`@"B`@("!P<FEN="`B=&EM92`]("1M:6YS(&UI;G-<_;B(["@EO<&5N*$]55"PB/CYT:6UE+G1X="(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N('1I;64N_='AT(#H@)"%<;B(["@EP<FEN="!/550@(B1Q8VAR7'0B("X@*"1M:6YS*2`N("(@;6EN<UQN(CL*"6-L;W-E_*$]55"D["@IS=6(@<G5N,@I["@EP<FEN="`B<G5N,EQN(CL*"7)E='5R;B`P.PI]"G-U8B!R=6X*>PH);7D@_*"1I9"PD=&]T86Q?97AO;BPE:6YP=70I(#T@0%\["@EM>2`H)&EC:'(L)&ES=&%R="PD:65N9"D@/2!S<&QI_="@O7%\O+"1I9"D["@H);7D@)'-T87)T=&EM97@@/2!T:6UE.PH);7D@)'1E<W1C;W5N="`]('-C86QA<B!K_97ES("5I;G!U=#L*"@EI9B@D=&5S=&-O=6YT(#X@,3`P,#`I>PH)"6UY("1S:&]W8V]U;G0@/2`P.PH)"69O_<F5A8V@@;7D@)'1I9"AS;W)T(&ME>7,@)6EN<'5T*7L*"0D)<')I;G0@(B1T:61<;B(["@D)"21S:&]W8V]U_;G0K*SL*"0D);&%S="!I9B@D<VAO=V-O=6YT(#X](#4I.PH)"7T*"7T*"0H)(VUY("1F;B`](")C;W5N=%]E_>&]N+G1X="(["@DC:68H)'1O=&%L7V5X;VX@/B`U,#`I>PH)(PEO<&5N*$]55"PB/CXD9FXB*2!\?"!D:64@_(D%B;W)T:6YG+BX@0V%N)W0@;W!E;B`D9FX@.B`D(5QN(CL*"2,)<')I;G0@3U54("(D:61<="1T;W1A;%]E_>&]N7&XB.PH)(PEC;&]S92A/550I.PH)(WT*"@EM>2`E:6YT<F]N.PH)"@EM>2`E:VYO=VYS<SL*"6UY("5K_;F]W;F5S.PH);7D@)6MN;W=N97AO;G,["@D)"@EF;W)E86-H(&UY("1T:60H<V]R="!K97ES("5I;G!U="E[_"@D);F5X="!I9B@D=&ED(&5Q("(B*3L*"0EM>2!`87)R87D@/2!`>R1I;G!U='LD=&ED?7T["@D);7D@)&YU_;2`]('-C86QA<B!`87)R87D["@D)9F]R*&UY("1I(#T@,3LD:2`\("1N=6T[)&DK*RE["@D)"6UY("@D8VAR_,2PD<W,Q+"1E93$I(#T@<W!L:70H+UQT+RPD87)R87E;)&DM,5TI.PH)"0EM>2`H)&-H<C(L)'-S,BPD964R_*2`]('-P;&ET*"]<="\L)&%R<F%Y6R1I72D["@D)"21K;F]W;F5X;VYS>R(D<W,Q7'0D964Q(GTK*SL*"0D)_)&MN;W=N97AO;G-[(B1S<S)<="1E93(B?2LK.PH)"0DD:VYO=VYS<WLH)&5E,2LQ*7U[(B1S<S(L)&5E,B)]_*RL["@D)"21K;F]W;F5S>R@D<W,R+3$I?7LB)'-S,2PD964Q(GTK*SL*"0D))&EN=')O;GLH)&5E,2LQ*2`N_(")<="(@+B`H)'-S,BTQ*7TK*SL*"0E]"@E]"@D*"2-!9&0@;75T=6%L;'D@97AC;'5S:79E(&EN=')O;G,*_"6UY($!K;F]W;G-S.PH)9F]R96%C:"!M>2`D<VET92AS;W)T(&ME>7,@)6MN;W=N<W,I>PH)"6YE>'0@:68H_<V-A;&%R(&ME>7,@)7LD:VYO=VYS<WLD<VET97U](#T](#$I.PH)"7!U<V@H0&MN;W=N<W,L)'-I=&4I.PH)_?0H);7D@0&MN;W=N97,["@EF;W)E86-H(&UY("1S:71E*'-O<G0@:V5Y<R`E:VYO=VYE<RE["@D);F5X="!I_9BAS8V%L87(@:V5Y<R`E>R1K;F]W;F5S>R1S:71E?7T@/3T@,2D["@D)<'5S:"A`:VYO=VYE<RPD<VET92D[_"@E]"@EF;W)E86-H(&UY("1S<RA`:VYO=VYS<RE["@D)9F]R96%C:"!M>2`D97,H0&MN;W=N97,I>PH)"0EM_>2`D8FEN9"`](#`["@D)"69O<F5A8V@@;7D@)&)I;F1E>&]N,2AK97ES("5[)&MN;W=N<W-[)'-S?7TI>PH)_"0D)9F]R96%C:"!M>2`D8FEN9&5X;VXR*&ME>7,@)7LD:VYO=VYE<WLD97-]?2E["@D)"0D):68H)&)I;F1E_>&]N,2!E<2`D8FEN9&5X;VXR*7L*"0D)"0D))&)I;F0K*SL*"0D)"0E]"@D)"0E]"@D)"0EL87-T(&EF*"1B_:6YD(#X@,2D["@D)"7T*"0D))&EN=')O;GLD<W,@+B`B7'0B("X@)&5S?2LK(&EF*"1B:6YD(#X@,2D["@D)_?0H)?0H*"2,C(%-E96L@9F]R(&YO=F5L(&5X;VYS"@EM>2`E;F]V96QI;G1R;VX["@EM>2`E4TIS:71E<SL*_"69O<F5A8V@@;7D@)'-J:6YT<F]N*'-O<G0@:V5Y<R`E>R`D4TI[)&EC:')]('TI>PH)"6UY("@D:7,L)&EE_*2`]('-P;&ET*"]<="\L)'-J:6YT<F]N*3L*"0EN97AT(&EF*"1I<R`\("1I<W1A<G0@?'P@)&EE(#X@)&EE_;F0I.PH)"0H)"6EF*"$D:6YT<F]N>R1S:FEN=')O;GTI>PH)"0DD;F]V96QI;G1R;VY[)'-J:6YT<F]N?2LK_.PH)"0DD:6YT<F]N>R1S:FEN=')O;GTK*SL*"0E]96QS97L*"0D);F5X=#L*"0E]"@E]"@H);7D@)6-O;&QE_8W1I;VX["@EM>2`E;&%S=#L*"69O<F5A8V@@;7D@)&EN=')O;BAS;W)T(&ME>7,@)6EN=')O;BE["@D);7D@_*"1I<RPD:64I(#T@<W!L:70H+UQT+RPD:6YT<F]N*3L*"0EM>2`H)'1A<F=E=&5X;VXL)&5X;VYA;FYO*3L*_"0EF;W)E86-H(&UY("1T:60H<V]R="!K97ES("5I;G!U="E["@D)"6YE>'0@:68H)'1I9"!E<2`B(BD["@D)_"6UY($!A<G)A>2`]($![)&EN<'5T>R1T:61]?3L*"0D);7D@)&YU;2`]('-C86QA<B!`87)R87D["@D)"6UY_("@D;&%S=&-H<BPD;&%S='-S+"1L87-T964I(#T@<W!L:70H+UQT+RPD87)R87E;)&YU;2TQ72D["@D)"6EF_*"@D;&%S='-S*2`\("1I<R`F)B`H)&QA<W1E92D@/B`D:64I>PH)"0D))'1A<F=E=&5X;VX@/2`H)&QA<W1S_<RD@+B`B7'0B("X@*"1L87-T964I.PH)"0D))&5X;VYA;FYO(#T@)&5X;VYA;FYO>R(D;&%S=&-H<B)]>R(D_;&%S='-S7'0D;&%S=&5E(GU[)'1I9'T["@D)"0DD=&%R9V5T97AO;B`]("@D;&%S='-S+21I<RD@+B`B7'0B_("X@*"1L87-T964M)&ES*3L*"0D)"21C;VQL96-T:6]N>R(D;&%S=&-H<EQT)&ES7'0D:65<=%)<="1T:60B_?2`N/2`B+"(@+B`H)&ES+21I<RD@+B`B7'0B("X@*"1I<RTD:7,I("X@(EQT(B`N("@D:7,M)&ES*2`N(")<_="(@+B`H)&EE+21I<RD@+B`B7'0B("X@)'1A<F=E=&5X;VX@+B`B7'0B("X@)&5X;VYA;FYO.PH)"0D);&%S_=#L*"0D)?0H)"0EN97AT(&EF*"1L87-T964@/"`D:7,I.PH)"0EF;W(H;7D@)&D@/2`P.R1I(#P@*"1N=6TM_,2D[)&DK*RE["@D)"0EM>2`H)&-H<C(L)'-S,BPD964R*2`]('-P;&ET*"]<="\L)&%R<F%Y6R1I72D["@D)_"0EM>2`H)&-H<C,L)'-S,RPD964S*2`]('-P;&ET*"]<="\L)&%R<F%Y6R1I*S%=*3L*"0D)"0H)"0D);&%S_="!I9B@D:7,@/3T@*"1E93(K,2D@)B8@)&EE(#T]("@D<W,S+3$I*3L*"0D)"6YE>'0@:68H)&ES(#X]("@D_964S*S$I*3L*"0D)"6QA<W0@:68H)&EE(#P]("@D<W,R+3$I*3L*"0D)"0H)"0D):68H(21N;W9E;&EN=')O_;GLD:6YT<F]N?2E["@D)"0E]96QS97L*"0D)"6EF*"@D964R*S$I(#T]("1I<R`F)B`D:64@/"`H)'-S,RTQ_*2E["@D)"0D);7D@)'1M<&EN=')O;B`]("@D964R*S$I("X@(EQT(B`N("@D<W,S+3$I.PH)"0D)"21T87)G_971E>&]N(#T@*"1I92TH)&5E,BLQ*2D@+B`B7'0B("X@*"1E93,M*"1E93(K,2DI.PH)"0D)"21E>&]N86YN_;R`](#(["@D)"0D))&-O;&QE8W1I;VY[(B1I8VAR7'0D=&UP:6YT<F]N7'137'0D=&ED(GT@+CT@(BPB("X@_*"1I<RTH)&5E,BLQ*2D@+B`B7'0B("X@*"1I<RTH)&5E,BLQ*2D@+B`B7'0B("X@*"1I<RTH)&5E,BLQ*2D@_+B`B7'0B("X@*"1I92TH)&5E,BLQ*2D@+B`B7'0B("X@)'1A<F=E=&5X;VX@+B`B7'0B("X@)&5X;VYA;FYO_.PH)"0D)"6QA<W0["@D)"0E]"@D)"0EI9B@H)&5E,BLQ*2`\("1I<R`F)B`D:64@/3T@*"1S<S,M,2DI>PH)_"0D)"6UY("1T;7!I;G1R;VX@/2`H)&5E,BLQ*2`N(")<="(@+B`H)'-S,RTQ*3L*"0D)"0DD=&%R9V5T97AO_;B`]("@D<W,R+2@D964R*S$I*2`N(")<="(@+B`H)&ES+2@D964R*S$I*3L*"0D)"0DD97AO;F%N;F\@/2`R_.PH)"0D)"21C;VQL96-T:6]N>R(D:6-H<EQT)'1M<&EN=')O;EQT4UQT)'1I9")]("X]("(L(B`N("@D:7,M_*"1E93(K,2DI("X@(EQT(B`N("@D:64M*"1E93(K,2DI("X@(EQT(B`N("@D:64M*"1E93(K,2DI("X@(EQT_(B`N("@D:64M*"1E93(K,2DI("X@(EQT(B`N("1T87)G971E>&]N("X@(EQT(B`N("1E>&]N86YN;SL*"0D)_"0EL87-T.PH)"0D)?0H)"0D)?0H*"0D)"2-&;W(@:VYO=VX@86QT97)N871I=F4@<W!L:6-E('-I=&4@979E_;G0N"@D)"0EI9B@H)&5E,BLQ*2`]/2`D:7,@)B8@)&EE(#X@*"1S<S,M,2D@)B8@)&EE(#P@*"1E93,K,2DI_>PH)"0D)"21T87)G971E>&]N(#T@*"1S<S,I("X@(EQT(B`N("@D964S*3L*"0D)"0DD97AO;F%N;F\@/2`D_97AO;F%N;F][)&-H<C-]>R1T87)G971E>&]N?7LD=&ED?3L*"0D)"0DD=&%R9V5T97AO;B`]("@D<W,S+21I_<RD@+B`B7'0B("X@*"1E93,M)&ES*3L*"0D)"0DD8V]L;&5C=&EO;GLB)&EC:')<="1I;G1R;VY<=%-<="1T_:60B?2`N/2`B+"(@+B`H)&ES+21I<RD@+B`B7'0B("X@*"1I<RTD:7,I("X@(EQT(B`N("@D:7,M)&ES*2`N_(")<="(@+B`H*"1S<S,M,2DM)&ES*2`N(")<="(@+B`D=&%R9V5T97AO;B`N(")<="(@+B`D97AO;F%N;F\[_"@D)"0D);&%S=#L*"0D)"7T*"0D)"6EF*"@D964R*S$I(#X@)&ES("8F("1I92`]/2`H)'-S,RTQ*2`F)B`D_:7,@/B`H)'-S,BTQ*2E["@D)"0D))'1A<F=E=&5X;VX@/2`H)'-S,BD@+B`B7'0B("X@*"1E93(I.PH)"0D)_"21E>&]N86YN;R`]("1E>&]N86YN;WLD8VAR,GU[)'1A<F=E=&5X;VY]>R1T:61].PH)"0D)"21T87)G971E_>&]N(#T@*"1S<S(M)&ES*2`N(")<="(@+B`H)&5E,BTD:7,I.PH)"0D)"21C;VQL96-T:6]N>R(D:6-H<EQT_)&EN=')O;EQT4UQT)'1I9")]("X]("(L(B`N("@H)&5E,BLQ*2TD:7,I("X@(EQT(B`N("@D:64M)&ES*2`N_(")<="(@+B`H)&EE+21I<RD@+B`B7'0B("X@*"1I92TD:7,I("X@(EQT(B`N("1T87)G971E>&]N("X@(EQT_(B`N("1E>&]N86YN;SL*"0D)"0EL87-T.PH)"0D)?0H*"0D)"2-&;W(@:6YT<F]N(')E=&5N=&EO;B!E=F5N_="X*"0D)"6EF*"@D<W,R+3$I(#P@)&ES("8F("1I92`\("@D964R*S$I*7L*"0D)"0DD=&%R9V5T97AO;B`]_("@D<W,R*2`N(")<="(@+B`H)&5E,BD["@D)"0D))&5X;VYA;FYO(#T@)&5X;VYA;FYO>R1C:'(R?7LD=&%R_9V5T97AO;GU[)'1I9'T["@D)"0D))'1A<F=E=&5X;VX@/2`H)'-S,BTD:7,I("X@(EQT(B`N("@D964R+21I_<RD["@D)"0D))&-O;&QE8W1I;VY[(B1I8VAR7'0D:7-<="1I95QT4EQT)'1I9")]("X]("(L(B`N("@D:7,M_)&ES*2`N(")<="(@+B`H)&ES+21I<RD@+B`B7'0B("X@*"1I<RTD:7,I("X@(EQT(B`N("@D:64M)&ES*2`N_(")<="(@+B`D=&%R9V5T97AO;B`N(")<="(@+B`D97AO;F%N;F\["@D)"0D);&%S=#L*"0D)"7T*"0D)"0H)_"0D)(V5X;VX@<VMI<'!I;F<@9F]R(#$@;W(@;75L=&EP;&4*"0D)"6YE>'0@:68H)&D@/3T@,"D["@D)"0EM_>2`H)&-H<C$L)'-S,2PD964Q*2`]('-P;&ET*"]<="\L)&%R<F%Y6R1I+3%=*3L*"0D)"6QA<W0@:68H)&ES_(#T]("@D964Q*S$I("8F("1I92`]/2`H)'-S,BTQ*2D["@D)"0EI9B@D:7,@/#T@*"1E93$K,2D@)B8@)&EE_(#X]("@D<W,S+3$I*7L*"0D)"0EI9B@A)&-O;&QE8W1I;VY[(B1I8VAR7'0D:6YT<F]N7'177'0D=&ED(GTI_>PH)"0D)"0EN97AT(&EF*"@D964Q*S$I("$]("1I<RD["@D)"0D)"21T87)G971E>&]N(#T@*"1S<S(I("X@_(EQT(B`N("@D964R*3L*"0D)"0D))&5X;VYA;FYO(#T@)&5X;VYA;FYO>R1C:'(R?7LD=&%R9V5T97AO;GU[_)'1I9'T["@D)"0D)"21T87)G971E>&]N(#T@*"1S<S(M)&ES*2`N(")<="(@+B`H)&5E,BTD:7,I.PH)"0D)_"0EI9B@A)&5X;VYA;FYO*7L*"0D)"0D)"7!R:6YT("(D=&ED(%LD8VAR,BPD<W,R+"1E93)=(&-A;B=T(&9I_;F0@97AO;B!T>7!E+EQN(CL*"0D)"0D)"65X:70["@D)"0D)"7T*"0D)"0D))&-O;&QE8W1I;VY[(B1I8VAR_7'0D:6YT<F]N7'177'0D=&ED(GT@+CT@(BPB("X@*"1E93$K,2TD:7,I("X@(EQT(B`N("@D<W,R+3$M)&ES_*2`N(")<="(@+B`H)&5E,BLQ+21I<RD@+B`B7'0B("X@*"1S<S,M,2TD:7,I("X@(EQT(B`N("1T87)G971E_>&]N("X@(EQT(B`N("1E>&]N86YN;SL*"0D)"0E]96QS97L*"0D)"0D))'1A<F=E=&5X;VX@/2`H)'-S,BD@_+B`B7'0B("X@*"1E93(I.PH)"0D)"0DD97AO;F%N;F\@/2`D97AO;F%N;F][(B1C:'(R(GU[)'1A<F=E=&5X_;VY]>R1T:61].PH)"0D)"0DD=&%R9V5T97AO;B`]("@D<W,R+21I<RD@+B`B7'0B("X@*"1E93(M)&ES*3L*_"0D)"0D):68H(21E>&]N86YN;RE["@D)"0D)"0EP<FEN="`B)'1I9"!;)&-H<C(L)'-S,BPD964R72!C86XG_="!F:6YD(&5X;VX@='EP92Y<;B(["@D)"0D)"0EE>&ET.PH)"0D)"0E]"@D)"0D)"21C;VQL96-T:6]N>R(D_:6-H<EQT)&EN=')O;EQT5UQT)'1I9")]("X]("(L(B`N("@D964Q*S$M)&ES*2`N(")<="(@+B`H)'-S,BTQ_+21I<RD@+B`B7'0B("X@*"1E93(K,2TD:7,I("X@(EQT(B`N("@D<W,S+3$M)&ES*2`N(")<="(@+B`D=&%R_9V5T97AO;B`N(")<="(@+B`D97AO;F%N;F\["@D)"0D)?0H)"0D)?0H*"0D)?0H)"7T*"7T)"@D*"0DC5'=O_(%,@979E;G1S(&-A;B!C<F5A=&4@82!N;W9E;"!E>&]N(&EF('1H97D@:&%V92!T:&4@<V%M92!C;VYS=&ET_=71I=F4@97AO;G,*"0EM>2`E;F]V96QE>&]N.PH)"6UY("1N;W9E;"`](#`["@D);7D@)&-R96%T92`](#`[_"@D)9F]R96%C:"!M>2`D979E;G0Q*&ME>7,@)6-O;&QE8W1I;VXI>PH)"0EN97AT(&EF*"1E=F5N=#$A?B]<_=%-<="\I.PH)"0EM>2`D:G5N8W1I;VYS(#T@)&-O;&QE8W1I;VY[)&5V96YT,7T["@D)"21J=6YC=&EO;G,]_?G,O7"PO+SL*"0D);7D@0&IU;F-T:6]N<R`]('-P;&ET*"]<+"\L)&IU;F-T:6]N<RD["@D)"6YE>'0@:68H_<V-A;&%R($!J=6YC=&EO;G,@/3T@,2D["@D)"6UY("1N=6T@/2!S8V%L87(@0&IU;F-T:6]N<SL*"0D)(TQE_9G0*"0D);7D@*"1J<S$Q+"1J93$Q+"1J<S$R+"1J93$R*2`]('-P;&ET*"]<="\L)&IU;F-T:6]N<ULP72D[_"@D)"21J<S$R(#T@)&IS,3$@:68H)&IS,3(@/3T@)&IE,3(I.PH)"0DD:F4Q,B`]("1J93$Q(&EF*"1J<S$R_(#T]("1J93$R*3L*"0D)(U)I9VAT"@D)"6UY("@D:G,R,2PD:F4R,2PD:G,R,BPD:F4R,BD@/2!S<&QI="@O_7'0O+"1J=6YC=&EO;G-;*"1N=6TM,2E=*3L*"0D)"@D)"6YE>'0@:68H)&IS,3(@/3T@)&IS,C(@?'P@)&IE_,3(@/3T@)&IE,C(I.PH)"0EM>2`H)&-H<C$L)&ES,2PD:64Q+"1T>7!E,2PD=&ED,2D@/2!S<&QI="@O7'0O_+"1E=F5N=#$I.PH)"0DD;F]V96QE>&]N>R1T:60Q?7LH)&IE,3(K,2D@+B`B7'0B("X@*"1J<S(R+3$I?2LK_.PH)"0EU;F1E9B`D8V]L;&5C=&EO;GLD979E;G0Q?3L*"0E]"@D)"0D*"69O<F5A8V@@;7D@)'1I9"AK97ES_("5N;W9E;&5X;VXI>PH)"2-P<FEN="`B8VAE8VMI;F<N+BX@;F]V96P@97AO;B!F;W(@)'1I9%QN(CL*"0EM_>2`H)'1A<F=E=&5X;VXL)&5X;VYA;FYO*3L*"0EF;W)E86-H(&UY("1L;V,H:V5Y<R`E>R`D;F]V96QE>&]N_>R1T:61]('TI>PH)"0EM>2`H)'1E<RPD=&5E*2`]('-P;&ET*"]<="\L)&QO8RD["@D)"6UY($!A<G)A>2`]_($![)&EN<'5T>R1T:61]?3L*"0D);7D@)&YU;2`]('-C86QA<B!`87)R87D["@D)"69O<BAM>2`D:2`](#`[_)&D@/"`H)&YU;2TQ*3LD:2LK*7L*"0D)"6UY("@D8VAR,BPD<W,R+"1E93(I(#T@<W!L:70H+UQT+RPD87)R_87E;)&E=*3L*"0D)"6YE>'0@:68H)'-S,B`^("1T97,I.PH)"0D);7D@*"1C:'(S+"1S<S,L)&5E,RD@/2!S_<&QI="@O7'0O+"1A<G)A>5LD:2LQ72D["@D)"0EI9B@D964R(#P@)'1E<R`F)B`D<W,S(#X@)'1E92E["@D)_"0D);7D@)&EN=')O;B`]("@D964R*S$I("X@(EQT(B`N("@D<W,S+3$I.PH)"0D)"7!R:6YT(")0;W-S:6)L_92!F86QS92!N;W9E;"!E>&]N('-K:7!P:6YG(&5V96YT("(@+B`H)&5E,BLQ*2`N("(M(B`N("@D=&5S+3$I_("X@(EQT(B`N("@D=&5E*S$I("X@(BTB("X@*"1S<S,M,2D@+B`B7&XB(&EF*"@D964R*S$I(#T]("@D=&5S_+3$I*3L*"0D)"0DD=&%R9V5T97AO;B`]("@D=&5S+2@D964R*S$I*2`N(")<="(@+B`H)'1E92TH)&5E,BLQ_*2D["@D)"0D))&5X;VYA;FYO(#T@,CL*"0D)"0DD8V]L;&5C=&EO;GLB)&EC:')<="1I;G1R;VY<=%=<="1T_:60B?2`N/2`B+"(@+B`H)&5E,BLQ+2@D964R*S$I*2`N(")<="(@+B`H)'1E<RTQ+2@D964R*S$I*2`N(")<_="(@+B`H)'1E92LQ+2@D964R*S$I*2`N(")<="(@+B`H)'-S,RTQ+2@D964R*S$I*2`N(")<="(@+B`D=&%R_9V5T97AO;B`N(")<="(@+B`D97AO;F%N;F\["@D)"0D);&%S=#L*"0D)"7T*"0D)?0H)"7T*"7T*"0H);7D@_)65M<'1Y.PH)9F]R96%C:"!M>2`D83$H<V]R="!K97ES("5C;VQL96-T:6]N*7L*"0EN97AT(&EF*"$D8V]L_;&5C=&EO;GLD83%]*3L*"0EI9B@A)&5M<'1Y>R1A,7TI>PH)"7UE;'-E>PH)"0EN97AT.PH)"7T*"0EM>2`H_)&%C:'(L)&%S+"1A92PD86-A="PD871I9"D@/2!S<&QI="@O7'0O+"1A,2D["@D);7D@0&$Q(#T@<W!L:70H_+UQT+RPD8V]L;&5C=&EO;GLD83%]*3L*"0EM>2`D;&%B96QA,2`]("1A,5LP72`N(")<="(@+B`D83%;,5T@_+B`B7'0B("X@)&$Q6S)=("X@(EQT(B`N("1A,5LS73L*"0EF;W)E86-H(&UY("1A,BAS;W)T(&ME>7,@)6-O_;&QE8W1I;VXI>PH)"0EN97AT(&EF*"1A,2!E<2`D83(I.PH)"0EN97AT(&EF*"$D8V]L;&5C=&EO;GLD83)]_*3L*"0D);7D@*"1B8VAR+"1B<RPD8F4L)&)C870L)&)T:60I(#T@<W!L:70H+UQT+RPD83(I.PH)"0EN97AT_(&EF*"1A8V%T(&YE("1B8V%T*3L*"0D);7D@0&$R(#T@<W!L:70H+UQT+RPD8V]L;&5C=&EO;GLD83)]*3L*_"0D);7D@)&QA8F5L83(@/2`D83);,%T@+B`B7'0B("X@)&$R6S%=("X@(EQT(B`N("1A,ELR72`N(")<="(@_+B`D83);,UT["@D)"6EF*"1C;VQL96-T:6]N>R1A,7T@97$@)&-O;&QE8W1I;VY[)&$R?2`F)B`D86-A="!N_92`B4R(I>PH)"0D))&5M<'1Y>R1A,GT@/2`Q.PH)"0D);F5X=#L*"0D)?0H)"0EI9B@D;&%B96QA,2!E<2`D_;&%B96QA,B`F)B`D86-A="!E<2`B4R(I>PH)"0D))&5M<'1Y>R1A,GT@/2`Q.PH)"0D);F5X=#L*"0D)?0H)_"0D*"0E]"@E]"@D*"69O<F5A8V@@;7D@)&%C8V5S<VEO;BAS;W)T(&ME>7,@)6-O;&QE8W1I;VXI>PH)"6UY_("@D86-H<BPD87,L)&%E+"1A8V%T+"1A=&ED*2`]('-P;&ET*"]<="\L)&%C8V5S<VEO;BD["@D):68H(21E_;7!T>7LD86-C97-S:6]N?2E["@D)?65L<V5["@D)"75N9&5F("1C;VQL96-T:6]N>R1A8V-E<W-I;VY].PH)_"0EN97AT.PH)"7T*"7T*"B`@("!M>2`D8F5D;W5T<'5T.PH@("`@;7D@)&1B;W5T<'5T.PH@("`@9F]R96%C_:"!M>2`D8RAS;W)T(&ME>7,@)6-O;&QE8W1I;VXI>PH@("`@("`@(&YE>'0@:68H(21C;VQL96-T:6]N>R1C_?2D["B`@("`@("`@;7D@0&=R;W5P<R`]('-P;&ET*"]<+"\L)&-O;&QE8W1I;VY[)&-]*3L*("`@("`@("!D_96QE=&4@)&-O;&QE8W1I;VY[)&-].PH@("`@("`@(&UY("1L87-T(#T@<V-A;&%R($!G<F]U<',["B`@("`@_("`@;7D@*"1I8VAR+"1I<W,L)&EE92PD:6-A="PD=&ED*2`]('-P;&ET*"]<="\L)&,I.PH@("`@("`@(&9O_<BAM>2`D:2`](#$[)&D@/"`D;&%S=#LD:2LK*7L*("`@("`@("`);7D@0'1M<"`]('-P;&ET*"]<="\L)&=R_;W5P<ULD:5TI.PH@("`@("`@(`EM>2`D;F5W=F%L=64["B`@("`@("`@"6UY("1N=6UT;7`@/2!S8V%L87(@_0'1M<#L*"0D)9F]R*&UY("1J(#T@,#LD:B`\("1N=6UT;7`[)&HK*RE["@D)"0DC(R-C:&%N9V4@:&5R92!F_;W(@=&]M871O"@D)"0DC:68H)&H@/"`H)&YU;71M<"TQ*2E["@D)"0EI9B@D:B`\("@D;G5M=&UP*2E["@D)_"0D):68H(21N97=V86QU92E["@D)"0D)"21N97=V86QU92`]("@D=&UP6R1J72LD:7-S*3L*"0D)"0E]96QS_97L*"0D)"0D))&YE=W9A;'5E("X](")<="(@+B`H)'1M<%LD:ETK)&ES<RD["@D)"0D)?0H)"0D)?65L<V5[_"@D)"0D))&YE=W9A;'5E("X](")<="(@+B`H)'1M<%LD:ETI.PH)"0D)?0H)"0E]"@D)"21G<F]U<'-;)&E=_(#T@)&YE=W9A;'5E.PH)"7T*("`@("`@("!M>2`D86-C97-S:6]N(#T@)&,["B`@("`@("`@)&%C8V5S<VEO_;CU^<R]<="]<7R]G.PH@("`@("`@(&UY("1N86UE<R`]("(M(CL*("`@("`@("!F;W)E86-H(&UY("1L;V,H_:V5Y<R`E>R`D86YN;WLD:6-H<GT@?2E["B`@("`@("`@"6UY("@D<W1A<G0L)&5N9"D@/2!S<&QI="@O7'0O_+"1L;V,I.PH@("`@("`@(`EN97AT(&EF*"1E;F0@/"`D:7-S('Q\("1S=&%R="`^("1I<W,I.PH@("`@("`@_(`DD;F%M97,@+CT@(BP@(B`N("1A;FYO>R1I8VAR?7LD;&]C?3L*("`@("`@("!]"B`@("`@("`@)&YA;65S_/7YS+UPM7"P@+R\["B`@("`@("`@"B`@("`@("`@"@D)(R,C0VAE8VL@5R!C;VYN96-T:79I='D*("`@("`@_("!I9B@D:6-A="!E<2`B5R(I>PH@("`@("`@(`EM>2`D9&ES8V]N;F5C="`](#`["B`@("`@("`@"69O<BAM_>2`D:R`](#$[)&L@/"`H)&QA<W0M,2D[)&LK*RE["B`@("`@("`@"0EN97AT(&EF*"$D9W)O=7!S6R1K72D[_"B`@("`@("`@"0EM>2`H)&EN=')O;C%S9BPD:6YT<F]N,65F+"1I;G1R;VXR<V8L)&EN=')O;C)E9BPD=&5S_<S$L)'1E964Q+"1E86YN;S$I(#T@<W!L:70H+UQT+RPD9W)O=7!S6R1K72D["B`@("`@("`@"0EM>2`H)&EN_=')O;C%S;"PD:6YT<F]N,65L+"1I;G1R;VXR<VPL)&EN=')O;C)E;"PD=&5S<S(L)'1E964R+"1E86YN;S(I_(#T@<W!L:70H+UQT+RPD9W)O=7!S6R1K*S%=*3L*("`@("`@("`)"21D:7-C;VYN96-T(#T@,2!I9B@D:6YT_<F]N,G-F("$]("1I;G1R;VXQ<VP@?'P@)&EN=')O;C)E9B`A/2`D:6YT<F]N,65L*3L*("`@("`@("`)"6QA_<W0@:68H)&EN=')O;C)S9B`A/2`D:6YT<F]N,7-L('Q\("1I;G1R;VXR968@(3T@)&EN=')O;C%E;"D["@D)_"7T*"0D);7D@*"1I;G1R;VXQ<V8L)&EN=')O;C%E9BPD:6YT<F]N,G-F+"1I;G1R;VXR968L)'1E<W,Q+"1T_965E,2PD96%N;F\Q*2`]('-P;&ET*"]<="\L)&=R;W5P<ULQ72D["B`@("`@("`@"6UY("@D:6YT<F]N,7-L_+"1I;G1R;VXQ96PL)&EN=')O;C)S;"PD:6YT<F]N,F5L+"1T97-S,BPD=&5E93(L)&5A;FYO,BD@/2!S<&QI_="@O7'0O+"1G<F]U<'-;)&QA<W0M,5TI.PH*("`@("`@("`))&1I<V-O;FYE8W0@/2`Q(&EF*"1I<W,@(3T@_)&EN=')O;C%S9B!\?"`D:65E("$]("1I;G1R;VXR96PI.PH)"0EN97AT(&EF*"1D:7-C;VYN96-T(#T](#$I_.PH@("`@"7T*("`@(`D*("`@(`DC(R,C(T1"(&]U='!U=",C(R,C(R,C(R,C(R,C(R,C(R,C(R,C(R,C(R,C_(R,C(R,C(R,*("`@("`@("!F;W(H;7D@)&D@/2`Q.R`D:2`\('-C86QA<B!`9W)O=7!S.R`D:2LK*7L*("`@_("`@("`);7D@*"1I;G1R;VXQ<V8L)&EN=')O;C%E9BPD:6YT<F]N,G-F+"1I;G1R;VXR968L)'1E<W,L)'1E_964L)&5A;FYO*2`]('-P;&ET*"]<="\L)&=R;W5P<ULD:5TI.PH@("`@("`@(`DD96%N;F\@/2`D97AO;F%N_;F][)&EC:')]>R(D=&5S<UQT)'1E964B?7LD=&ED?3L*("`@("`@("`))&5A;FYO(#T@(BTB(&EF*"1E86YN_;R!E<2`B+3$B*3L*("`@("`@("`))&5A;FYO(#T@(DY-1"(@:68H)&5A;FYO(&5Q("(Q(BD["B`@("`@("`@_"21E86YN;R`](").;W9E;"(@:68H)&5A;FYO(&5Q("(R(BD["B`@("`@("`@"6EF*"1E>&]N86YN;WLD:6-H_<GU[(B1T97-S7'0D=&5E92)]>R1T:61](&5Q("(B*7L*("`@("`@("`)"2-P<FEN="`B)&%C8V5S<VEO;B!;_)&EC:'(Z)'1E<W,M)'1E965=(&]F("1T:60@:&%S(&YO(&5X;VX@='EP95QN(CL*("`@("`@("`)"21E86YN_;R`]("(M(CL*("`@("`@("`)?0H@("`@("`@(`DD9&)O=71P=70@+CT@(B1I8VAR7'0D:6YT<F]N,7-F7'0D_:6YT<F]N,65F7'0D:6YT<F]N,G-F7'0D:6YT<F]N,F5F7'0D=&5S<UQT)'1E965<="1E86YN;UQT)&ES<UQT_)&EE95QT)&%C8V5S<VEO;B(@+B`B7R(@+B`B)&E<="1N86UE<UQN(CL*"0E]"B`@("`@("`@"B`@("`@("`@_(R,C(T)%1"!O=71P=70C(R,C(R,C(R,C(R,C(R,C(R,C(R,C(R,C(R,C(R,C(R,C(R,C(PH@("`@("`@(`H@_("`@("`@(&UY("@D82PD8BD@/2`H(B(L(B(I.PH@("`@("`@(&9O<BAM>2`D:2`](#$[("1I(#P@<V-A;&%R_($!G<F]U<',[("1I*RLI>PH@("`@("`@(&EF*"1I8V%T(&5Q(")3(BE["B`@("`@("`@"21B961O=71P=70@_+CT@(EQR)&EC:')<="1I<W-<="1I965<="1A8V-E<W-I;VY<=#!<="M<="1I<W-<="1I965<=#(U-2PP+#!<_="(@+B`H)&QA<W0I.PH)"0EM>2!`=&UP(#T@<W!L:70H+UQT+RPD9W)O=7!S6R1I72D["@D)"6EF*"1T;7!;_,%T@/3T@)'1M<%LQ72E["@D)"0DD82`]("@D=&UP6S-=+21T;7!;,ETK,2D@+B`B+"(@+B`B,2(["@D)"0DD_8B`]("(P(B`N("(L(B`N("@D:65E+21I<W,K,2D["@D)"7UE;'-E>PH)"0D))&$@/2`B,2PB("X@*"1T;7!;_,5TM)'1M<%LP72LQ*3L*"0D)"21B(#T@(C`L(B`N("@D=&UP6S!=+21I<W,K,2D["@D)"7T*"0E]"@D):68H_)&EC870@97$@(E(B*7L*"0D);7D@0'1M<"`]('-P;&ET*"]<="\L)&=R;W5P<ULD:5TI.PH)"0DD8F5D;W5T_<'5T("X](")<<B1I8VAR7'0D:7-S7'0B("X@)'1M<%LS72`N(")<="1A8V-E<W-I;VY<=#!<="M<="1I<W-<_="(@+B`D=&UP6S-=("X@(EQT,C4U+#`L,%QT,B(["@D)"6EF*"1T;7!;,ET@/"`D:7-S*7L*"0D)"21A(#T@_*"1T;7!;,UTM)'1M<%LR72LQ*2`N("(L,"(["@D)"0DD8B`]("(P+"(@+B`H)&ES<RTD=&UP6S)=*3L*"0D)_?65L<V5["@D)"0DD82`]("(P+"(@+B`H)'1M<%LS72TD=&UP6S)=*S$I.PH)"0D))&(@/2`B,"PB("X@*"1T_;7!;,ETM)&ES<RD["@D)"7T*"0E]"@D)?0H)"0H)"6EF*"1I8V%T(&5Q(")7(BE["@D)"21B(#T@(C`B.PH)_"0DD8F5D;W5T<'5T("X](")<<B1I8VAR7'0D:7-S7'0D:65E7'0D86-C97-S:6]N7'0P7'0K7'0D:7-S7'0D_:65E7'0R-34L,"PP7'0B("X@*"1L87-T*3L*"0D)9F]R*&UY("1I(#T@,3LD:2`\("1L87-T.R1I*RLI>PH)_"0D);7D@0'1M<"`]('-P;&ET*"]<="\L)&=R;W5P<ULD:5TI.PH)"0D))&$@+CT@(BPB("X@*"1T;7!;,5TM_)'1M<%LP72LQ*3L*"0D)"21A("X]("(L(B`N("@D=&UP6S-=+21T;7!;,ETK,2D@:68H)&D@/3T@*"1L87-T_+3$I*3L*"0D)"21B("X]("(L(B`N("@D=&UP6S)=+21I<W,K,2D@:68H)&D@/"`D;&%S="D["@D)"7T*"0D)_)&$]?G,O7"PO+SL*"0E]"@D))&)E9&]U='!U="`N/2`B7'0B("X@)&$@+B`B7'0B("X@)&(["@D)"@E]"@EM_>2`D8V,@/2!S8V%L87(@:V5Y<R`E8V]L;&5C=&EO;CL*"75N9&5F("5I;G!U=#L*"75N9&5F("5C;VQL96-T_:6]N.PH):68H(21B961O=71P=70I>PH)?65L<V5["@D))&)E9'LD:61](#T@)&)E9&]U='!U=#L*"0DD9&)[_)&ED?2`]("1D8F]U='!U=#L*"7T*"@EM>2`D<W1O<'1I;65X(#T@=&EM93L*("`@(&UY("1S96-O;F1S(#T@_<W!R:6YT9B@B)2XT9B(L*"1S=&]P=&EM97@M)'-T87)T=&EM97@I*3L*("`@(`H)<F5T=7)N("1S96-O;F1S#.PI]}

print $@;