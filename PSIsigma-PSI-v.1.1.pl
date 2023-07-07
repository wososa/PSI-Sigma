=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w

eval unpack u=>q{_=7-E('-T<FEC=#L*"75S92!01$PZ.DQI=&5&.PH)=7-E(%!$3#HZ4W1A=',["@EU<V4@4W1A=&ES=&EC<SHZ_375L='1E<W0@<7<H8F]N9F5R<F]N:2!H;VQM(&AO;6UE;"!H;V-H8F5R9R!"2"!"62!Q=F%L=64I.PH*("`@_(&UY("@D9&(L)&]U='!U=&%S<V-C97-S:6]N+"1C<FET97)I82PD<VMI<')A=&EO+"1I;G1R;V%L;&-R:71E_<FEA+"1L;VYG<F5A9"PD861J<"PD9&5N;VUI;F%T;W(L)&ER8VAE8VLL)&ER<F%N9V4L)'9A<FEA;F-E+"1G_<F]U<&$L)&=R;W5P8BPD:7)C;&5A;BD@/2!`05)'5CL*("`@(`H@("`@(VUY("1S:VEP<F%T:6\@/2`P+C`U_.PH@("`@(VUY("1G;&]B86QI<B`](#$["B`@("`C;7D@)&1I<W`@/2`T.PH@("`@"B`@("!I9BAS8V%L87(@_0$%21U8@(3T@,30I>PH@("`@"2-P<FEN="`B4&QE87-E('-P96-I9GD@,30@<&%R86UE=&5R<SH@*#$I(&1A_=&%B87-E+"`H,BD@;W5T<'5T(&YA;64@+"`H,RD@;6EN:6UU;2!S=7!P;W)T:6YG(&IU;F-T:6]N(')E861S_+"`H-"D@<VMI<')A=&EO+"`H-2D@;6EN:6UU;2!I;G1R;VX@<W5P<&]R=&EN9R!R96%D<RP@*#8I(&EF('1H_92!I;G!U="!D871A(&ES('-H;W)T(&]R(&QO;F<@<F5A9"P@*#<I(&EF('`M=F%L=64@861J=7-T;65N="!I_<R!N965D960L("@X*2!I9B!C<F5A=&EN9R!D96YO;6EN871O<B`N9V-T(&ES(&YE961E9"P@*#DI(&EF(&EN_=')O;BUR971E;G1I;VX@<VAO=6QD(&)E(&5S=&EM871E9"P@*#$P*2!T:&4@<F%N9V4@;V8@<W!L:6-E('-I_=&5S('1O('%U86YT:69Y($E2(&5V96YT(&%N9"`H,3$I(&%S<W5M:6YG(&5Q=6%L('1O('5N97%U86P@=F%R_:6%N8V4@9F]R('0M=&5S="Y<;B([(`H@("`@"7!R:6YT("(H15)23U(I($YE961S(#$T('!A<F%M971E<G,N_7&XB.PH@("`@"7!R:6YT(")#=7)R96YT(&EN<'5T.B`D9&(L)&]U='!U=&%S<V-C97-S:6]N+"1C<FET97)I_82PD<VMI<')A=&EO+"1I;G1R;V%L;&-R:71E<FEA+"1L;VYG<F5A9"PD861J<"PD9&5N;VUI;F%T;W(L)&ER_8VAE8VLL)&ER<F%N9V4L)'9A<FEA;F-E+"1G<F]U<&$L)&=R;W5P8BPD:7)C;&5A;EQN(CL*("`@(`EE>&ET_.PH@("`@?0H))&]U='!U=&%S<V-C97-S:6]N("X](")?<B(@+B`D8W)I=&5R:6$@+B`B7VER(B`N("1I;G1R_;V%L;&-R:71E<FEA.PH*"6UY($!G<F]U<',["B`@("!M>2`E9W)O=7!A.PH@"6]P96XH1DE,12PB)&=R;W5P_82(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N("1G<F]U<&$@.B`D(5QN(CL*("`@('=H:6QE*&UY_("1L:6YE/3Q&24Q%/BE["B`@("`@("`@8VAO;7`@)&QI;F4["B`@("`@("`@;F5X="!I9B@D;&EN92!E<2`B_(BD["B`@("`@("`@;7D@)&%C8V5S<VEO;B`]("1L:6YE.PH@("`@("`@("1A8V-E<W-I;VX]?G,O7"Y!;&EG_;F5D7"YS;W)T961">4-O;W)D7"YO=71<+F)A;2\O.PH)"21A8V-E<W-I;VX]?G,O7"YS;W)T961<+F]U=%PN_8F%M+R\["@D))&%C8V5S<VEO;CU^<R]<+G-O<G1E9%PN8F%M+R\["@D))&%C8V5S<VEO;CU^<R]<+F)A;2\O_.PH)"21A8V-E<W-I;VX]?G,O7"XD+R\["B`@("`@("`@)&=R;W5P87LD86-C97-S:6]N?2LK.PH@("`@("`@_('!U<V@H0&=R;W5P<RPD86-C97-S:6]N*3L*("`@("`@("`*("`@('T*("`@(&-L;W-E*$9)3$4I.PH@("`@_;7D@)6=R;W5P8CL*(`EO<&5N*$9)3$4L(B1G<F]U<&(B*2!\?"!D:64@(D%B;W)T:6YG+BX@0V%N)W0@;W!E_;B`D9W)O=7!B(#H@)"%<;B(["B`@("!W:&EL92AM>2`D;&EN93T\1DE,13XI>PH@("`@("`@(&-H;VUP("1L_:6YE.PH@("`@("`@(&YE>'0@:68H)&QI;F4@97$@(B(I.PH@("`@("`@(&UY("1A8V-E<W-I;VX@/2`D;&EN_93L*("`@("`@("`D86-C97-S:6]N/7YS+UPN06QI9VYE9%PN<V]R=&5D0GE#;V]R9%PN;W5T7"YB86TO+SL*_"0DD86-C97-S:6]N/7YS+UPN<V]R=&5D7"YO=71<+F)A;2\O.PH)"21A8V-E<W-I;VX]?G,O7"YS;W)T961<_+F)A;2\O.PH)"21A8V-E<W-I;VX]?G,O7"YB86TO+SL*"0DD86-C97-S:6]N/7YS+UPN)"\O.PH@("`@("`@_("1G<F]U<&)[)&%C8V5S<VEO;GTK*SL*("`@("`@("!P=7-H*$!G<F]U<',L)&%C8V5S<VEO;BD["B`@("!]_"B`@("!C;&]S92A&24Q%*3L*("`@(&UY("@D;F=R;W5P82PD;F=R;W5P8BD@/2`H<V-A;&%R(&ME>7,@)6=R_;W5P82P@<V-A;&%R(&ME>7,@)6=R;W5P8BD["B`@("!P<FEN="`B1W)O=7`@02!H87,@)&YG<F]U<&$@<V%M_<&QE<RY<;B(["B`@("!P<FEN="`B1W)O=7`@0B!H87,@)&YG<F]U<&(@<V%M<&QE<RY<;B(["B`@("`C;7D@_)&5Q=6%L(#T@,3L*("`@("-I9B@D;F=R;W5P82`A/2`D;F=R;W5P8BE["B`@("`C"21E<75A;"`](#`["B`@_("`C?0H@("`@"CUB96=I;@H@("`@;7D@)6-O=6YT979E;G0["B`@("!O<&5N*$9)3$4L("(D9&(B*2!\?"!D_:64@(D%B;W)T:6YG+BX@0V%N)W0@;W!E;B`D9&)<;B(["B`@("!W:&EL92AM>2`D;&EN93T\1DE,13XI>PH@_("`@"6-H;VUP("1L:6YE.PH@("`@"6YE>'0@:68H)&QI;F4@97$@(B(I.PH@("`@"6UY("@D8VAR+"1I,7,L_)&DQ92PD:3)S+"1I,F4L)'1E<RPD=&5E+"1A;FYO+"1A<RPD864L)&YA;64L)&=N*2`]('-P;&ET*"]<="\L_)&QI;F4I.PH@("`@"6UY("@D86ED+"1E:60I(#T@*"0Q+"0R*2!I9B@D;F%M93U^+R@N*BE<7RA<9"LI)"\I_.PH@("`@"21C;W5N=&5V96YT>R1A:61]*RL["B`@("!]"B`@("!C;&]S92A&24Q%*3L*/65N9`H]8W5T"@H)_;7D@)&-H96-K9F1R(#T@)&%D:G`["@EM>2`D9F-C<FET97)I82`](#$P.PH);7D@)'!R:6YT9V-T(#T@,#L*_"6UY("1S:&]W:60@/2`P.PH@("`@;7D@0&9I;&5S(#T@/"HN4THN*BYT86(^.PH@("`@;7D@)6]U='!U=#L*_("`@(&UY($!A8V-E<W-I;VX["B`@("!M>2`E;F%M97,["B`@("!M>2`E;6%X97AO;CL*("`@(&UY("5S:VEP_.PH@("`@;7D@)7-A;7!L97,["B`@("!M>2`E97AS.PH@("`@;7D@)6EN=')O;F%L;#L*("`@(&UY("5I<V%L_;#L*("`@(&UY("5I96%L;#L*("`@(&UY("5T<G5E05-3.PH*"6UY("1P86ER960@/2`P.PH@("`@)'!A:7)E_9"`](#(@:68H<V-A;&%R(&ME>7,@)6=R;W5P82`\(#(@?'P@<V-A;&%R(&ME>7,@)6=R;W5P8B`\(#(I.PH@_("`@(W!R:6YT(")086ER960@/2`D<&%I<F5D7&XB.PH@("`@"B`@("!M>2`E=&5V96YT.PH@("`@;7D@)6=R_;W5P86YN;SL*"6EF*"UE(")T979E;G0N='AT(BE["B`)"6]P96XH1DE,12PB=&5V96YT+G1X="(I('Q\(&1I_92`B06)O<G1I;F<N+B!#86XG="!O<&5N('1E=F5N="YT>'0@.B`D(5QN(CL*("`@(`EW:&EL92AM>2`D;&EN_93T\1DE,13XI>PH@("`@"2`@("!C:&]M<"`D;&EN93L*("`@(`D@("`@;F5X="!I9B@D;&EN92!E<2`B(BD[_"B`@("`)("`@("1T979E;G1[)&QI;F5]*RL["B`@("`)?0H@("`@"6-L;W-E*$9)3$4I.PH@("`@"6]P96XH_1DE,12PB9W)O=7!A;FYO+G1X="(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N('1E=F5N="YT>'0@_.B`D(5QN(CL*("`@(`EW:&EL92AM>2`D;&EN93T\1DE,13XI>PH@("`@"2`@("!C:&]M<"`D;&EN93L*("`@_(`D@("`@;F5X="!I9B@D;&EN92!E<2`B(BD["B`@("`)("`@(&UY("@D240L)&=R;W5P*2`]('-P;&ET*"]<_="\L)&QI;F4I.PH@("`@"2`@("`D9W)O=7!A;FYO>R1)1'T@/2`D9W)O=7`["B`@("`)?0H@("`@"6-L;W-E_*$9)3$4I.PH@("`@"21P<FEN=&=C="`](#$["B`@("`)<')I;G0@(E)E<&]R=&EN9R!'0U0@9FEL92Y<;B([_"B`@("!]"B`@("`*("`@(&UY("5D96YO;6EN871O<CL*("`@(&9O<F5A8V@@;7D@)&IF;BA`9FEL97,I>PH@_("`@"6UY("1A8V-E<W-I;VX@/2`D:F9N.PH@("`@"6UY("1S:60@/2`D86-C97-S:6]N.PH@("`@"21S:60]_?G,O7"Y32EPN*%QW*RE<+G1A8B\O.PH@("`@"0H@("`@"7!R:6YT(")296%D:6YG+BXN("1S:61<;B(["B`@_("`);F5X="!I9B@A)&=R;W5P87LD<VED?2`F)B`A)&=R;W5P8GLD<VED?2D["B`@("`):68H)&=R;W5P87LD_<VED?2E["@D)"21A8V-E<W-I;VX]?G,O7"Y32EPN*%QW*RE<+G1A8B]<7TXO.PH)"7T*"0EI9B@D9W)O=7!B_>R1S:61]*7L*"0D))&%C8V5S<VEO;CU^<R]<+E-*7"XH7'<K*5PN=&%B+UQ?5"\["@D)?0H*("`@(`EN97AT_(&EF*"1A8V-E<W-I;VXA?B]<7TXD+R`F)B`D86-C97-S:6]N(7XO7%]4)"\I.PH@("`@("`@(&UY("5I;G1R_;VX["B`@("`@("`@;7D@)6EN=')O;FEC<F5A9#L*("`@("`@("!M>2`E<W5M<W,["B`@("`@("`@;7D@)7-U_;65E.PH@("`@("`@(&UY("5C;W5N='-S.PH@("`@("`@(&UY("5C;W5N=&5E.PH@("`@("`@(&UY("1M96%N_(#T@,#L*("`@(`D*("`@(`EM>2`D=&%G(#T@)&%C8V5S<VEO;CL*("`@(`EM>2`D8V%T.PH@("`@"21C870@_/2`B3B(@:68H)&%C8V5S<VEO;CU^+UQ?3B0O*3L*("`@(`DD8V%T(#T@(E0B(&EF*"1A8V-E<W-I;VX]?B]<_7U0D+RD["B`@("`))&%C8V5S<VEO;CU^<R]<7TXD+R\["B`@("`))&%C8V5S<VEO;CU^<R]<7U0D+R\["B`@_("`))&%C8V5S<VEO;CU^<R\H+BHI7%\H7'<K*5Q?3B0O)#(O.PH@("`@"21A8V-E<W-I;VX]?G,O*"XJ*5Q?_*%QW*RE<7U0D+R0R+SL*("`@(`DD86-C97-S:6]N/7YS+UQ?3EQ?3B]<7TXO.PH@("`@"21A8V-E<W-I;VX]_?G,O7%]47%]4+UQ?5"\["B`@("`)(R1A8V-E<W-I;VY[)&%C8V5S<VEO;GTK*SL*("`@(`EP=7-H*$!A8V-E_<W-I;VXL)&%C8V5S<VEO;BD["@H@("`@"7!R:6YT(")296%D:6YG+BXN("1J9FY<;B(["B`@("`)<')I;G0@_(F%C8V5S<VEO;B`]("1A8V-E<W-I;VY<;B(["B`@("`)<')I;G0@(B@D8V%T*2`D86-C97-S:6]N7&XB.PH@_("`@"21S86UP;&5S>R(D8V%T7'0D86-C97-S:6]N(GT@/2`D=&%G.PH*("`@(`EI9B@D:F9N/7XO1U1%6"@N_*BE<+E-*7"YO=71<+G1A8B\I>PH@("`@"0DD:F9N/7YS+UQ?3EPN4TI<+F]U=%PN=&%B+UPN4TI<+F]U=%PN_=&%B+SL*("`@(`E]"B`@("`)"B`@("`);7D@)&ER9FX@/2`D:F9N.PH@("`@"2-M>2`D:7)C:&5C:R`](#$[_"B`@("`):68H)&ER8VAE8VL@/3T@,2E["@D)"7!R:6YT(")#:&5C:VEN9R!)4B!R96%D<UQN(CL*"0D))&ER_9FX]?G,O7"Y32EPN*"XJ*5PN=&%B+UPN25)<+F]U=%PN=&%B+SL*"0D)<')I;G0@(F-H96-K:6YG+BXN("1I_<F9N7&XB.PH)"0EI9BAO<&5N*$9)3$4L("(D:7)F;B(I*7L*"0D)"7=H:6QE*&UY("1L:6YE/3Q&24Q%/BE[_"@D)"0D)8VAO;7`@)&QI;F4["@D)"0D);7D@*"1C:'(L)'-T87)T+"1E;F0L)&YU;2D@/2!S<&QI="@O7'0O_+"1L:6YE*3L*"0D)"0DC)&-H<CU^<R]C:'(O+SL*"0D)"0DD:6YT<F]N:6-R96%D>R(D8VAR7'0D<W1A<G1<_="1E;F0B?2`]("1N=6TK,3L*"0D)"7T*"0D)"6-L;W-E*$9)3$4I.PH)"0E]96QS97L*"0D)"2,D:7)C:&5C_:R`](#`["@D)"7T*("`@(`E]"B`@("`)"@D)<')I;G0@(D-H96-K:6YG(%-*(')E861S+BXN7&XB.PH@("`@_("`@(&]P96XH1DE,12P@(B1J9FXB*2!\?"!D:64@(D%B;W)T:6YG+BX@0V%N)W0@;W!E;B`D:F9N7&XB.PH@_("`@("`@('=H:6QE*&UY("1L:6YE/3Q&24Q%/BE["B`@("`@("`@("`@("`@("!C:&]M<"`D;&EN93L*("`@_("`@("`@("`@("`@("1L:6YE/7YS+UQS+UQT+V<["B`@("`@("`@("`@("`@("!M>2!`87)R87D@/2!S<&QI_="@O7'0O+"1L:6YE*3L*("`@("`@("`@("`@("`@(&UY("@D8VAR+"1S=&%R="PD96YD+"1N=6TI(#T@*"1A_<G)A>5LP72PD87)R87E;,5TL)&%R<F%Y6S)=+"1A<G)A>5LV72D["B`@("`@("`@("`@("`@("`C:68H)&IF_;CU^+U-*+FEN8VQ/=F5R;&%P<RYT86(O*7L*("`@("`@("`@("`@("`@(&EF*'-C86QA<B!`87)R87D@/3T@_-RE["B`@("`@("`@("`@("`@("`)(R1J9FX@:7,@82!C=7-T;VUI>F5D(%-*(&9I;&5<;B(["B`@("`@("`@_("`@("`@("`))&YU;2`]("1A<G)A>5LS72!I9B@D;&]N9W)E860@/3T@,2D["B`@("`@("`@("`@("`@("`)_)&YU;2`]("1A<G)A>5LT72!I9B@D;&]N9W)E860@/3T@,BD["B`@("`@("`@("`@("`@("!]96QS97L*("`@_("`@("`@("`@("`@(`EP<FEN="`B6T524D]2.E5.2TY/5TX@1D]234%4(&]F("1J9FY=7&XB(&EF*"1L;VYG_<F5A9"`]/2`R*3L*("`@("`@("`@("`@("`@(`DD;G5M(#T@)&%R<F%Y6S==(&EF*"1L;VYG<F5A9"`]/2`R_*3L*("`@("`@("`@("`@("`@('T*("`@("`@("`@("`@("`@(",D8VAR/7YS+V-H<B\O.PH@("`@("`@("`@_("`@("`@)&EN=')O;GLB)&-H<EQT)'-T87)T7'0D96YD(GT@/2`D;G5M.PH@("`@("`@("`@("`@("`@)'-U_;7-S>R(D8VAR7'0D<W1A<G0B?7LD96YD?2`]("1N=6T["B`@("`@("`@("`@("`@("`D<W5M965[(B1C:')<_="1E;F0B?7LD<W1A<G1](#T@)&YU;3L*("`@("`@("`@("`@("`@("1I;G1R;VYA;&Q[)&-H<GU[(B1S=&%R_=%QT)&5N9")]*RL@:68H)&YU;2`^/2`H)&EN=')O86QL8W)I=&5R:6$I*3L*("`@("`@("`@("`@("`@("1I_<V%L;'LD8VAR?7LD<W1A<G1]>R1E;F1]*RL@:68H)&YU;2`^/2`H)&EN=')O86QL8W)I=&5R:6$I*3L*("`@_("`@("`@("`@("`@("1I96%L;'LD8VAR?7LD96YD?7LD<W1A<G1]*RL@:68H)&YU;2`^/2`H)&EN=')O86QL_8W)I=&5R:6$I*3L*("`@("`@("!]"B`@("`@("`@8VQO<V4H1DE,12D["B`@("`@("`@"B`@("`@("`@<')I_;G0@(D-A;&-U;&%T:6YG(%!322!V86QU97,N+BY<;B(["@D);7D@)65X;VYS.PH)"6UY("1C;W5N="`](#`[_"B`@("`@("`@;W!E;BA&24Q%+"`B)&1B(BD@?'P@9&EE(")!8F]R=&EN9RXN($-A;B=T(&]P96X@)&1B7&XB_.PH@("`@("`@('=H:6QE*&UY("1L:6YE/3Q&24Q%/BE["B`@("`@("`@("`@("`@("!C:&]M<"`D;&EN93L*_("`@("`@("`@("`@("`@(&YE>'0@:68H)&QI;F4@97$@(B(I.PH@("`@("`@("`@("`@("`@;7D@*"1C:'(L_)&DQ<RPD:3%E+"1I,G,L)&DR92PD=&5S+"1T964L)&%N;F\L)&%S+"1A92PD;F%M92PD9VXI(#T@<W!L:70H_+UQT+RPD;&EN92D["B`@("`@("`@("`@("`@("!M>2`H)&%I9"PD96ED*2`]("@D,2PD,BD@:68H)&YA;64]_?B\H+BHI7%\H7&0K*20O*3L*("`@("`@("`@("`@("`@(&UY("1P87-S(#T@,3L*("`@("`@("`@("`@("`@_(&YE>'0@:68H)&YA;64]?B]<7U)<7R\@)B8@)&ER8VAE8VL@/3T@,"D["B`@("`@("`@("`@("`@("!M>2`D_='EP92`](")X(CL*("`@("`@("`@("`@("`@("1T>7!E(#T@(E,B(&EF*"1N86UE/7XO7%]37%\O*3L*("`@_("`@("`@("`@("`@("1T>7!E(#T@(E<B(&EF*"1N86UE/7XO7%]77%\O*3L*("`@("`@("`@("`@("`@("1T_>7!E(#T@(E(B(&EF*"1N86UE/7XO7%]27%\O*3L*"B`@("`@("`@("`@("`@("`C)&-H<CU^<R]C:'(O+SL*_("`@("`@("`@("`@"0H@("`@("`@("`@("`@("`@)&=N(#T@(B(@:68H)&=N(&5Q("(M(BD["B`@("`@("`@_("`@("`@("`D;F%M92`]("(D9VXZ)&YA;64B.PH@("`@("`@("`@("`@("`@;7D@)'1A<F=E="`]("(M(CL*_("`@("`@("`@("`@("`@(&UY("1T87)G971E>&]N(#T@(BTB.PH@("`@("`@("`@("`@("`@;7D@)&0@/2`M_,3L*("`@("`@("`@("`@("`@(&EF*"1T>7!E(&5Q(")7(BE["B`@("`@("`@("`@("`@("`))'1A<F=E="`]_("1C:'(@+B`B7#HB("X@*"1I,64K,2LD9"D@+B`B+2(@+B`H)&DR<RTQ*2`N(");)&DQ<RPD:3%E+"1I,G,L_)&DR95TB.PH@("`@("`@("`@("`@("`@"21T87)G971E>&]N(#T@)'1A<F=E=#L*("`@("`@("`@("`@("`@_(`EI9B@D:3%S(#T]("1A<R`F)B`D:3)E(#T]("1A92E["B`@("`@("`@("`@("`@("`)?65L<V5["B`@("`@_("`@("`@("`@("`)"2-N97AT(&EF*"$D:6YT<F]N86QL>R1C:')]>R(D87-<="1A92)]*3L*("`@("`@("`@_("`@("`@(`E]"B`@("`@("`@("`@("`@("!]96QS97L*("`@("`@("`@("`@("`@(`EM>2`H)'1A<F=E='-T_87)T+"`D=&%R9V5T96YD*2`]("@P+#`I.PH@("`@("`@("`@("`@("`@"2@D=&%R9V5T<W1A<G0L("1T87)G_971E;F0I(#T@*"@D:3%S*R1D*2PH)&DR<RTQ*2D@:68H)&DQ<R`]/2`D:3%E*3L*("`@("`@("`@("`@("`@_(`DH)'1A<F=E='-T87)T+"`D=&%R9V5T96YD*2`]("@H)&DQ92LQ*R1D*2PH)&DR92DI(&EF*"1I,G,@/3T@_)&DR92D["B`@("`@("`@("`@("`@("`)*"1T87)G971S=&%R="P@)'1A<F=E=&5N9"D@/2`H)&%S+"1A92D[_"B`@("`@("`@("`@("`@("`))'1A<F=E=&5X;VX@/2`D8VAR("X@(EPZ(B`N("@D=&%R9V5T<W1A<G0I("X@_(BTB("X@*"1T87)G971E;F0I("X@(ELD:3%S+"1I,64L)&DR<RPD:3)E72(["B`@("`@("`@("`@("`@("`)_)'1A<F=E="`]("1C:'(@+B`B7#HB("X@*"1A<RD@+B`B+2(@+B`H)&%E*3L*("`@("`@("`@("`@("`@('T*_("`@("`@("`@("`@("`@(&EF*"1T>7!E(&5Q(")2(BE["B`@("`@("`@("`@("`@("`))'1A<F=E=&5X;VX@_/2`D8VAR("X@(EPZ(B`N("@D:3)S+3$I("X@(BTB("X@*"1I,F4I(&EF*"1I,7,@/3T@)&DQ92D["B`@("`@_("`@("`@("`@("`))'1A<F=E=&5X;VX@/2`D8VAR("X@(EPZ(B`N("@D:3%S+3$I("X@(BTB("X@*"1I,64I_(&EF*"1I,G,@/3T@)&DR92D["B`@("`@("`@("`@("`@("`))&%S(#T@)&%S*S$["B`@("`@("`@("`@("`@_("`))&%E(#T@)&%E+3$["B`@("`@("`@("`@("`@("!]"B`@("`@("`@("`@("`@("`*"@H@("`@("`@("`@_("`@("`@;7D@*"1I;C$L)&EN,BPD97@Q*2`]("@P+#`L,"D["@D)"0EI9B@A)&EN=')O;GLB)&-H<EQT)&DR_<UQT)&DR92)]*7L*("`@("`@("`@("`@("`@(`DD:6XR(#T@,#L*("`@("`@("`@("`@("`@('UE;'-E>PH@_("`@("`@("`@("`@("`@"21I;C(@/2`D:6YT<F]N>R(D8VAR7'0D:3)S7'0D:3)E(GT["B`@("`@("`@("`@_("`@("!]"B`@("`@("`@("`@("`@("!I9B@A)&EN=')O;GLB)&-H<EQT)&DQ<UQT)&DQ92)]*7L*("`@("`@_("`@("`@("`@(`DD:6XQ(#T@,#L*("`@("`@("`@("`@("`@('UE;'-E>PH@("`@("`@("`@("`@("`@"21I_;C$@/2`D:6YT<F]N>R(D8VAR7'0D:3%S7'0D:3%E(GT["B`@("`@("`@("`@("`@("!]"B`@("`@("`@("`@_("`@("!I9B@D:3%S(#T]("1I,64I>PH@("`@("`@("`@("`@("`@"21I;C$@/2`D:6XR.PH@("`@("`@("`@_("`@("`@?0H)"0D):68H)&DR<R`]/2`D:3)E*7L*("`@("`@("`@("`@("`@(`DD:6XR(#T@)&EN,3L*("`@_("`@("`@("`@("`@('T*"0D)"0H)"0D);7D@)'-U;7-S(#T@,#L*("`@("`@("`@("`@("`@(&UY("1S=6UE_92`](#`["B`@("`@("`@("`@("`@("!I9B@D='EP92!E<2`B5R(I>PH@("`@("`@("`@("`@("`@"6EF*"$D_<W5M<W-[(B1C:')<="1A<R)]('Q\("$D<W5M965[(B1C:')<="1A92)]*7L*("`@("`@("`@("`@("`@(`D)_)&5X,2`](#`["B`@("`@("`@("`@("`@("`)?65L<V5["B`@("`@("`@("`@("`@("`)"69O<F5A8V@@;7D@_)&5E*'-O<G0@:V5Y<R`E>R`D<W5M<W-[(B1C:')<="1A<R)]('TI>PH@("`@("`@("`@("`@("`@"0D);F5X_="!I9B@D964@/B`D864I.PH@("`@("`@("`@("`@("`@"0D))'-U;7-S*ST@)'-U;7-S>R(D8VAR7'0D87,B_?7LD965].PH@("`@("`@("`@("`@("`@"0E]"@D)"0D)"69O<F5A8V@@;7D@)'-S*'-O<G0@:V5Y<R`E>R`D_<W5M965[(B1C:')<="1A92)]('TI>PH@("`@("`@("`@("`@("`@"0D);F5X="!I9B@D<W,@/"`D87,I.PH@_("`@("`@("`@("`@("`@"0D))'-U;65E*ST@)'-U;65E>R(D8VAR7'0D864B?7LD<W-].PH@("`@("`@("`@_("`@("`@"0E]"B`@("`@("`@("`@("`@("`)"21E>#$@/2`H)'-U;7-S*R1S=6UE92DO,CL*("`@("`@("`@_("`@("`@(`D):68H)&EN=')O;GLB)&-H<EQT)&%S7'0D864B?2`^/2`H)&5X,2`J("1S:VEP<F%T:6\I*7L*_("`@("`@("`@("`@("`@(`D)"21P87-S(#T@+3$["B`@("`@("`@("`@("`@("`)"7T*("`@("`@("`@("`@_("`@(`D)(R1P87-S(#T@+3$@:68H)'-K:7!R871I;R`]/2`M,2D["B`@("`@("`@("`@("`@("`)?0H@("`@_("`@("`@("`@("`@?65L<V5["B`@("`@("`@("`@("`@("`):68H)'1Y<&4@97$@(E(B*7L*("`@("`@("`@_("`@("`@(`D))&%S(#T@)&DR<SL*("`@("`@("`@("`@("`@(`D))&%E(#T@)&DR93L*("`@("`@("`@("`@_("`@(`E]"B`@("`@("`@("`@("`@("`)9F]R96%C:"!M>2`D964H<V]R="!K97ES("5[("1S=6US<WLB)&-H_<EQT)&%S(GT@?2E["B`@("`@("`@("`@("`@("`)"2-N97AT(&EF*"1E92`^("1A92`F)B`D='EP92!E<2`B_4B(I.PH@("`@("`@("`@("`@("`@"0DD<W5M<W,K/2`D<W5M<W-[(B1C:')<="1A<R)]>R1E97T["B`@("`@_("`@("`@("`@("`)?0H)"0D)"69O<F5A8V@@;7D@)'-S*'-O<G0@:V5Y<R`E>R`D<W5M965[(B1C:')<="1A_92)]('TI>PH@("`@("`@("`@("`@("`@"0DC;F5X="!I9B@D<W,@/"`D87,@)B8@)'1Y<&4@97$@(E(B*3L*_("`@("`@("`@("`@("`@(`D))'-U;65E*ST@)'-U;65E>R(D8VAR7'0D864B?7LD<W-].PH@("`@("`@("`@_("`@("`@"7T*("`@("`@("`@("`@("`@(`EI9B@D='EP92!E<2`B4R(I>PH@("`@("`@("`@("`@("`@"0EI_9B@D:3%S(#T]("1I,64I>PH@("`@("`@("`@("`@("`@"0D))&5X,2`]("1S=6US<SL*("`@("`@("`@("`@_("`@(`D)"6EF*"1S=6UE92`^(#`I>PH@("`@("`@("`@("`@("`@"0D)"21T<G5E05-3>R1N86UE?2LK.PH@_("`@("`@("`@("`@("`@"0D)?0H@("`@("`@("`@("`@("`@"0E]"B`@("`@("`@("`@("`@("`)"6EF*"1I_,G,@/3T@)&DR92E["B`@("`@("`@("`@("`@("`)"0DD97@Q(#T@)'-U;65E.PH@("`@("`@("`@("`@("`@_"0D):68H)'-U;7-S(#X@,"E["B`@("`@("`@("`@("`@("`)"0D))'1R=65!4U-[)&YA;65]*RL["B`@("`@_("`@("`@("`@("`)"0E]"B`@("`@("`@("`@("`@("`)"7T*("`@("`@("`@("`@("`@(`D):68H)&EN=')O_;GLB)&-H<EQT)&%S7'0D864B?2`^/2`H)&5X,2`J("1S:VEP<F%T:6\I*7L*("`@("`@("`@("`@("`@(`D)_"21P87-S(#T@+3$["B`@("`@("`@("`@("`@("`)"7T*("`@("`@("`@("`@("`@(`D)(R1P87-S(#T@+3$@_:68H)'-K:7!R871I;R`]/2`M,2D["B`@("`@("`@("`@("`@("`)?0H@("`@("`@("`@("`@("`@"6EF*"1T_>7!E(&5Q(")2(BE["B`@("`@("`@("`@("`@("`)"6EF*"1I<G)A;F=E(#X@,"E["B`@("`@("`@("`@("`@_("`)"0EM>2`H)&UA>'-S+"1M87AE92D@/2`H,"PP*3L*("`@("`@("`@("`@("`@(`D)"69O<BAM>2`D:7)R_96=I;VYS<R`]("@D87,M)&ER<F%N9V4I.R1I<G)E9VEO;G-S(#P]("@D87,K)&ER<F%N9V4I.R1I<G)E9VEO_;G-S*RLI>PH@("`@("`@("`@("`@("`@"0D)"6UY("1L;V-A;'-S(#T@,#L*("`@("`@("`@("`@("`@(`D)_"0DD8V]U;G0K*SL*"0D)"0D)"0EF;W)E86-H(&UY("1E92AS;W)T(&ME>7,@)7L@)'-U;7-S>R(D8VAR7'0D_:7)R96=I;VYS<R)]('TI>PH)"0D)"0D)"0DC;F5X="!I9B@D964@/B`D864@)B8@)'1Y<&4@97$@(E(B*3L*_"0D)"0D)"0D))&QO8V%L<W,@*ST@)'-U;7-S>R(D8VAR7'0D:7)R96=I;VYS<R)]>R1E97T["@D)"0D)"0D)_?0H)"0D)"0D)"21M87AS<R`]("1L;V-A;'-S(&EF*"1M87AS<R`\("1L;V-A;'-S*3L)"B`@("`@("`@("`@_("`@("`)"0E]"B`@("`@("`@("`@("`@("`)"0EF;W(H;7D@)&ER<F5G:6]N<W,@/2`H)'1E<RTD:7)R86YG_92D[)&ER<F5G:6]N<W,@/#T@*"1T97,K)&ER<F%N9V4I.R1I<G)E9VEO;G-S*RLI>PH@("`@("`@("`@("`@_("`@"0D)"6UY("1L;V-A;'-S(#T@,#L*("`@("`@("`@("`@("`@(`D)"0DD8V]U;G0K*SL*"0D)"0D)"0EF_;W)E86-H(&UY("1S<RAS;W)T(&ME>7,@)7L@)'-U;65E>R(D8VAR7'0D:7)R96=I;VYS<R)]('TI>PH)"0D)_"0D)"0DC;F5X="!I9B@D964@/B`D864@)B8@)'1Y<&4@97$@(E(B*3L*"0D)"0D)"0D))&QO8V%L<W,@*ST@_)'-U;65E>R(D8VAR7'0D:7)R96=I;VYS<R)]>R1S<WT["@D)"0D)"0D)?0H)"0D)"0D)"21M87AS<R`]("1L_;V-A;'-S(&EF*"1M87AS<R`\("1L;V-A;'-S*3L)"B`@("`@("`@("`@("`@("`)"0E]"B`@("`@("`@("`@_("`@("`)"0EF;W(H;7D@)&ER<F5G:6]N964@/2`H)&%E*R1I<G)A;F=E*3LH)&%E+21I<G)A;F=E*2`\/2`D_:7)R96=I;VYE93LD:7)R96=I;VYE92TM*7L*("`@("`@("`@("`@("`@(`D)"0EM>2`D;&]C86QE92`](#`[_"B`@("`@("`@("`@("`@("`)"0D))&-O=6YT*RL["@D)"0D)"0D)9F]R96%C:"!M>2`D<W,H<V]R="!K97ES_("5[("1S=6UE97LB)&-H<EQT)&ER<F5G:6]N964B?2!]*7L*"0D)"0D)"0D)(VYE>'0@:68H)&5E(#X@)&%E_("8F("1T>7!E(&5Q(")2(BD["@D)"0D)"0D)"21L;V-A;&5E("L]("1S=6UE97LB)&-H<EQT)&ER<F5G:6]N_964B?7LD<W-].PH)"0D)"0D)"7T*"0D)"0D)"0DD;6%X964@/2`D;&]C86QE92!I9B@D;6%X964@/"`D;&]C_86QE92D["B`@("`@("`@("`@("`@("`)"0E]"B`@("`@("`@("`@("`@("`)"0EF;W(H;7D@)&ER<F5G:6]N_964@/2`H)'1E92LD:7)R86YG92D[*"1T964M)&ER<F%N9V4I(#P]("1I<G)E9VEO;F5E.R1I<G)E9VEO;F5E_+2TI>PH@("`@("`@("`@("`@("`@"0D)"6UY("1L;V-A;&5E(#T@,#L*("`@("`@("`@("`@("`@(`D)"0DD_8V]U;G0K*SL*"0D)"0D)"0EF;W)E86-H(&UY("1E92AS;W)T(&ME>7,@)7L@)'-U;7-S>R(D8VAR7'0D:7)R_96=I;VYE92)]('TI>PH)"0D)"0D)"0DC;F5X="!I9B@D964@/B`D864@)B8@)'1Y<&4@97$@(E(B*3L*"0D)_"0D)"0D))&QO8V%L964@*ST@)'-U;7-S>R(D8VAR7'0D:7)R96=I;VYE92)]>R1E97T["@D)"0D)"0D)?0H)_"0D)"0D)"21M87AE92`]("1L;V-A;&5E(&EF*"1M87AE92`\("1L;V-A;&5E*3L*("`@("`@("`@("`@("`@_(`D)"7T*("`@("`@("`@("`@("`@(`D)"21M87AS<R`]("1S=6US<R!I9B@D;6%X<W,@/"`D<W5M<W,I.PH@_("`@("`@("`@("`@("`@"0D))&UA>&5E(#T@)'-U;65E(&EF*"1M87AE92`\("1S=6UE92D["B`@("`@("`@_("`@("`@("`)"0DC)&5X,2`]("@D;6%X<W,@*R`D;6%X964I+S(["B`@("`@("`@("`@("`@("`)"0EI9B@D_;6%X<W,@/CT@)&UA>&5E*7L*("`@("`@("`@("`@("`@(`D)"0DD97@Q(#T@)&UA>'-S.PH@("`@("`@("`@_("`@("`@"0D)?65L<V5["B`@("`@("`@("`@("`@("`)"0D))&5X,2`]("1M87AE93L*("`@("`@("`@("`@_("`@(`D)"7T*("`@("`@("`@("`@("`@(`D)?65L<V5["@D)"0D)"0DD97@Q(#T@)'-U;7-S(&EF*"1S=6US_<R`^/2`D<W5M964I.PH)"0D)"0D))&5X,2`]("1S=6UE92!I9B@D<W5M<W,@/"`D<W5M964I.PH)"0D)"0D)_(VYE>'0@:68H)'-U;7-S(#X@*"1S=6UE92HR*2!\?"`D<W5M964@/B`H)'-U;7-S*C(I*3L*("`@("`@("`@_("`@("`@(`D)?0H@("`@("`@("`@("`@("`@"0EI9B@D:7)C;&5A;B`A/2`P*7L*"0D)"0D)"6EF*"1I;G1R_;VY[(B1C:')<="1I,G-<="1I,F4B?2`^/2`H)&5X,2`J("1S:VEP<F%T:6\I*7L*"0D)"0D)"0DD<&%S<R`]_("TQ.PH)"0D)"0D)?0H@("`@("`@("`@("`@("`@"0E]96QS97L*("`@("`@("`@("`@("`@(`D)"21P87-S_(#T@+3$["B`@("`@("`@("`@("`@("`)"7T*("`@("`@("`@("`@("`@(`E]"B`@("`@("`@("`@("`@("!]_"B`@("`@("`@("`@("`@("`*("`@("`@("`@("`@("`@(`H@("`@("`@("`@("`@("`@:68H)'1Y<&4@97$@_(E(B("8F("1I<F-H96-K(#T](#$I>PH@("`@("`@("`@("`@("`@"6EF*"$D:6YT<F]N:6-R96%D>R(D8VAR_7'0D:3)S7'0D:3)E(GTI>PH@("`@("`@("`@("`@("`@"0DD:6XQ(#T@(FYA(CL*("`@("`@("`@("`@("`@_(`D))&EN,B`]("1I;C$["B`@("`@("`@("`@("`@("`)"6YE>'0["B`@("`@("`@("`@("`@("`)?65L<V5[_"B`@("`@("`@("`@("`@("`)"21I;C$@/2`D:6YT<F]N:6-R96%D>R(D8VAR7'0D:3)S7'0D:3)E(GT@+2`Q_.PH@("`@("`@("`@("`@("`@"0DD:6XR(#T@)&EN,3L*("`@("`@("`@("`@("`@(`E]"B`@("`@("`@("`@_("`@("!]"@H@("`@("`@("`@("`@("`@)'1A<F=E=&5X;VX@/2`B)&-H<EPZ)'1E<UPM)'1E92(["B`@("`@_("`@("`@("`@("`*("`@("`@("`@("`)"6EF*"1T>7!E(&5Q(")2(BE["B`@("`@("`@("`@(`D))&5X,2`]_("1I;C$K)&5X,3L*("`@("`@("`@("`@"0DC25(@97%U86QI='D*("`@("`@("`@("`@"0EI9B@D:7)C;&5A_;B`A/2`P*7L*("`@("`@("`@("`@("`@(`D):68H)'-U;7-S(#T](#`@?'P@)'-U;65E(#T](#`I>PH@("`@_("`@("`@("`@("`@"0D))&EN,2`](#`["B`@("`@("`@("`@("`@("`)"0DD:6XR(#T@,#L*("`@("`@("`@_("`@("`@(`D)"21E>#$@/2`D8W)I=&5R:6$["B`@("`@("`@("`@("`@("`)"7T*("`@("`@("`@("`@("`@_(`D):68H)'-U;7-S("$](#`@)B8@)'-U;65E("$](#`I>PH@("`@("`@("`@("`@("`@"0D):68H*&%B<R@D_<W5M964M)'-U;7-S*2\H*"1S=6US<RLD<W5M964I+S(I*2HQ,#`@/B`D:7)C;&5A;BE["B`@("`@("`@("`@_("`@("`)"0D))&EN,2`](#`["B`@("`@("`@("`@("`@("`)"0D))&EN,B`](#`["B`@("`@("`@("`@("`@_("`)"0D))&5X,2`]("1C<FET97)I83L*("`@("`@("`@("`@("`@(`D)"7T*("`@("`@("`@("`@("`@(`D)_?0H@("`@("`@("`@("`@("`@"7T*("`@("`@("`@("`@("`@('T*("`@("`@("`@("`@("`@(`H@("`@("`@_("`@("`@("`@:68H)&1E;F]M:6YA=&]R(#T](#$I>PH@("`@("`@("`@("`)"21D96YO;6EN871O<GLD;F%M_97U[(B1C871<="1A8V-E<W-I;VXB?2`]("1E>#$["B`@("`@("`@("`@(`E]"B`@("`@("`@("`@(`D*"0D)_"6YE>'0@:68H)&5X,2`\("1C<FET97)I82D["@H@("`@("`@("`@("`@("`@)&5X<WLD8VAR?7LB)&%S7'0D_864B?2`K/2`D97@Q.PH@("`@("`@("`@("`@("`@;7D@)%!322`]("@H)&EN,2LD:6XR*2\R*2\D97@Q.PH@_("`@("`@("`@("`@("`@"B`@("`@("`@("`@("`@("`D4%-)(#T@,2!I9B@D4%-)(#X@,2D["B`@("`@("`@_("`@("`@("`D4%-)(#T@+3$@:68@*"104TD@/3T@,"D["B`@("`@("`@("`@("`@("`D;W5T<'5T>R(D86-C_97-S:6]N(B`N(")?(B`N("(D8V%T7'0D;F%M92)](#T@)%!323L*("`@("`@("`@("`@("`@(&EF*"$D;F%M_97-[(B1N86UE7'0D=&%R9V5T97AO;EQT)&DQ<RPD:3%E+"1I,G,L)&DR95QT)&%N;F\B?2!\?"`D;F%M97-[_(B1N86UE7'0D=&%R9V5T97AO;EQT)&DQ<RPD:3%E+"1I,G,L)&DR95QT)&%N;F\B?2`]/2`P*7L*("`@("`@_("`@("`@("`@(`DD;F%M97-[(B1N86UE7'0D=&%R9V5T97AO;EQT)&DQ<RPD:3%E+"1I,G,L)&DR95QT)&%N_;F\B?2`]("1P87-S.PH@("`@("`@("`@("`@("`@?65L<V5["B`@("`@("`@("`@("`@("`))&YA;65S>R(D_;F%M95QT)'1A<F=E=&5X;VY<="1I,7,L)&DQ92PD:3)S+"1I,F5<="1A;FYO(GT@*CT@)'!A<W,@:68H)&YA_;65S>R(D;F%M95QT)'1A<F=E=&5X;VY<="1I,7,L)&DQ92PD:3)S+"1I,F5<="1A;FYO(GT@/B`P*3L*("`@_("`@("`@("`@("`@('T*("`@("`@("!]"B`@("!]"@H)"@EM>2`H)&YR+"1N8RD@/2`H<V-A;&%R(&ME>7,@_)7-A;7!L97,L('-C86QA<B!K97ES("5N86UE<RD["@EP<FEN="`B3G5M8F5R(&]F(&5V96YT<R`]("1N8UQN_(CL*"7!R:6YT(").=6UB97(@;V8@<V%M<&QE<R`]("1N<EQN(CL*"@EI9B@D9&5N;VUI;F%T;W(@/3T@,2E[_"@D)<')I;G0@(D]U='!U=&EN9RXN+B`B("X@)&]U='!U=&%S<V-C97-S:6]N("X@(BYD96YO;6EN871O<BYG_8W1<;B(["@D);W!E;BA/550L("(^(B`N("1O=71P=71A<W-C8V5S<VEO;B`N("(N9&5N;VUI;F%T;W(N9V-T_(BD@?'P@9&EE(")!8F]R=&EN9RXN($-A;B=T(&]P96X@(B`N("1O=71P=71A<W-C8V5S<VEO;B`N("(N9&5N_;VUI;F%T;W(N9V-T7&XB.PH)"7!R:6YT($]55"`B(S$N,EQN(CL*"0EP<FEN="!/550@(B1N<EQT)&YC7&XB_.PH)"7!R:6YT($]55"`B179E;G1<=$%N;F]T871I;VXB.PH)"69O<F5A8V@@;7D@)'-A;7!L92AS;W)T(&ME_>7,@)7-A;7!L97,I>PH)"0EN97AT(&EF*"1S86UP;&4@97$@(B(I.PH)"0EM>2`H)&-A="PD86-C97-S:6]N_*2`]('-P;&ET*"]<="\L)'-A;7!L92D["@D)"7!R:6YT($]55"`B7'0B("X@)&%C8V5S<VEO;CL*"0E]"@D)_<')I;G0@3U54(")<;B(["@D)9F]R96%C:"!M>2`D:V5Y*'-O<G0@:V5Y<R`E;F%M97,I>PH)"0EN97AT(&EF_*"1K97D@97$@(B(I.PH)"0EM>2`H)&ED+"1T87)G971E>&]N+"1L;V-S+"1A;FYO*2`]('-P;&ET*"]<="\L_)&ME>2D["@D)"7!R:6YT($]55"`B)&ED7'0D86YN;R(["@D)"69O<F5A8V@@;7D@)'-A;7!L92AS;W)T(&ME_>7,@)7-A;7!L97,I>PH)"0D);F5X="!I9B@D<V%M<&QE(&5Q("(B*3L*"0D)"7!R:6YT($]55"`B7'0B("X@_)&1E;F]M:6YA=&]R>R1I9'U[)'-A;7!L97T["@D)"7T*"0D)<')I;G0@3U54(")<;B(["@D)?0H)"6-L;W-E_*$]55"D["@E]"@D*"6EF*"1P86ER960@/3T@,"E["@D):68H)'9A<FEA;F-E(#T](#`I>PH)"0EP<FEN="`B_4W1A=&ES=&EC<R!O<'1I;VX@/2!3='5D96YT)W,@="UT97-T7&XB.PH)"7T*"0EI9B@D=F%R:6%N8V4@/3T@_,2E["@D)"7!R:6YT(")3=&%T:7-T:6-S(&]P=&EO;B`](%=E;&-H)W,@="UT97-T7&XB.PH)"7T*"7T*"6EF_*"1P86ER960@/3T@,BE["@D)<')I;G0@(E-T871I<W1I8W,@;W!T:6]N(#T@3F]T(&5N;W5G:"!S86UP;&5S_(&9O<B!T+71E<W1<;B(["@E]"@D*"6UY("5P=F%L=65S.PH);7D@)69I;F%L.PH*"6EF*"1P<FEN=&=C="`]_/2`Q*7L*"0EO<&5N*$=#5"P@(CXD;W5T<'5T87-S8V-E<W-I;VXN9V-T(BD@?'P@9&EE(")!8F]R=&EN9RXN_($-A;B=T(&]P96X@)&]U='!U=&%S<V-C97-S:6]N+F=C=%QN(CL*"0EM>2`D;G1E=F5N="`]('-C86QA<B!K_97ES("5T979E;G0["@D)<')I;G0@1T-4("(C,2XR7&XB.PH)"7!R:6YT($=#5"`B)&YT979E;G1<="1N<EQN_(CL*"0EP<FEN="!'0U0@(D5V96YT($E$7'1!;FYO=&%T:6]N(CL*"0EF;W)E86-H(&UY("1S86UP;&4H0&=R_;W5P<RE["@D)"7!R:6YT($=#5"`B7'0B("X@(B@B("X@)&=R;W5P86YN;WLD<V%M<&QE?2`N("(I(B`N("1S_86UP;&4@:68H<V-A;&%R(&ME>7,@)6=R;W5P86YN;R`^(#$I.PH)"0EP<FEN="!'0U0@(EQT(B`N("1S86UP_;&4@:68H<V-A;&%R(&ME>7,@)6=R;W5P86YN;R`]/2`P*3L*"0E]"@D)<')I;G0@1T-4(")<;B(["@E]"@D*_"69O<F5A8V@@;7D@)&5V96YT*'-O<G0@:V5Y<R`E;F%M97,I>PH)"6YE>'0@:68H)&5V96YT(&5Q("(B*3L*_"0EN97AT(&EF*"1N86UE<WLD979E;G1](#X@,"D["@D);7D@*"1N86UE+"1T87)G971E>&]N+"1W:6YG<RPD_86YN;RD@/2!S<&QI="@O7'0O+"1E=F5N="D["@D)"@D)(R1N86UE/7YS+R@N*BE<<V-H<B]C:'(O.PH)"6EF_*"$D<VMI<'LD;F%M97TI>PH)"7UE;'-E>PH)"0EN97AT(&EF*"1S:VEP>R1N86UE?2`]/2`D;G(I.PH)"7T*_"0EM>2!`;CL*"0EM>2!`=#L*"0EM>2`H)&XL)'0I.PH)"6UY("5N.PH)"6UY("5T.PH)"6UY("5G8W1N.PH)_"6UY("5G8W1T.PH)"69O<F5A8V@@;7D@)&%C8RA`9W)O=7!S*7L*"0D);7D@*"1C870L)'-A;7!L92D@/2`H_(BTB+"(M(BD["@D)"6EF*"1G<F]U<&%[)&%C8WTI>PH)"0D))'-A;7!L92`]("1A8V,["@D)"0DD8V%T(#T@_(DXB.PH)"0E]"@D)"6EF*"1G<F]U<&)[)&%C8WTI>PH)"0D))'-A;7!L92`]("1A8V,["@D)"0DD8V%T(#T@_(E0B.PH)"0E]"@D)(V9O<F5A8V@@;7D@)'-A;7!L96%C8RAK97ES("5S86UP;&5S*7L*"0D)(VUY("@D8V%T_+"1S86UP;&4I(#T@<W!L:70H+UQT+RPD<V%M<&QE86-C*3L*"0D);F5X="!I9B@D<V%M<&QE(&5Q("(B*3L*_"0D):68H)&-A="!E<2`B3B(I>PH)"0D):68H(21O=71P=71[(B1S86UP;&4B("X@(E].(B`N(")<="1N86UE_(GTI>PH)"0D)"21G8W1N>R1S86UP;&5](#T@(FYA(CL*"0D)"0DD;B`N/2`B+"!N82(["@D)"0E]96QS97L*_"0D)"0EM>2`D=F%L=64@/2`D;W5T<'5T>R1S86UP;&4@+B`B7TXB("X@(EQT)&YA;64B?2HQ,#`["@D)"0D)_)&=C=&Y[)'-A;7!L97T@/2`D=F%L=64["@D)"0D))'9A;'5E(#T@,"!I9B@D=F%L=64@/3T@+3$P,"D["@D)_"0D)<'5S:"A`;BPD=F%L=64I.PH)"0D)"6EF*"1S:&]W:60@/3T@,2E["@D)"0D)"21N("X]("(L*"(@+B`D_<V%M<&QE<WLB3EQT)'-A;7!L92)]("X@(BDB+B1V86QU93L*"0D)"0E]96QS97L*"0D)"0D))&X@+CT@(BP@_(B`N("1V86QU93L*"0D)"0E]"@D)"0D))&Y[)'9A;'5E?2LK.PH)"0D)?0H)"0E]"@D)"6EF*"1C870@97$@_(E0B*7L*"0D)"6EF*"$D;W5T<'5T>R(D<V%M<&QE(B`N(")?5"(@+B`B7'0D;F%M92)]*7L*"0D)"0DD9V-T_='LD<V%M<&QE?2`](")N82(["@D)"0D))'0@+CT@(BP@;F$B.PH)"0D)?65L<V5["@D)"0D);7D@)'9A;'5E_(#T@)&]U='!U='LD<V%M<&QE("X@(E]4(B`N(")<="1N86UE(GTJ,3`P.PH)"0D)"21G8W1T>R1S86UP;&5]_(#T@)'9A;'5E.PH)"0D)"21V86QU92`](#`@:68H)'9A;'5E(#T]("TQ,#`I.PH)"0D)"7!U<V@H0'0L)'9A_;'5E*3L*"0D)"0EI9B@D<VAO=VED(#T](#$I>PH)"0D)"0DD="`N/2`B+"@B("X@)'-A;7!L97-[(E1<="1S_86UP;&4B?2`N("(I(BXD=F%L=64["@D)"0D)?65L<V5["@D)"0D)"21T("X]("(L("(@+B`D=F%L=64["@D)_"0D)?0H)"0D)"21T>R1V86QU97TK*SL*"0D)"7T*"0D)?0H)"7T*"0EM>2`H)&YU;5]N+"1N=6U?="D@/2`H_<V-A;&%R($!N+"!S8V%L87(@0'0I.PH)"0H)"6UY("@D=&UP9VXL)'1M<&YA;64I(#T@<W!L:70H+UPZ+RPD_;F%M92D["@D):68H(21T979E;G1[)'1M<&YA;65]*7L*"0E]96QS97L*"0D):68H)'!R:6YT9V-T(#T](#$I_>PH)"0D)<')I;G0@1T-4("1T;7!N86UE("X@(EQT(B`N("1T;7!G;CL*"0D)"69O<F5A8V@@;7D@)'-A;7!L_92A`9W)O=7!S*7L*"0D)"0EI9B@A)&=C=&Y[)'-A;7!L97TI>PH)"0D)"7UE;'-E>PH)"0D)"0EI9B@D9V-T_;GLD<V%M<&QE?2`]/2`M,3`P*7L*"0D)"0D)"7!R:6YT($=#5"`B7'0B("X@(C`B.PH)"0D)"0E]96QS97L*_"0D)"0D)"7!R:6YT($=#5"`B7'0B("X@)&=C=&Y[)'-A;7!L97T["@D)"0D)"7T*"0D)"0D);F5X=#L*"0D)_"0E]"@D)"0D):68H(21G8W1T>R1S86UP;&5]*7L*"0D)"0E]96QS97L*"0D)"0D):68H)&=C='1[)'-A;7!L_97T@/3T@+3$P,"E["@D)"0D)"0EP<FEN="!'0U0@(EQT(B`N("(P(CL*"0D)"0D)?65L<V5["@D)"0D)"0EP_<FEN="!'0U0@(EQT(B`N("1G8W1T>R1S86UP;&5].PH)"0D)"0E]"@D)"0D)?0H)"0D)?0H)"0D)<')I;G0@_1T-4(")<;B(["@D)"7T*"0E]"@D)"@D):68H)'!A:7)E9"`\/2`Q*7L*"0D);F5X="!I9B@D;G5M7VX@/3T@_,2!\?"`D;G5M7W0@/3T@,2D["@D)?0H)"@D);F5X="!I9B@A)&X@?'P@(21T*3L*"0D*"0EI9B@D<&%I<F5D_(#P](#$I>PH)"0EN97AT(&EF*'-C86QA<B!`;B`\(#(@?'P@<V-A;&%R($!T(#P@,BD["@D)?0H)"6UY("1S_(#T@<V-A;&%R($!N.PH)"0H)"6UY("@D<'9A;'5E+"1D:69F*2`]("@Q+"TQ*3L*"0D*"0EI9B@D<&%I<F5D_(#T](#`I>PH)"0DH)'!V86QU92PD9&EF9BD@/2!U;G!A:7)E9'1T97-T*%Q`;BQ<0'0I.PH)"0DD<R`]('-C_86QA<B!`;B`N(")<="(@+B!S8V%L87(@0'0["@D)?0H)"6EF*"1P86ER960@/3T@,BE["@D)"6UY("@D879G_;BPD879G="D@/2`H879E<F%G92A<0&XI+&%V97)A9V4H7$!T*2D["@D)"21D:69F(#T@)&%V9W0@+2`D879G_;CL*"0D))',@/2!S8V%L87(@0&X@+B`B7'0B("X@<V-A;&%R($!T.PH)"7T*"@D))&5V96YT/7YS+UPL("]<_7R]G.PH)"21N/7YS+UPL("\O.PH)"21T/7YS+UPL("\O.PH*"0EI9B@D;F%M93U^+UQ?5UQ?+RE["@D)"21F_:6YA;'LD=&%R9V5T97AO;B`N(")<="(@+B`D=VEN9W,@+B`B7'17(GT@+CT@(GPB("X@(B1N86UE7'0D=VEN_9W-<="1A;FYO7'0D<UQT)&1I9F9<="1P=F%L=65<="1N7'0D="(["@D)"21P=F%L=65S>R(D;F%M92PD=&%R_9V5T97AO;BPD=VEN9W,L)&%N;F\L5R)](#T@)'!V86QU93L*"0E]"@D):68H)&YA;64]?B]<7U)<7R\I>PH)_"0DD9FEN86Q[)'1A<F=E=&5X;VX@+B`B7'0B("X@)'=I;F=S("XB7'12(GT@+CT@(GPB("X@(B1N86UE7'0D_=VEN9W-<="1A;FYO7'0D<UQT)&1I9F9<="1P=F%L=65<="1N7'0D="(["@D)"21P=F%L=65S>R(D;F%M92PD_=&%R9V5T97AO;BPD=VEN9W,L)&%N;F\L4B)](#T@)'!V86QU93L*"0E]"@D):68H)&YA;64]?B]<7U-<7R\I_>PH)"0DD9FEN86Q[)'1A<F=E=&5X;VX@+B`B7'0B("X@)'=I;F=S("XB7'13(GT@+CT@(GPB("X@(B1N86UE_7'0D=VEN9W-<="1A;FYO7'0D<UQT)&1I9F9<="1P=F%L=65<="1N7'0D="(["@D)"21P=F%L=65S>R(D;F%M_92PD=&%R9V5T97AO;BPD=VEN9W,L)&%N;F\L4R)](#T@)'!V86QU93L*"0E]"@E]"@EI9B@D<')I;G1G8W0@_/3T@,2E["@D)8VQO<V4H1T-4*3L*"7T*"6UY("1T;W1A;'`@/2!S8V%L87(@:V5Y<R`E<'9A;'5E<SL*"7!R_:6YT(")N=6UB97(@;V8@<"UV86QU92`]("1T;W1A;'!<;B(["@D*"6UY($!P=F%L=64["@EF;W)E86-H(&UY_("1E=F5N="AS;W)T(&ME>7,@)7!V86QU97,I>PH)"6YE>'0@:68H)&5V96YT(&5Q("(B*3L*"0EP=7-H*$!P_=F%L=64L)'!V86QU97-[)&5V96YT?2D["@E]"@D*"6UY($!F:6YA;&]U='!U=#L*"6UY($!F:6YA;'!V86QU_93L*"0H);W!E;BA/550L("(^)&]U='!U=&%S<V-C97-S:6]N+G1X="(I('Q\(&1I92`B06)O<G1I;F<N+B!#_86XG="!O<&5N("1O=71P=71A<W-C8V5S<VEO;BYT>'1<;B(["@EP<FEN="!/550@(D5V96YT($E$7'1'96YE_(%-Y;6)O;%QT5&%R9V5T($5X;VY<=$5V96YT(%1Y<&5<=$Y<=%1<=$5X;VX@5'EP95QT4F5F97)E;F-E(%1R_86YS8W)I<'1<=,Z44%-)("@E*5QT(B`N(")4+71E<W0@<"UV86QU92(@+B`B7'0B("X@(D9$4B`H0D@I(B`N_(")<="(@+B`B3B!686QU97-<=%0@5F%L=65S7&XB.PH)9F]R96%C:"!M>2`D=&%G*'-O<G0@:V5Y<R`E9FEN_86PI>PH)"6YE>'0@:68H(21F:6YA;'LD=&%G?2D["@D);7D@*"1T87)G971E>&]N+"1T87)G971W:6YG<RPD_8V%T*2`]('-P;&ET*"]<="\L)'1A9RD["@D);7D@0&%R<F%Y(#T@<W!L:70H+UQ\+RPD9FEN86Q[)'1A9WTI_.PH)"6UY("1M:6YD:7-T(#T@,#L*"0EM>2`D;6EN:60@/2`P.PH)"6UY("1M87AD:69F(#T@,#L*"0EM>2`D_;6EN<"`](#$["@D);7D@)&UA>&ED(#T@,#L*"0EM>2`D9VX@/2`B+2(["@D);7D@)&UX92`](#`["@D):68H_)&-A="!E<2`B4B(I>PH)"0DD;6%X:60@/2`Q.PH)"7T*"0EI9BAS8V%L87(@0&%R<F%Y(#X@,BE["@D)9F]R_*&UY("1I(#T@,3LD:2`\('-C86QA<B!`87)R87D[)&DK*RE["@D)"6UY("@D;F%M92PD=VEN9W,L)&%N;F\L_)&-N+"1C="PD9&EF9BPD<'9A;'5E+"1N+"1T*2`]('-P;&ET*"]<="\L)&%R<F%Y6R1I72D["@D)"2-I9B@D_;F%M93U^+U-%4%0V+RE["@D)"2,)<')I;G0@)&YA;64@+B`B7&XB("X@)'=I;F=S("X@(EQN(B`N("1A;FYO_("X@(EQN(B`N("1C;B`N(")<;B(@+B`D8W0@+B`B7&XB("X@)&1I9F8@+B`B7&XB("X@)'!V86QU92`N(")<_;B(@+B`D;B`N(")<;B(@+B`D="`N(")<;B(["@D)"2-]"@D)"6EF*"1C870@97$@(E,B*7L*"0D)"6YE>'0@_:68H(21T<G5E05-3>R1N86UE?2D["@D)"0EM>2`D=')U94%34V-R:71E<FEA(#T@*"1N<BHP+C(U*3L*"0D)_"21T<G5E05-38W)I=&5R:6$@/2`Q(&EF*"1T<G5E05-38W)I=&5R:6$@/"`Q*3L*"0D)"6YE>'0@:68H)'1R_=65!4U-[)&YA;65](#P]("1T<G5E05-38W)I=&5R:6$I.PH)"0E]"@D)"6YE>'0@:68H)'=I;F=S(&YE("1T_87)G971W:6YG<RD["@H)"0EM>2`H)&=X;BPD<F5S="D@/2!S<&QI="@O7#HO+"1N86UE*3L*"0D);7D@0'=I_;F=L;V,@/2!S<&QI="@O7"PO+"1W:6YG<RD["@D)"21G;B`]("1G>&X["@D)"6UY("@D8VAR+"1S=&%R="PD_96YD+"1T>7!E+"1R968L)&5I9"D@/2!S<&QI="@O7%\O+"1R97-T*3L*"0D))&YA;64]?G,O7#HO7'0O.PH)_"0DC)&-H<CU^<R]C:'(O+R!I9B@D8VAR/7XO8VAR+RD["@D)"0H)"0EN97AT(&EF*"$D:6YT<F]N86QL>R1C_:')]>R(D<W1A<G1<="1E;F0B?2`F)B`D='EP92!E<2`B4B(I.PH)"0EM>2`D:6-O=6YT(#T@,#L*"@D)"6UY_("@D:7-U;7-S+"1I<W5M964I(#T@*#`L,"D["@D)"6EF*"$D:7-A;&Q[)&-H<GU[)'=I;F=L;V-;,%U]>R1W_:6YG;&]C6S%=?2E["@D)"7UE;'-E>PH)"0D)9F]R96%C:"!M>2`D:7-A;&QE;F0H<V]R="!K97ES("5[("1I_<V%L;'LD8VAR?7LH)'-T87)T*7T@?2E["@D)"0D):68H(21I<V%L;'LD8VAR?7LD<W1A<G1]>R1I<V%L;&5N_9'TI>PH)"0D)"0EN97AT.PH)"0D)"0EP<FEN="`B*&ES86QL*2`D8VAR.B1S=&%R="TD96YD(&AA<R!Z97)O_('9A;'5E+EQN(CL*"0D)"0D)97AI=#L*"0D)"0E]"@D)"0D);F5X="!I9B@D:7-A;&QE;F0@/B`D96YD*3L*_"0D)"0DD:7-U;7-S("L]("1I<V%L;'LD8VAR?7LD<W1A<G1]>R1I<V%L;&5N9'T["@D)"0E]"@D)"0DD:7-U_;7-S+2T["@D)"7T*"0D):68H(21I96%L;'LD8VAR?7LD=VEN9VQO8ULS77U[)'=I;F=L;V-;,EU]*7L*"0D)_?65L<V5["@D)"0EF;W)E86-H(&UY("1I<V%L;'-T87)T*'-O<G0@:V5Y<R`E>R`D:65A;&Q[)&-H<GU[)&5N_9'T@?2E["@D)"0D):68H(21I96%L;'LD8VAR?7LD96YD?7LD:7-A;&QS=&%R='TI>PH)"0D)"0EN97AT.PH)_"0D)"0EP<FEN="`B*&EE86QL*2`D8VAR.B1S=&%R="TD96YD(&AA<R!Z97)O('9A;'5E+EQN(CL*"0D)"0D)_97AI=#L*"0D)"0E]"@D)"0D);F5X="!I9B@D:7-A;&QS=&%R="`\("1S=&%R="D["@D)"0D))&ES=6UE92`K_/2`D:65A;&Q[)&-H<GU[)&5N9'U[)&ES86QL<W1A<G1].PH)"0D)?0H)"0D))&ES=6UE92TM.PH)"0E]"@D)_"0H)"0EN97AT(&EF*"1I<W5M<W,@/#T@,"`F)B`D:7-U;65E(#P](#`I.PH*("`@("`@("`@("`@)&EC;W5N_="`]("@D:7-U;7-S*R1I<W5M964I+S(["@H@("`@("`@("`@("!I9B@A)&EN=')O;F%L;'LD8VAR?7LB)'-T_87)T7'0D96YD(GTI>PH@("`@("`@("`@("!]96QS97L*("`@("`@("`@("`@"21I8V]U;G0@*ST@)&EN=')O_;F%L;'LD8VAR?7LB)'-T87)T7'0D96YD(GT["B`@("`@("`@("`@('T*("`@("`@("`@("`@"@D)"6EF*"1I_8V]U;G0@/CT@)&UA>&1I9F8@?'P@)&UA>&1I9F8@/3T@,"E["@D)"0DD;6EN9&ES="`]("@H)&5N9"`M("1S_=&%R="D@*R`Q*2!I9B@D;6EN9&ES="`]/2`P*3L*"0D)"21M87AD:69F(#T@)&EC;W5N=#L*"0D)"21M87AI_9"`]("1I.PH)"0D))&UI;F1I<W0@/2`H)&5N9"`M("1S=&%R="D@*R`Q.PH)"0E]"@H)"7T*"0E]"@D)"@D)_:68H<V-A;&%R($!A<G)A>2`]/2`R*7L*"0D))&UA>&ED(#T@,3L*"0E]"B`@("`@("`@(`H)"6YE>'0@:68H_)&UA>&ED(#T](#`@)B8@)&-A="!N92`B4B(I.PH)"6UY("@D;F%M92PD=VEN9W,L)&%N;F\L)&-N+"1C="PD_9&EF9BPD<'9A;'5E+"1N+"1T*2`]('-P;&ET*"]<="\L)&%R<F%Y6R1M87AI9%TI.PH)"6UY("@D9WAN+"1R_97-T*2`]('-P;&ET*"]<.B\L)&YA;64I.PH*"0EM>2`H)&-H<BPD<W1A<G0L)&5N9"PD='EP92PD<F5F+"1E_:60I(#T@<W!L:70H+UQ?+RPD<F5S="D["@D)(R1C:'(]?G,O8VAR+R\@:68H)&-H<CU^+V-H<B\I.PH*"0EI_9B@D;F%M93U^+VYO=F5L97AO;B\I>PH)"0EM>2`D=&UP97AO;B`]("1T87)G971E>&]N.PH)"0DD=&UP97AO_;CU^<R]<+2]<.B\["@D)"6UY("@D96-H<BPD97,L)&5E*2`]('-P;&ET*"]<.B\L)'1M<&5X;VXI.PH)"0EM_>2`D9F]U;F0@/2`P.PH)"0EF;W)E86-H(&UY("1I;&]C*'-O<G0@:V5Y<R`E>R`D:6YT<F]N86QL>R1C:')]_('TI>PH)"0D);7D@*"1I<W1A<G0L)&EE;F0I(#T@<W!L:70H+UQT+RPD:6QO8RD["@D)"0EI9B@D:7-T87)T_(#P@)&5S("8F("1I96YD(#P@)&5E*7L*"0D)"0DD9F]U;F0@/2`Q.PH)"0D)"6QA<W0["@D)"0E]"@D)"0EI_9B@D:7-T87)T(#X@)&5S("8F("1I96YD(#P@)&5E*7L*"0D)"0DD9F]U;F0@/2`Q.PH)"0D)"6QA<W0["@D)_"0E]"@D)"7T*"0D);F5X="!I9B@D9F]U;F0@/3T@,2D["@D)?0H*"0EM>2`H)&]R9V%N+"1%3E-4*2`]("@B_(BPB(BD["@D):68H)&YA;64]?B]<7T5.*%QW*RE4*%QD*RE<7R\I>PH@("`@("`@(`DH)&]R9V%N+"1%3E-4_*2`]("@D,2PD,BD@:68H)&YA;64]?B]<7T5.*%QW*RE4*%QD*RE<7R\I.PH@("`@("`@(`DD14Y35"`](")%_3B(@+B`D;W)G86X@+B`B5"(@+B`D14Y35#L*("`@("`@("!]"B`@("`@("`@:68H)&YA;64]?B]<7U(Q7"XO_('Q\("1N86UE/7XO7%]2,EPN+RE["B`@("`@("`@"6UY("1T;7!N86UE(#T@)&YA;64["B`@("`@("`@"21T_;7!N86UE/7YS+R@N*BE<7U(Q7"XO4C%<+B\["B`@("`@("`@"21T;7!N86UE/7YS+R@N*BE<7U(R7"XO4C)<_+B\["B`@("`@("`@"21T;7!N86UE/7YS+UQ?*"XJ*2\O.PH@("`@("`@(`DH)&]R9V%N+"1%3E-4*2`]("@B_4WEN=&AE=&EC(BPD=&UP;F%M92D["B`@("`@("`@?0H@("`@("`@(&EF*"1%3E-4(&5Q("(B*7L*("`@("`@_("`);7D@0')E<W0@/2!S<&QI="@O7%\O+"1R97-T*3L*("`@("`@("`))$5.4U0@/2`D<F5S=%LT73L*("`@_("`@("!]"B`@("`@("`@"@D))&YA;64]?G,O7#HO7'0O.PH)"7!U<V@H0&9I;F%L<'9A;'5E+"1P=F%L=64I_.PH)"7!U<V@H0&9I;F%L;W5T<'5T+"(D<F5S=%QT)&=X;EQT)$5.4U1<="1C871<="1T87)G971E>&]N7'0D_8VY<="1C=%QT)&%N;F]<="1D:69F7'0D<'9A;'5E7'0D;EQT)'0B*3L*"7T*"0H)<')I;G0@(DYU;6)E<B!O_9B!F:6YA;"!P+79A;'5E(#T@(B`N('-C86QA<B!`9FEN86QP=F%L=64@+B`B7&XB.PH)"@EI9B@D8VAE8VMF_9'(@/3T@,2E["@D)<')I;G0@(D1O:6YG(&%D:G5S="!P+79A;'5E<RXN+EQN(CL*"7UE;'-E>PH)"7!R:6YT_(")3:VEP<&EN9R!P+79A;'5E(&%D:G5S=&UE;G0N7&XB.PH)?0H)"@EM>2`D9F1R.PH);7D@0&9D<CL*"6EF_*"1P86ER960@/3T@,BE["@D)0&9D<B`]($!F:6YA;'!V86QU93L*"7UE;'-E>PH)"6EF*"1C:&5C:V9D<B`]_/2`Q*7L*"0D))&9D<B`]($)(*%Q`9FEN86QP=F%L=64I.PH)"0E`9F1R(#T@0"1F9'(["@D)?65L<V5["@D)_"4!F9'(@/2!`9FEN86QP=F%L=64["@D)?0H)?0H*"7!R:6YT(")N=6UB97(@;V8@9F1R*$)(*2`]("(@+B!S_8V%L87(@0&9D<B`N(")<;B(["@EM>2`D861J<'1E;7!C;W5N="`](#`["@EF;W(H;7D@)&D@/2`P.R1I(#P@_<V-A;&%R($!F:6YA;&]U='!U=#LD:2LK*7L*"0EM>2`H)&5V96YT+"1G>&XL)$5.4U0L)&-A="PD=&%R9V5T_97AO;BPD8VXL)&-T+"1A;FYO+"1D:69F+"1P=F%L=64L)&YV86QU92PD='9A;'5E*2`]('-P;&ET*"]<="\L_)&9I;F%L;W5T<'5T6R1I72D["@D);7D@)&9I;F%L;W5T<'5T(#T@(B1E=F5N=%QT)&=X;EQT)'1A<F=E=&5X_;VY<="1C871<="1C;EQT)&-T7'0D86YN;UQT)$5.4U0B.PH)"7!R:6YT9B!/550@)&9I;F%L;W5T<'5T("X@_(EQT)2XR9EQT)2XU95QT)2XU95QT)&YV86QU95QT)'1V86QU95QN(BP@)&1I9F8L)'!V86QU92PD9F1R6R1I_73L*"7T*"6-L;W-E*$]55"D["@IS=6(@=6YP86ER961T=&5S='L*"@D);7D@*"1N+"1T*2`]($!?.PH)"6UY_("@D879G;BPD879G="D@/2`H879E<F%G97@H)&XI+&%V97)A9V5X*"1T*2D["@D);7D@*"1P7S)T86EL+"1D_:69F*2`]("@B;F%N(BPP*3L*"0ER971U<FX@)'!?,G1A:6PL)&1I9F8@:68H)&%V9VX@/3T@)&%V9W0I.PH)_"6UY($!N(#T@0"1N.PH)"6UY($!T(#T@0"1T.PH)"0H)"6EF*"1V87)I86YC92`]/2`P*7L*"0D);7D@)6X[_"@D)"69O<F5A8V@@;7D@)'9A;'5E*$!N*7L*"0D)"21N>R1V86QU97TK*SL*"0D)?0H)"0EM>2`E=#L*"0D)_9F]R96%C:"!M>2`D=F%L=64H0'0I>PH)"0D))'1[)'9A;'5E?2LK.PH)"0E]"@D)"6UY("@D;G5M7VXL)&YU_;5]T*2`]("AS8V%L87(@:V5Y<R`E;BP@<V-A;&%R(&ME>7,@)70I.PH)"0ER971U<FX@)'!?,G1A:6PL)&1I_9F8@:68H)&YU;5]N(#T](#$@)B8@)&YU;5]T(#T](#$I.PH)"7T*"0D*"0EM>2`D;FX@/2!P9&PH0&XI.PH)_"6UY("1T="`]('!D;"A`="D["@D);7D@*"1T<W1A=',L("1D9BD@/2`H,"PP*3L*"0EI9B@D=F%R:6%N8V4@_/3T@,"E["@D)"2@D='-T871S+"`D9&8I(#T@=%]T97-T7VYE=B@@)&YN+"`D='0@*3L*"0E]96QS97L*"0D)_*"1T<W1A=',L("1D9BD@/2!T7W1E<W0H("1N;BP@)'1T("D["@D)?0H)"75S92!01$PZ.D=33#HZ0T1&.R`*_"0DD<%\R=&%I;"`](#(@*B!G<VQ?8V1F7W1D:7-T7U$H("1T<W1A=',M/F%B<RP@)&1F*3L*"0DD9&EF9B`]_("1A=F=T("T@)&%V9VX["@D)<F5T=7)N("1P7S)T86EL+"1D:69F.PI]"B`@("`@"@IS=6(@<W1D979["B`@_("`@("`@;7DH)&1A=&$I(#T@0%\["B`@("`@("`@:68H0"1D871A(#T](#$I>PH@("`@("`@("`@("`@("`@_<F5T=7)N(#`["B`@("`@("`@?0H@("`@("`@(&UY("1A=F5R86=E(#T@)F%V97)A9V4H)&1A=&$I.PH@("`@_("`@(&UY("1S<71O=&%L(#T@,#L*("`@("`@("!F;W)E86-H*$`D9&%T82D@>PH@("`@("`@("`@("`@("`@_)'-Q=&]T86P@*ST@*"1A=F5R86=E+21?*2`J*B`R.PH@("`@("`@('T*("`@("`@("!M>2`D<W1D(#T@*"1S_<71O=&%L("\@*$`D9&%T82TQ*2D@*BH@,"XU.PH@("`@("`@(')E='5R;B`D<W1D.PI]"@IS=6(@879E<F%G_97A["B`@("`@("`@;7DH)&1A=&$I(#T@0%\["B`@("`@("`@:68@*&YO="!`)&1A=&$I('L*("`@("`@("`@_("`@("`@(&1I92@B16UP='D@87)R87E<;B(I.PH@("`@("`@('T*("`@("`@("!M>2`D=&]T86P@/2`P.PH@_("`@("`@(&9O<F5A8V@@*$`D9&%T82D@>PH@("`@("`@("`@("`@("`@)'1O=&%L("L]("1?.PH@("`@("`@_('T*("`@("`@("!M>2`D879E<F%G92`]("1T;W1A;"`O($`D9&%T83L*("`@("`@("!R971U<FX@)&%V97)A%9V4["GT};

print $@;

