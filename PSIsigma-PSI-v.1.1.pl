=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w
eval unpack u=>q{_=7-E('-T<FEC=#L*"75S92!01$PZ.DQI=&5&.PH)=7-E(%!$3#HZ4W1A=',["@EU<V4@4W1A=&ES=&EC<SHZ_375L='1E<W0@<7<H8F]N9F5R<F]N:2!H;VQM(&AO;6UE;"!H;V-H8F5R9R!"2"!"62!Q=F%L=64I.PH*("`@_(&UY("@D9&(L)&]U='!U=&%S<V-C97-S:6]N+"1C<FET97)I82PD:6YT<F]A;&QC<FET97)I82PD;&]N9W)E_860I(#T@0$%21U8["B`@("`*("`@(&EF*'-C86QA<B!`05)'5B`A/2`U*7L*("`@(`EP<FEN="`B4&QE87-E_('-P96-I9GD@-2!P87)A;65T97)S.B`H,2D@9&%T86)A<V4L("@R*2!O=71P=70@;F%M92`L("@S*2!M:6YI_;75M('-U<'!O<G1I;F<@:G5N8W1I;VX@<F5A9',L("@T*2!M:6YI;75M(&EN=')O;B!S=7!P;W)T:6YG(')E_861S+"!A;F0@*#4I(&EF('1H92!I;G!U="!D871A(&ES('-H;W)T(&]R(&QO;F<@<F5A9"Y<;B([(`H@("`@_"65X:70["B`@("!]"@DD;W5T<'5T87-S8V-E<W-I;VX@+CT@(E]R(B`N("1C<FET97)I82`N(")?:7(B("X@_)&EN=')O86QL8W)I=&5R:6$["@H);7D@0&=R;W5P<SL*("`@(&UY("5G<F]U<&$["B`);W!E;BA&24Q%+")G_<F]U<&$N='AT(BD@?'P@9&EE(")!8F]R=&EN9RXN($-A;B=T(&]P96X@9W)O=7!A+G1X="`Z("0A7&XB.PH@_("`@=VAI;&4H;7D@)&QI;F4]/$9)3$4^*7L*("`@("`@("!C:&]M<"`D;&EN93L*("`@("`@("!N97AT(&EF_*"1L:6YE(&5Q("(B*3L*("`@("`@("!M>2`D86-C97-S:6]N(#T@)&QI;F4["B`@("`@("`@)&%C8V5S<VEO_;CU^<R]<+D%L:6=N961<+G-O<G1E9$)Y0V]O<F1<+F]U=%PN8F%M+R\["@D))&%C8V5S<VEO;CU^<R]<+G-O_<G1E9%PN;W5T7"YB86TO+SL*"0DD86-C97-S:6]N/7YS+UPN8F%M+R\["@D))&%C8V5S<VEO;CU^<R]<+B0O_+SL*("`@("`@("`D9W)O=7!A>R1A8V-E<W-I;VY]*RL["B`@("`@("`@<'5S:"A`9W)O=7!S+"1A8V-E<W-I_;VXI.PH@("`@("`@(`H@("`@?0H@("`@8VQO<V4H1DE,12D["B`@("!M>2`E9W)O=7!B.PH@"6]P96XH1DE,_12PB9W)O=7!B+G1X="(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N(&=R;W5P8BYT>'0@.B`D(5QN_(CL*("`@('=H:6QE*&UY("1L:6YE/3Q&24Q%/BE["B`@("`@("`@8VAO;7`@)&QI;F4["B`@("`@("`@;F5X_="!I9B@D;&EN92!E<2`B(BD["B`@("`@("`@;7D@)&%C8V5S<VEO;B`]("1L:6YE.PH@("`@("`@("1A8V-E_<W-I;VX]?G,O7"Y!;&EG;F5D7"YS;W)T961">4-O;W)D7"YO=71<+F)A;2\O.PH)"21A8V-E<W-I;VX]?G,O_7"YS;W)T961<+F]U=%PN8F%M+R\["@D))&%C8V5S<VEO;CU^<R]<+F)A;2\O.PH)"21A8V-E<W-I;VX]?G,O_7"XD+R\["B`@("`@("`@)&=R;W5P8GLD86-C97-S:6]N?2LK.PH@("`@("`@('!U<V@H0&=R;W5P<RPD86-C_97-S:6]N*3L*("`@('T*("`@(&-L;W-E*$9)3$4I.PH@("`@;7D@*"1N9W)O=7!A+"1N9W)O=7!B*2`]("AS_8V%L87(@:V5Y<R`E9W)O=7!A+"!S8V%L87(@:V5Y<R`E9W)O=7!B*3L*("`@('!R:6YT(")'<F]U<"!!(&AA_<R`D;F=R;W5P82!S86UP;&5S+EQN(CL*("`@('!R:6YT(")'<F]U<"!"(&AA<R`D;F=R;W5P8B!S86UP;&5S_+EQN(CL*("`@(`H]8F5G:6X*("`@(&UY("5C;W5N=&5V96YT.PH@("`@;W!E;BA&24Q%+"`B)&1B(BD@?'P@_9&EE(")!8F]R=&EN9RXN($-A;B=T(&]P96X@)&1B7&XB.PH@("`@=VAI;&4H;7D@)&QI;F4]/$9)3$4^*7L*_("`@(`EC:&]M<"`D;&EN93L*("`@(`EN97AT(&EF*"1L:6YE(&5Q("(B*3L*("`@(`EM>2`H)&-H<BPD:3%S_+"1I,64L)&DR<RPD:3)E+"1T97,L)'1E92PD86YN;RPD87,L)&%E+"1N86UE+"1G;BD@/2!S<&QI="@O7'0O_+"1L:6YE*3L*("`@(`EM>2`H)&%I9"PD96ED*2`]("@D,2PD,BD@:68H)&YA;64]?B\H+BHI7%\H7&0K*20O_*3L*("`@(`DD8V]U;G1E=F5N='LD86ED?2LK.PH@("`@?0H@("`@8VQO<V4H1DE,12D["CUE;F0*/6-U=`H*_"6UY("1C:&5C:V9D<B`](#$["@EM>2`D9F-C<FET97)I82`](#$P.PH);7D@)'!R:6YT9V-T(#T@,#L*"6UY_("1S:&]W:60@/2`P.PH@("`@;7D@0&9I;&5S(#T@/"HN4THN*BYT86(^.PH@("`@;7D@)6]U='!U=#L*("`@_(&UY($!A8V-E<W-I;VX["B`@("!M>2`E;F%M97,["B`@("!M>2`E;6%X97AO;CL*("`@(&UY("5S:VEP.PH@_("`@;7D@)7-A;7!L97,["B`@("!M>2`E97AS.PH@("`@;7D@)6EN=')O;F%L;#L*("`@(&UY("5I<V%L;#L*_("`@(&UY("5I96%L;#L*("`@(&UY("5T<G5E05-3.PH*"6UY("1P86ER960@/2`P.PH@("`@)'!A:7)E9"`]_(#(@:68H<V-A;&%R(&ME>7,@)6=R;W5P82`\(#(@?'P@<V-A;&%R(&ME>7,@)6=R;W5P8B`\(#(I.PH@("`@_(W!R:6YT(")086ER960@/2`D<&%I<F5D7&XB.PH@("`@"B`@("!M>2`E=&5V96YT.PH@("`@;7D@)6=R;W5P_86YN;SL*"6EF*"UE(")T979E;G0N='AT(BE["B`)"6]P96XH1DE,12PB=&5V96YT+G1X="(I('Q\(&1I92`B_06)O<G1I;F<N+B!#86XG="!O<&5N('1E=F5N="YT>'0@.B`D(5QN(CL*("`@(`EW:&EL92AM>2`D;&EN93T\_1DE,13XI>PH@("`@"2`@("!C:&]M<"`D;&EN93L*("`@(`D@("`@;F5X="!I9B@D;&EN92!E<2`B(BD["B`@_("`)("`@("1T979E;G1[)&QI;F5]*RL["B`@("`)?0H@("`@"6-L;W-E*$9)3$4I.PH@("`@"6]P96XH1DE,_12PB9W)O=7!A;FYO+G1X="(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N('1E=F5N="YT>'0@.B`D_(5QN(CL*("`@(`EW:&EL92AM>2`D;&EN93T\1DE,13XI>PH@("`@"2`@("!C:&]M<"`D;&EN93L*("`@(`D@_("`@;F5X="!I9B@D;&EN92!E<2`B(BD["B`@("`)("`@(&UY("@D240L)&=R;W5P*2`]('-P;&ET*"]<="\L_)&QI;F4I.PH@("`@"2`@("`D9W)O=7!A;FYO>R1)1'T@/2`D9W)O=7`["B`@("`)?0H@("`@"6-L;W-E*$9)_3$4I.PH@("`@"21P<FEN=&=C="`](#$["B`@("`)<')I;G0@(E)E<&]R=&EN9R!'0U0@9FEL92Y<;B(["B`@_("!]"B`@("`@(`H@("`@9F]R96%C:"!M>2`D:F9N*$!F:6QE<RE["B`@("`);7D@)&%C8V5S<VEO;B`]("1J_9FX["B`@("`);7D@)'-I9"`]("1A8V-E<W-I;VX["B`@("`))'-I9#U^<R]<+E-*7"XH7'<K*5PN=&%B+R\[_"B`@("`)"B`@("`)<')I;G0@(E)E861I;F<N+BX@)'-I9%QN(CL*("`@(`EN97AT(&EF*"$D9W)O=7!A>R1S_:61]("8F("$D9W)O=7!B>R1S:61]*3L*("`@(`EI9B@D9W)O=7!A>R1S:61]*7L*"0D))&%C8V5S<VEO;CU^_<R]<+E-*7"XH7'<K*5PN=&%B+UQ?3B\["@D)?0H)"6EF*"1G<F]U<&)[)'-I9'TI>PH)"0DD86-C97-S:6]N_/7YS+UPN4TI<+BA<=RLI7"YT86(O7%]4+SL*"0E]"@H@("`@"6YE>'0@:68H)&%C8V5S<VEO;B%^+UQ?3B0O_("8F("1A8V-E<W-I;VXA?B]<7U0D+RD["B`@("`@("`@;7D@)6EN=')O;CL*("`@("`@("!M>2`E:6YT<F]N_:6-R96%D.PH@("`@("`@(&UY("5S=6US<SL*("`@("`@("!M>2`E<W5M964["B`@("`@("`@;7D@)6-O=6YT_<W,["B`@("`@("`@;7D@)6-O=6YT964["B`@("`@("`@;7D@)&UE86X@/2`P.PH@("`@"0H@("`@"6UY("1T_86<@/2`D86-C97-S:6]N.PH@("`@"6UY("1C870["B`@("`))&-A="`](").(B!I9B@D86-C97-S:6]N/7XO_7%].)"\I.PH@("`@"21C870@/2`B5"(@:68H)&%C8V5S<VEO;CU^+UQ?5"0O*3L*("`@(`DD86-C97-S:6]N_/7YS+UQ?3B0O+SL*("`@(`DD86-C97-S:6]N/7YS+UQ?5"0O+SL*("`@(`DD86-C97-S:6]N/7YS+R@N*BE<_7RA<=RLI7%].)"\D,B\["B`@("`))&%C8V5S<VEO;CU^<R\H+BHI7%\H7'<K*5Q?5"0O)#(O.PH@("`@"21A_8V-E<W-I;VX]?G,O7%].7%].+UQ?3B\["B`@("`))&%C8V5S<VEO;CU^<R]<7U1<7U0O7%]4+SL*("`@(`DC_)&%C8V5S<VEO;GLD86-C97-S:6]N?2LK.PH@("`@"7!U<V@H0&%C8V5S<VEO;BPD86-C97-S:6]N*3L*"B`@_("`)<')I;G0@(E)E861I;F<N+BX@)&IF;EQN(CL*("`@(`EP<FEN="`B86-C97-S:6]N(#T@)&%C8V5S<VEO_;EQN(CL*("`@(`EP<FEN="`B*"1C870I("1A8V-E<W-I;VY<;B(["B`@("`))'-A;7!L97-[(B1C871<="1A_8V-E<W-I;VXB?2`]("1T86<["@H@("`@"6EF*"1J9FX]?B]'5$58*"XJ*5PN4TI<+F]U=%PN=&%B+RE["B`@_("`)"21J9FX]?G,O7%].7"Y32EPN;W5T7"YT86(O7"Y32EPN;W5T7"YT86(O.PH@("`@"7T*("`@(`D*("`@_(`EM>2`D:7)F;B`]("1J9FX["B`@("`);7D@)&ER8VAE8VL@/2`Q.PH@("`@"7!R:6YT(")#:&5C:VEN9R!)_4B!R96%D<UQN(B!I9B@D:7)C:&5C:R`]/2`Q*3L*("`@(`DD:7)F;CU^<R]<+E-*7"XH+BHI7"YT86(O7"Y)_4EPN;W5T7"YT86(O.PH@("`@"7!R:6YT(")C:&5C:VEN9RXN+B`D:7)F;EQN(CL*("`@(`EI9BAO<&5N*$9)_3$4L("(D:7)F;B(I*7L*"0D)=VAI;&4H;7D@)&QI;F4]/$9)3$4^*7L*("`@("`@("`@("`@"6-H;VUP("1L_:6YE.PH@("`@("`@("`@("`);7D@*"1C:'(L)'-T87)T+"1E;F0L)&YU;2D@/2!S<&QI="@O7'0O+"1L:6YE_*3L*("`@("`@("`@("`@"21C:'(]?G,O8VAR+R\["B`@("`)"0DD:6YT<F]N:6-R96%D>R(D8VAR7'0D<W1A_<G1<="1E;F0B?2`]("1N=6T["B`@("`)"7T*("`@(`D)8VQO<V4H1DE,12D["B`@("`)?65L<V5["B`@("`)_"21I<F-H96-K(#T@,#L*("`@(`E]"B`@("`)"@D)<')I;G0@(D-H96-K:6YG(%-*(')E861S+BXN7&XB.PH@_("`@("`@(&]P96XH1DE,12P@(B1J9FXB*2!\?"!D:64@(D%B;W)T:6YG+BX@0V%N)W0@;W!E;B`D:F9N7&XB_.PH@("`@("`@('=H:6QE*&UY("1L:6YE/3Q&24Q%/BE["B`@("`@("`@("`@("`@("!C:&]M<"`D;&EN93L*_("`@("`@("`@("`@("`@("1L:6YE/7YS+UQS+UQT+V<["B`@("`@("`@("`@("`@("!M>2!`87)R87D@/2!S_<&QI="@O7'0O+"1L:6YE*3L*("`@("`@("`@("`@("`@(&UY("@D8VAR+"1S=&%R="PD96YD+"1N=6TI(#T@_*"1A<G)A>5LP72PD87)R87E;,5TL)&%R<F%Y6S)=+"1A<G)A>5LV72D["B`@("`@("`@("`@("`@("`C:68H_)&IF;CU^+U-*+FEN8VQ/=F5R;&%P<RYT86(O*7L*("`@("`@("`@("`@("`@(&EF*'-C86QA<B!`87)R87D@_/3T@-RE["B`@("`@("`@("`@("`@("`)(W!R:6YT("(D:F9N(&ES(&$@8W5S=&]M:7IE9"!32B!F:6QE7&XB_.PH@("`@("`@("`@("`@("`@"21N=6T@/2`D87)R87E;,UT@:68H)&QO;F=R96%D(#T](#$I.PH@("`@("`@_("`@("`@("`@"21N=6T@/2`D87)R87E;,UT@:68H)&QO;F=R96%D(#T](#(I.PH@("`@("`@("`@("`@("`@_?65L<V5["B`@("`@("`@("`@("`@("`))&YU;2`]("1A<G)A>5LW72!I9B@D;&]N9W)E860@/3T@,BD["B`@_("`@("`@("`@("`@("!]"B`@("`@("`@("`@("`@("`D8VAR/7YS+V-H<B\O.PH@("`@("`@("`@("`@("`@_)&EN=')O;GLB)&-H<EQT)'-T87)T7'0D96YD(GT@/2`D;G5M.PH@("`@("`@("`@("`@("`@)'-U;7-S>R(D_8VAR7'0D<W1A<G0B?7LD96YD?2`]("1N=6T["B`@("`@("`@("`@("`@("`D<W5M965[(B1C:')<="1E;F0B_?7LD<W1A<G1](#T@)&YU;3L*("`@("`@("`@("`@("`@("1I;G1R;VYA;&Q[)&-H<GU[(B1S=&%R=%QT)&5N_9")]*RL@:68H)&YU;2`^/2`H)&EN=')O86QL8W)I=&5R:6$I*3L*("`@("`@("`@("`@("`@("1I<V%L;'LD_8VAR?7LD<W1A<G1]>R1E;F1]*RL@:68H)&YU;2`^/2`H)&EN=')O86QL8W)I=&5R:6$I*3L*("`@("`@("`@_("`@("`@("1I96%L;'LD8VAR?7LD96YD?7LD<W1A<G1]*RL@:68H)&YU;2`^/2`H)&EN=')O86QL8W)I=&5R_:6$I*3L*("`@("`@("!]"B`@("`@("`@8VQO<V4H1DE,12D["B`@("`@("`@"@D);7D@)65X;VYS.PH)"6UY_("1C;W5N="`](#`["B`@("`@("`@;W!E;BA&24Q%+"`B)&1B(BD@?'P@9&EE(")!8F]R=&EN9RXN($-A;B=T_(&]P96X@)&1B7&XB.PH@("`@("`@('=H:6QE*&UY("1L:6YE/3Q&24Q%/BE["B`@("`@("`@("`@("`@("!C_:&]M<"`D;&EN93L*("`@("`@("`@("`@("`@(&YE>'0@:68H)&QI;F4@97$@(B(I.PH@("`@("`@("`@("`@_("`@;7D@*"1C:'(L)&DQ<RPD:3%E+"1I,G,L)&DR92PD=&5S+"1T964L)&%N;F\L)&%S+"1A92PD;F%M92PD_9VXI(#T@<W!L:70H+UQT+RPD;&EN92D["B`@("`@("`@("`@("`@("!M>2`H)&%I9"PD96ED*2`]("@D,2PD_,BD@:68H)&YA;64]?B\H+BHI7%\H7&0K*20O*3L*("`@("`@("`@("`@("`@(&YE>'0@:68H)&YA;64]?B]<_7U)<7R\@)B8@)&ER8VAE8VL@/3T@,"D["B`@("`@("`@("`@("`@("!M>2`D='EP92`](")X(CL*("`@("`@_("`@("`@("`@("1T>7!E(#T@(E,B(&EF*"1N86UE/7XO7%]37%\O*3L*("`@("`@("`@("`@("`@("1T>7!E_(#T@(E<B(&EF*"1N86UE/7XO7%]77%\O*3L*("`@("`@("`@("`@("`@("1T>7!E(#T@(E(B(&EF*"1N86UE_/7XO7%]27%\O*3L*"B`@("`@("`@("`@("`@("`D8VAR/7YS+V-H<B\O.PH@("`@("`@("`@("`@("`@"B`@_("`@("`@("`@("`@("`D9VX@/2`B(B!I9B@D9VX@97$@(BTB*3L*("`@("`@("`@("`@("`@("1N86UE(#T@_(B1G;CHD;F%M92(["B`@("`@("`@("`@("`@("!M>2`D=&%R9V5T(#T@(BTB.PH@("`@("`@("`@("`@("`@_;7D@)'1A<F=E=&5X;VX@/2`B+2(["B`@("`@("`@("`@("`@("!M>2`D9"`]("TQ.PH@("`@("`@("`@("`@_("`@:68H)'1Y<&4@97$@(E<B*7L*("`@("`@("`@("`@("`@(`DD=&%R9V5T(#T@)&-H<B`N(")<.B(@+B`H_)&DQ92LQ*R1D*2`N("(M(B`N("@D:3)S+3$I("X@(ELD:3%S+"1I,64L)&DR<RPD:3)E72(["B`@("`@("`@_("`@("`@("`))'1A<F=E=&5X;VX@/2`D=&%R9V5T.PH@("`@("`@("`@("`@("`@"6EF*"1I,7,@/3T@)&%S_("8F("1I,F4@/3T@)&%E*7L*("`@("`@("`@("`@("`@(`E]96QS97L*("`@("`@("`@("`@("`@(`D);F5X_="!I9B@A)&EN=')O;F%L;'LD8VAR?7LB)&DQ<UQT)&DR92)]*3L*("`@("`@("`@("`@("`@(`E]"B`@("`@_("`@("`@("`@("!]96QS97L*("`@("`@("`@("`@("`@(`EM>2`H)'1A<F=E='-T87)T+"`D=&%R9V5T96YD_*2`]("@P+#`I.PH@("`@("`@("`@("`@("`@"2@D=&%R9V5T<W1A<G0L("1T87)G971E;F0I(#T@*"@D:3%S_*R1D*2PH)&DR<RTQ*2D@:68H)&DQ<R`]/2`D:3%E*3L*("`@("`@("`@("`@("`@(`DH)'1A<F=E='-T87)T_+"`D=&%R9V5T96YD*2`]("@H)&DQ92LQ*R1D*2PH)&DR92DI(&EF*"1I,G,@/3T@)&DR92D["B`@("`@("`@_("`@("`@("`)*"1T87)G971S=&%R="P@)'1A<F=E=&5N9"D@/2`H)&%S+"1A92D["B`@("`@("`@("`@("`@_("`))'1A<F=E=&5X;VX@/2`D8VAR("X@(EPZ(B`N("@D=&%R9V5T<W1A<G0I("X@(BTB("X@*"1T87)G971E_;F0I("X@(ELD:3%S+"1I,64L)&DR<RPD:3)E72(["B`@("`@("`@("`@("`@("`))'1A<F=E="`]("1C:'(@_+B`B7#HB("X@*"1A<RD@+B`B+2(@+B`H)&%E*3L*("`@("`@("`@("`@("`@('T*("`@("`@("`@("`@("`@_(&EF*"1T>7!E(&5Q(")2(BE["B`@("`@("`@("`@("`@("`))'1A<F=E=&5X;VX@/2`D8VAR("X@(EPZ(B`N_("@D:3)S+3$I("X@(BTB("X@*"1I,F4I(&EF*"1I,7,@/3T@)&DQ92D["B`@("`@("`@("`@("`@("`))'1A_<F=E=&5X;VX@/2`D8VAR("X@(EPZ(B`N("@D:3%S+3$I("X@(BTB("X@*"1I,64I(&EF*"1I,G,@/3T@)&DR_92D["B`@("`@("`@("`@("`@("`))&%S(#T@)&%S*S$["B`@("`@("`@("`@("`@("`))&%E(#T@)&%E+3$[_"B`@("`@("`@("`@("`@("!]"@H@("`@("`@("`@("`@("`@;7D@*"1I;C$L)&EN,BPD97@Q*2`]("@P+#`L_,"D["@D)"0EI9B@A)&EN=')O;GLB)&-H<EQT)&DR<UQT)&DR92)]*7L*("`@("`@("`@("`@("`@(`DD:6XR_(#T@,#L*("`@("`@("`@("`@("`@('UE;'-E>PH@("`@("`@("`@("`@("`@"21I;C(@/2`D:6YT<F]N>R(D_8VAR7'0D:3)S7'0D:3)E(GT["B`@("`@("`@("`@("`@("!]"B`@("`@("`@("`@("`@("!I9B@A)&EN=')O_;GLB)&-H<EQT)&DQ<UQT)&DQ92)]*7L*("`@("`@("`@("`@("`@(`DD:6XQ(#T@,#L*("`@("`@("`@("`@_("`@('UE;'-E>PH@("`@("`@("`@("`@("`@"21I;C$@/2`D:6YT<F]N>R(D8VAR7'0D:3%S7'0D:3%E(GT[_"B`@("`@("`@("`@("`@("!]"B`@("`@("`@("`@("`@("!I9B@D:3%S(#T]("1I,64I>PH@("`@("`@("`@_("`@("`@"21I;C$@/2`D:6XR.PH@("`@("`@("`@("`@("`@?0H)"0D):68H)&DR<R`]/2`D:3)E*7L*("`@_("`@("`@("`@("`@(`DD:6XR(#T@)&EN,3L*("`@("`@("`@("`@("`@('T*"0D)"0H)"0D);7D@)'-U;7-S_(#T@,#L*("`@("`@("`@("`@("`@(&UY("1S=6UE92`](#`["B`@("`@("`@("`@("`@("!I9B@D='EP92!E_<2`B5R(I>PH@("`@("`@("`@("`@("`@"6EF*"$D<W5M<W-[(B1C:')<="1A<R)]('Q\("$D<W5M965[(B1C_:')<="1A92)]*7L*("`@("`@("`@("`@("`@(`D))&5X,2`](#`["B`@("`@("`@("`@("`@("`)?65L<V5[_"B`@("`@("`@("`@("`@("`)"69O<F5A8V@@;7D@)&5E*&ME>7,@)7L@)'-U;7-S>R(D8VAR7'0D87,B?2!]_*7L*("`@("`@("`@("`@("`@(`D)"6YE>'0@:68H)&5E(#X@)&%E*3L*("`@("`@("`@("`@("`@(`D)"21S_=6US<RL]("1S=6US<WLB)&-H<EQT)&%S(GU[)&5E?3L*("`@("`@("`@("`@("`@(`D)?0H)"0D)"0EF;W)E_86-H(&UY("1S<RAK97ES("5[("1S=6UE97LB)&-H<EQT)&%E(GT@?2E["B`@("`@("`@("`@("`@("`)"0EN_97AT(&EF*"1S<R`\("1A<RD["B`@("`@("`@("`@("`@("`)"0DD<W5M964K/2`D<W5M965[(B1C:')<="1A_92)]>R1S<WT["B`@("`@("`@("`@("`@("`)"7T*("`@("`@("`@("`@("`@(`D))&5X,2`]("@D<W5M<W,K_)'-U;65E*2\R.PH@("`@("`@("`@("`@("`@"7T*("`@("`@("`@("`@("`@('UE;'-E>PH@("`@("`@("`@_("`@("`@"6EF*"1T>7!E(&5Q(")2(BE["B`@("`@("`@("`@("`@("`)"21A<R`]("1I,G,["B`@("`@("`@_("`@("`@("`)"21A92`]("1I,F4["B`@("`@("`@("`@("`@("`)?0H@("`@("`@("`@("`@("`@"69O<F5A_8V@@;7D@)&5E*&ME>7,@)7L@)'-U;7-S>R(D8VAR7'0D87,B?2!]*7L*("`@("`@("`@("`@("`@(`D);F5X_="!I9B@D964@/B`D864@)B8@)'1Y<&4@97$@(E(B*3L*("`@("`@("`@("`@("`@(`D))'-U;7-S*ST@)'-U_;7-S>R(D8VAR7'0D87,B?7LD965].PH@("`@("`@("`@("`@("`@"7T*"0D)"0EF;W)E86-H(&UY("1S<RAK_97ES("5[("1S=6UE97LB)&-H<EQT)&%E(GT@?2E["B`@("`@("`@("`@("`@("`)"6YE>'0@:68H)'-S(#P@_)&%S("8F("1T>7!E(&5Q(")2(BD["B`@("`@("`@("`@("`@("`)"21S=6UE92L]("1S=6UE97LB)&-H<EQT_)&%E(GU[)'-S?3L*("`@("`@("`@("`@("`@(`E]"B`@("`@("`@("`@("`@("`):68H)'1Y<&4@97$@(E,B_*7L*("`@("`@("`@("`@("`@(`D):68H)&DQ<R`]/2`D:3%E*7L*("`@("`@("`@("`@("`@(`D)"21E>#$@_/2`D<W5M<W,["B`@("`@("`@("`@("`@("`)"0EI9B@D<W5M964@/B`P*7L*("`@("`@("`@("`@("`@(`D)_"0DD=')U94%34WLD;F%M97TK*SL*("`@("`@("`@("`@("`@(`D)"7T*("`@("`@("`@("`@("`@(`D)?0H@_("`@("`@("`@("`@("`@"0EI9B@D:3)S(#T]("1I,F4I>PH@("`@("`@("`@("`@("`@"0D))&5X,2`]("1S_=6UE93L*("`@("`@("`@("`@("`@(`D)"6EF*"1S=6US<R`^(#`I>PH@("`@("`@("`@("`@("`@"0D)"21T_<G5E05-3>R1N86UE?2LK.PH@("`@("`@("`@("`@("`@"0D)?0H@("`@("`@("`@("`@("`@"0E]"B`@("`@_("`@("`@("`@("`)?0H@("`@("`@("`@("`@("`@"6EF*"1T>7!E(&5Q(")2(BE["B`@("`@("`@("`@(`D)_"21E>#$@/2`D<W5M<W,@:68H)'-U;7-S(#X]("1S=6UE92D["B`@("`@("`@("`@(`D)"21E>#$@/2`D<W5M_964@:68H)'-U;7-S(#P@)'-U;65E*3L*("`@("`@("`@("`@("`@(`D);F5X="!I9B@D<W5M<W,@/B`H)'-U_;65E*C(I('Q\("1S=6UE92`^("@D<W5M<W,J,BDI.PH@("`@("`@("`@("`@("`@"7T*("`@("`@("`@("`@_("`@('T*("`@("`@("`@("`@("`@(`H@("`@("`@("`@("`@("`@:68H)'1Y<&4@97$@(E(B("8F("1I<F-H_96-K(#T](#$I>PH@("`@("`@("`@("`@("`@"6EF*"$D:6YT<F]N:6-R96%D>R(D8VAR7'0D:3)S7'0D:3)E_(GTI>PH@("`@("`@("`@("`@("`@"0DD:6XQ(#T@,#L*("`@("`@("`@("`@("`@(`D))&EN,B`]("1I;C$[_"B`@("`@("`@("`@("`@("`)?65L<V5["B`@("`@("`@("`@("`@("`)"21I;C$@/2`D:6YT<F]N:6-R96%D_>R(D8VAR7'0D:3)S7'0D:3)E(GT["B`@("`@("`@("`@("`@("`)"21I;C(@/2`D:6XQ.PH@("`@("`@("`@_("`@("`@"7T*("`@("`@("`@("`@("`@('T*"B`@("`@("`@("`@("`@("`D=&%R9V5T97AO;B`]("(D8VAR_7#HD=&5S7"TD=&5E(CL*("`@("`@("`@("`@("`@(`H@("`@("`@("`@(`D):68H)'1Y<&4@97$@(E(B*7L*_("`@("`@("`@("`@"0DD97@Q(#T@)&EN,2LD97@Q.PH@("`@("`@("`@("`@("`@"6YE>'0@:68H)'-U;7-S_(#T](#`@?'P@)'-U;65E(#T](#`I.PH@("`@("`@("`@("`@("`@?0H@("`@("`@("`@("`@("`@"@D)"0EN_97AT(&EF*"1E>#$@/"`D8W)I=&5R:6$I.PH*("`@("`@("`@("`@("`@("1E>'-[)&-H<GU[(B1A<UQT)&%E_(GT@*ST@)&5X,3L*("`@("`@("`@("`@("`@(&UY("104TD@/2`H*"1I;C$K)&EN,BDO,BDO)&5X,3L*("`@_("`@("`@("`@"B`@("`@("`@("`@("`@("`D4%-)(#T@,2!I9B@D4%-)(#X@,2D["B`@("`@("`@("`@("`@_("`D4%-)(#T@+3$@:68@*"104TD@/3T@,"D["B`@("`@("`@("`@("`@("`D;W5T<'5T>R(D86-C97-S:6]N_(B`N(")?(B`N("(D8V%T7'0D;F%M92)](#T@)%!323L*("`@("`@("`@("`@("`@("1N86UE<WLB)&YA;65<_="1T87)G971E>&]N7'0D:3%S+"1I,64L)&DR<RPD:3)E7'0D86YN;R)]*RL["@H@("`@("`@('T*("`@('T*_"@EM>2`H)&YR+"1N8RD@/2`H<V-A;&%R(&ME>7,@)7-A;7!L97,L('-C86QA<B!K97ES("5N86UE<RD["@EP_<FEN="`B3G5M8F5R(&]F(&5V96YT<R`]("1N8UQN(CL*"7!R:6YT(").=6UB97(@;V8@<V%M<&QE<R`]("1N_<EQN(CL*"6EF*"1P86ER960@/3T@,"E["@D)<')I;G0@(E-T871I<W1I8W,@;W!T:6]N(#T@5'=O('-A;7!L_92!T+71E<W1<;B(["@E]"@EI9B@D<&%I<F5D(#T](#(I>PH)"7!R:6YT(")3=&%T:7-T:6-S(&]P=&EO;B`]_($YO="!E;F]U9V@@<V%M<&QE<R!F;W(@="UT97-T7&XB.PH)?0H)"@EM>2`E<'9A;'5E<SL*"6UY("5F:6YA_;#L*"@EI9B@D<')I;G1G8W0@/3T@,2E["@D);W!E;BA'0U0L("(^)&]U='!U=&%S<V-C97-S:6]N+F=C="(I_('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N("1O=71P=71A<W-C8V5S<VEO;BYG8W1<;B(["@D);7D@_)&YT979E;G0@/2!S8V%L87(@:V5Y<R`E=&5V96YT.PH)"7!R:6YT($=#5"`B(S$N,EQN(CL*"0EP<FEN="!'_0U0@(B1N=&5V96YT7'0D;G)<;B(["@D)<')I;G0@1T-4(")%=F5N="!)1%QT06YN;W1A=&EO;B(["@D)9F]R_96%C:"!M>2`D<V%M<&QE*$!G<F]U<',I>PH)"0EP<FEN="!'0U0@(EQT(B`N("(H(B`N("1G<F]U<&%N;F][_)'-A;7!L97T@+B`B*2(@+B`D<V%M<&QE(&EF*'-C86QA<B!K97ES("5G<F]U<&%N;F\@/B`Q*3L*"0D)<')I_;G0@1T-4(")<="(@+B`D<V%M<&QE(&EF*'-C86QA<B!K97ES("5G<F]U<&%N;F\@/3T@,"D["@D)?0H)"7!R_:6YT($=#5"`B7&XB.PH)?0H)9F]R96%C:"!M>2`D979E;G0H:V5Y<R`E;F%M97,I>PH)"6YE>'0@:68H)&5V_96YT(&5Q("(B*3L*"0D*"0EM>2`H)&YA;64L)'1A<F=E=&5X;VXL)'=I;F=S+"1A;FYO*2`]('-P;&ET*"]<_="\L)&5V96YT*3L*"0DD;F%M93U^<R\H+BHI7'-C:'(O8VAR+SL*"0EI9B@A)'-K:7![)&YA;65]*7L*"0E]_96QS97L*"0D);F5X="!I9B@D<VMI<'LD;F%M97T@/3T@)&YR*3L*"0E]"@D);7D@0&X["@D);7D@0'0["@D)_;7D@*"1N+"1T*3L*"0EM>2`E;CL*"0EM>2`E=#L*"0EM>2`E9V-T;CL*"0EM>2`E9V-T=#L*"0EF;W)E86-H_(&UY("1S86UP;&5A8V,H:V5Y<R`E<V%M<&QE<RE["@D)"6UY("@D8V%T+"1S86UP;&4I(#T@<W!L:70H+UQT_+RPD<V%M<&QE86-C*3L*"0D);F5X="!I9B@D<V%M<&QE(&5Q("(B*3L*"0D):68H)&-A="!E<2`B3B(I>PH)_"0EI9B@A)&]U='!U='LB)'-A;7!L92(@+B`B7TXB("X@(EQT)&YA;64B?2E["@D)"0DD9V-T;GLD<V%M<&QE_?2`]("(M(CL*"0D)?65L<V5["@D)"0EM>2`D=F%L=64@/2`D;W5T<'5T>R1S86UP;&4@+B`B7TXB("X@(EQT_)&YA;64B?2HQ,#`["@D)"0DD9V-T;GLD<V%M<&QE?2`]("1V86QU93L*"0D)"21V86QU92`](#`@:68H)'9A_;'5E(#T]("TQ,#`I.PH)"0D)<'5S:"A`;BPD=F%L=64I.PH)"0D):68H)'-H;W=I9"`]/2`Q*7L*"0D)"0DD_;BX](BPH(B`N("1S86UP;&5S>R).7'0D<V%M<&QE(GT@+B`B*2(N)'9A;'5E.PH)"0D)?65L<V5["@D)"0D)_)&XN/2(L("(N)'9A;'5E.PH)"0D)?0H)"0D))&Y[)'9A;'5E?2LK.PH)"0E]"@D)"7T*"0D):68H)&-A="!E_<2`B5"(I>PH)"0EI9B@A)&]U='!U='LB)'-A;7!L92(@+B`B7U0B("X@(EQT)&YA;64B?2E["@D)"0DD9V-T_='LD<V%M<&QE?2`]("(M(CL*"0D)?65L<V5["@D)"0EM>2`D=F%L=64@/2`D;W5T<'5T>R1S86UP;&4@+B`B_7U0B("X@(EQT)&YA;64B?2HQ,#`["@D)"0DD9V-T='LD<V%M<&QE?2`]("1V86QU93L*"0D)"21V86QU92`]_(#`@:68H)'9A;'5E(#T]("TQ,#`I.PH)"0D)<'5S:"A`="PD=F%L=64I.PH)"0D):68H)'-H;W=I9"`]/2`Q_*7L*"0D)"0DD="X](BPH(B`N("1S86UP;&5S>R)47'0D<V%M<&QE(GT@+B`B*2(N)'9A;'5E.PH)"0D)?65L_<V5["@D)"0D))'0N/2(L("(N)'9A;'5E.PH)"0D)?0H)"0D))'1[)'9A;'5E?2LK.PH)"0E]"@D)"7T*"0E]_"@D);7D@*"1N=6U?;BPD;G5M7W0I(#T@*'-C86QA<B!K97ES("5N+"!S8V%L87(@:V5Y<R`E="D["@D)"@D)_;7D@*"1T;7!G;BPD=&UP;F%M92D@/2!S<&QI="@O7#HO+"1N86UE*3L*"0EI9B@A)'1E=F5N='LD=&UP;F%M_97TI>PH)"7UE;'-E>PH)"0EI9B@D<')I;G1G8W0@/3T@,2E["@D)"0EP<FEN="!'0U0@)'1M<&YA;64@+B`B_7'0B("X@)'1M<&=N.PH)"0D)9F]R96%C:"!M>2`D<V%M<&QE*$!G<F]U<',I>PH)"0D)"6EF*"$D9V-T;GLD_<V%M<&QE?2E["@D)"0D)?65L<V5["@D)"0D)"6EF*"1G8W1N>R1S86UP;&5](#T]("TQ,#`I>PH)"0D)"0D)_<')I;G0@1T-4(")<="(@+B`B,"(["@D)"0D)"7UE;'-E>PH)"0D)"0D)<')I;G0@1T-4(")<="(@+B`D9V-T_;GLD<V%M<&QE?3L*"0D)"0D)?0H)"0D)"0EN97AT.PH)"0D)"7T*"0D)"0EI9B@A)&=C='1[)'-A;7!L97TI_>PH)"0D)"7UE;'-E>PH)"0D)"0EI9B@D9V-T='LD<V%M<&QE?2`]/2`M,3`P*7L*"0D)"0D)"7!R:6YT($=#_5"`B7'0B("X@(C`B.PH)"0D)"0E]96QS97L*"0D)"0D)"7!R:6YT($=#5"`B7'0B("X@)&=C='1[)'-A;7!L_97T["@D)"0D)"7T*"0D)"0E]"@D)"0E]"@D)"0EP<FEN="!'0U0@(EQN(CL*"0D)?0H)"7T*"0D*"0EI9B@D_<&%I<F5D(#P](#$I>PH)"0EN97AT(&EF*"1N=6U?;B`]/2`Q('Q\("1N=6U?="`]/2`Q*3L*"0E]"@D*"0EN_97AT(&EF*"$D;B!\?"`A)'0I.PH)"0H)"6EF*"1P86ER960@/#T@,2E["@D)"6YE>'0@:68H<V-A;&%R($!N_(#P@,B!\?"!S8V%L87(@0'0@/"`R*3L*"0E]"@D);7D@)',@/2!S8V%L87(@0&X["@D)"@D);7D@*"1P=F%L_=64L)&1I9F8I(#T@*#$L+3$I.PH)"0H)"6EF*"1P86ER960@/3T@,"E["@D)"2@D<'9A;'5E+"1D:69F*2`]_('5N<&%I<F5D='1E<W0H7$!N+%Q`="D["@D)"21S(#T@<V-A;&%R($!N("X@(EQT(B`N('-C86QA<B!`=#L*_"0E]"@D):68H)'!A:7)E9"`]/2`R*7L*"0D);7D@*"1A=F=N+"1A=F=T*2`]("AA=F5R86=E*%Q`;BDL879E_<F%G92A<0'0I*3L*"0D))&1I9F8@/2`D879G="`M("1A=F=N.PH)"0DD<R`]('-C86QA<B!`;B`N(")<="(@_+B!S8V%L87(@0'0["@D)?0H*"@D))&5V96YT/7YS+UPL("]<7R]G.PH)"21N/7YS+UPL("\O.PH)"21T/7YS_+UPL("\O.PH*"0EI9B@D;F%M93U^+UQ?5UQ?+RE["@D)"21F:6YA;'LD=&%R9V5T97AO;B`N(")<="(@+B`D_=VEN9W,@+B`B7'17(GT@+CT@(GPB("X@(B1N86UE7'0D=VEN9W-<="1A;FYO7'0D<UQT)&1I9F9<="1P=F%L_=65<="1N7'0D="(["@D)"21P=F%L=65S>R(D;F%M92PD=&%R9V5T97AO;BPD=VEN9W,L)&%N;F\L5R)](#T@_)'!V86QU93L*"0E]"@D):68H)&YA;64]?B]<7U)<7R\I>PH)"0DD9FEN86Q[)'1A<F=E=&5X;VX@+B`B7'0B_("X@)'=I;F=S("XB7'12(GT@+CT@(GPB("X@(B1N86UE7'0D=VEN9W-<="1A;FYO7'0D<UQT)&1I9F9<="1P_=F%L=65<="1N7'0D="(["@D)"21P=F%L=65S>R(D;F%M92PD=&%R9V5T97AO;BPD=VEN9W,L)&%N;F\L4B)]_(#T@)'!V86QU93L*"0E]"@D):68H)&YA;64]?B]<7U-<7R\I>PH)"0DD9FEN86Q[)'1A<F=E=&5X;VX@+B`B_7'0B("X@)'=I;F=S("XB7'13(GT@+CT@(GPB("X@(B1N86UE7'0D=VEN9W-<="1A;FYO7'0D<UQT)&1I9F9<_="1P=F%L=65<="1N7'0D="(["@D)"21P=F%L=65S>R(D;F%M92PD=&%R9V5T97AO;BPD=VEN9W,L)&%N;F\L_4R)](#T@)'!V86QU93L*"0E]"@E]"@EI9B@D<')I;G1G8W0@/3T@,2E["@D)8VQO<V4H1T-4*3L*"7T*"6UY_("1T;W1A;'`@/2!S8V%L87(@:V5Y<R`E<'9A;'5E<SL*"7!R:6YT(")N=6UB97(@;V8@<"UV86QU92`]("1T_;W1A;'!<;B(["@D*"6UY($!P=F%L=64["@EF;W)E86-H(&UY("1E=F5N="AS;W)T(&ME>7,@)7!V86QU97,I_>PH)"6YE>'0@:68H)&5V96YT(&5Q("(B*3L*"0EP=7-H*$!P=F%L=64L)'!V86QU97-[)&5V96YT?2D["@E]_"@D*"6UY($!F:6YA;&]U='!U=#L*"6UY($!F:6YA;'!V86QU93L*"0H);W!E;BA/550L("(^)&]U='!U=&%S_<V-C97-S:6]N+G1X="(I('Q\(&1I92`B06)O<G1I;F<N+B!#86XG="!O<&5N("1O=71P=71A<W-C8V5S<VEO_;BYT>'1<;B(["@EP<FEN="!/550@(D5V96YT($E$7'1'96YE(%-Y;6)O;%QT5&%R9V5T($5X;VY<=$5V96YT_(%1Y<&5<=$Y<=%1<=$5X;VX@5'EP95QT4F5F97)E;F-E(%1R86YS8W)I<'1<=,Z44%-)("@E*5QT(B`N(")4_+71E<W0@<"UV86QU92(@+B`B7'0B("X@(D9$4B`H0D@I(B`N(")<="(@+B`B3B!686QU97-<=%0@5F%L=65S_7&XB.PH)9F]R96%C:"!M>2`D=&%G*'-O<G0@:V5Y<R`E9FEN86PI>PH)"6YE>'0@:68H(21F:6YA;'LD=&%G_?2D["@D);7D@*"1T87)G971E>&]N+"1T87)G971W:6YG<RPD8V%T*2`]('-P;&ET*"]<="\L)'1A9RD["@D)_;7D@0&%R<F%Y(#T@<W!L:70H+UQ\+RPD9FEN86Q[)'1A9WTI.PH)"6UY("1M:6YD:7-T(#T@,#L*"0EM>2`D_;6EN:60@/2`P.PH)"6UY("1M87AD:69F(#T@,#L*"0EM>2`D;6EN<"`](#$["@D);7D@)&UA>&ED(#T@,#L*_"0EM>2`D9VX@/2`B+2(["@D);7D@)&UX92`](#`["@D):68H)&-A="!E<2`B4B(I>PH)"0DD;6%X:60@/2`Q_.PH)"7T*"0EI9BAS8V%L87(@0&%R<F%Y(#X@,BE["@D)9F]R*&UY("1I(#T@,3LD:2`\('-C86QA<B!`87)R_87D[)&DK*RE["@D)"6UY("@D;F%M92PD=VEN9W,L)&%N;F\L)&-N+"1C="PD9&EF9BPD<'9A;'5E+"1N+"1T_*2`]('-P;&ET*"]<="\L)&%R<F%Y6R1I72D["@D)"6EF*"1C870@97$@(E,B*7L*"0D)"6YE>'0@:68H(21T_<G5E05-3>R1N86UE?2D["@D)"0EM>2`D=')U94%34V-R:71E<FEA(#T@*"1N<BHP+C(U*3L*"0D)"21T<G5E_05-38W)I=&5R:6$@/2`Q(&EF*"1T<G5E05-38W)I=&5R:6$@/"`Q*3L*"0D)"6YE>'0@:68H)'1R=65!4U-[_)&YA;65](#P]("1T<G5E05-38W)I=&5R:6$I.PH)"0E]"@D)"6YE>'0@:68H)'=I;F=S(&YE("1T87)G971W_:6YG<RD["@H)"0EM>2`H)&=X;BPD<F5S="D@/2!S<&QI="@O7#HO+"1N86UE*3L*"0D);7D@0'=I;F=L;V,@_/2!S<&QI="@O7"PO+"1W:6YG<RD["@D)"21G;B`]("1G>&X["@D)"6UY("@D8VAR+"1S=&%R="PD96YD+"1T_>7!E+"1R968L)&5I9"D@/2!S<&QI="@O7%\O+"1R97-T*3L*"0D))&YA;64]?G,O7#HO7'0O.PH)"0DD8VAR_/7YS+V-H<B\O(&EF*"1C:'(]?B]C:'(O*3L*"0D)"@D)"6YE>'0@:68H(21I;G1R;VYA;&Q[)&-H<GU[(B1S_=&%R=%QT)&5N9")]("8F("1T>7!E(&5Q(")2(BD["@D)"6UY("1I8V]U;G0@/2`P.PH*"0D);7D@*"1I<W5M_<W,L)&ES=6UE92D@/2`H,"PP*3L*"0D):68H(21I<V%L;'LD8VAR?7LD=VEN9VQO8ULP77U[)'=I;F=L;V-;_,5U]*7L*"0D)?65L<V5["@D)"0EF;W)E86-H(&UY("1I<V%L;&5N9"AK97ES("5[("1I<V%L;'LD8VAR?7LH_)'-T87)T*7T@?2E["@D)"0D):68H(21I<V%L;'LD8VAR?7LD<W1A<G1]>R1I<V%L;&5N9'TI>PH)"0D)"0EN_97AT.PH)"0D)"0EP<FEN="`B*&ES86QL*2`D8VAR.B1S=&%R="TD96YD(&AA<R!Z97)O('9A;'5E+EQN(CL*_"0D)"0D)97AI=#L*"0D)"0E]"@D)"0D);F5X="!I9B@D:7-A;&QE;F0@/B`D96YD*3L*"0D)"0DD:7-U;7-S_("L]("1I<V%L;'LD8VAR?7LD<W1A<G1]>R1I<V%L;&5N9'T["@D)"0E]"@D)"0DD:7-U;7-S+2T["@D)"7T*_"0D):68H(21I96%L;'LD8VAR?7LD=VEN9VQO8ULS77U[)'=I;F=L;V-;,EU]*7L*"0D)?65L<V5["@D)"0EF_;W)E86-H(&UY("1I<V%L;'-T87)T*&ME>7,@)7L@)&EE86QL>R1C:')]>R1E;F1]('TI>PH)"0D)"6EF*"$D_:65A;&Q[)&-H<GU[)&5N9'U[)&ES86QL<W1A<G1]*7L*"0D)"0D);F5X=#L*"0D)"0D)<')I;G0@(BAI96%L_;"D@)&-H<CHD<W1A<G0M)&5N9"!H87,@>F5R;R!V86QU92Y<;B(["@D)"0D)"65X:70["@D)"0D)?0H)"0D)_"6YE>'0@:68H)&ES86QL<W1A<G0@/"`D<W1A<G0I.PH)"0D)"21I<W5M964@*ST@)&EE86QL>R1C:')]>R1E_;F1]>R1I<V%L;'-T87)T?3L*"0D)"7T*"0D)"21I<W5M964M+3L*"0D)?0H)"0D*"0D);F5X="!I9B@D:7-U_;7-S(#P](#`@)B8@)&ES=6UE92`\/2`P*3L*"B`@("`@("`@("`@("1I8V]U;G0@/2`H)&ES=6US<RLD:7-U_;65E*2\R.PH*("`@("`@("`@("`@:68H(21I;G1R;VYA;&Q[)&-H<GU[(B1S=&%R=%QT)&5N9")]*7L*("`@_("`@("`@("`@?65L<V5["B`@("`@("`@("`@(`DD:6-O=6YT("L]("1I;G1R;VYA;&Q[)&-H<GU[(B1S=&%R_=%QT)&5N9")].PH@("`@("`@("`@("!]"B`@("`@("`@("`@(`H)"0EI9B@D:6-O=6YT(#X]("1M87AD:69F_('Q\("1M87AD:69F(#T](#`I>PH)"0D))&UI;F1I<W0@/2`H*"1E;F0@+2`D<W1A<G0I("L@,2D@:68H)&UI_;F1I<W0@/3T@,"D["@D)"0DD;6%X9&EF9B`]("1I8V]U;G0["@D)"0DD;6%X:60@/2`D:3L*"0D)"21M:6YD_:7-T(#T@*"1E;F0@+2`D<W1A<G0I("L@,3L*"0D)?0H*"0E]"@D)?0H)"0H)"6EF*'-C86QA<B!`87)R87D@_/3T@,BE["@D)"21M87AI9"`](#$["@D)?0H@("`@("`@("`*"0EN97AT(&EF*"1M87AI9"`]/2`P("8F("1C_870@;F4@(E(B*3L*"0EM>2`H)&YA;64L)'=I;F=S+"1A;FYO+"1C;BPD8W0L)&1I9F8L)'!V86QU92PD;BPD_="D@/2!S<&QI="@O7'0O+"1A<G)A>5LD;6%X:61=*3L*"0EM>2`H)&=X;BPD<F5S="D@/2!S<&QI="@O7#HO_+"1N86UE*3L*"@D);7D@*"1C:'(L)'-T87)T+"1E;F0L)'1Y<&4L)')E9BPD96ED*2`]('-P;&ET*"]<7R\L_)')E<W0I.PH)"21C:'(]?G,O8VAR+R\@:68H)&-H<CU^+V-H<B\I.PH*"0EI9B@D;F%M93U^+VYO=F5L97AO_;B\I>PH)"0EM>2`D=&UP97AO;B`]("1T87)G971E>&]N.PH)"0DD=&UP97AO;CU^<R]<+2]<.B\["@D)"6UY_("@D96-H<BPD97,L)&5E*2`]('-P;&ET*"]<.B\L)'1M<&5X;VXI.PH)"0EM>2`D9F]U;F0@/2`P.PH)"0EF_;W)E86-H(&UY("1I;&]C*&ME>7,@)7L@)&EN=')O;F%L;'LD8VAR?2!]*7L*"0D)"6UY("@D:7-T87)T+"1I_96YD*2`]('-P;&ET*"]<="\L)&EL;V,I.PH)"0D):68H)&ES=&%R="`\("1E<R`F)B`D:65N9"`\("1E92E[_"@D)"0D))&9O=6YD(#T@,3L*"0D)"0EL87-T.PH)"0D)?0H)"0D):68H)&ES=&%R="`^("1E<R`F)B`D:65N_9"`\("1E92E["@D)"0D))&9O=6YD(#T@,3L*"0D)"0EL87-T.PH)"0D)?0H)"0E]"@D)"6YE>'0@:68H)&9O_=6YD(#T](#$I.PH)"7T*"@D);7D@*"1O<F=A;BPD14Y35"D@/2`H(B(L(B(I.PH)"6EF*"1N86UE/7XO7%]%_3BA<=RLI5"A<9"LI7%\O*7L*("`@("`@("`)*"1O<F=A;BPD14Y35"D@/2`H)#$L)#(I(&EF*"1N86UE/7XO_7%]%3BA<=RLI5"A<9"LI7%\O*3L*("`@("`@("`))$5.4U0@/2`B14XB("X@)&]R9V%N("X@(E0B("X@)$5._4U0["B`@("`@("`@?0H@("`@("`@(&EF*"1N86UE/7XO7%]2,5PN+R!\?"`D;F%M93U^+UQ?4C)<+B\I>PH@_("`@("`@(`EM>2`D=&UP;F%M92`]("1N86UE.PH@("`@("`@(`DD=&UP;F%M93U^<R\H+BHI7%]2,5PN+U(Q_7"XO.PH@("`@("`@(`DD=&UP;F%M93U^<R\H+BHI7%]2,EPN+U(R7"XO.PH@("`@("`@(`DD=&UP;F%M93U^_<R]<7R@N*BDO+SL*("`@("`@("`)*"1O<F=A;BPD14Y35"D@/2`H(E-Y;G1H971I8R(L)'1M<&YA;64I.PH@_("`@("`@('T*("`@("`@("!I9B@D14Y35"!E<2`B(BE["B`@("`@("`@"6UY($!R97-T(#T@<W!L:70H+UQ?_+RPD<F5S="D["B`@("`@("`@"21%3E-4(#T@)')E<W1;-%T["B`@("`@("`@?0H@("`@("`@(`H)"21N86UE_/7YS+UPZ+UQT+SL*"0EP=7-H*$!F:6YA;'!V86QU92PD<'9A;'5E*3L*"0EP=7-H*$!F:6YA;&]U='!U="PB_)')E<W1<="1G>&Y<="1%3E-47'0D8V%T7'0D=&%R9V5T97AO;EQT)&-N7'0D8W1<="1A;FYO7'0D9&EF9EQT_)'!V86QU95QT)&Y<="1T(BD["@E]"@D*"7!R:6YT(").=6UB97(@;V8@9FEN86P@<"UV86QU92`]("(@+B!S_8V%L87(@0&9I;F%L<'9A;'5E("X@(EQN(CL*"7!R:6YT(")$;VEN9R!A9&IU<W0@<"UV86QU97,N+BY<;B([_"@D*"6UY("1F9'(["@EM>2!`9F1R.PH):68H)'!A:7)E9"`]/2`R*7L*"0E`9F1R(#T@0&9I;F%L<'9A;'5E_.PH)?65L<V5["@D):68H)&-H96-K9F1R(#T](#$I>PH)"0DD9F1R(#T@0D@H7$!F:6YA;'!V86QU92D["@D)_"4!F9'(@/2!`)&9D<CL*"0E]96QS97L*"0D)0&9D<B`]($!F:6YA;'!V86QU93L*"0E]"@E]"@H)<')I;G0@_(FYU;6)E<B!O9B!F9'(H0D@I(#T@(B`N('-C86QA<B!`9F1R("X@(EQN(CL*"6UY("1A9&IP=&5M<&-O=6YT_(#T@,#L*"69O<BAM>2`D:2`](#`[)&D@/"!S8V%L87(@0&9I;F%L;W5T<'5T.R1I*RLI>PH)"6UY("@D979E_;G0L)&=X;BPD14Y35"PD8V%T+"1T87)G971E>&]N+"1C;BPD8W0L)&%N;F\L)&1I9F8L)'!V86QU92PD;G9A_;'5E+"1T=F%L=64I(#T@<W!L:70H+UQT+RPD9FEN86QO=71P=71;)&E=*3L*"0EM>2`D9FEN86QO=71P=70@_/2`B)&5V96YT7'0D9WAN7'0D=&%R9V5T97AO;EQT)&-A=%QT)&-N7'0D8W1<="1A;FYO7'0D14Y35"(["@D)_<')I;G1F($]55"`D9FEN86QO=71P=70@+B`B7'0E+C)F7'0E+C5E7'0E+C5E7'0D;G9A;'5E7'0D='9A;'5E_7&XB+"`D9&EF9BPD<'9A;'5E+"1F9');)&E=.PH)?0H)8VQO<V4H3U54*3L*"G-U8B!U;G!A:7)E9'1T97-T_>PH*"0EM>2`H)&XL)'0I(#T@0%\["@D);7D@0&X@/2!`)&X["@D);7D@0'0@/2!`)'0["@D);7D@)&YN(#T@_<&1L*$!N*3L*"0EM>2`D='0@/2!P9&PH0'0I.PH)"6UY("@D='-T871S+"`D9&8I(#T@=%]T97-T*"`D;FXL_("1T="`I.PH)"75S92!01$PZ.D=33#HZ0T1&.R`*"0EM>2`D<%\R=&%I;"`](#(@*B!G<VQ?8V1F7W1D:7-T_7U$H("1T<W1A=',M/F%B<RP@)&1F*3L*"0EM>2`H)&%V9VXL)&%V9W0I(#T@*&%V97)A9V5X*"1N*2QA=F5R_86=E>"@D="DI.PH)"6UY("1D:69F(#T@)&%V9W0@+2`D879G;CL*"0ER971U<FX@)'!?,G1A:6PL)&1I9F8[_"GT*("`@("`*"G-U8B!S=&1E=GL*("`@("`@("!M>2@D9&%T82D@/2!`7SL*("`@("`@("!I9BA`)&1A=&$@_/3T@,2E["B`@("`@("`@("`@("`@("!R971U<FX@,#L*("`@("`@("!]"B`@("`@("`@;7D@)&%V97)A9V4@_/2`F879E<F%G92@D9&%T82D["B`@("`@("`@;7D@)'-Q=&]T86P@/2`P.PH@("`@("`@(&9O<F5A8V@H0"1D_871A*2!["B`@("`@("`@("`@("`@("`D<W%T;W1A;"`K/2`H)&%V97)A9V4M)%\I("HJ(#(["B`@("`@("`@_?0H@("`@("`@(&UY("1S=&0@/2`H)'-Q=&]T86P@+R`H0"1D871A+3$I*2`J*B`P+C4["B`@("`@("`@<F5T_=7)N("1S=&0["GT*"G-U8B!A=F5R86=E>'L*("`@("`@("!M>2@D9&%T82D@/2!`7SL*("`@("`@("!I9B`H_;F]T($`D9&%T82D@>PH@("`@("`@("`@("`@("`@9&EE*")%;7!T>2!A<G)A>5QN(BD["B`@("`@("`@?0H@_("`@("`@(&UY("1T;W1A;"`](#`["B`@("`@("`@9F]R96%C:"`H0"1D871A*2!["B`@("`@("`@("`@("`@_("`D=&]T86P@*ST@)%\["B`@("`@("`@?0H@("`@("`@(&UY("1A=F5R86=E(#T@)'1O=&%L("\@0"1D871A<.PH@("`@("`@(')E='5R;B`D879E<F%G93L*?0}
