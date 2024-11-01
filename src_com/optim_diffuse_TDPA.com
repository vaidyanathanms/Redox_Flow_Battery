! To optimize the structure
%mem=2GB
%chk=optim_gasphase_3.chk
#p B3LYP/6-31+G(d,p) opt=(maxcycle=200,maxstep=20)  scrf=(pcm,solvent=diethylether) freq=noraman

 TDPA

0 1
C	-5.278772	-4.534404	-0.845251
C	-3.919082	-5.070719	-1.322916
N	-2.824797	-4.936221	-0.36724
C	-2.511811	-6.079249	0.483904
C	-3.135349	-6.042959	1.888935
C	-2.135413	-3.73117	-0.23118
C	-1.074124	-3.578964	0.694166
C	-0.36785	-2.384138	0.804368
C	-0.666517	-1.279122	-0.004715
C	-1.718314	-1.415289	-0.920282
C	-2.446673	-2.598644	-1.021413
N	0.065200	-0.061603	0.10992
C	-0.629242	1.174944	0.176624
C	-1.803453	1.303563	0.94117
C	-2.487223	2.511738	1.010513
C	-2.036556	3.658866	0.323165
C	-0.86035	3.523108	-0.437402
C	-0.174298	2.307824	-0.510279
N	-2.735627	4.886761	0.455111
C	-4.181009	4.849339	0.178484
C	-4.55655	4.786235	-1.310823
C	-2.035599	6.095909	0.017169
C	-2.641886	7.386832	0.573747
C	1.485210	-0.093374	0.126414
C	2.219359	0.730550	0.999437
C	3.608330	0.692523	1.025166
C	4.341961	-0.172787	0.185049
C	3.600288	-0.998835	-0.679673
C	2.203490	-0.958327	-0.708258
N	5.757272	-0.216574	0.272153
C	6.460208	1.072703	0.170602
C	6.558509	1.649351	-1.250986
C	6.423793	-1.356068	-0.36116
C	7.859071	-1.571475	0.126008
H	-5.634951	-5.086798	0.029872
H	-6.026894	-4.638968	-1.640088
H	-5.217127	-3.476829	-0.572176
H	-4.015957	-6.135972	-1.55887
H	-3.637024	-4.581414	-2.26315
H	-2.856518	-6.977701	-0.039073
H	-1.422337	-6.180259	0.564902
H	-4.22845	-6.053644	1.835575
H	-2.835006	-5.144396	2.435902
H	-2.814901	-6.917678	2.467341
H	-0.795983	-4.391933	1.353121
H	0.435248	-2.309678	1.53141
H	-1.982849	-0.576209	-1.556676
H	-3.25860	-2.631589	-1.736987
H	-2.172436	0.447897	1.497778
H	-3.365743	2.575454	1.644696
H	-0.47439	4.361193	-1.006088
H	0.723978	2.240794	-1.11591
H	-4.641569	5.726312	0.643271
H	-4.603021	3.983595	0.692493
H	-4.217526	5.673657	-1.855841
H	-5.645153	4.724501	-1.422359
H	-4.114966	3.905348	-1.788389
H	-1.972818	6.172494	-1.083939
H	-1.007126	6.015296	0.381833
H	-3.619822	7.617343	0.139904
H	-1.978483	8.226608	0.341184
H	-2.752278	7.327013	1.661543
H	1.690813	1.396058	1.67469
H	4.131155	1.316765	1.742744
H	4.104372	-1.67253	-1.362798
H	1.666952	-1.60639	-1.394086
H	7.461432	0.951233	0.594752
H	5.949154	1.790532	0.814974
H	7.129867	0.995729	-1.918717
H	7.062450	2.622406	-1.228521
H	5.563317	1.791468	-1.685043
H	6.424129	-1.285056	-1.46441
H	5.838979	-2.245377	-0.107536
H	8.545094	-0.791126	-0.217643
H	8.232502	-2.524873	-0.26255
H	7.897255	-1.607576	1.219774
