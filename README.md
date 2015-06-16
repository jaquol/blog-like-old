# blog-like
A place to share how I solved questions that I could not quickly find in Google and/or took me some time to figure out.

- 01. How to interpret structural variants in BEDPE format?

## 01. How to interpret structural variants in BEDPE format?

[Lumpy](https://github.com/arq5x/lumpy-sv) is a tool for the discovery of structural variants using information from paired-end sequencing data. The breakpoints of the identified structural variants are reported in the [VCF](http://www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41) or BEDPE format. Although the BEDPE format is described [here](http://bedtools.readthedocs.org/en/latest/content/general-usage.html), at first glance, it did not seem to me very straightforward to understand the genomic coordinates changes underlying each variant and how the fields in the BEDPE file compared to those in the VCF.

The following VCF record (last row, all the above is the VCF header) corresponds to a duplication identified by Lumpy.
```

##fileformat=VCFv4.2									
##source=LUMPY									
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">									
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">									
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">									
##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">									
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">									
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">									
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">									
##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description="Confidence interval (95%) around POS for imprecise variants">									
##INFO=<ID=CIEND95,Number=2,Type=Integer,Description="Confidence interval (95%) around END for imprecise variants">									
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">									
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">									
##INFO=<ID=SECONDARY,Number=0,Type=Flag,Description="Secondary breakend in a multi-line variants">									
##INFO=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">									
##INFO=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">									
##INFO=<ID=SR,Number=.,Type=Integer,Description="Number of split reads supporting the variant across all samples">									
##INFO=<ID=EV,Number=.,Type=String,Description="Type of LUMPY evidence contributing to the variant call">									
##INFO=<ID=PRPOS,Number=.,Type=String,Description="LUMPY probability curve of the POS breakend">									
##INFO=<ID=PREND,Number=.,Type=String,Description="LUMPY probability curve of the END breakend">									
##ALT=<ID=DEL,Description="Deletion">									
##ALT=<ID=DUP,Description="Duplication">									
##ALT=<ID=INV,Description="Inversion">									
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">									
##ALT=<ID=INS,Description="Insertion of novel sequence">									
##ALT=<ID=CNV,Description="Copy number variable region">									
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">									
##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">									
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant">									
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting the variant">									
##FORMAT=<ID=BD,Number=1,Type=Integer,Description="Amount of BED evidence supporting the variant">									
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	3128_T47D
chr1	1584528	11	N	<DUP>	.	.	SVTYPE=DUP;STRANDS=-+:6;SVLEN=66541;END=1651069;CIPOS=-449,10;CIEND=-10,503;CIPOS95=-374,0;CIEND95=0,366;IMPRECISE;SU=6;PE=6;SR=0;PRPOS=2.44491e-07,3.26458e-07,4.33697e-07,5.7608e-07,7.57136e-07,9.87322e-07,1.28796e-06,1.66553e-06,2.14862e-06,2.75607e-06,3.52645e-06,4.48654e-06,5.68684e-06,7.14253e-06,8.95608e-06,1.11311e-05,1.38146e-05,1.70081e-05,2.09074e-05,2.54614e-05,3.0947e-05,3.7325e-05,4.47552e-05,5.33515e-05,6.32491e-05,7.45332e-05,8.7438e-05,0.00010173,0.000117821,0.000135629,0.000155258,0.000177125,0.000200949,0.000227194,0.000255827,0.000286183,0.000319177,0.000353501,0.000390203,0.000428746,0.000468895,0.000511378,0.000555476,0.000600756,0.000648001,0.000695831,0.000745204,0.000795631,0.000845517,0.000896997,0.00094752,0.000998877,0.00104988,0.00110053,0.00115056,0.00120038,0.00124937,0.00129849,0.00134617,0.00139206,0.00143762,0.0014814,0.00152386,0.00156521,0.00160519,0.00164404,0.00168131,0.00171717,0.00175139,0.00178475,0.00181657,0.00184663,0.00187563,0.00190327,0.00192975,0.00195479,0.00197842,0.00200183,0.00202369,0.00204549,0.00206583,0.00208582,0.00210445,0.00212224,0.00213857,0.00215452,0.00216896,0.00218263,0.00219584,0.00220845,0.00221975,0.00223061,0.00224049,0.00224956,0.00225841,0.00226617,0.00227317,0.00227931,0.00228505,0.00228996,0.00229495,0.00229948,0.00230347,0.00230755,0.00231109,0.0023146,0.00231767,0.00232087,0.0023235,0.00232593,0.00232843,0.00233037,0.00233266,0.00233471,0.0023365,0.00233822,0.00233974,0.00234121,0.00234261,0.00234391,0.00234524,0.00234656,0.00234787,0.00234901,0.00235006,0.0023513,0.0023522,0.00235324,0.00235429,0.00235519,0.00235614,0.00235695,0.00235786,0.00235864,0.00235957,0.00236033,0.00236117,0.00236198,0.00236288,0.00236365,0.00236441,0.00236514,0.00236581,0.00236671,0.0023675,0.00236851,0.0023693,0.00237023,0.00237106,0.00237196,0.00237283,0.00237371,0.00237471,0.00237556,0.00237647,0.00237733,0.00237819,0.00237897,0.00237973,0.00238056,0.00238145,0.00238222,0.00238294,0.00238394,0.0023848,0.00238545,0.00238618,0.0023871,0.00238796,0.0023889,0.00238976,0.00239053,0.00239142,0.00239222,0.0023932,0.00239412,0.00239504,0.00239599,0.00239704,0.00239815,0.00239902,0.00239984,0.00240073,0.00240184,0.00240288,0.00240385,0.00240489,0.00240597,0.00240695,0.00240787,0.00240872,0.00240953,0.00241046,0.00241129,0.00241228,0.00241313,0.00241411,0.00241502,0.00241603,0.00241712,0.00241811,0.00241934,0.00242034,0.00242141,0.00242256,0.00242369,0.00242488,0.00242598,0.00242723,0.00242845,0.00242958,0.00243082,0.00243182,0.00243282,0.00243389,0.00243485,0.00243593,0.00243699,0.00243818,0.00243921,0.00244038,0.00244162,0.00244278,0.00244386,0.002445,0.00244621,0.00244734,0.00244858,0.00244991,0.0024512,0.00245246,0.00245379,0.00245524,0.00245665,0.00245774,0.00245913,0.00246048,0.00246179,0.00246295,0.00246419,0.00246564,0.00246697,0.00246825,0.00246954,0.0024709,0.00247228,0.00247353,0.002475,0.00247639,0.00247812,0.0024795,0.00248075,0.00248215,0.00248359,0.00248517,0.00248655,0.00248793,0.00248927,0.00249072,0.00249223,0.00249351,0.00249512,0.00249648,0.00249801,0.00249976,0.00250133,0.00250293,0.00250464,0.00250618,0.00250758,0.00250932,0.00251113,0.00251256,0.00251406,0.00251567,0.00251749,0.00251919,0.00252051,0.00252197,0.00252351,0.00252518,0.00252675,0.00252834,0.00252988,0.00253136,0.00253307,0.00253469,0.00253652,0.0025383,0.0025403,0.00254212,0.00254387,0.00254575,0.00254754,0.00254942,0.0025514,0.00255316,0.00255534,0.0025573,0.00255935,0.00256123,0.00256313,0.00256486,0.00256668,0.00256838,0.00257033,0.00257235,0.00257408,0.00257612,0.00257812,0.00258004,0.00258208,0.00258394,0.00258613,0.00258805,0.00259035,0.00259238,0.00259449,0.00259664,0.00259868,0.00260091,0.0026033,0.00260538,0.00260742,0.00260951,0.00261187,0.00261382,0.00261579,0.00261808,0.00262029,0.00262253,0.00262481,0.00262695,0.00262929,0.00263173,0.00263396,0.0026362,0.00263851,0.00264069,0.00264344,0.0026458,0.00264821,0.00265024,0.00265247,0.00265488,0.00265735,0.00265965,0.00266185,0.00266402,0.00266663,0.00266919,0.00267144,0.00267383,0.00267621,0.002679,0.00268127,0.00268393,0.00268632,0.00268891,0.00269145,0.00269389,0.00269661,0.00269912,0.00270165,0.00270404,0.0027067,0.00270939,0.00271195,0.00271455,0.00271722,0.00271969,0.00272252,0.0027254,0.002728,0.00273082,0.00273369,0.00273648,0.00273952,0.00274227,0.00274511,0.00274816,0.00275103,0.00275394,0.00275687,0.00275977,0.00276241,0.00276495,0.00276767,0.00277064,0.00277351,0.00277641,0.00277909,0.00278202,0.00278452,0.00278729,0.00278998,0.00279275,0.0027958,0.00279878,0.0028016,0.00280472,0.00280771,0.0028109,0.00281415,0.00281717,0.00282038,0.00282346,0.00282626,0.00282919,0.00283207,0.00283518,0.00283843,0.00284122,0.0028439,0.00284708,0.00284994,0.00285294,0.00285587,0.00285848,0.00286126,0.0028643,0.00286701,0.00287003,0.00287284,0.00287589,0.00287886,0.00288176,0.00288469,0.00288762,0.0028909,0.00289392,0.00289664,0.00289945,0.00290244,0.00290543,0.00290824,0.00291089,0.00291354,0.00291634,0.00291902,0.00292172,0.00292413,0.0029271,0.00292984,0.002933,0.0029355,0.00293802,0.0029405,0.00294333,0.00294613,0.00117372,0.000467584,0.000186264,7.42093e-05,2.95642e-05,1.17781e-05,4.69226e-06,1.86919e-06,7.44696e-07,2.96676e-07;PREND=3.02026e-15,4.78471e-14,7.57973e-13,1.20067e-11,1.90178e-10,3.01282e-09,4.77324e-08,7.56151e-07,1.19779e-05,0.000189743,0.00300582,0.0030039,0.00300127,0.00299837,0.00299541,0.00299255,0.00299047,0.00298793,0.00298555,0.00298267,0.00298033,0.0029782,0.00297579,0.00297295,0.00297017,0.00296769,0.00296562,0.0029632,0.00296072,0.00295816,0.00295557,0.00295308,0.00295062,0.0029479,0.0029447,0.00294211,0.00293963,0.00293734,0.00293449,0.00293173,0.00292904,0.0029266,0.00292425,0.00292122,0.00291848,0.00291567,0.00291333,0.00291077,0.00290832,0.00290537,0.00290246,0.00289957,0.00289704,0.00289385,0.00289114,0.0028884,0.00288549,0.00288284,0.00287995,0.00287726,0.00287463,0.00287144,0.00286885,0.00286569,0.0028629,0.00286073,0.00285801,0.00285472,0.00285196,0.00284879,0.00284594,0.00284346,0.00284056,0.00283777,0.00283502,0.0028323,0.00282939,0.00282645,0.00282353,0.00282053,0.00281778,0.00281457,0.00281189,0.00280891,0.00280615,0.00280355,0.00280057,0.00279736,0.00279468,0.00279185,0.00278916,0.00278699,0.0027844,0.00278164,0.00277904,0.00277575,0.00277311,0.00277044,0.00276776,0.00276461,0.00276186,0.00275906,0.00275677,0.00275435,0.00275168,0.00274942,0.00274617,0.00274363,0.00274107,0.00273819,0.0027355,0.00273297,0.00273059,0.00272816,0.00272563,0.00272344,0.00272085,0.00271833,0.00271609,0.00271369,0.00271133,0.00270894,0.0027064,0.00270412,0.00270171,0.0026994,0.00269703,0.00269418,0.00269178,0.00268936,0.00268728,0.00268506,0.00268312,0.00268036,0.00267787,0.00267564,0.00267347,0.00267118,0.00266877,0.00266653,0.00266441,0.00266271,0.00266054,0.00265825,0.0026559,0.00265383,0.00265161,0.00264968,0.00264732,0.00264539,0.00264295,0.00264098,0.00263869,0.00263636,0.00263434,0.00263215,0.00263023,0.00262822,0.00262628,0.00262422,0.00262236,0.00262036,0.00261831,0.00261628,0.00261451,0.00261257,0.00261059,0.00260871,0.00260687,0.00260501,0.00260267,0.00260094,0.00259868,0.00259649,0.00259478,0.00259308,0.00259143,0.00258955,0.00258721,0.00258551,0.00258347,0.00258174,0.00258005,0.00257832,0.00257656,0.00257497,0.00257335,0.00257165,0.00256984,0.00256803,0.00256624,0.00256436,0.00256254,0.0025607,0.00255925,0.00255792,0.00255616,0.00255439,0.0025528,0.00255149,0.00254971,0.00254832,0.00254703,0.00254555,0.00254428,0.00254283,0.00254166,0.00254023,0.00253862,0.00253728,0.002536,0.00253436,0.00253252,0.00253096,0.00252943,0.00252819,0.00252656,0.00252508,0.00252379,0.00252243,0.00252084,0.00251975,0.00251831,0.00251666,0.00251521,0.00251391,0.00251267,0.00251144,0.00251024,0.00250897,0.00250758,0.00250646,0.00250525,0.00250396,0.00250248,0.00250084,0.00249962,0.00249836,0.00249712,0.00249567,0.00249435,0.00249341,0.00249217,0.00249093,0.00248988,0.00248849,0.00248716,0.00248589,0.00248488,0.00248398,0.00248264,0.00248155,0.00248047,0.00247944,0.00247824,0.00247724,0.00247609,0.00247492,0.00247375,0.00247262,0.00247147,0.0024706,0.00246888,0.0024678,0.00246675,0.00246572,0.00246487,0.00246414,0.00246316,0.0024619,0.00246083,0.00245995,0.00245894,0.002458,0.0024569,0.00245599,0.00245492,0.0024538,0.00245275,0.0024515,0.00245035,0.0024494,0.00244835,0.00244731,0.00244625,0.00244542,0.00244473,0.0024439,0.00244295,0.00244172,0.00244076,0.00243982,0.00243907,0.00243812,0.00243725,0.00243639,0.00243555,0.00243464,0.0024338,0.00243277,0.00243145,0.00243048,0.00242965,0.0024285,0.00242736,0.00242654,0.00242575,0.00242483,0.00242391,0.00242298,0.00242193,0.00242078,0.00241981,0.00241862,0.00241771,0.00241663,0.0024155,0.0024146,0.00241367,0.0024125,0.00241116,0.00241008,0.00240902,0.00240729,0.00240548,0.00240419,0.00240253,0.00240081,0.00239908,0.00239741,0.00239517,0.0023929,0.00239052,0.00238798,0.00238505,0.00238189,0.00237904,0.00237565,0.00237201,0.00236829,0.00236374,0.00235923,0.00235403,0.00234858,0.00234288,0.00233677,0.00233038,0.00232438,0.0023171,0.00230981,0.00230198,0.00229329,0.00228476,0.00227551,0.00226594,0.00225598,0.00224564,0.00223392,0.00222138,0.00220881,0.00219607,0.00218362,0.00217038,0.00215639,0.00214182,0.00212751,0.00211102,0.00209334,0.00207535,0.00205704,0.00203711,0.00201595,0.00199504,0.00197257,0.00194934,0.0019241,0.00189856,0.0018717,0.00184359,0.00181463,0.00178446,0.00175283,0.0017206,0.0016859,0.00165096,0.00161563,0.00157919,0.00154094,0.00150137,0.00146059,0.00141945,0.00137825,0.0013367,0.00129363,0.00124908,0.00120422,0.00115974,0.00111459,0.00106967,0.00102485,0.000979928,0.000935386,0.000891343,0.000847616,0.000805325,0.000763099,0.000721605,0.000680876,0.000641258,0.000602415,0.000564516,0.000529212,0.000494813,0.000461726,0.000429703,0.000398541,0.000369522,0.000342268,0.000316742,0.000292816,0.000270237,0.000249087,0.000228998,0.000210083,0.000192775,0.000176194,0.000160682,0.000146541,0.000133475,0.000121171,0.000110027,9.98983e-05,9.04025e-05,8.1898e-05,7.37485e-05,6.64247e-05,5.991e-05,5.34662e-05,4.79254e-05,4.28124e-05,3.80863e-05,3.37871e-05,2.98488e-05,2.64158e-05,2.31778e-05,2.02341e-05,1.77006e-05,1.53314e-05,1.32501e-05,1.13287e-05,9.72139e-06,8.17381e-06,6.84904e-06,5.73097e-06,4.74769e-06,3.91083e-06,3.20076e-06,2.6024e-06,2.08961e-06,1.69223e-06,1.36759e-06,1.08247e-06,8.47588e-07,6.69622e-07,5.17384e-07,4.00954e-07,3.11072e-07,2.34849e-07,1.76165e-07,1.31673e-07,9.73364e-08,7.19841e-08,5.32892e-08,3.90585e-08,2.83402e-08,2.02418e-08,1.4579e-08,1.06446e-08,7.59531e-09,5.37474e-09,3.79538e-09,2.64724e-09,1.84573e-09,1.29566e-09,9.10203e-10,6.14793e-10,4.26199e-10,2.87287e-10,1.936e-10,1.30732e-10,9.00388e-11,6.2542e-11,4.32799e-11,2.90205e-11,1.95416e-11,1.3365e-11,8.69465e-12,5.75237e-12,3.71099e-12,2.54394e-12,1.73591e-12,1.14175e-12,7.46542e-13,4.68409e-13,3.08728e-13,1.94589e-13,1.17053e-13,7.58227e-14,4.59278e-14,2.89681e-14,1.97988e-14,1.19182e-14,7.29947e-15,4.37438e-15	GT:SU:PE:SR	./.:6:6:0

```

Essentially, we have:

- chromosome and position where the variant starts (#CHROM and POS, respectively)
- base found in the reference genome sequence to which reads were mapped (REF)
- structural variant found in the data (ALT); the use of '<>' symbols (e.g. <DUP>) indicates that it is an imprecise variant (in the sense that it was defined with a give probability).
- INFO: string of metainformation about the variant. For instance:
  - SVLEN = length of the structural variant (bp)
  - END = end of the structural variant. It should correspond to POS + SVLEN; e.g. 1,584,528 + 66,541 = 1,651,069 

And these to the same duplication in BEDPE format:

```
chr1	1584078	1584538	chr1	1651058	1651572	11	.	-	+	DUP	.	SVTYPE=DUP;SVLEN=66541;END=1651069;STRANDS=-+:6;IMPRECISE;CIPOS=-449,10;CIEND=-10,503;CIPOS95=-374,0;CIEND95=0,366;SU=6;PE=6;SR=0;PRPOS=2.44491e-07,3.26458e-07,4.33697e-07,5.7608e-07,7.57136e-07,9.87322e-07,1.28796e-06,1.66553e-06,2.14862e-06,2.75607e-06,3.52645e-06,4.48654e-06,5.68684e-06,7.14253e-06,8.95608e-06,1.11311e-05,1.38146e-05,1.70081e-05,2.09074e-05,2.54614e-05,3.0947e-05,3.7325e-05,4.47552e-05,5.33515e-05,6.32491e-05,7.45332e-05,8.7438e-05,0.00010173,0.000117821,0.000135629,0.000155258,0.000177125,0.000200949,0.000227194,0.000255827,0.000286183,0.000319177,0.000353501,0.000390203,0.000428746,0.000468895,0.000511378,0.000555476,0.000600756,0.000648001,0.000695831,0.000745204,0.000795631,0.000845517,0.000896997,0.00094752,0.000998877,0.00104988,0.00110053,0.00115056,0.00120038,0.00124937,0.00129849,0.00134617,0.00139206,0.00143762,0.0014814,0.00152386,0.00156521,0.00160519,0.00164404,0.00168131,0.00171717,0.00175139,0.00178475,0.00181657,0.00184663,0.00187563,0.00190327,0.00192975,0.00195479,0.00197842,0.00200183,0.00202369,0.00204549,0.00206583,0.00208582,0.00210445,0.00212224,0.00213857,0.00215452,0.00216896,0.00218263,0.00219584,0.00220845,0.00221975,0.00223061,0.00224049,0.00224956,0.00225841,0.00226617,0.00227317,0.00227931,0.00228505,0.00228996,0.00229495,0.00229948,0.00230347,0.00230755,0.00231109,0.0023146,0.00231767,0.00232087,0.0023235,0.00232593,0.00232843,0.00233037,0.00233266,0.00233471,0.0023365,0.00233822,0.00233974,0.00234121,0.00234261,0.00234391,0.00234524,0.00234656,0.00234787,0.00234901,0.00235006,0.0023513,0.0023522,0.00235324,0.00235429,0.00235519,0.00235614,0.00235695,0.00235786,0.00235864,0.00235957,0.00236033,0.00236117,0.00236198,0.00236288,0.00236365,0.00236441,0.00236514,0.00236581,0.00236671,0.0023675,0.00236851,0.0023693,0.00237023,0.00237106,0.00237196,0.00237283,0.00237371,0.00237471,0.00237556,0.00237647,0.00237733,0.00237819,0.00237897,0.00237973,0.00238056,0.00238145,0.00238222,0.00238294,0.00238394,0.0023848,0.00238545,0.00238618,0.0023871,0.00238796,0.0023889,0.00238976,0.00239053,0.00239142,0.00239222,0.0023932,0.00239412,0.00239504,0.00239599,0.00239704,0.00239815,0.00239902,0.00239984,0.00240073,0.00240184,0.00240288,0.00240385,0.00240489,0.00240597,0.00240695,0.00240787,0.00240872,0.00240953,0.00241046,0.00241129,0.00241228,0.00241313,0.00241411,0.00241502,0.00241603,0.00241712,0.00241811,0.00241934,0.00242034,0.00242141,0.00242256,0.00242369,0.00242488,0.00242598,0.00242723,0.00242845,0.00242958,0.00243082,0.00243182,0.00243282,0.00243389,0.00243485,0.00243593,0.00243699,0.00243818,0.00243921,0.00244038,0.00244162,0.00244278,0.00244386,0.002445,0.00244621,0.00244734,0.00244858,0.00244991,0.0024512,0.00245246,0.00245379,0.00245524,0.00245665,0.00245774,0.00245913,0.00246048,0.00246179,0.00246295,0.00246419,0.00246564,0.00246697,0.00246825,0.00246954,0.0024709,0.00247228,0.00247353,0.002475,0.00247639,0.00247812,0.0024795,0.00248075,0.00248215,0.00248359,0.00248517,0.00248655,0.00248793,0.00248927,0.00249072,0.00249223,0.00249351,0.00249512,0.00249648,0.00249801,0.00249976,0.00250133,0.00250293,0.00250464,0.00250618,0.00250758,0.00250932,0.00251113,0.00251256,0.00251406,0.00251567,0.00251749,0.00251919,0.00252051,0.00252197,0.00252351,0.00252518,0.00252675,0.00252834,0.00252988,0.00253136,0.00253307,0.00253469,0.00253652,0.0025383,0.0025403,0.00254212,0.00254387,0.00254575,0.00254754,0.00254942,0.0025514,0.00255316,0.00255534,0.0025573,0.00255935,0.00256123,0.00256313,0.00256486,0.00256668,0.00256838,0.00257033,0.00257235,0.00257408,0.00257612,0.00257812,0.00258004,0.00258208,0.00258394,0.00258613,0.00258805,0.00259035,0.00259238,0.00259449,0.00259664,0.00259868,0.00260091,0.0026033,0.00260538,0.00260742,0.00260951,0.00261187,0.00261382,0.00261579,0.00261808,0.00262029,0.00262253,0.00262481,0.00262695,0.00262929,0.00263173,0.00263396,0.0026362,0.00263851,0.00264069,0.00264344,0.0026458,0.00264821,0.00265024,0.00265247,0.00265488,0.00265735,0.00265965,0.00266185,0.00266402,0.00266663,0.00266919,0.00267144,0.00267383,0.00267621,0.002679,0.00268127,0.00268393,0.00268632,0.00268891,0.00269145,0.00269389,0.00269661,0.00269912,0.00270165,0.00270404,0.0027067,0.00270939,0.00271195,0.00271455,0.00271722,0.00271969,0.00272252,0.0027254,0.002728,0.00273082,0.00273369,0.00273648,0.00273952,0.00274227,0.00274511,0.00274816,0.00275103,0.00275394,0.00275687,0.00275977,0.00276241,0.00276495,0.00276767,0.00277064,0.00277351,0.00277641,0.00277909,0.00278202,0.00278452,0.00278729,0.00278998,0.00279275,0.0027958,0.00279878,0.0028016,0.00280472,0.00280771,0.0028109,0.00281415,0.00281717,0.00282038,0.00282346,0.00282626,0.00282919,0.00283207,0.00283518,0.00283843,0.00284122,0.0028439,0.00284708,0.00284994,0.00285294,0.00285587,0.00285848,0.00286126,0.0028643,0.00286701,0.00287003,0.00287284,0.00287589,0.00287886,0.00288176,0.00288469,0.00288762,0.0028909,0.00289392,0.00289664,0.00289945,0.00290244,0.00290543,0.00290824,0.00291089,0.00291354,0.00291634,0.00291902,0.00292172,0.00292413,0.0029271,0.00292984,0.002933,0.0029355,0.00293802,0.0029405,0.00294333,0.00294613,0.00117372,0.000467584,0.000186264,7.42093e-05,2.95642e-05,1.17781e-05,4.69226e-06,1.86919e-06,7.44696e-07,2.96676e-07;PREND=3.02026e-15,4.78471e-14,7.57973e-13,1.20067e-11,1.90178e-10,3.01282e-09,4.77324e-08,7.56151e-07,1.19779e-05,0.000189743,0.00300582,0.0030039,0.00300127,0.00299837,0.00299541,0.00299255,0.00299047,0.00298793,0.00298555,0.00298267,0.00298033,0.0029782,0.00297579,0.00297295,0.00297017,0.00296769,0.00296562,0.0029632,0.00296072,0.00295816,0.00295557,0.00295308,0.00295062,0.0029479,0.0029447,0.00294211,0.00293963,0.00293734,0.00293449,0.00293173,0.00292904,0.0029266,0.00292425,0.00292122,0.00291848,0.00291567,0.00291333,0.00291077,0.00290832,0.00290537,0.00290246,0.00289957,0.00289704,0.00289385,0.00289114,0.0028884,0.00288549,0.00288284,0.00287995,0.00287726,0.00287463,0.00287144,0.00286885,0.00286569,0.0028629,0.00286073,0.00285801,0.00285472,0.00285196,0.00284879,0.00284594,0.00284346,0.00284056,0.00283777,0.00283502,0.0028323,0.00282939,0.00282645,0.00282353,0.00282053,0.00281778,0.00281457,0.00281189,0.00280891,0.00280615,0.00280355,0.00280057,0.00279736,0.00279468,0.00279185,0.00278916,0.00278699,0.0027844,0.00278164,0.00277904,0.00277575,0.00277311,0.00277044,0.00276776,0.00276461,0.00276186,0.00275906,0.00275677,0.00275435,0.00275168,0.00274942,0.00274617,0.00274363,0.00274107,0.00273819,0.0027355,0.00273297,0.00273059,0.00272816,0.00272563,0.00272344,0.00272085,0.00271833,0.00271609,0.00271369,0.00271133,0.00270894,0.0027064,0.00270412,0.00270171,0.0026994,0.00269703,0.00269418,0.00269178,0.00268936,0.00268728,0.00268506,0.00268312,0.00268036,0.00267787,0.00267564,0.00267347,0.00267118,0.00266877,0.00266653,0.00266441,0.00266271,0.00266054,0.00265825,0.0026559,0.00265383,0.00265161,0.00264968,0.00264732,0.00264539,0.00264295,0.00264098,0.00263869,0.00263636,0.00263434,0.00263215,0.00263023,0.00262822,0.00262628,0.00262422,0.00262236,0.00262036,0.00261831,0.00261628,0.00261451,0.00261257,0.00261059,0.00260871,0.00260687,0.00260501,0.00260267,0.00260094,0.00259868,0.00259649,0.00259478,0.00259308,0.00259143,0.00258955,0.00258721,0.00258551,0.00258347,0.00258174,0.00258005,0.00257832,0.00257656,0.00257497,0.00257335,0.00257165,0.00256984,0.00256803,0.00256624,0.00256436,0.00256254,0.0025607,0.00255925,0.00255792,0.00255616,0.00255439,0.0025528,0.00255149,0.00254971,0.00254832,0.00254703,0.00254555,0.00254428,0.00254283,0.00254166,0.00254023,0.00253862,0.00253728,0.002536,0.00253436,0.00253252,0.00253096,0.00252943,0.00252819,0.00252656,0.00252508,0.00252379,0.00252243,0.00252084,0.00251975,0.00251831,0.00251666,0.00251521,0.00251391,0.00251267,0.00251144,0.00251024,0.00250897,0.00250758,0.00250646,0.00250525,0.00250396,0.00250248,0.00250084,0.00249962,0.00249836,0.00249712,0.00249567,0.00249435,0.00249341,0.00249217,0.00249093,0.00248988,0.00248849,0.00248716,0.00248589,0.00248488,0.00248398,0.00248264,0.00248155,0.00248047,0.00247944,0.00247824,0.00247724,0.00247609,0.00247492,0.00247375,0.00247262,0.00247147,0.0024706,0.00246888,0.0024678,0.00246675,0.00246572,0.00246487,0.00246414,0.00246316,0.0024619,0.00246083,0.00245995,0.00245894,0.002458,0.0024569,0.00245599,0.00245492,0.0024538,0.00245275,0.0024515,0.00245035,0.0024494,0.00244835,0.00244731,0.00244625,0.00244542,0.00244473,0.0024439,0.00244295,0.00244172,0.00244076,0.00243982,0.00243907,0.00243812,0.00243725,0.00243639,0.00243555,0.00243464,0.0024338,0.00243277,0.00243145,0.00243048,0.00242965,0.0024285,0.00242736,0.00242654,0.00242575,0.00242483,0.00242391,0.00242298,0.00242193,0.00242078,0.00241981,0.00241862,0.00241771,0.00241663,0.0024155,0.0024146,0.00241367,0.0024125,0.00241116,0.00241008,0.00240902,0.00240729,0.00240548,0.00240419,0.00240253,0.00240081,0.00239908,0.00239741,0.00239517,0.0023929,0.00239052,0.00238798,0.00238505,0.00238189,0.00237904,0.00237565,0.00237201,0.00236829,0.00236374,0.00235923,0.00235403,0.00234858,0.00234288,0.00233677,0.00233038,0.00232438,0.0023171,0.00230981,0.00230198,0.00229329,0.00228476,0.00227551,0.00226594,0.00225598,0.00224564,0.00223392,0.00222138,0.00220881,0.00219607,0.00218362,0.00217038,0.00215639,0.00214182,0.00212751,0.00211102,0.00209334,0.00207535,0.00205704,0.00203711,0.00201595,0.00199504,0.00197257,0.00194934,0.0019241,0.00189856,0.0018717,0.00184359,0.00181463,0.00178446,0.00175283,0.0017206,0.0016859,0.00165096,0.00161563,0.00157919,0.00154094,0.00150137,0.00146059,0.00141945,0.00137825,0.0013367,0.00129363,0.00124908,0.00120422,0.00115974,0.00111459,0.00106967,0.00102485,0.000979928,0.000935386,0.000891343,0.000847616,0.000805325,0.000763099,0.000721605,0.000680876,0.000641258,0.000602415,0.000564516,0.000529212,0.000494813,0.000461726,0.000429703,0.000398541,0.000369522,0.000342268,0.000316742,0.000292816,0.000270237,0.000249087,0.000228998,0.000210083,0.000192775,0.000176194,0.000160682,0.000146541,0.000133475,0.000121171,0.000110027,9.98983e-05,9.04025e-05,8.1898e-05,7.37485e-05,6.64247e-05,5.991e-05,5.34662e-05,4.79254e-05,4.28124e-05,3.80863e-05,3.37871e-05,2.98488e-05,2.64158e-05,2.31778e-05,2.02341e-05,1.77006e-05,1.53314e-05,1.32501e-05,1.13287e-05,9.72139e-06,8.17381e-06,6.84904e-06,5.73097e-06,4.74769e-06,3.91083e-06,3.20076e-06,2.6024e-06,2.08961e-06,1.69223e-06,1.36759e-06,1.08247e-06,8.47588e-07,6.69622e-07,5.17384e-07,4.00954e-07,3.11072e-07,2.34849e-07,1.76165e-07,1.31673e-07,9.73364e-08,7.19841e-08,5.32892e-08,3.90585e-08,2.83402e-08,2.02418e-08,1.4579e-08,1.06446e-08,7.59531e-09,5.37474e-09,3.79538e-09,2.64724e-09,1.84573e-09,1.29566e-09,9.10203e-10,6.14793e-10,4.26199e-10,2.87287e-10,1.936e-10,1.30732e-10,9.00388e-11,6.2542e-11,4.32799e-11,2.90205e-11,1.95416e-11,1.3365e-11,8.69465e-12,5.75237e-12,3.71099e-12,2.54394e-12,1.73591e-12,1.14175e-12,7.46542e-13,4.68409e-13,3.08728e-13,1.94589e-13,1.17053e-13,7.58227e-14,4.59278e-14,2.89681e-14,1.97988e-14,1.19182e-14,7.29947e-15,4.37438e-15	GT:SU:PE:SR	./.:6:6:0

```

To understand how the structural variant is encoded in the BEDPE format we need to keep in mind that Lumpy assigns a probability distribution to the breakpoint, reflecting the relative uncertainty in the definition of the start and end breakpoints. Specifically, the metainformation encoded in the INFO fields contains the confidence interval for the given breakpoints.

In the example above, Lumpy identifies a duplication in the coordinates:

```
chr1  1584528 1651069
```

But because there is probability distribution underlying these start/end breakpoints Lumpy expresses them in confidence intervals in the BEDPE format. The CIPOS and CIEND each consists in two comma-separated integers representing the offset relative to the POS and END, respectively. Thus:

```
<chrom> POS-CIPOS[first integer]-1  POS+CIPOS[second integer] END-CIEND[first integer]-1  END+CIEND[second integer] ...
````

*Note that '-1' is added to convert from VCF (1-based) to BED (0-based) coordinates.*

Which results in the original BEDPE entry:

```
chr1	1584078	1584538	chr1	1651058	1651572	11	.	-	+	DUP	.	SVTYPE=DUP;SVLEN=66541;END=1651069;STRANDS=-+:6;IMPRECISE;CIPOS=-449,10;CIEND=-10,503;CIPOS95=-374,0;CIEND95=0,366; ...
```
