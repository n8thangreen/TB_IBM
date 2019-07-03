/*=============================================================================
PLOT MODEL RESULTS

This ...

*/

/*
* FILE * plot = fopen ("plots.tex","w");               //Write to tex file.
* fprintf(plot,"\\input /z/local/texlib/pictex.tex\n"  //Write typesetting headers.
* "\\input /z/local/texlib/color.tex\n"
* "\\input /z/local/texlib/ten.tex\n"
* "\\nopagenumbers\n"
* "$$ \\beginpicture\n"                   //Print plot with all rob for
*                                         //comparison of overall rates, not
*                                         //stratified by age/sex.
* "\\setcoordinatesystem units < 40pt, .9pt>\n"
* "\\setplotarea x from 0 to 10, y from 0 to 300\n"
* "\\plotheading {Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 300 by 50 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
* "\\green {"
* "\\plot\n"                             //Print UK born target rates
* "0   4.882405657 1   4.841432095  2  4.999251035 3  4.551971732  "
* "4   4.198331034 5   4.363563123  6  4.409015903 7	4.254926727  "
* "8	  4.346247257 9	  4.481558323  10	 4.642186697 / }%%UK-born\n"
* "\\red {"
* "\\plot\n"                             //Print non-UK born target rates
* "0	  73.41779446 1  79.78174887 2  78.52894338 3  79.32366578 "
* "4	  76.73495741 5  79.78629177 6  86.62330678 7  82.48493586 "
* "8	  78.7280715  9  79.7119641  10  83.28745895 / }%%non-UK-born\n"
* "\\blue {"
* "\\plot\n"                              //Print SSA target rates
* "0  149.0582631  1  184.4386793 2  176.1974421  3  210.5059923 "
* "4  210.9035237 5  222.3185632  6  221.260557   7  208.0021695 "
* "8  180.9518076 9  168.6456347  10  159.327162   / }%%SSA-born\n"
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
* for(r=(1+SSAV); r>=0; r--)                  //Print rates by year and rob
* { fprintf(plot, "\\plot\n");                //for model results (rates not
*   for(y=(1999-(int)t0); y<RT; y++)          //stratified by age/sex).
*   { tot=0; tot2=0;
*     for(s=0; s<2; s++)
*     for(a=0; a<4; a++)
*     { tot+=(repc[a][s][r][0][y]+repc[a][s][r][1][y]);
*       tot2+=N2[a][s][r][y]; }
*   fprintf(plot,"%d  %f  ",y+(int)t0-1999,100000*tot/tot2);
*   }
*   fprintf(plot,"/\n");
* }
* fprintf(plot,"\\endpicture $$\n");
*
* if(SSA)                                    //Print male SSA targets by age/sex
* fprintf(plot,"$$ \\beginpicture\n"
* "\\setcoordinatesystem units < 33pt, .9pt>\n"
* "\\setplotarea x from 0 to 10, y from 0 to 350\n"
* "\\plotheading {SSA (male) Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 350 by 50 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
* "\\black{"
* "\\plot\n"                                //Plot SSA, male 0-14
* "0 81.30759035   1 57.01610906    2 129.6436327   3 81.41575912   4 123.8523177 "
* "5 96.4971874    6 100.9776334    7 70.1857364    8 73.87644849 "
* "9 98.33944666   10  65.35784865        / }%%SSA-born male 0-14\n  "
* "\\blue{"
* "\\plot\n"                                //Plot SSA, male 15-44
* "0 182.384953    1 278.653995    2  224.0329861    3 289.8213705 "
* "4 309.4255671   5 300.5422008   6  296.3128659    7 302.7763847 "
* "8 263.5493666   9 249.8301211   10  253.7371514  / }%%SSA-born male 15-44\n  "
* "\\green{"
* "\\plot\n"                                //Plot SSA, male 45-64
* "0 89.19127906   1 75.72285512    2 122.5953662    3 108.9171898 "
* "4 120.298722    5 127.9306008    6 109.5978314    7 89.45759394 "
* "8 103.160141    9 93.28360838    10 98.08401075   / }%%SSA-born male 45-64\n  "
* "\\red{"
* "\\plot\n"                                //Plot SSA, male 65+
* "0 114.629405    1 100.0240003    2 75.02807027    3 97.895443  "
* "4 166.638709    5 101.6333684    6 148.3237679    7 135.4717205 "
* "8 86.4362854    9 68.80688022    10 101.7862427   / }%%SSA-born male 65+\n         "
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
* a=0;                                 //Print model results for SSA male,
* { fprintf(plot,"\\black{\\plot\n");  //age 0-14
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][SSA][0][y]+repc[a][M][SSA][1][y])/N2[a][M][SSA][y]);
*   fprintf(plot," / }%%SSA-born male 0-14\n  "); }
* a=1;                                 //Print model results for SSA male,
* { fprintf(plot,"\\blue{\\plot\n");  //age 15-44
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][SSA][0][y]+repc[a][M][SSA][1][y])/N2[a][M][SSA][y]);
*   fprintf(plot," / }%%SSA-born male 15-44\n  "); }
* a=2;                                 //Print model results for SSA male,
* { fprintf(plot,"\\green{\\plot\n");  //age 45-64
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][SSA][0][y]+repc[a][M][SSA][1][y])/N2[a][M][SSA][y]);
*   fprintf(plot," / }%%SSA-born male 45-64\n  "); }
* a=3;                                 //Print model results for SSA male,
* { fprintf(plot,"\\red{\\plot\n");  //age 65+
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][SSA][0][y]+repc[a][M][SSA][1][y])/N2[a][M][SSA][y]);
*   fprintf(plot," / }%%SSA-born male 65+\n  "); }
*
* fprintf(plot,"\\endpicture $$\n\n");
*
* fprintf(plot,"$$ \\beginpicture\n"      //Print female SSA targets by age.
* "\\setcoordinatesystem units < 40pt, .9pt>\n"
* "\\setplotarea x from 0 to 10, y from 0 to 300\n"
* "\\plotheading {SSA (female) Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 300 by 50 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
* "\\black{"
* "\\plot\n"                              //Plot SSA, female 0-14
* "0 82.7811354    1 127.5292519    2 156.4169798     3 210.1566924 "
* "4 132.287454    5 104.7914454    6 170.8107238     7 118.251692  "
* "8 121.387803    9 80.58741831    10 94.82097604     / }%%SSA-born female\n"
*
* "\\blue{"
* "\\plot\n"                               //Plot SSA, female 15-44
* "0 168.2151267   1 198.9762592    2 182.1682152     3 265.0961369 "
* "4 246.292152    5 290.4649611    6 292.0925175     7 255.9139333 "
* "8 228.765335    9 236.1763474    10 205.670382    / }%%SSA-born female\n"
*
* "\\green{"
* "\\plot\n"                               //Plot SSA, female 45-64
* "0 107.010524    1 116.7803685    2 106.5195623    3 85.98164027 "
* "4 110.309601    5 132.1772013    6 125.2363129    7 119.4989706 "
* "8 91.9206762    9 88.3215587     10 79.6751643    / }%%SSA-born female\n"
*
* "\\red{"
* "\\plot\n"                               //Plot SSA, female 65+
* "0 147.261002 1 138.7758718       2 82.6503574        3 162.6130578          "
* "4 99.7959936 5 156.9537087       6 158.675179        7 137.4043861          "
* "8 146.449757 9 153.6576178          10 165.7163552          / }%%SSA-born female\n"
*
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
* a=0;                                  //Print model results for SSA female,
* { fprintf(plot,"\\black{\\plot\n");   //age 0-14
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][SSA][0][y]+repc[a][F][SSA][1][y])/N2[a][F][SSA][y]);
*   fprintf(plot," / }%%SSA-born female\n"); }
*
* a=1;                                 //Print model results for SSA female,
* { fprintf(plot,"\\blue{\\plot\n");   //age 15-44
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][SSA][0][y]+repc[a][F][SSA][1][y])/N2[a][F][SSA][y]);
*   fprintf(plot," / }%%SSA-born female\n"); }
*
* a=2;                                 //Print model results for SSA female,
* { fprintf(plot,"\\green{\\plot\n");  //age 45-64
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][SSA][0][y]+repc[a][F][SSA][1][y])/N2[a][F][SSA][y]);
*   fprintf(plot," / }%%SSA-born female\n"); }
*
* a=3;                                 //Print model results for SSA female,
* { fprintf(plot,"\\red{\\plot\n");  //age 65+
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][SSA][0][y]+repc[a][F][SSA][1][y])/N2[a][F][SSA][y]);
*   fprintf(plot," / }%%SSA-born female\n"); }
*
* fprintf(plot,"\\endpicture $$\n\n");
*
* fprintf(plot,"$$ \\beginpicture\n"      //Print male ONUK-born targets by age.
* "\\setcoordinatesystem units < 40pt, 2pt>\n"
* "\\setplotarea x from 0 to 10, y from 0 to 150\n"
* "\\plotheading {ONUK (male) Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 150 by 50 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
*
* "\\black {"
* "\\plot\n"                              //Plot ONUK, male 0-14
* "0   20.519719    1 25.566153    2  21.948115    3 24.35223981 "
* "4   19.65760638   5      27.24035011  "
* "6   23.25884889  "
* "7   16.45524567  "
* "8   18.96954226  "
* "9   26.76212875  "
* "10   14.44366651  / }%%ONUK-born male\n"
*
* "\\blue {"
* "\\plot\n"                              //Plot ONUK, male 15-44
* "0 96.94117507   "
* "1 104.8781725   "
* "2 110.4507574   "
* "3 110.7537371   "
* "4 116.6851735   "
* "5 124.7103596   "
* "6 132.7951884   "
* "7 112.7432602   "
* "8 109.6465015   "
* "9 107.8022318   "
* "10 112.2679753    / }%%ONUK-born male\n"
*
* "\\green {"
* "\\plot\n"                              //Plot ONUK, male 45-64
* "0 69.65455563   "
* "1 71.30366807   "
* "2 69.65196017   "
* "3 61.0478407    "
* "4 74.44788532   "
* "5 68.55880151   "
* "6 77.70914784   "
* "7 79.40596595   "
* "8 79.41794348   "
* "9 85.56817894   "
* "10 87.70151445   / }%%ONUK-born male\n"
*
* "\\red {"
* "\\plot\n"                              //Plot ONUK, male 65+
* "0 83.06314326 "
* "1 88.31642929 "
* "2 80.0098383   "
* "3 94.72864732 "
* "4 81.87619637 "
* "5 94.98006197 "
* "6 104.4316362 "
* "7 95.38988861 "
* "8 97.90801625 "
* "9 95.7266136   "
* "10 111.0048541 / }%%ONUK-born male\n"
*
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
*
* a=0;                                 //Print model results for ONUK male,
* { fprintf(plot,"\\black {\\plot\n");    //age 0-14
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][NUK][0][y]+repc[a][M][NUK][1][y])/N2[a][M][NUK][y]);
*   fprintf(plot," / }%%ONUK-born male\n"); }
*
* a=1;                                 //Print model results for ONUK male,
* { fprintf(plot,"\\blue {\\plot\n");    //age 15-44
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][NUK][0][y]+repc[a][M][NUK][1][y])/N2[a][M][NUK][y]);
*   fprintf(plot," / }%%ONUK-born male\n"); }
*
* a=2;                                 //Print model results for ONUK male,
* { fprintf(plot,"\\green {\\plot\n");    //age 45-64
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][NUK][0][y]+repc[a][M][NUK][1][y])/N2[a][M][NUK][y]);
*   fprintf(plot," / }%%ONUK-born male\n"); }
*
* a=3;                                 //Print model results for ONUK male,
* { fprintf(plot,"\\red {\\plot\n");    //age 65+
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][NUK][0][y]+repc[a][M][NUK][1][y])/N2[a][M][NUK][y]);
*   fprintf(plot," / }%%ONUK-born male\n"); }
*
* fprintf(plot,"\\endpicture $$\n\n");
*
* fprintf(plot,"$$ \\beginpicture\n"      //Print female ONUK-born targets.
* "\\setcoordinatesystem units < 40pt, 2pt>\n"
* "\\setplotarea x from 0 to 10, y from 0 to 150\n"
* "\\plotheading {ONUK (female) Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 150 by 50 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
*
* "\\black {"
* "\\plot\n"                              //Plot ONUK, female 0-14
* "0 27.73995734 "
* "1 30.39899016 "
* "2 25.76278076 "
* "3 31.33697557 "
* "4 16.24994779 "
* "5 28.82177021 "
* "6 22.62440091 "
* "7 23.00380747 "
* "8 28.74556095 "
* "9 25.06660064 "
* "10 19.85357244 / }%%ONUK-born female\n"
*
* "\\blue {"
* "\\plot\n"                              //Plot ONUK, female 15-44
* "0 77.18337502 "
* "1 74.22282112 "
* "2 73.50583426 "
* "3 82.63738795 "
* "4 80.60739072 "
* "5 81.29275612 "
* "6 88.44092526 "
* "7 86.31178136 "
* "8 78.67907355 "
* "9 84.59571177 "
* "10 86.32410676 / }%%ONUK-born female\n"
*
* "\\green {"
* "\\plot\n"                              //Plot ONUK, female 45-64
* "0 58.67618753 "
* "1 66.4596999 "
* "2 66.68204267 "
* "3 64.96463655 "
* "4 57.7791004 "
* "5 56.26930911 "
* "6 61.58206859 "
* "7 65.01659931 "
* "8 62.26166914 "
* "9 66.08390209 "
* "10 77.21027931 / }%%ONUK-born female\n"
*
* "\\red {"
* "\\plot\n"                              //Plot ONUK, female 65+
* "0 63.67875963 "
* "1 73.80980762 "
* "2 66.69174421 "
* "3 73.21850893 "
* "4 63.21143319 "
* "5 71.02642027 "
* "6 79.35773464 "
* "7 73.52549031 "
* "8 73.92518981 "
* "9 67.04379303 "
* "10 67.75765216 / }%%ONUK-born female\n"
*
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
* a=0;                                 //Print model results for ONUK female,
* { fprintf(plot,"\\black {\\plot\n");   //age 0-14
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][NUK][0][y]+repc[a][F][NUK][1][y])/N2[a][F][NUK][y]);
*   fprintf(plot," / }%%ONUK-born female\n"); }
*
* a=1;                                 //Print model results for ONUK female,
* { fprintf(plot,"\\blue {\\plot\n");   //age 15-44
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][NUK][0][y]+repc[a][F][NUK][1][y])/N2[a][F][NUK][y]);
*   fprintf(plot," / }%%ONUK-born female\n"); }
*
* a=2;                                 //Print model results for ONUK female,
* { fprintf(plot,"\\green {\\plot\n");   //age 45-64
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][NUK][0][y]+repc[a][F][NUK][1][y])/N2[a][F][NUK][y]);
*   fprintf(plot," / }%%ONUK-born female\n"); }
*
* a=3;                                 //Print model results for ONUK female,
* { fprintf(plot,"\\red {\\plot\n");   //age 65+
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][NUK][0][y]+repc[a][F][NUK][1][y])/N2[a][F][NUK][y]);
*   fprintf(plot," / }%%ONUK-born female\n"); }
*
* fprintf(plot,"\\endpicture $$\n\n");
*
* fprintf(plot,"$$ \\beginpicture\n"      //Print male UK-born targets.
* "\\setcoordinatesystem units < 40pt, 9pt>\n"
* "\\setplotarea x from 0 to 10, y from 0 to 25\n"
* "\\plotheading {UK (male) Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 25 by 5 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
*
*"\\black {"
* "\\plot\n"                              //Plot UK, male 0-14
* "0 3.023932701 "
* "1 2.434645646 "
* "2 2.91520634 "
* "3 2.859970201 "
* "4 2.289011313 "
* "5 2.834470504 "
* "6 3.075712648 "
* "7 2.165160104 "
* "8 3.488689857 "
* "9 3.026708721 "
* "10 3.048902211 / }%%UK-born male\n"
*
* "\\blue {"
* "\\plot\n"                              //Plot UK, male 15-44
* "0 4.481943171 "
* "1 4.903896279 "
* "2 5.204785661 "
* "3 4.704959903 "
* "4 4.584048133 "
* "5 4.853866888 "
* "6 5.412551468 "
* "7 5.19964456 "
* "8 5.218619994 "
* "9 5.454976673 "
* "10 5.701486305 / }%%UK-born male\n"
*
* "\\green {"
* "\\plot\n"                              //Plot UK, male 45-64
* "0 5.955786744 "
* "1 6.501474101 "
* "2 5.521479266 "
* "3 5.876003713 "
* "4 5.220587734 "
* "5 5.413838509 "
* "6 5.288704578 "
* "7 4.842168807 "
* "8 4.557295579 "
* "9 4.641536439 "
* "10 5.588852857 / }%%UK-born male\n"
*
* "\\red {"
* "\\plot\n"                              //Plot UK, male 65+
* "0 13.26327417 "
* "1 12.15852607 "
* "2 12.2248443 "
* "3 10.57645431 "
* "4 9.9171165 "
* "5 9.493676244 "
* "6 9.036515876 "
* "7 8.659189327 "
* "8 7.719622881 "
* "9 7.40689583 "
* "10 7.709898619 / }%%UK-born male\n"
*
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
* a=0;                                 //Print model results for UK male,
* { fprintf(plot,"\\black {\\plot\n");   //age 0-14
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][UK][0][y]+repc[a][M][UK][1][y])/N2[a][M][UK][y]);
*   fprintf(plot," / }%%UK-born male\n"); }
*
* a=1;                                 //Print model results for UK male,
* { fprintf(plot,"\\blue {\\plot\n");   //age 15-44
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][UK][0][y]+repc[a][M][UK][1][y])/N2[a][M][UK][y]);
*   fprintf(plot," / }%%UK-born male\n"); }
*
* a=2;                                 //Print model results for UK male,
* { fprintf(plot,"\\green {\\plot\n");   //age 45-64
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][UK][0][y]+repc[a][M][UK][1][y])/N2[a][M][UK][y]);
*   fprintf(plot," / }%%UK-born male\n"); }
*
* a=3;                                 //Print model results for UK male,
* { fprintf(plot,"\\red {\\plot\n");   //age 65+
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][UK][0][y]+repc[a][M][UK][1][y])/N2[a][M][UK][y]);
*   fprintf(plot," / }%%UK-born male\n"); }
*
* fprintf(plot,"\\endpicture $$\n\n");
*
* fprintf(plot,"$$ \\beginpicture\n"      //Print female UK-born targets.
* "\\setcoordinatesystem units < 40pt, 9pt>\n"
* "\\setplotarea x from 0 to 10, y from 0 to 25\n"
* "\\plotheading {UK (female) Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 25 by 5 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
*
* "\\black {"
* "\\plot\n"                              //Plot UK, female 0-14
* "0 3.02260629 "
* "1 2.957328801 "
* "2 3.226992269 "
* "3 2.767763268 "
* "4 2.152409205 "
* "5 3.463091521 "
* "6 2.902636838 "
* "7 3.107509714 "
* "8 3.613310622 "
* "9 4.182510838 "
* "10 3.296668796 / }%%UK-born female\n"
*
* "\\blue {"
* "\\plot\n"                              //Plot UK, female 15-44
* "0 3.976018679 "
* "1 3.981311386 "
* "2 4.453794826 "
* "3 3.97786918 "
* "4 3.801788355 "
* "5 3.701147955 "
* "6 3.984906183 "
* "7 4.145043247 "
* "8 4.427543892 "
* "9 4.789110856 "
* "10 4.677290326 / }%%UK-born female\n"
*
* "\\green {"
* "\\plot\n"                              //Plot UK, female 45-64
* "0 3.152642642 "
* "1 3.197683351 "
* "2 2.850284366 "
* "3 2.458880309 "
* "4 2.453128651 "
* "5 2.308908671 "
* "6 2.054682702 "
* "7 2.491064001 "
* "8 2.159942763 "
* "9 2.356540014 "
* "10 2.513674773 / }%%UK-born female\n"
*
* "\\red {"
* "\\plot\n"                              //Plot UK, female 65+
* "0 6.843814026 "
* "1 6.184312144 "
* "2 6.987058766 "
* "3 5.946829075 "
* "4 5.237953214 "
* "5 4.873922856 "
* "6 4.629780019 "
* "7 4.107429197 "
* "8 4.07669363 "
* "9 4.059207216 "
* "10 4.616098263 / }%%UK-born female\n"
*
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
* a=0;                                   //Print model results for UK female,
* { fprintf(plot,"\\black {\\plot\n");   //age 0-14
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][UK][0][y]+repc[a][F][UK][1][y])/N2[a][F][UK][y]);
*   fprintf(plot," / }%%UK-born female\n"); }
*
* a=1;                                   //Print model results for UK female,
* { fprintf(plot,"\\blue {\\plot\n");   //age 15-44
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][UK][0][y]+repc[a][F][UK][1][y])/N2[a][F][UK][y]);
*   fprintf(plot," / }%%UK-born female\n"); }
*
* a=2;                                   //Print model results for UK female,
* { fprintf(plot,"\\green {\\plot\n");   //age 45-64
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][UK][0][y]+repc[a][F][UK][1][y])/N2[a][F][UK][y]);
*   fprintf(plot," / }%%UK-born female\n"); }
*
* a=3;                                   //Print model results for UK female,
* { fprintf(plot,"\\red {\\plot\n");   //age 65+
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][F][UK][0][y]+repc[a][F][UK][1][y])/N2[a][F][UK][y]);
*   fprintf(plot," / }%%UK-born female\n"); }
*
* fprintf(plot,"\\endpicture $$\n\n");
*
*
* if(SSA)                                    //Print SSA targets by age/sex
* fprintf(plot,"$$ \\beginpicture\n"         //with male and female combined.
* "\\setcoordinatesystem units < 33pt, .9pt>\n"
* "\\setplotarea x from 0 to 10, y from 0 to 350\n"
* "\\plotheading {SSA Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 350 by 50 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
* "\\black{"
* "\\plot\n"                                //Plot SSA 0-14
* "0  82.03774654 "
* "1  94.08170948 "
* "2  144.7273558 "
* "3  134.2011174 "
* "4  122.7118467 "
* "5  95.78528893 "
* "6  128.3448033 "
* "7  90.84580758 "
* "8  91.38075925 "
* "9  83.59879999 "
* "10 76.18733086        / }%%SSA-born 0-14\n  "
*
* "\\blue{"
* "\\plot\n"                                //Plot SSA 15-44
* "0  174.8934842 "
* "1  238.1059085 "
* "2  205.3501223 "
* "3  271.576088 "
* "4  263.5173978 "
* "5  280.9410356 "
* "6  281.7366477 "
* "7  270.234899 "
* "8  234.6065855 "
* "9  229.1450422 "
* "10 216.5974111 / }%%SSA-born 15-44\n        "
*
* "\\green{"
* "\\plot\n"                                //Plot SSA 45-64
* "0  97.78510654 "
* "1  97.55037528 "
* "2  115.9485722 "
* "3  94.50657628 "
* "4  110.0519747 "
* "5  123.8334078 "
* "6  112.6557587 "
* "7  101.4149688 "
* "8  93.33417879 "
* "9  85.6770944 "
* "10 83.47263249 / }%%SSA-born45-64\n        "
*
* "\\red{"
* "\\plot\n"                                //Plot SSA 65+
* "0  130.0698757 "
* "1  123.4689494 "
* "2  80.16079334 "
* "3  129.9506796 "
* "4  124.8235816 "
* "5  122.6193504 "
* "6  147.0251662 "
* "7  132.9420172 "
* "8  112.7581354 "
* "9  109.0314517 "
* "10 122.439279 / }%%SSA-born 65+\n        "
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
* a=0;                                 //Print model results for SSA
* { fprintf(plot,"\\black{\\plot\n");  //age 0-14
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][SSA][0][y]+repc[a][M][SSA][1][y]+repc[a][F][SSA][0][y]+repc[a][F][SSA][1][y])/(N2[a][M][SSA][y]+N2[a][F][SSA][y]));
*   fprintf(plot," / }%%SSA-born 0-14\n  "); }
* a=1;                                 //Print model results for SSA
* { fprintf(plot,"\\blue{\\plot\n");  //age 15-44
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][SSA][0][y]+repc[a][M][SSA][1][y]+repc[a][F][SSA][0][y]+repc[a][F][SSA][1][y])/(N2[a][M][SSA][y]+N2[a][F][SSA][y]));
*   fprintf(plot," / }%%SSA-born 15-44\n  "); }
* a=2;                                 //Print model results for SSA
* { fprintf(plot,"\\green{\\plot\n");  //age 45-64
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][SSA][0][y]+repc[a][M][SSA][1][y]+repc[a][F][SSA][0][y]+repc[a][F][SSA][1][y])/(N2[a][M][SSA][y]+N2[a][F][SSA][y]));
*   fprintf(plot," / }%%SSA-born 45-64\n  "); }
* a=3;                                 //Print model results for SSA
* { fprintf(plot,"\\red{\\plot\n");  //age 65+
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][SSA][0][y]+repc[a][M][SSA][1][y]+repc[a][F][SSA][0][y]+repc[a][F][SSA][1][y])/(N2[a][M][SSA][y]+N2[a][F][SSA][y]));
*   fprintf(plot," / }%%SSA-born 65+\n  "); }
*
* fprintf(plot,"\\endpicture $$\n\n");
*
* fprintf(plot,"$$ \\beginpicture\n"         //Print ONUK-born targets by age
* "\\setcoordinatesystem units < 40pt, 2pt>\n" //(male & female combined).
* "\\setplotarea x from 0 to 10, y from 0 to 150\n"
* "\\plotheading {ONUK Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 150 by 50 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
*
* "\\black {"
* "\\plot\n"                              //Plot ONUK 0-14
* "0  24.02098065 "
* "1  28.46693141 "
* "2  24.00232762 "
* "3  27.17560391 "
* "4  17.27289609 "
* "5  26.63434804 "
* "6  21.97843145 "
* "7  19.12031748 "
* "8  23.08720562 "
* "9  24.49023784 "
* "10 16.15446939 / }%%ONUK-born\n"
*
* "\\blue {"
* "\\plot\n"                              //Plot ONUK 15-44
* "0  86.52713278 "
* "1  90.86627454 "
* "2  92.4359259 "
* "3  94.39820126 "
* "4  94.25449233 "
* "5  97.43242588 "
* "6  105.429993 "
* "7  97.10657 "
* "8  90.39583233 "
* "9  90.90346552 "
* "10 94.47990594 / }%%ONUK-born\n"
*
* "\\green {"
* "\\plot\n"                              //Plot ONUK 45-64
* "0  63.56315986 "
* "1  70.08882761 "
* "2  69.06216761 "
* "3  61.84219086 "
* "4  62.59499762 "
* "5  58.93462331 "
* "6  66.19250694 "
* "7  69.79914131 "
* "8  67.16253875 "
* "9  70.97218255 "
* "10 78.13041805 / }%%ONUK-born\n"
*
* "\\red {"
* "\\plot\n"                              //Plot ONUK 65+
* "0  72.36852793 "
* "1  82.07087066 "
* "2  73.74532326 "
* "3  81.23397247 "
* "4  68.42519524 "
* "5  78.10144215 "
* "6  86.56882958 "
* "7  80.86337041 "
* "8  80.76897784 "
* "9  74.70486862 "
* "10 82.03072393 / }%%ONUK-born\n"
*
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
* a=0;                                 //Print model results for ONUK
* { fprintf(plot,"\\black {\\plot\n");    //age 0-14
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][NUK][0][y]+repc[a][M][NUK][1][y]+repc[a][F][NUK][0][y]+repc[a][F][NUK][1][y])/(N2[a][M][NUK][y]+N2[a][F][NUK][y]));
*   fprintf(plot," / }%%ONUK-born\n"); }
*
* a=1;                                 //Print model results for ONUK
* { fprintf(plot,"\\blue {\\plot\n");    //age 15-44
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][NUK][0][y]+repc[a][M][NUK][1][y]+repc[a][F][NUK][0][y]+repc[a][F][NUK][1][y])/(N2[a][M][NUK][y]+N2[a][F][NUK][y]));
*   fprintf(plot," / }%%ONUK-born\n"); }
*
* a=2;                                 //Print model results for ONUK
* { fprintf(plot,"\\green {\\plot\n");    //age 45-64
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][NUK][0][y]+repc[a][M][NUK][1][y]+repc[a][F][NUK][0][y]+repc[a][F][NUK][1][y])/(N2[a][M][NUK][y]+N2[a][F][NUK][y]));
*   fprintf(plot," / }%%ONUK-born\n"); }
*
* a=3;                                 //Print model results for ONUK
* { fprintf(plot,"\\red {\\plot\n");    //age 65+
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][NUK][0][y]+repc[a][M][NUK][1][y]+repc[a][F][NUK][0][y]+repc[a][F][NUK][1][y])/(N2[a][M][NUK][y]+N2[a][F][NUK][y]));
*   fprintf(plot," / }%%ONUK-born\n"); }
*
* fprintf(plot,"\\endpicture $$\n\n");
*
*
* fprintf(plot,"$$ \\beginpicture\n"      //Print UK-born targets.
* "\\setcoordinatesystem units < 40pt, 9pt>\n"
* "\\setplotarea x from 0 to 10, y from 0 to 25\n"
* "\\plotheading {UK Notification Rates, Obs versus Model} \n"
* "\\axis bottom label {Year}   ticks withvalues {1999} {2001} {2003} {2005} {2007} {2009} /"
* "\nat 0 2 4 6 8 10 / /\n"
* "\\axis left   label {\\stack {C,a,s,e,s,,,p,e,r,,,1,0,0,0,0,0}}  ticks numbered from 0 to 25 by 5 /\n"
* "\\axis right /\n"
* "\\axis top   /\n"
*
*"\\black {"
* "\\plot\n"                              //Plot UK 0-14
* "0  3.0232877    1 2.689615229   2 3.066922319    3 2.814877657 "
* "4  2.22224511   5 3.142075021   6 2.991190523    7 2.624804578 "
* "8  3.549264801  9 3.590365428   10 3.170301759 / }%%UK-born\n"
*
* "\\blue {"
* "\\plot\n"                              //Plot UK 15-44
* "0  4.228559458   1 4.442424061   2 4.828536249   3 4.341597683 "
* "4  4.193767608   5 4.280249544   6 4.702626549   7 4.673069059 "
* "8  4.823973987   9 5.123673878   10 5.194507641 / }%%UK-born\n"
*
* "\\green {"
* "\\plot\n"                              //Plot UK 45-64
* "0  4.551111707   1 4.838700256   2 4.174979113   3 4.154885534 "
* "4  3.82816526    5 3.850275705   6 3.658435322   7 3.656021656 "
* "8  3.346534496   9 3.485718073  10 4.035165492 / }%%UK-born\n"
*
* "\\red {"
* "\\plot\n"                              //Plot UK 65+
* "0  9.552168139   1 8.716099461   2 9.222565877   3 7.935611695 "
* "4  7.263341107   5 6.87774158    6 6.558932691   7 6.11420074  "
* "8 5.690732185    9 5.555474313   10  6.001019774 / }%%UK-born\n"
*
* "\\linethickness=.5pt\n"
* "\\setdashes\n");
*
* a=0;                                 //Print model results for UK 0-14
* { fprintf(plot,"\\black {\\plot\n");
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][UK][0][y]+repc[a][M][UK][1][y]+repc[a][F][UK][0][y]+repc[a][F][UK][1][y])/(N2[a][M][UK][y]+N2[a][F][UK][y]));
*   fprintf(plot," / }%%UK-born\n"); }
*
* a=1;                                 //Print model results for UK 15-44
* { fprintf(plot,"\\blue {\\plot\n");
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][UK][0][y]+repc[a][M][UK][1][y]+repc[a][F][UK][0][y]+repc[a][F][UK][1][y])/(N2[a][M][UK][y]+N2[a][F][UK][y]));
*   fprintf(plot," / }%%UK-born\n"); }
*
* a=2;                                 //Print model results for UK 45-64
* { fprintf(plot,"\\green {\\plot\n");
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][UK][0][y]+repc[a][M][UK][1][y]+repc[a][F][UK][0][y]+repc[a][F][UK][1][y])/(N2[a][M][UK][y]+N2[a][F][UK][y]));
*   fprintf(plot," / }%%UK-born\n"); }
*
* a=3;                                 //Print model results for UK 65+
* { fprintf(plot,"\\red {\\plot\n");
*   for(y=(1999-(int)t0); y<RT; y++)
*     fprintf(plot,"%d  %f\t",y+(int)t0-1999,100000*(repc[a][M][UK][0][y]+repc[a][M][UK][1][y]+repc[a][F][UK][0][y]+repc[a][F][UK][1][y])/(N2[a][M][UK][y]+N2[a][F][UK][y]));
*   fprintf(plot," / }%%UK-born\n"); }
*
* fprintf(plot,"\\endpicture $$\n\n");
*
* fprintf(plot,"\\end \n");                  //Close the tex file.
* fclose(plot);
*
*/

