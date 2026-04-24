# Scale continuous traits

This function standardizes continuous traits. It can be useful before
computing functional space. You will have to choose which standardized
method to use based on your data. For this function to work, there must
be no NA in your `sp_tr` data frame.

## Usage

``` r
tr.cont.scale(sp_tr, std_method = "scale_center", stop_if_NA = TRUE)
```

## Arguments

- sp_tr:

  a data frame of traits values (columns) for each species (rows). Note
  that species names **must be** specified in the row names and traits
  must be **continuous**.

- std_method:

  a character string referring to the standardization method. Possible
  values: `range` (standardize by the range), `center` (use the center
  transformation: \\x' = x - mean(x)\\), `scale` (use the scale
  transformation: \\x' = \frac{x}{sd(x)}\\), or `scale_center` (use the
  scale-center transformation: \\x' = \frac{x - mean(x)}{sd(x)}\\).
  Default is `scale_center`.

- stop_if_NA:

  a logical value to stop or not the process if the `sp_tr` data frame
  contains NA. Functional measures are sensitive to missing traits. For
  further explanations, see the Note section. Default is `TRUE`.

## Value

A data frame of standardized trait values (columns) for each species
(rows).

## Author

Camille Magneville and Sebastien Villeger

## Examples

``` r
load(system.file('extdata', 'sp_tr_cestes_df', package = 'mFD'))

mFD::tr.cont.scale(sp_tr = sp_tr, std_method = 'scale_center', 
stop_if_NA = TRUE)
#>                                    logM         OgSf         OgSh         OgPo
#> Achirus_lineatus            -0.58299801 -0.681449298 -0.263909738  2.554334103
#> Anchoa_mitchilli            -1.97956543  1.130469237  1.743827528  0.317005423
#> Archosargus_probatocephalus -0.12376213 -0.595576382 -0.613948846 -0.988102973
#> Archosargus_rhomboidalis     0.49201000 -0.818845965 -0.806198305 -0.751637341
#> Ariopsis_felis               0.28738128  0.185867157 -1.051044313 -0.419675971
#> Bagre_marinus               -0.59902883  0.829914029 -1.040162268  0.230604519
#> Bairdiella_chrysoura        -0.12847708 -0.003053259 -0.006368011  0.148751031
#> Bairdiella_ronchus           0.76641994 -0.432417841 -0.610321498 -0.428770803
#> Cathorops_melanopus         -0.12093316 -0.149037217 -0.840658113 -0.756184757
#> Cetengraulis_edentulus      -0.35385156  0.297501948  2.073916221 -0.788016669
#> Chaetodipterus_faber         1.08420739 -0.896131589 -0.588557409 -0.383296643
#> Chilomycterus_schoepfii      2.20636488 -0.878957006 -1.509903869 -0.064977522
#> Chloroscombrus_chrysurus    -1.57313696 -0.209148258  1.268644905  1.335626611
#> Citharichthys_spilopterus   -0.64334933 -0.767322215 -0.996634089  2.554334103
#> Conodon_nobilis              0.14498987 -0.217735550  0.019023426  0.385216663
#> Cynoscion_arenarius         -0.43023371  0.580882572  0.323720681  1.008212657
#> Cynoscion_nebulosus          0.95407485  1.061770904 -0.289301176  1.012760074
#> Cynoscion_nothus            -0.40288702  0.340438406  0.359994164  0.962738497
#> Dasyatis_sabina              2.64674097 -0.853195131 -1.132659649 -1.993081913
#> Diapterus_auratus           -0.55282234 -0.801671381  0.525038510  0.030518214
#> Diapterus_rhombeus          -0.54716441 -0.878957006 -0.097051718 -0.110451682
#> Dorosoma_petenense           0.06012082  0.177279865  0.934928864  0.812673769
#> Eucinostomus_argenteus      -1.68912467 -0.664274715  0.554057296  0.139656199
#> Eucinostomus_gula           -0.88664062 -0.810258673  0.347298445  0.180582943
#> Eucinostomus_melanopterus   -1.00262833 -0.586989090  0.675573463  0.048707878
#> Gobionellus_oceanicus       -0.01626133  1.722992359  0.330975378  0.212414855
#> Menticirrhus_americanus      0.34301766 -0.372306800 -0.329202007 -0.969913309
#> Menticirrhus_saxatilis       0.23740284 -0.475354299 -0.621203543 -0.856227909
#> Micropogonias_undulatus      0.10727029 -0.372306800 -0.282046480 -0.919891733
#> Opsanus_beta                 1.17944933  2.263991733 -0.898695686  0.762652193
#> Orthopristis_chrysoptera    -0.15676676  0.091406949  0.129657549 -0.601572612
#> Peprilus_paru               -0.44626453 -0.810258673  0.156862661  0.203320023
#> Polydactylus_octonemus       0.41751383  0.220216323 -0.287487502 -0.992650389
#> Prionotus_carolinus         -0.73387632  1.765928818 -0.532333510 -1.074503878
#> Prionotus_scitulus          -0.06718276  1.207754861 -0.523265140 -1.169999614
#> Selene_setapinnis           -1.03563296 -0.501116174  2.284302420 -0.292348323
#> Selene_vomer                -1.25912147 -0.561227216  3.007958400 -0.165020675
#> Sphoeroides_nephelus         0.63911636 -0.535465341 -0.545029229 -0.178662923
#> Sphoeroides_pachygaster      0.55519030 -0.501116174 -0.289301176 -0.242326747
#> Sphoeroides_testudineus      1.45951721 -0.896131589 -1.067367380  0.003233718
#> Stellifer_lanceolatus       -0.71973148 -0.020227843 -0.240331975  0.230604519
#> Symphurus_plagiusa          -0.49718597 -0.836020548 -0.564979645  2.554334103
#> Synodus_foetens              0.66834903  3.758180477  1.009289504  0.144203615
#> Trichiurus_lepturus         -0.25578066  1.250691320  1.560646441  0.307910591
#> Urobatis_jamaicensis         2.55527099 -0.758734923 -1.277753580 -1.993081913
#>                                     EySz       GrLg        GtLg         EyPo
#> Achirus_lineatus            -1.753090654 -0.8300723  0.30349122  2.179709071
#> Anchoa_mitchilli             1.127439829  2.1200914 -0.72954982 -0.905437898
#> Archosargus_probatocephalus -0.584825969 -0.5590599  1.06741279 -0.401135028
#> Archosargus_rhomboidalis    -0.477809356 -0.5822895  1.89460104 -0.438216121
#> Ariopsis_felis              -0.754268938  0.1610588  0.58299775 -0.438216121
#> Bagre_marinus                0.084027859  0.2694638  0.65476294 -0.690367556
#> Bairdiella_chrysoura         0.797471941  0.7805158 -0.41132786 -0.504962089
#> Bairdiella_ronchus           0.547766512  0.4243280 -0.31878853 -0.282475529
#> Cathorops_melanopus          0.182126420  0.6024219  0.50462261 -0.504962089
#> Cetengraulis_edentulus       1.296882799  4.0558944  3.15238044 -1.224335301
#> Chaetodipterus_faber        -0.968302163 -0.6674649  3.74255470 -0.638454026
#> Chilomycterus_schoepfii     -0.442137152 -0.8300723  1.94275926  0.355319277
#> Chloroscombrus_chrysurus    -0.388628846  0.4630441 -0.80981352 -0.237978217
#> Citharichthys_spilopterus    0.877734400 -0.7139241 -0.46137463  2.179709071
#> Conodon_nobilis              0.547766512  0.5017601 -0.22813776  0.318238184
#> Cynoscion_arenarius          0.868816349  1.0437850 -0.46326319 -0.252810654
#> Cynoscion_nebulosus          0.164290318  0.7882590 -0.45382040  0.206994903
#> Cynoscion_nothus             0.485340155  0.9895825 -0.53125127 -0.193480905
#> Dasyatis_sabina              1.609014585 -0.8300723 -1.06476881  2.179709071
#> Diapterus_auratus            0.529930410 -0.6752081 -0.43682339 -1.202086645
#> Diapterus_rhombeus           0.717209482 -0.6752081 -0.13087705 -1.172421770
#> Dorosoma_petenense           1.207702288  0.5637058  0.31954396 -0.720032431
#> Eucinostomus_argenteus       1.145275931 -0.6752081 -0.21869498 -0.957351429
#> Eucinostomus_gula            0.645865074 -0.6984377 -0.24607906 -0.979600085
#> Eucinostomus_melanopterus    0.904488553 -0.7061809 -0.43682339 -1.305913706
#> Gobionellus_oceanicus       -0.647252326 -0.8223291  0.09574987  0.993114083
#> Menticirrhus_americanus     -0.245940030 -0.4506549 -0.48498161 -0.008075438
#> Menticirrhus_saxatilis      -0.406464948 -0.4506549 -0.45382040  0.118000279
#> Micropogonias_undulatus      0.075109808 -0.2648179 -0.44343334 -0.401135028
#> Opsanus_beta                -0.433219101 -0.5590599 -0.19508801  1.193351987
#> Orthopristis_chrysoptera    -1.110990980 -0.3732228 -0.29990296  0.162497591
#> Peprilus_paru               -0.852367500 -0.4429117 -0.17336959 -1.291081269
#> Polydactylus_octonemus       1.065013472  2.1123482 -0.63701050 -0.920270335
#> Prionotus_carolinus         -0.549153765  0.5404762 -0.46609603  0.548140963
#> Prionotus_scitulus          -0.834531398  0.5095034 -0.36694675  0.518476088
#> Selene_setapinnis           -1.735254552 -0.1796425 -0.86174886 -0.912854117
#> Selene_vomer                -1.690664297 -0.2880475 -0.58129805 -0.564291839
#> Sphoeroides_nephelus        -1.530139378 -0.8300723  0.01265334  0.740962648
#> Sphoeroides_pachygaster     -1.075318775 -0.8300723  0.08725136  0.370151714
#> Sphoeroides_testudineus     -1.494467174 -0.7216674  0.96448639  0.963449208
#> Stellifer_lanceolatus        0.003765399  1.2760813 -0.37355670 -0.245394436
#> Symphurus_plagiusa          -1.191253439 -0.8300723 -0.39621940  2.179709071
#> Synodus_foetens              0.717209482 -0.6364920 -0.81170208  0.140248935
#> Trichiurus_lepturus          1.412817462 -0.2493314 -1.06193597 -0.134151156
#> Urobatis_jamaicensis         2.153015698 -0.8300723 -0.77676376  2.179709071
#>                                    BdSh        BdSf        PfPo         PfSh
#> Achirus_lineatus            -1.17042446  0.17755278 -2.42187116 -2.042954308
#> Anchoa_mitchilli             0.50227657  3.46001748  1.08394608 -0.420805564
#> Archosargus_probatocephalus  0.39550842 -0.71122088  0.59866478  0.796405014
#> Archosargus_rhomboidalis     0.40135145 -0.86207612  0.47167528  2.071121604
#> Ariopsis_felis              -0.70404418 -0.74938907  1.23361228  0.041638612
#> Bagre_marinus               -0.70670010  0.03033261  0.58052343  0.157848677
#> Bairdiella_chrysoura        -0.05015565 -0.60398643  0.37643316  0.223740982
#> Bairdiella_ronchus          -0.16860987 -0.90569692  0.21769628  0.029658193
#> Cathorops_melanopus         -0.60683735 -0.35498440  1.10208744 -0.377676056
#> Cetengraulis_edentulus       0.24093364 -0.26229021  1.38781381  0.209963500
#> Chaetodipterus_faber         0.91766311 -0.42405066  0.51702868 -0.912002747
#> Chilomycterus_schoepfii     -0.81559299 -1.13652361 -0.75286632 -1.648798520
#> Chloroscombrus_chrysurus     1.34420453  0.61921332 -0.20409027  1.412198555
#> Citharichthys_spilopterus   -1.17945460  0.55559966 -2.42187116 -2.042954308
#> Conodon_nobilis             -0.17232816 -0.20412916  0.39003918  0.420818876
#> Cynoscion_arenarius         -0.02040930 -0.33317400  0.22223162  0.652639985
#> Cynoscion_nebulosus         -0.23713271 -1.21467753  0.26758502  0.050024905
#> Cynoscion_nothus            -0.07671489 -0.44404352  0.29026171  0.044633717
#> Dasyatis_sabina             -1.19698370 -1.13470608 -2.42187116  0.006296376
#> Diapterus_auratus            0.66588150  0.48835094  0.62587682  1.497858551
#> Diapterus_rhombeus           0.62285553  0.45200028  0.59866478  1.353494501
#> Dorosoma_petenense           0.56389401  0.95000433  1.21093558  0.262677344
#> Eucinostomus_argenteus       0.32804795  0.94455173  0.52609935  0.861698298
#> Eucinostomus_gula            0.37585458  0.16119499  0.46260460  1.491868342
#> Eucinostomus_melanopterus    0.07148567  1.04269852  0.42632189  0.921001372
#> Gobionellus_oceanicus       -0.43792058 -0.40042273 -0.81182573 -0.463935073
#> Menticirrhus_americanus     -0.28387698 -0.25320255  0.26304968 -0.262065012
#> Menticirrhus_saxatilis      -0.28812646 -0.36952466  0.29479705  0.242310632
#> Micropogonias_undulatus     -0.10061821 -0.20594669  0.24037298  0.759265715
#> Opsanus_beta                -0.80390693 -0.78755727 -0.94788590 -1.370852797
#> Orthopristis_chrysoptera     0.12407297 -1.17469181  0.34468578 -0.719717020
#> Peprilus_paru                1.54233647  0.39747428 -0.04988873  0.115917210
#> Polydactylus_octonemus       0.07892226 -0.25683762  0.73472496  0.952150462
#> Prionotus_carolinus         -0.68067205  1.32078108  0.25851434 -0.626868773
#> Prionotus_scitulus          -0.70510655  0.51379640  0.36736248 -0.033838028
#> Selene_setapinnis            3.46097608  1.76244161 -0.31293841  0.425611044
#> Selene_vomer                 3.45141475  2.50036003 -0.21769628  0.573569219
#> Sphoeroides_nephelus        -0.70085707  0.27933464  0.06803009 -1.310351681
#> Sphoeroides_pachygaster     -0.60311906  1.25898495 -0.10884814 -1.449324542
#> Sphoeroides_testudineus     -0.63233423 -0.79846247 -0.35829180 -1.449324542
#> Stellifer_lanceolatus       -0.23182086  0.19572812  0.20862561  0.707749913
#> Symphurus_plagiusa          -1.14492759 -0.20776422 -2.42187116 -2.042954308
#> Synodus_foetens             -0.67057953 -0.80209753 -0.25397900 -0.034437049
#> Trichiurus_lepturus          0.48687221 -1.23103533  0.75740166  0.390867828
#> Urobatis_jamaicensis        -1.18529764 -1.28192626 -2.42187116  0.535830899
#>                                    CpHt       CfSh        FsRt         FsSf
#> Achirus_lineatus            -0.91103742 -0.9176609 -1.07107750 -0.635307836
#> Anchoa_mitchilli             0.04457675  1.1875436 -0.48969182  0.628536554
#> Archosargus_probatocephalus  0.10795933  0.5445614  0.68922913 -0.335794028
#> Archosargus_rhomboidalis     0.33223612  0.8062507  0.15167809 -0.337525322
#> Ariopsis_felis               0.81004321  2.0364599 -0.17016041 -0.211140883
#> Bagre_marinus                1.23700493  1.5040886 -0.35011312 -0.185171478
#> Bairdiella_chrysoura        -0.32457449 -0.4680231 -0.42740050 -0.185171478
#> Bairdiella_ronchus          -0.27651518 -0.4176636 -0.22091630 -0.228453820
#> Cathorops_melanopus          0.76198389  1.6686561 -0.31550683 -0.439671650
#> Cetengraulis_edentulus       0.19850586  1.2046298 -0.56467212 -0.289914746
#> Chaetodipterus_faber         0.17203841  0.6066114 -0.29589660 -0.892404948
#> Chilomycterus_schoepfii      0.02228838 -0.6065115  2.12885081 -0.856913428
#> Chloroscombrus_chrysurus     2.48167149  1.0013935 -0.04673131 -0.057921392
#> Citharichthys_spilopterus   -1.00158395 -0.9176609 -1.07107750 -0.817959320
#> Conodon_nobilis             -0.48477220 -0.4896057 -0.02481399 -0.344450497
#> Cynoscion_arenarius         -0.37124078 -0.6397847 -0.64311305 -0.437940356
#> Cynoscion_nebulosus         -0.34477333 -0.5462601 -0.52545166 -0.282989571
#> Cynoscion_nothus            -0.42487219 -0.6226985 -0.59235715 -0.503729516
#> Dasyatis_sabina             -1.69182713 -1.6074054 -1.07107750  3.418516327
#> Diapterus_auratus            0.39910125  0.8718978 -0.38933358 -0.224991233
#> Diapterus_rhombeus           0.32318147  0.7891644 -0.16323915 -0.199021827
#> Dorosoma_petenense           0.23193842  0.3943824 -0.58658944 -0.049264924
#> Eucinostomus_argenteus       0.06547211  0.3835911 -0.18977064 -0.110725849
#> Eucinostomus_gula            0.28139076  0.8107470 -0.29935723 -0.270004868
#> Eucinostomus_melanopterus    0.32736054  0.6911434 -0.33857769 -0.269139222
#> Gobionellus_oceanicus       -1.03501651 -1.3376227 -0.20130607  2.180641343
#> Menticirrhus_americanus     -0.59203502 -0.6964391  0.39392212  0.009599062
#> Menticirrhus_saxatilis      -0.56835362 -0.7018348 -0.03058171 -0.083025150
#> Micropogonias_undulatus     -0.31412682 -0.7773739 -0.44816428 -0.161799013
#> Opsanus_beta                -0.72158621 -0.7710790  1.98927210 -0.432746475
#> Orthopristis_chrysoptera     0.37193729  1.1021124  0.30971348 -0.021564225
#> Peprilus_paru                1.16178166  0.7684811  0.46890242 -0.340122263
#> Polydactylus_octonemus      -0.07034769  0.7352079 -0.47469577  0.538509282
#> Prionotus_carolinus         -0.11004886 -0.6757558  3.58462209  2.076763721
#> Prionotus_scitulus          -0.13651631 -0.7476978  3.17511432  1.829188725
#> Selene_setapinnis            2.76166923  1.0076884  0.19205209 -0.052727511
#> Selene_vomer                 2.84594716  1.5724336 -0.02596753 -0.268273575
#> Sphoeroides_nephelus         0.10308374 -0.4626274  0.58887089 -0.381673311
#> Sphoeroides_pachygaster     -0.22427679 -0.5525550  0.81150469 -0.489879166
#> Sphoeroides_testudineus     -0.35731055 -0.6595688  0.31086702 -0.344450497
#> Stellifer_lanceolatus       -0.46805592 -0.7917623 -0.14939664 -0.488147873
#> Symphurus_plagiusa          -1.69182713 -1.6074054 -1.07107750 -1.637726880
#> Synodus_foetens              0.46318034  0.5427629 -0.40432964 -0.176515009
#> Trichiurus_lepturus         -1.69182713 -1.6074054 -1.07107750 -0.866435543
#> Urobatis_jamaicensis        -1.69182713 -1.6074054 -1.07107750  3.228939669
```
