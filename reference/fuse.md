# Compute FUSE (Functionally Unique, Specialized and Endangered)

This index takes into account species functional uniqueness (also called
Functional Originality), species specialisation and species IUCN status.

## Usage

``` r
fuse(sp_dist, sp_faxes_coord, nb_NN = 5, GE, standGE = FALSE)
```

## Arguments

- sp_dist:

  a dist object provided by
  [`funct.dist`](https://cmlmagneville.github.io/mFD/reference/funct.dist.md),
  [`daisy`](https://rdrr.io/pkg/cluster/man/daisy.html) or
  [`dist.ktab`](https://adeverse.github.io/ade4/reference/dist.ktab.html).

- sp_faxes_coord:

  a data frame with the coordinates of the species on a multidimensional
  space based on a selected number of axes derived from a Principal
  Coordinate Analysis (PCoA). The species are in rows and the PCOA axes
  are in column.

- nb_NN:

  a numerical value giving the number of nearest neighbor to consider.
  Default: `nb_NN = 5`.

- GE:

  a numerical vector giving the IUCN status rank (DD = NA, LC = 0, NT =
  1, VU = 2, EN = 3, CR = 4) or the IUCN extinction probability
  associated with each status. See Mooers *et al.* (2008) for further
  information. For example, DD = NA, LC = 0, NT = 0.1, VU = 0.4, EN =
  0.666, and CR = 0.999).

- standGE:

  a logical value to standardize the GE values.

## Value

A data frame with species in rows and the following metrics in columns:

- FUSE: functionally unique, specialized and endangered (see Pimiento
  *et al.* (2020);

- FUn_std: functional uniqueness standardized between 0 and 1 (see
  Mouillot *et al.* (2013);

- FSp_std: functional specialization standardized between 0 and 1 (see
  Mouillot *et al.* (2013);

## References

Mouillot *et al.* (2013) Rare species support vulnerable functions in
high-diversity ecosystems. *PLoS Biology*, **11**, e1001569.  
Pimiento *et al.* (2020) Functional diversity of marine megafauna in the
Anthropocene. *Science Advances*, **6**, eaay7650.  
Violle *et al.* (2007) Let the concept of trait be functional! *Oikos*,
**116**, 882-892.

## Author

Fabien Leprieur and Camille Albouy

## Examples

``` r
# Load species traits data:
sp_tr <- read.csv(system.file('extdata', 'data_traits_MMA_ursus.csv', 
  package = 'mFD'), dec = ',', sep = ';', header = TRUE, row.names = 1,
  na.strings='NA')

# Trait compilation and ordination:
dimorphism      <- ordered(sp_tr$dimorphism)
breeding_site   <- ordered(sp_tr$breeding_site)
social_behavior <- ordered(sp_tr$social_behavior)
weight_max      <- log(sp_tr$adult_weight_max)
social_group    <- log(sp_tr$social_group_mean)
 
# Trait Matrix construction:
sp_tr_end <- data.frame(
  main_diet = sp_tr$main_diet, 
  foraging_water_depth = sp_tr$foraging_water_depth,
  foraging_location = sp_tr$foraging_location, 
  fasting_strategy = sp_tr$fasting_strategy,
  female_sexual_maturity = sp_tr$female_sexual_maturity, 
  weaning = sp_tr$weaning,
  gestation = sp_tr$gestation, inter_litter = sp_tr$inter_litter,
  breeding_site = sp_tr$breeding_site, 
  social_group = sp_tr$social_group_mean,
  social_behavior = sp_tr$social_behavior, 
  weight_max = sp_tr$adult_weight_max,
  dimorphism = sp_tr$dimorphism)
  
rownames(sp_tr_end) <- rownames(sp_tr)
 
# Function weigthing vector:
v <- c(0.25, 0.25, 0.25, 0.25, 0.20, 0.20, 0.20, 0.20, 0.20, 0.5, 0.5, 0.5, 
  0.5)
   
# Gower distance calculation:
sp_tr_end$main_diet <- as.factor(sp_tr_end$main_diet)
sp_tr_end$foraging_water_depth <- as.factor(sp_tr_end$foraging_water_depth)
sp_tr_end$foraging_location <- as.factor(sp_tr_end$foraging_location)
sp_tr_end$breeding_site <- as.factor(sp_tr_end$breeding_site)
sp_tr_end$social_behavior <- as.factor(sp_tr_end$social_behavior)

sp_dist_tr <- cluster::daisy(sp_tr_end, metric = c('gower'), 
  type = list(symm = c(4)), weights = v)
  
# Principal coordinate analyses
Pcoa <- ade4::dudi.pco(ade4::quasieuclid(sp_dist_tr), scann = FALSE, 
                       nf = 40)
#> Warning: Zero distance(s)

sp_faxes_coord <- Pcoa$li[1:40]
 
# FUSE calculation:
 FUSE_res <- mFD::fuse(
    sp_dist        = sp_dist_tr, 
    sp_faxes_coord = as.matrix(sp_faxes_coord), 
    nb_NN          = 5,  
    GE             = sp_tr$IUCN_num,
    standGE        = TRUE)
 FUSE_res
#>                                   FUSE      FUn_std    FSp_std
#> Arctocephalus_australis     0.00000000 0.4066626897 0.44574291
#> Arctocephalus_forsteri      0.00000000 0.5170579368 0.51871477
#> Arctocephalus_galapagoensis 0.65704590 0.5592652556 0.47871747
#> Arctocephalus_gazella       0.00000000 0.6208798411 0.56465339
#> Arctocephalus_philippii     0.00000000 0.2779525182 0.44247109
#> Arctocephalus_pusillus      0.00000000 0.2442786896 0.37884040
#> Arctocephalus_townsendi     0.00000000 0.5852680899 0.64653619
#> Arctocephalus_tropicalis    0.00000000 0.2297455947 0.32036900
#> Balaena_mysticetus          0.00000000 0.7452467855 0.57306402
#> Balaenoptera_acutorostrata  0.00000000 0.4460157412 0.41489069
#> Balaenoptera_bonaerensis    0.13238904 0.2184186166 0.32978308
#> Balaenoptera_borealis       0.41103021 0.2139234148 0.39976514
#> Balaenoptera_edeni          0.00000000 0.4661991166 0.27554077
#> Balaenoptera_musculus       0.80941909 0.6292520413 0.70171703
#> Balaenoptera_omurai                 NA 0.2386976599 0.41958574
#> Balaenoptera_physalus       0.56665878 0.3536749038 0.52385972
#> Berardius_arnuxii                   NA 0.3691785951 0.24980784
#> Berardius_bairdii                   NA 0.3306752012 0.23813008
#> Callorhinus_ursinus         0.48846635 0.5012945401 0.60635396
#> Caperea_marginata           0.00000000 0.1832726206 0.40111351
#> Cephalorhynchus_commersonii 0.00000000 0.1835307609 0.17480251
#> Cephalorhynchus_eutropia    0.01817766 0.0069987566 0.06626082
#> Cephalorhynchus_heavisidii  0.01830039 0.0083171263 0.06542231
#> Cephalorhynchus_hectori     0.01749032 0.0070959999 0.01634259
#> Cystophora_cristata         0.45529572 0.5865158316 0.43824460
#> Delphinapterus_leucas       0.00000000 0.6739629047 0.32847428
#> Delphinus_capensis                  NA 0.0264311644 0.16472623
#> Delphinus_delphis           0.00000000 0.3984163792 0.16127757
#> Dugong_dugon                0.38031579 0.5379056252 0.30543867
#> Enhydra_lutris              0.71394550 0.7016950229 0.45056295
#> Erignathus_barbatus         0.00000000 0.5531416419 0.35107768
#> Eschrichtius_robustus       0.00000000 0.8271900919 0.61300225
#> Eubalaena_australis         0.00000000 0.3144018773 0.78268639
#> Eubalaena_glacialis         0.53448288 0.2810689133 0.54593451
#> Eubalaena_japonica          0.64071030 0.4661757614 0.54157645
#> Eumetopias_jubatus          0.23895925 0.4516548620 0.56433150
#> Feresa_attenuata            0.00000000 0.2955099743 0.27261560
#> Globicephala_macrorhynchus  0.00000000 0.5255589486 0.35104810
#> Globicephala_melas          0.00000000 0.4990595861 0.35763988
#> Grampus_griseus             0.00000000 0.5233269607 0.29251186
#> Halichoerus_grypus          0.00000000 0.1032026841 0.35586820
#> Histriophoca_fasciata       0.00000000 0.6521823871 0.28824480
#> Hydrurga_leptonyx           0.00000000 0.6622353793 0.44212015
#> Hyperoodon_ampullatus               NA 0.2547878402 0.12270144
#> Hyperoodon_planifrons       0.00000000 0.2557952131 0.12249240
#> Indopacetus_pacificus               NA 0.2850555343 0.28554866
#> Kogia_breviceps                     NA 0.2839209537 0.25751272
#> Kogia_sima                          NA 0.1456685909 0.25492545
#> Lagenodelphis_hosei         0.00000000 0.7241772659 0.51940242
#> Lagenorhynchus_acutus       0.00000000 0.4544319444 0.19736356
#> Lagenorhynchus_albirostris  0.00000000 0.0907709847 0.14280471
#> Lagenorhynchus_australis    0.00000000 0.0010393047 0.21058380
#> Lagenorhynchus_cruciger     0.00000000 0.0020371106 0.15621825
#> Lagenorhynchus_obliquidens  0.00000000 0.3327438642 0.20430708
#> Lagenorhynchus_obscurus     0.00000000 0.4240208007 0.17545748
#> Leptonychotes_weddellii     0.00000000 0.5404210014 0.55524903
#> Lissodelphis_borealis       0.00000000 0.5812804661 0.50809556
#> Lissodelphis_peronii        0.00000000 0.2947422923 0.43157141
#> Lobodon_carcinophaga        0.00000000 0.5833807958 0.41597696
#> Lontra_felina               0.59408663 0.5749777422 0.35413972
#> Megaptera_novaeangliae      0.00000000 0.3564239828 0.45460052
#> Mesoplodon_bidens                   NA 0.2401439269 0.23803164
#> Mesoplodon_bowdoini                 NA 0.0644576583 0.20799830
#> Mesoplodon_carlhubbsi               NA 0.0008237359 0.07159363
#> Mesoplodon_densirostris             NA 0.2058523276 0.19967333
#> Mesoplodon_europaeus                NA 0.1358313979 0.09000617
#> Mesoplodon_ginkgodens               NA 0.0013549553 0.07202895
#> Mesoplodon_grayi                    NA 0.0544755693 0.20125592
#> Mesoplodon_hectori                  NA 0.0533639848 0.20213100
#> Mesoplodon_layardii                 NA 0.0156247025 0.09936896
#> Mesoplodon_mirus                    NA 0.2615344649 0.02711371
#> Mesoplodon_perrini                  NA 0.0007873510 0.07105382
#> Mesoplodon_peruvianus               NA 0.0022081471 0.07148019
#> Mesoplodon_stejnegeri               NA 0.0019316220 0.07240253
#> Mesoplodon_traversii                NA 0.0000000000 0.09237153
#> Mirounga_angustirostris     0.00000000 1.0000000000 1.00000000
#> Mirounga_leonina            0.00000000 0.6948598959 0.62780483
#> Monachus_monachus           0.25265231 0.1977684465 0.16152114
#> Monachus_schauinslandi      0.69523452 0.5791027790 0.52972758
#> Monodon_monoceros           0.00000000 0.7032597866 0.45028422
#> Neophoca_cinerea            0.36345229 0.1434988631 0.39804367
#> Neophocaena_phocaenoides    0.32176325 0.5420666711 0.17076625
#> Odobenus_rosmarus                   NA 0.5803602059 0.45366482
#> Ommatophoca_rossii          0.00000000 0.6659562055 0.34130521
#> Orcaella_brevirostris       0.04133240 0.0101282373 0.04578858
#> Orcaella_heinsohni          0.05624410 0.0292679517 0.08519702
#> Orcinus_orca                        NA 0.7536273607 0.39332028
#> Otaria_flavescens           0.00000000 0.0808815948 0.34394151
#> Pagophilus_groenlandicus    0.00000000 0.7570041744 0.46023398
#> Peponocephala_electra       0.00000000 0.5314927079 0.48063075
#> Phoca_largha                0.00000000 0.5982829480 0.33477382
#> Phoca_vitulina              0.00000000 0.4543455128 0.37496303
#> Phocarctos_hookeri          0.29111541 0.0782092491 0.35171865
#> Phocoena_dioptrica          0.00000000 0.2513166389 0.41532085
#> Phocoena_phocoena           0.00000000 0.0736635094 0.00000000
#> Phocoena_sinus              0.30989600 0.1864775088 0.14901742
#> Phocoena_spinipinnis        0.11003937 0.2651030237 0.18774230
#> Phocoenoides_dalli          0.00000000 0.2916871646 0.15983045
#> Physeter_macrocephalus      0.54891552 0.7564270292 0.51249076
#> Pontoporia_blainvillei      0.27464151 0.3710816347 0.22018302
#> Pseudorca_crassidens        0.17203888 0.4744533574 0.24712980
#> Pusa_caspica                0.36513400 0.1941091818 0.34349375
#> Pusa_hispida                0.00000000 0.1823213934 0.15240645
#> Sotalia_guianensis          0.01946782 0.0546637568 0.02364730
#> Sousa_chinensis             0.20619793 0.5623262797 0.31006940
#> Sousa_teuszii               0.38826398 0.2726857890 0.15850980
#> Stenella_attenuata          0.00000000 0.4500796407 0.27718723
#> Stenella_clymene            0.00000000 0.3747462296 0.40334210
#> Stenella_coeruleoalba       0.00000000 0.2311412497 0.18363765
#> Stenella_frontalis          0.00000000 0.3490478318 0.23845778
#> Stenella_longirostris       0.00000000 0.5133002403 0.18285982
#> Steno_bredanensis           0.00000000 0.1334810508 0.21805988
#> Tasmacetus_shepherdi                NA 0.3474343349 0.17676343
#> Trichechus_manatus          0.33814860 0.5037922863 0.24035976
#> Trichechus_senegalensis     0.36824108 0.5315342856 0.28350122
#> Tursiops_aduncus            0.04686502 0.0517896426 0.13834135
#> Tursiops_truncatus          0.00000000 0.3052280508 0.19082721
#> Ursus_maritimus             0.67482557 0.9090386502 0.70012282
#> Zalophus_californianus      0.00000000 0.1204667495 0.62677747
#> Zalophus_wollebaeki         0.32738815 0.0749554707 0.41799919
#> Ziphius_cavirostris         0.00000000 0.1383689569 0.09353647
    
 FUSE_res2 <- mFD::fuse(
    sp_dist        = sp_dist_tr, 
    sp_faxes_coord = as.matrix(sp_faxes_coord), 
    nb_NN          = 5,
    GE             = sp_tr$IUCN_50,
    standGE        = TRUE)
 FUSE_res2
#>                                     FUSE      FUn_std    FSp_std
#> Arctocephalus_australis     0.0000000000 0.4066626897 0.44574291
#> Arctocephalus_forsteri      0.0000000000 0.5170579368 0.51871477
#> Arctocephalus_galapagoensis 0.4051935747 0.5592652556 0.47871747
#> Arctocephalus_gazella       0.0000000000 0.6208798411 0.56465339
#> Arctocephalus_philippii     0.0000000000 0.2779525182 0.44247109
#> Arctocephalus_pusillus      0.0000000000 0.2442786896 0.37884040
#> Arctocephalus_townsendi     0.0000000000 0.5852680899 0.64653619
#> Arctocephalus_tropicalis    0.0000000000 0.2297455947 0.32036900
#> Balaena_mysticetus          0.0000000000 0.7452467855 0.57306402
#> Balaenoptera_acutorostrata  0.0000000000 0.4460157412 0.41489069
#> Balaenoptera_bonaerensis    0.0022311864 0.2184186166 0.32978308
#> Balaenoptera_borealis       0.2482137302 0.2139234148 0.39976514
#> Balaenoptera_edeni          0.0000000000 0.4661991166 0.27554077
#> Balaenoptera_musculus       0.5062324442 0.6292520413 0.70171703
#> Balaenoptera_omurai                   NA 0.2386976599 0.41958574
#> Balaenoptera_physalus       0.3468953473 0.3536749038 0.52385972
#> Berardius_arnuxii                     NA 0.3691785951 0.24980784
#> Berardius_bairdii                     NA 0.3306752012 0.23813008
#> Callorhinus_ursinus         0.0562359308 0.5012945401 0.60635396
#> Caperea_marginata           0.0000000000 0.1832726206 0.40111351
#> Cephalorhynchus_commersonii 0.0000000000 0.1835307609 0.17480251
#> Cephalorhynchus_eutropia    0.0002983037 0.0069987566 0.06626082
#> Cephalorhynchus_heavisidii  0.0003002586 0.0083171263 0.06542231
#> Cephalorhynchus_hectori     0.0101183562 0.0070959999 0.01634259
#> Cystophora_cristata         0.0520745374 0.5865158316 0.43824460
#> Delphinapterus_leucas       0.0000000000 0.6739629047 0.32847428
#> Delphinus_capensis                    NA 0.0264311644 0.16472623
#> Delphinus_delphis           0.0000000000 0.3984163792 0.16127757
#> Dugong_dugon                0.0429309730 0.5379056252 0.30543867
#> Enhydra_lutris              0.4434975661 0.7016950229 0.45056295
#> Erignathus_barbatus         0.0000000000 0.5531416419 0.35107768
#> Eschrichtius_robustus       0.0000000000 0.8271900919 0.61300225
#> Eubalaena_australis         0.0000000000 0.3144018773 0.78268639
#> Eubalaena_glacialis         0.3270161133 0.2810689133 0.54593451
#> Eubalaena_japonica          0.3945009065 0.4661757614 0.54157645
#> Eumetopias_jubatus          0.0041331511 0.4516548620 0.56433150
#> Feresa_attenuata            0.0000000000 0.2955099743 0.27261560
#> Globicephala_macrorhynchus  0.0000000000 0.5255589486 0.35104810
#> Globicephala_melas          0.0000000000 0.4990595861 0.35763988
#> Grampus_griseus             0.0000000000 0.5233269607 0.29251186
#> Halichoerus_grypus          0.0000000000 0.1032026841 0.35586820
#> Histriophoca_fasciata       0.0000000000 0.6521823871 0.28824480
#> Hydrurga_leptonyx           0.0000000000 0.6622353793 0.44212015
#> Hyperoodon_ampullatus                 NA 0.2547878402 0.12270144
#> Hyperoodon_planifrons       0.0000000000 0.2557952131 0.12249240
#> Indopacetus_pacificus                 NA 0.2850555343 0.28554866
#> Kogia_breviceps                       NA 0.2839209537 0.25751272
#> Kogia_sima                            NA 0.1456685909 0.25492545
#> Lagenodelphis_hosei         0.0000000000 0.7241772659 0.51940242
#> Lagenorhynchus_acutus       0.0000000000 0.4544319444 0.19736356
#> Lagenorhynchus_albirostris  0.0000000000 0.0907709847 0.14280471
#> Lagenorhynchus_australis    0.0000000000 0.0010393047 0.21058380
#> Lagenorhynchus_cruciger     0.0000000000 0.0020371106 0.15621825
#> Lagenorhynchus_obliquidens  0.0000000000 0.3327438642 0.20430708
#> Lagenorhynchus_obscurus     0.0000000000 0.4240208007 0.17545748
#> Leptonychotes_weddellii     0.0000000000 0.5404210014 0.55524903
#> Lissodelphis_borealis       0.0000000000 0.5812804661 0.50809556
#> Lissodelphis_peronii        0.0000000000 0.2947422923 0.43157141
#> Lobodon_carcinophaga        0.0000000000 0.5833807958 0.41597696
#> Lontra_felina               0.3649494018 0.5749777422 0.35413972
#> Megaptera_novaeangliae      0.0000000000 0.3564239828 0.45460052
#> Mesoplodon_bidens                     NA 0.2401439269 0.23803164
#> Mesoplodon_bowdoini                   NA 0.0644576583 0.20799830
#> Mesoplodon_carlhubbsi                 NA 0.0008237359 0.07159363
#> Mesoplodon_densirostris               NA 0.2058523276 0.19967333
#> Mesoplodon_europaeus                  NA 0.1358313979 0.09000617
#> Mesoplodon_ginkgodens                 NA 0.0013549553 0.07202895
#> Mesoplodon_grayi                      NA 0.0544755693 0.20125592
#> Mesoplodon_hectori                    NA 0.0533639848 0.20213100
#> Mesoplodon_layardii                   NA 0.0156247025 0.09936896
#> Mesoplodon_mirus                      NA 0.2615344649 0.02711371
#> Mesoplodon_perrini                    NA 0.0007873510 0.07105382
#> Mesoplodon_peruvianus                 NA 0.0022081471 0.07148019
#> Mesoplodon_stejnegeri                 NA 0.0019316220 0.07240253
#> Mesoplodon_traversii                  NA 0.0000000000 0.09237153
#> Mirounga_angustirostris     0.0000000000 1.0000000000 1.00000000
#> Mirounga_leonina            0.0000000000 0.6948598959 0.62780483
#> Monachus_monachus           0.1497520468 0.1977684465 0.16152114
#> Monachus_schauinslandi      0.4302127336 0.5791027790 0.52972758
#> Monodon_monoceros           0.0000000000 0.7032597866 0.45028422
#> Neophoca_cinerea            0.2192750379 0.1434988631 0.39804367
#> Neophocaena_phocaenoides    0.0362881475 0.5420666711 0.17076625
#> Odobenus_rosmarus                     NA 0.5803602059 0.45366482
#> Ommatophoca_rossii          0.0000000000 0.6659562055 0.34130521
#> Orcaella_brevirostris       0.0240062354 0.0101282373 0.04578858
#> Orcaella_heinsohni          0.0058839286 0.0292679517 0.08519702
#> Orcinus_orca                          NA 0.7536273607 0.39332028
#> Otaria_flavescens           0.0000000000 0.0808815948 0.34394151
#> Pagophilus_groenlandicus    0.0000000000 0.7570041744 0.46023398
#> Peponocephala_electra       0.0000000000 0.5314927079 0.48063075
#> Phoca_largha                0.0000000000 0.5982829480 0.33477382
#> Phoca_vitulina              0.0000000000 0.4543455128 0.37496303
#> Phocarctos_hookeri          0.1750436555 0.0782092491 0.35171865
#> Phocoena_dioptrica          0.0000000000 0.2513166389 0.41532085
#> Phocoena_phocoena           0.0000000000 0.0736635094 0.00000000
#> Phocoena_sinus              0.3098959997 0.1864775088 0.14901742
#> Phocoena_spinipinnis        0.0018432814 0.2651030237 0.18774230
#> Phocoenoides_dalli          0.0000000000 0.2916871646 0.15983045
#> Physeter_macrocephalus      0.0642642678 0.7564270292 0.51249076
#> Pontoporia_blainvillei      0.0302045496 0.3710816347 0.22018302
#> Pseudorca_crassidens        0.0029361868 0.4744533574 0.24712980
#> Pusa_caspica                0.2193438717 0.1941091818 0.34349375
#> Pusa_hispida                0.0000000000 0.1823213934 0.15240645
#> Sotalia_guianensis          0.0003188826 0.0546637568 0.02364730
#> Sousa_chinensis             0.0035493076 0.5623262797 0.31006940
#> Sousa_teuszii               0.3882639834 0.2726857890 0.15850980
#> Stenella_attenuata          0.0000000000 0.4500796407 0.27718723
#> Stenella_clymene            0.0000000000 0.3747462296 0.40334210
#> Stenella_coeruleoalba       0.0000000000 0.2311412497 0.18363765
#> Stenella_frontalis          0.0000000000 0.3490478318 0.23845778
#> Stenella_longirostris       0.0000000000 0.5133002403 0.18285982
#> Steno_bredanensis           0.0000000000 0.1334810508 0.21805988
#> Tasmacetus_shepherdi                  NA 0.3474343349 0.17676343
#> Trichechus_manatus          0.0379151532 0.5037922863 0.24035976
#> Trichechus_senegalensis     0.0414988099 0.5315342856 0.28350122
#> Tursiops_aduncus            0.0007741038 0.0517896426 0.13834135
#> Tursiops_truncatus          0.0000000000 0.3052280508 0.19082721
#> Ursus_maritimus             0.0811703393 0.9090386502 0.70012282
#> Zalophus_californianus      0.0000000000 0.1204667495 0.62677747
#> Zalophus_wollebaeki         0.1982794545 0.0749554707 0.41799919
#> Ziphius_cavirostris         0.0000000000 0.1383689569 0.09353647
    
 FUSE_res3 <- mFD::fuse(
    sp_dist        = sp_dist_tr, 
    sp_faxes_coord = as.matrix(sp_faxes_coord), 
    nb_NN          = 5, 
    GE             = sp_tr$IUCN_100,
    standGE        = TRUE)
 FUSE_res3
#>                                     FUSE      FUn_std    FSp_std
#> Arctocephalus_australis     0.0000000000 0.4066626897 0.44574291
#> Arctocephalus_forsteri      0.0000000000 0.5170579368 0.51871477
#> Arctocephalus_galapagoensis 0.5992815838 0.5592652556 0.47871747
#> Arctocephalus_gazella       0.0000000000 0.6208798411 0.56465339
#> Arctocephalus_philippii     0.0000000000 0.2779525182 0.44247109
#> Arctocephalus_pusillus      0.0000000000 0.2442786896 0.37884040
#> Arctocephalus_townsendi     0.0000000000 0.5852680899 0.64653619
#> Arctocephalus_tropicalis    0.0000000000 0.2297455947 0.32036900
#> Balaena_mysticetus          0.0000000000 0.7452467855 0.57306402
#> Balaenoptera_acutorostrata  0.0000000000 0.4460157412 0.41489069
#> Balaenoptera_bonaerensis    0.0054747614 0.2184186166 0.32978308
#> Balaenoptera_borealis       0.3731209005 0.2139234148 0.39976514
#> Balaenoptera_edeni          0.0000000000 0.4661991166 0.27554077
#> Balaenoptera_musculus       0.7405518753 0.6292520413 0.70171703
#> Balaenoptera_omurai                   NA 0.2386976599 0.41958574
#> Balaenoptera_physalus       0.5159872370 0.3536749038 0.52385972
#> Berardius_arnuxii                     NA 0.3691785951 0.24980784
#> Berardius_bairdii                     NA 0.3306752012 0.23813008
#> Callorhinus_ursinus         0.1087457494 0.5012945401 0.60635396
#> Caperea_marginata           0.0000000000 0.1832726206 0.40111351
#> Cephalorhynchus_commersonii 0.0000000000 0.1835307609 0.17480251
#> Cephalorhynchus_eutropia    0.0007324479 0.0069987566 0.06626082
#> Cephalorhynchus_heavisidii  0.0007372515 0.0083171263 0.06542231
#> Cephalorhynchus_hectori     0.0157191188 0.0070959999 0.01634259
#> Cystophora_cristata         0.1007823971 0.5865158316 0.43824460
#> Delphinapterus_leucas       0.0000000000 0.6739629047 0.32847428
#> Delphinus_capensis                    NA 0.0264311644 0.16472623
#> Delphinus_delphis           0.0000000000 0.3984163792 0.16127757
#> Dugong_dugon                0.0832219750 0.5379056252 0.30543867
#> Enhydra_lutris              0.6522096715 0.7016950229 0.45056295
#> Erignathus_barbatus         0.0000000000 0.5531416419 0.35107768
#> Eschrichtius_robustus       0.0000000000 0.8271900919 0.61300225
#> Eubalaena_australis         0.0000000000 0.3144018773 0.78268639
#> Eubalaena_glacialis         0.4866156416 0.2810689133 0.54593451
#> Eubalaena_japonica          0.5841787823 0.4661757614 0.54157645
#> Eumetopias_jubatus          0.0101348519 0.4516548620 0.56433150
#> Feresa_attenuata            0.0000000000 0.2955099743 0.27261560
#> Globicephala_macrorhynchus  0.0000000000 0.5255589486 0.35104810
#> Globicephala_melas          0.0000000000 0.4990595861 0.35763988
#> Grampus_griseus             0.0000000000 0.5233269607 0.29251186
#> Halichoerus_grypus          0.0000000000 0.1032026841 0.35586820
#> Histriophoca_fasciata       0.0000000000 0.6521823871 0.28824480
#> Hydrurga_leptonyx           0.0000000000 0.6622353793 0.44212015
#> Hyperoodon_ampullatus                 NA 0.2547878402 0.12270144
#> Hyperoodon_planifrons       0.0000000000 0.2557952131 0.12249240
#> Indopacetus_pacificus                 NA 0.2850555343 0.28554866
#> Kogia_breviceps                       NA 0.2839209537 0.25751272
#> Kogia_sima                            NA 0.1456685909 0.25492545
#> Lagenodelphis_hosei         0.0000000000 0.7241772659 0.51940242
#> Lagenorhynchus_acutus       0.0000000000 0.4544319444 0.19736356
#> Lagenorhynchus_albirostris  0.0000000000 0.0907709847 0.14280471
#> Lagenorhynchus_australis    0.0000000000 0.0010393047 0.21058380
#> Lagenorhynchus_cruciger     0.0000000000 0.0020371106 0.15621825
#> Lagenorhynchus_obliquidens  0.0000000000 0.3327438642 0.20430708
#> Lagenorhynchus_obscurus     0.0000000000 0.4240208007 0.17545748
#> Leptonychotes_weddellii     0.0000000000 0.5404210014 0.55524903
#> Lissodelphis_borealis       0.0000000000 0.5812804661 0.50809556
#> Lissodelphis_peronii        0.0000000000 0.2947422923 0.43157141
#> Lobodon_carcinophaga        0.0000000000 0.5833807958 0.41597696
#> Lontra_felina               0.5413762905 0.5749777422 0.35413972
#> Megaptera_novaeangliae      0.0000000000 0.3564239828 0.45460052
#> Mesoplodon_bidens                     NA 0.2401439269 0.23803164
#> Mesoplodon_bowdoini                   NA 0.0644576583 0.20799830
#> Mesoplodon_carlhubbsi                 NA 0.0008237359 0.07159363
#> Mesoplodon_densirostris               NA 0.2058523276 0.19967333
#> Mesoplodon_europaeus                  NA 0.1358313979 0.09000617
#> Mesoplodon_ginkgodens                 NA 0.0013549553 0.07202895
#> Mesoplodon_grayi                      NA 0.0544755693 0.20125592
#> Mesoplodon_hectori                    NA 0.0533639848 0.20213100
#> Mesoplodon_layardii                   NA 0.0156247025 0.09936896
#> Mesoplodon_mirus                      NA 0.2615344649 0.02711371
#> Mesoplodon_perrini                    NA 0.0007873510 0.07105382
#> Mesoplodon_peruvianus                 NA 0.0022081471 0.07148019
#> Mesoplodon_stejnegeri                 NA 0.0019316220 0.07240253
#> Mesoplodon_traversii                  NA 0.0000000000 0.09237153
#> Mirounga_angustirostris     0.0000000000 1.0000000000 1.00000000
#> Mirounga_leonina            0.0000000000 0.6948598959 0.62780483
#> Monachus_monachus           0.2283726931 0.1977684465 0.16152114
#> Monachus_schauinslandi      0.6345952150 0.5791027790 0.52972758
#> Monodon_monoceros           0.0000000000 0.7032597866 0.45028422
#> Neophoca_cinerea            0.3298536309 0.1434988631 0.39804367
#> Neophocaena_phocaenoides    0.0703478821 0.5420666711 0.17076625
#> Odobenus_rosmarus                     NA 0.5803602059 0.45366482
#> Ommatophoca_rossii          0.0000000000 0.6659562055 0.34130521
#> Orcaella_brevirostris       0.0371819980 0.0101282373 0.04578858
#> Orcaella_heinsohni          0.0115106174 0.0292679517 0.08519702
#> Orcinus_orca                          NA 0.7536273607 0.39332028
#> Otaria_flavescens           0.0000000000 0.0808815948 0.34394151
#> Pagophilus_groenlandicus    0.0000000000 0.7570041744 0.46023398
#> Peponocephala_electra       0.0000000000 0.5314927079 0.48063075
#> Phoca_largha                0.0000000000 0.5982829480 0.33477382
#> Phoca_vitulina              0.0000000000 0.4543455128 0.37496303
#> Phocarctos_hookeri          0.2640009489 0.0782092491 0.35171865
#> Phocoena_dioptrica          0.0000000000 0.2513166389 0.41532085
#> Phocoena_phocoena           0.0000000000 0.0736635094 0.00000000
#> Phocoena_sinus              0.3098959997 0.1864775088 0.14901742
#> Phocoena_spinipinnis        0.0045236417 0.2651030237 0.18774230
#> Phocoenoides_dalli          0.0000000000 0.2916871646 0.15983045
#> Physeter_macrocephalus      0.1239917183 0.7564270292 0.51249076
#> Pontoporia_blainvillei      0.0587425276 0.3710816347 0.22018302
#> Pseudorca_crassidens        0.0072022892 0.4744533574 0.24712980
#> Pusa_caspica                0.3310632694 0.1941091818 0.34349375
#> Pusa_hispida                0.0000000000 0.1823213934 0.15240645
#> Sotalia_guianensis          0.0007830123 0.0546637568 0.02364730
#> Sousa_chinensis             0.0087042852 0.5623262797 0.31006940
#> Sousa_teuszii               0.3882639834 0.2726857890 0.15850980
#> Stenella_attenuata          0.0000000000 0.4500796407 0.27718723
#> Stenella_clymene            0.0000000000 0.3747462296 0.40334210
#> Stenella_coeruleoalba       0.0000000000 0.2311412497 0.18363765
#> Stenella_frontalis          0.0000000000 0.3490478318 0.23845778
#> Stenella_longirostris       0.0000000000 0.5133002403 0.18285982
#> Steno_bredanensis           0.0000000000 0.1334810508 0.21805988
#> Tasmacetus_shepherdi                  NA 0.3474343349 0.17676343
#> Trichechus_manatus          0.0735594963 0.5037922863 0.24035976
#> Trichechus_senegalensis     0.0804618871 0.5315342856 0.28350122
#> Tursiops_aduncus            0.0019004117 0.0517896426 0.13834135
#> Tursiops_truncatus          0.0000000000 0.3052280508 0.19082721
#> Ursus_maritimus             0.1560437697 0.9090386502 0.70012282
#> Zalophus_californianus      0.0000000000 0.1204667495 0.62677747
#> Zalophus_wollebaeki         0.2973796407 0.0749554707 0.41799919
#> Ziphius_cavirostris         0.0000000000 0.1383689569 0.09353647
```
