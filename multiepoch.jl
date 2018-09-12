


oifitsfiles = ["./DATA/2011Sep02.lam_And_prepped.oifits", "./DATA/2011Sep06.lam_And_prepped.oifits",
"./DATA/2011Sep10.lam_And_prepped.oifits","./DATA/2011Sep14.lam_And_prepped.oifits",
"./DATA/2011Sep19.lam_And_prepped.oifits","./DATA/2011Sep24.lam_And_prepped.oifits"];

nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles);
