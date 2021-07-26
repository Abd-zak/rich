#################################################################################################################
#################################################################################################################
cd Plots_files/
hadd hist_iphi_08_10.root  hist_iphi.P08.root hist_iphi.P09.root hist_iphi.P10.root
hadd hist_ephi_08_10.root hist_ephi.P08.root hist_ephi.P09.root hist_ephi.P10.root
hadd hist_K0_08_10.root hist_K0.P08.root hist_K0.P09.root hist_K0.P10.root
hadd hist_Lambda_08_10.root hist_Lambda.P08.root hist_Lambda.P09.root hist_Lambda.P10.root
#################################################################################################################
#################################################################################################################
hadd hist_iphi_07_10.root  hist_iphi.P07.root hist_iphi.P08.root hist_iphi.P09.root hist_iphi.P10.root
hadd hist_ephi_07_10.root hist_ephi.P07.root hist_ephi.P08.root hist_ephi.P09.root hist_ephi.P10.root
hadd hist_K0_07_10.root hist_K0.P07.root hist_K0.P08.root hist_K0.P09.root hist_K0.P10.root
hadd hist_Lambda_07_10.root hist_Lambda.P07.root hist_Lambda.P08.root hist_Lambda.P09.root hist_Lambda.P10.root
#################################################################################################################
#################################################################################################################
hadd hist_iphi_06_10.root  hist_iphi.P06.root hist_iphi.P07.root hist_iphi.P08.root hist_iphi.P09.root hist_iphi.P10.root
hadd hist_ephi_06_10.root hist_ephi.P06.root hist_ephi.P07.root hist_ephi.P08.root hist_ephi.P09.root hist_ephi.P10.root
hadd hist_K0_06_10.root hist_K0.P06.root hist_K0.P07.root hist_K0.P08.root hist_K0.P09.root hist_K0.P10.root
hadd hist_Lambda_06_10.root hist_Lambda.P06.root hist_Lambda.P07.root hist_Lambda.P08.root hist_Lambda.P09.root hist_Lambda.P10.root
#################################################################################################################
#################################################################################################################
hadd hist_iphi_05_10.root  hist_iphi.P05.root hist_iphi.P06.root  hist_iphi.P07.root hist_iphi.P08.root hist_iphi.P09.root hist_iphi.P10.root
hadd hist_ephi_05_10.root hist_ephi.P05.root hist_ephi.P06.root hist_ephi.P07.root hist_ephi.P08.root hist_ephi.P09.root hist_ephi.P10.root
hadd hist_K0_05_10.root hist_K0.P05.root hist_K0.P06.root hist_K0.P07.root hist_K0.P08.root hist_K0.P09.root hist_K0.P10.root
hadd hist_Lambda_05_10.root hist_Lambda.P05.root hist_Lambda.P06.root hist_Lambda.P07.root hist_Lambda.P08.root hist_Lambda.P09.root hist_Lambda.P10.root

#################################################################################################################
#################################################################################################################
hadd hist_iphi_07_09.root  hist_iphi.P07.root hist_iphi.P08.root hist_iphi.P09.root
hadd hist_ephi_07_09.root hist_ephi.P07.root hist_ephi.P08.root hist_ephi.P09.root
hadd hist_K0_07_09.root hist_K0.P07.root hist_K0.P08.root hist_K0.P09.root
hadd hist_Lambda_07_09.root hist_Lambda.P07.root hist_Lambda.P08.root hist_Lambda.P09.root
#################################################################################################################
#################################################################################################################




hadd hist_iphi_0708_10.root  hist_iphi.P07.root hist_iphi.P08.root hist_iphi.P10.root
hadd hist_ephi_0708_10.root hist_ephi.P07.root hist_ephi.P08.root hist_ephi.P09.root hist_ephi.P10.root
hadd hist_K0_0708_10.root hist_K0.P07.root hist_K0.P08.root hist_K0.P09.root hist_K0.P10.root
hadd hist_Lambda_0708_10.root hist_Lambda.P07.root hist_Lambda.P08.root hist_Lambda.P09.root hist_Lambda.P10.root
#################################################################################################################
#################################################################################################################
hadd hist_iphi_0709_10.root  hist_iphi.P07.root hist_iphi.P08.root hist_iphi.P09.root hist_iphi.P10.root
hadd hist_ephi_0709_10.root hist_ephi.P07.root hist_ephi.P08.root hist_ephi.P09.root hist_ephi.P10.root
hadd hist_K0_0709_10.root hist_K0.P07.root hist_K0.P08.root hist_K0.P09.root hist_K0.P10.root
hadd hist_Lambda_0709_10.root hist_Lambda.P07.root hist_Lambda.P08.root hist_Lambda.P09.root hist_Lambda.P10.root


#################################################################################################################
#################################################################################################################
hadd hist_iphi_09_10.root  hist_iphi.P09.root hist_iphi.P10.root
hadd hist_ephi_09_10.root hist_ephi.P09.root hist_ephi.P10.root
hadd hist_K0_09_10.root hist_K0.P09.root hist_K0.P10.root
hadd hist_Lambda_09_10.root hist_Lambda.P09.root hist_Lambda.P10.root

#################################################################################################################
#################################################################################################################
hadd hist_iphi_08_09.root  hist_iphi.P09.root hist_iphi.P08.root
hadd hist_ephi_08_09.root hist_ephi.P09.root hist_ephi.P08.root
hadd hist_K0_08_09.root hist_K0.P09.root hist_K0.P08.root
hadd hist_Lambda_08_09.root hist_Lambda.P09.root hist_Lambda.P08.root



#################################################################################################################
#################################################################################################################
mkdir files_fit

mv hist_*_*.root ./files_fit
cd ..
#################################################################################################################
#################################################################################################################
