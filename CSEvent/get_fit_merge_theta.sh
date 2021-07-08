# ------------------------------------------------------------------------------------
# Merge theta bin on the 0th bin and the 2th bin to the 1st
# ------------------------s------------------------------------------------------------
mkdir results_fit_theta_cols
mkdir results_fit_theta_cols/P08_10
mkdir results_fit_theta_cols/P07_10
mkdir results_fit_theta_cols/P06_10
mkdir results_fit_theta_cols/P05_10
mkdir results_fit_theta_cols/P0708_10
mkdir results_fit_theta_cols/P0709_10
sed -i 's/merged_h: true/merged_h: false/g' options_fit/options_fit*
sed -i 's/cols_theta_bin: false/cols_theta_bin: true/g' options_fit/options_fit*
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------<<<8->10>>>>------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi.root+hist_file_iphi: ./hist_iphi_08_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi.root+hist_file_ephi: ./hist_ephi_08_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0.root+hist_file_K0: ./hist_K0_08_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda.root+hist_file_Lam: ./hist_Lambda_08_10.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich.root+out_file: rich_08_10.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------
cp Plots_files/files_fit/hist_*_08_10.root . 
fit_table -f options_fit/options_fit.dat  fit |& tee log_fit.log
mv *.root *.txt *.pdf *.log ./results_fit_theta_cols/P08_10
# ------------------------------------------------------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi_08_10.root+hist_file_iphi: ./hist_iphi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi_08_10.root+hist_file_ephi: ./hist_ephi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0_08_10.root+hist_file_K0: ./hist_K0.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda_08_10.root+hist_file_Lam: ./hist_Lambda.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich_08_10.root+out_file: rich.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------<<<7->10>>>>------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi.root+hist_file_iphi: ./hist_iphi_07_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi.root+hist_file_ephi: ./hist_ephi_07_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0.root+hist_file_K0: ./hist_K0_07_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda.root+hist_file_Lam: ./hist_Lambda_07_10.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich.root+out_file: rich_07_10.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------  
cp Plots_files/files_fit/hist_*_07_10.root .
fit_table -f options_fit/options_fit.dat  fit |& tee log_fit.log
mv *.root *.txt *.pdf *.log ./results_fit_theta_cols/P07_10
# ------------------------------------------------------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi_07_10.root+hist_file_iphi: ./hist_iphi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi_07_10.root+hist_file_ephi: ./hist_ephi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0_07_10.root+hist_file_K0: ./hist_K0.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda_07_10.root+hist_file_Lam: ./hist_Lambda.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich_07_10.root+out_file: rich.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------<<<6->10>>>>------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi.root+hist_file_iphi: ./hist_iphi_06_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi.root+hist_file_ephi: ./hist_ephi_06_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0.root+hist_file_K0: ./hist_K0_06_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda.root+hist_file_Lam: ./hist_Lambda_06_10.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich.root+out_file: rich_06_10.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------ 
cp Plots_files/files_fit/hist_*_06_10.root . 
fit_table -f options_fit/options_fit.dat  fit |& tee log_fit.log
mv *.root *.txt *.pdf *.log ./results_fit_theta_cols/P06_10
# ------------------------------------------------------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi_06_10.root+hist_file_iphi: ./hist_iphi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi_06_10.root+hist_file_ephi: ./hist_ephi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0_06_10.root+hist_file_K0: ./hist_K0.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda_06_10.root+hist_file_Lam: ./hist_Lambda.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich_06_10.root+out_file: rich.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------<<<5->10>>>>------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi.root+hist_file_iphi: ./hist_iphi_05_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi.root+hist_file_ephi: ./hist_ephi_05_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0.root*+hist_file_K0: ./hist_K0_05_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda.root+hist_file_Lam: ./hist_Lambda_05_10.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich.root+out_file: rich_05_10.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------ 
cp Plots_files/files_fit/hist_*_05_10.root . 
fit_table -f options_fit/options_fit.dat  fit |& tee log_fit.log
mv *.root *.txt *.pdf *.log ./results_fit_theta_cols/P05_10
# ------------------------------------------------------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi_05_10.root+hist_file_iphi: ./hist_iphi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi_05_10.root+hist_file_ephi: ./hist_ephi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0_05_10.root+hist_file_K0: ./hist_K0.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda_05_10.root+hist_file_Lam: ./hist_Lambda.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich_05_10.root+out_file: rich.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------



# ------------------------------------<<<0709->10>>>>------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi.root+hist_file_iphi: ./hist_iphi_0709_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi.root+hist_file_ephi: ./hist_ephi_0709_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0.root+hist_file_K0: ./hist_K0_0709_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda.root+hist_file_Lam: ./hist_Lambda_0709_10.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich.root+out_file: rich_0709_10.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------ 
cp Plots_files/files_fit/hist_*_0709_10.root . 
fit_table -f options_fit/options_fit.dat  fit |& tee log_fit.log
mv *.root *.txt *.pdf *.log ./results_fit_theta_cols/P0709_10
# ------------------------------------------------------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi_0709_10.root+hist_file_iphi: ./hist_iphi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi_0709_10.root+hist_file_ephi: ./hist_ephi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0_0709_10.root+hist_file_K0: ./hist_K0.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda_0709_10.root+hist_file_Lam: ./hist_Lambda.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich_0709_10.root+out_file: rich.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------<<<070810>>>>------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi.root+hist_file_iphi: ./hist_iphi_0708_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi.root+hist_file_ephi: ./hist_ephi_0708_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0.root+hist_file_K0: ./hist_K0_0708_10.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda.root+hist_file_Lam: ./hist_Lambda_0708_10.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich.root+out_file: rich_0708_10.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------ 
cp Plots_files/files_fit/hist_*_0708_10.root . 
fit_table -f options_fit/options_fit.dat  fit |& tee log_fit.log
mv *.root *.txt *.pdf *.log ./results_fit_theta_cols/P0708_10
# ------------------------------------------------------------------------------------
sed -i 's+hist_file_iphi: ./hist_iphi_0708_10.root+hist_file_iphi: ./hist_iphi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_ephi: ./hist_ephi_0708_10.root+hist_file_ephi: ./hist_ephi.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_K0: ./hist_K0_0708_10.root+hist_file_K0: ./hist_K0.root+g' options_fit/options_fit.dat
sed -i 's+hist_file_Lam: ./hist_Lambda_0708_10.root+hist_file_Lam: ./hist_Lambda.root+g' options_fit/options_fit.dat
sed -i 's+out_file: rich_0708_10.root+out_file: rich.root+g' options_fit/options_fit.dat
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
