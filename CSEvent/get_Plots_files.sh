# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
sed -i 's/analysis: K0L/analysis: phi/' options_fit/options_fit*
fit_table -f options_fit/options_fit.P05.dat -v plots
fit_table -f options_fit/options_fit.P06.dat -v plots
fit_table -f options_fit/options_fit.P07.dat -v plots
fit_table -f options_fit/options_fit.P08.dat -v plots
fit_table -f options_fit/options_fit.P09.dat -v plots
fit_table -f options_fit/options_fit.P10.dat -v plots
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
sed -i 's/analysis: phi/analysis: K0L/' options_fit/options_fit*
fit_table -f options_fit/options_fit.P05.dat -v plots
fit_table -f options_fit/options_fit.P06.dat -v plots
fit_table -f options_fit/options_fit.P07.dat -v plots
fit_table -f options_fit/options_fit.P08.dat -v plots
fit_table -f options_fit/options_fit.P09.dat -v plots
fit_table -f options_fit/options_fit.P10.dat -v plots
mkdir Plots_files							
mv *.root *.txt *.pdf ./Plots_files					
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

