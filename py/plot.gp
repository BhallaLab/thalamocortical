# plot 'data/2009_05_06/Vm.plot' u ($0*1e-3):($1*1e3) w l, '../nrn/mydata/Vm.plot' u ($1):($2) w l
# plot 'data/2009_04_24/m.plot', 'data/2009_04_22/m.bak.plot'
# plot 'cal_xa.plot', 'cal_xa.plot.bak'
# plot 'cal_xb.plot', 'cal_xb.plot.bak'
# plot 'data/2009_04_25/Ca.plot' u ($0*1e-2):($1*1e3) w l, '../nrn/mydata/Ca.plot' u ($1):($2) w l
# plot 'data/2009_04_24/Ca.plot' u ($0*1e-2):($1)
# plot 'data/2009_04_25/m_kahp.plot' u ($0*1e-2):($1/9.42e-6), '../nrn/mydata/Vm.plot' u ($1):($3) 
# plot 'data/2009_04_25/m_kahp.plot' u ($0*1e-2):($1)
# plot 'beta.txt' u ($1*1e3):($2), '../nrn/mydata/Vm.plot' u ($2):($3)
# plot  '../nrn/mydata/Vm.plot' u ($2):($3)
# plot '~/src/sim/cortical/nrn/dat/spinstell_v_F.dat' w l, 'data/2009_05_03/Vm.plot' u ($0*1e-3):($1*1e3) w l

# plot '~/src/sim/cortical/nrn/dat/spinstell_v_F.dat' w l

# plot 'data/2009_05_06/Gk_NaF2.plot' u ($0 * 1e-3): ($1) w l, '../nrn/mydata/Vm.plot' u ($1): ($2) w l
# plot 'data/2009_05_06/Gk_NaF2.plot'
# plot 'data/2009_05_06/Vm.plot'

set terminal svg size 600, 800
set output '20111018_short_term_plasticity.svg'
set multiplot
set size 1,0.25
set origin 0, 0.75
set ytics 0.5
set xlabel 'Time (s)'
plot './_model_net_DeepBasket_0_comp_30_ampa_from_TCR.Pr.txt' u ($0*1e-4):($1) w l ti 'AMPA: Release probability'
set origin 0, 0.5
set ytics 5
set xlabel 'Time (s)'
plot './_model_net_DeepBasket_0_comp_30_ampa_from_TCR.F.txt' u ($0*1e-4):($1) w l ti 'AMPA.F'
set origin 0, 0.25
set ytics 0.5
set xlabel 'Time (s)'
plot './_model_net_DeepBasket_0_comp_30_ampa_from_TCR.D1.txt' u ($0*1e-4):($1) w l ti 'AMPA.D1'
set ytics 0.5
set origin 0, 0.0
set xlabel 'Time (s)'
plot './_model_net_DeepBasket_0_comp_30_ampa_from_TCR.D2.txt' u ($0*1e-4):($1) w l ti 'AMPA.D2'

