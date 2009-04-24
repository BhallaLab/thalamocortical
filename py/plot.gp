plot 'data/2009_04_24/Vm.plot' u ($0*1e-2):($1*1e3) w l, '../nrn/mydata/Vm.plot' u ($1):($2) w l
plot 'data/2009_04_24/m.plot', 'data/2009_04_22/m.bak.plot'
plot 'cal_xa.plot', 'cal_xa.plot.bak'
plot 'cal_xb.plot', 'cal_xb.plot.bak'
plot 'data/2009_04_25/Ca.plot' u ($0*1e-2):($1*1e3) w l, '../nrn/mydata/Ca.plot' u ($1):($2) w l
plot 'data/2009_04_24/Ca.plot' u ($0*1e-2):($1)
plot 'data/2009_04_25/m_kahp.plot' u ($0*1e-2):($1/9.42e-10) , '../nrn/mydata/Vm.plot' u ($1):($3) 
