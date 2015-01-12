python py_polar.py
python py_polar2.py
python py_polar3.py
python py_polar4.py
python py_polar_every.py
python py_polar_rotation.py
python py_polar_subplot.py

python py_timediffseries23.py #time series of diff hz at each position 23
python py_timediffseries12.py #time series of diff hz at each position 12

python py_diff_max.py>diff_max.d #extract max diff hz at each poisition
python py_diff_min.py>diff_min.d #extract min diff hz at each poisition

python py_diffmaxminplot.py #extract max min at each position (differential)

python py_abso_max.py>abso_max.d #extract max abso hz at each position
python py_abso_min.py>abso_min.d #extract min abso hz at each position

python py_absomaxminplot.py #extract max min at each position (absolute)

python py_lissajous_absolute.py #lissajous absolute
python py_lissajous_differential.py #lissajous differential
# python py_lissajous.py #absolute + differential

python py_model.py #sig,myu model

python py_diff.py # plot timeseries diff hz with and without crack at one 23 receiver position
python py_diff12.py # " at 12 receiver position

#residual
# py_residual.py # residual plus combine bad bad bad
python py_residual_inout.py #from nocrack and verXX to residual
python py_combine_residual.py #from residual to input format

#reverse time
python py_rev_residual.py # reverse residual time and combine

#combine after XXpattern2_12diff.d XXpattern2_23diff.d
python py_combine.py

#plot one time combine_pattern or combine_residual_abso
python py_pickup_one_time.py

#plot contour of combine_pattern or combine_residual_
python py_contour.py

#derivation of combine_pattern or comvine_residual_
python py_derivation.py

#matlab wiggle
# wigb_control_abso.m #wiggle absolute
# wigb_control_diff.m #wiggle differential
# wigb_control_diff_rev.m
# wigb_control_abso_rev.m
# wigb_control_12diff.m
# wigb_control_23diff.m