from bokeh.plotting import *
from bokeh.models import LinearAxis, Range1d
import datetime as dt
import pandas as pd
import numpy as np


# Read data

data=pd.read_csv('filtered_.csv')

uniq=np.unique(data['Instrument'])

for i in range(0,len(uniq)):
	prefix=uniq[i]
	name_file=prefix+'.html'
	output_file(name_file)
	subset_temp=data[data['Instrument']==prefix]
	x = [dt.datetime.strptime(d,'%d / %b %Y %H: %M') for d in subset_temp['Time']]
	p = figure(width=1500, height=550, x_axis_type="datetime",y_range=(-0.1, 1))
	p.scatter(x, subset_temp['OD'], color='navy', alpha=0.5)
	p.extra_y_ranges = {"foo": Range1d(start=10, end=50)}
	p.line(x, subset_temp['Temperature'], color='firebrick', alpha=0.5, y_range_name="foo")
	p.add_layout(LinearAxis(y_range_name="foo"), 'left')
	save(p)
	


