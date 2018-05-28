from bokeh.plotting import *
from bokeh.models import LinearAxis, Range1d
import datetime as dt
import pandas as pd
import numpy as np
import datetime as dt
import pandas as pd
import numpy as np


# Read data
data=pd.read_csv('bench_2.csv')
# Strain names
uniq=np.unique(data['Instrument'])
# Subset the data
i=6 #number to choose from
prefix=uniq[i]
name_file=prefix+'march.html'
output_file(name_file)
test=data[data['Instrument']==prefix]
#index dataframe by time
x = [dt.datetime.strptime(d,'%d / %b %Y %H: %M') for d in test['Time']]
test.index=x

#Subset the timeseries

lol=test.ix['2016-02-22':'2016-03-30']

#Count the drops

p = figure(width=1500, height=550, x_axis_type="datetime",y_range=(-0.1, 1))
p.scatter(lol.index,lol['OD'], color='navy', alpha=0.5)
p.extra_y_ranges = {"foo": Range1d(start=10, end=50)}
p.line(lol.index,lol['Temperature'], color='firebrick', alpha=0.5, y_range_name="foo")
p.add_layout(LinearAxis(y_range_name="foo"), 'left')
save(p)

