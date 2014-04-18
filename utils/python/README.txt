tccon.py: Library to read common tccon files

How to use it:

---
import tccon

# General TCCON file
tf = tccon.tccon_file('file.xxx')

# eof file
# Warning: some pre-2013 GGG versions produced broken eof files
# that will not be readable. This is not a tccon.py problem!
tf = tccon.eof_file('file.eof')

# col file
tf = tccon.col_file('file.col')

# gop file
tf = tccon.gop_file('file.gop')

# grl file
tf = tccon.grl_file('file.grl')

# map file
tf = tccon.map_file('file.map')

# oof file
tf = tccon.oof_file('file.oof')

# ray file
tf = tccon.ray_file('file.ray')

# tav file
tf = tccon.tav_file('file.tav')

# tsw file
tf = tccon.tsw_file('file.tsw')

# vav file
tf = tccon.vav_file('file.vav')

# vsw file
tf = tccon.vsw_file('file.vsw')

# file header as list (1 entry per line)
tf.header

# column names of the data block: names will be generic for tccon_file()
# other classes use real names from header
tf.columns

# complete data block
tf.data

# complete column from the data block as numpy vector
tf.data['colname']

# alternatively: 1 column (index 0) as numpy vector
tf.data[xxx.columns[0]]
---


Example: plot XCO2 over time from oof file

---
import tccon
oof = tccon.oof_file('file.oof')

import matplotlib.pyplot as plt
plt.plot(oof.data['day'] + oof.data['hour'], oof.data['xco2(ppm)'])
plt.show()
---

Note: tf.units is currently only defined for map files!