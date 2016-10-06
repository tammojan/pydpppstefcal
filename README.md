# Python-DPPP calibration demo

*Note: this code is not by any means production ready!*

This is a demonstration of how to call a python step from [DPPP](http://svn.astron.nl/LOFAR/). The data from DPPP is passed to python, where in this case stefcal is performed. The solutions are not yet saved.

The code is currently painfully slow, due to non-efficient reordering of the data.

**Usage:**: put the files `dpstefcal.py`, `stefcal.py` and `DPPP.parset` in a directory with a measurement set. Then just call `DPPP`.

Example output

```
MSReader
  input MS:       /data/scratch/dijkema/nikki/pydppp/tmp.MS
  band            0
  startchan:      0  (0)
  nchan:          40  (0)
  ncorrelations:  4
  nbaselines:     946
  ntimes:         105
  time interval:  4.00556
  DATA column:    DATA
  WEIGHT column:  WEIGHT_SPECTRUM
  autoweight:     false
PythonStep pystep. class=DPStefcal
MSUpdater msout.
  MS:             /data/scratch/dijkema/nikki/pydppp/tmp.MS
  datacolumn:     DATA
  weightcolumn    WEIGHT_SPECTRUM
  Compressed:     no

  flush:          0

Processing 105 time slots ...

0%....10....20....30....40....50....60....70....80....90....100%
Finishing processing ...

NaN/infinite data flagged in reader
===================================

Percentage of flagged visibilities detected per correlation:
  [0,0,0,0] out of 3973200 visibilities   [0%, 0%, 0%, 0%]
0 missing time slots were inserted
    Time spent in getting data:    0.476607322693
    Time spent in checking flags:  1.73473906517
    Time spent in reordering data: 679.8027215
    Time spent in actual stefcal:  11.7115616798



Total NDPPP time     694.62 real      692.83 user        0.67 system
    0.1% MSReader
   99.9% PythonStep pystep. class=DPStefcal
    0.0% MSUpdater msout.
```
