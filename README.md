# wisp_analysis

The WISP emission line detection and verification code for slitless 
grism spectra. 

This package detects emission line candidates in WISP 
spectra using a peak finder that performs a continuous wavelet transform
(CWT) on the data
(see SciPys's [find_peaks_cwt](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks_cwt.html)).
The code then iterates through the detected emission lines, allowing the 
reviewers to identify lines or reject spurious detections quickly and 
efficiently. 

This package merges line fitting code from Alaina Henry with the interactive
redshift cataloger from Nathaniel R. Ross. 


## Getting Started

### Prerequisites

* Python 2.7  (Python 3 is not supported)
* DS9
* XPA Messaging System

We recommend using [astroconda](https://astroconda.readthedocs.io/en/latest/) - 
the Conda channel maintained by STScI. Everything should also work if you are 
managing your own Python distribution or are using 
[Ureka](http://ssb.stsci.edu/ureka/).

`wisp_analysis` uses the XPA Messaging System to display the 2d grism stamps 
in DS9 with their wavelength solutions as well as to pan direct images based 
on source coordinates. You will specifically need the executables 
`xpans`, `xpaget`, and `xpaset`. 
You can download and install XPA from 
[GitHub](https://github.com/ericmandel/xpa). 

If you are using a Mac and run into trouble, try these [tips](http://staff.washington.edu/rowen/ds9andxpa.html#Installing). 

For more information about XPA, see the help pages [here](http://hea-www.harvard.edu/RD/xpa/).


### Setting up wisp_analysis
Either download and unzip the repo, or clone it:

```
git clone git@github.com:HSTWISP/wisp_analysis.git
```

Add wisp\_analysis to your $PYTHONPATH. 
For example, if you put wisp\_analysis in /Users/ahenry/software, then 
in your bashrc:
```
export $PYTHONPATH="/Users/ahenry/software:$PYTHONPATH"
``` 
Or in your cshrc:
```
setenv PYTHONPATH /Users/ahenry/software:${PYTHONPATH}
```

### Running wisp_analysis

* Go to the `Spectra` subdirectory of a WISP field. 
* Copy the `default.config` file to this directory
* Open a ds9 window
* Launch Python or ipython
* To load the package:
```
import wisp_analysis as wisp
```
* To run the peak finder and create a file of emission line candidates:
```
wisp.loop_field_cwt()
```
* If an emission line candidate file (`linelists/ParXXXlines.dat`) has already 
been created, run the interactive part of the program to inspect the emission 
lines:
```
wisp.measure_z_interactive()
```

