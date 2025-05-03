# AstroPyPanchanga
Calculate _pañcāṅga_ at given time and location using [astropy](https://www.astropy.org/). 

This set of codes also calucates _tithi_, _vāra_, _nakṣatra_, _yoga_ and _karaṇa_ that make up the _pañcāṅga_.

The included `plot_panchanga.py` script allows to visualise the positions of the _grahas_. Here is an example output of the `make_circle_plot` function :  
<p align="center">
<img src="https://github.com/arvindb95/AstroPyPanchanga/blob/main/panchanga_at_test_time-1.jpg" title="Sample circle plot" />
</p>

The code also can calculate the _jātakaṃ_ (in the south indian traditional format) at the given time. Use the `make_jatakam_plot` function for this :
<p align="center">
<img src="https://github.com/arvindb95/AstroPyPanchanga/blob/main/jatakam_at_test_time-1.jpg" title="Sample jatakam plot"  width="50%"/>
</p>

The code by default prints all information in _Devanāgarī_. But it supports transliteration to the following indian languages using the [indic_transliteration](https://github.com/indic-transliteration/indic_transliteration_py) package :
- Grantha
- Kannada
- Malayalam
- Tamil
- Telugu


## Dependancies

### Python Packages 
- The file `requirements.txt` lists the python dependencies that are required to run this code.
   
### External Requirements
- [XeTeX](https://tug.org/xetex/) - For rendering information in indian languages.  
- Noto serif fonts for all the above supported languages

_Edit the LaTeX preamble as you wish to change the typesetting (you might need to change the files that have the names of naks&#803;atras etc.)_

## To do
- Better comments on code for readability 
- Calculate _pañcāṅga_ for the whole day to get extent of each _aṅga_ and when it changes 
