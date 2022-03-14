# How to use the MS converter

To convert `*.wiff` and `*.wiff.scan` files to `mzML` file, use https://gnps-quickstart.ucsd.edu/conversion, upload there both files and click **Convert Uploaded Files**, then  extract the downloaded archive and convert the `mzML` file with `mzML2csv.py` converter. If necessary, change the resolution of m/z in the `mzML2csv.py`. Default value is 1.

The converter only converts MS data (from HPLC experiment or using only syringe). If `*.wiff` contains some HPLC UV data, they will not be converted.