# Time At Temperature
This python script calculates predicted time for a thermosetting resin to cure, given data from a differential scanning calorimeter (DSC).

## Motivation
Carbon fiber / glass fiber pultrusion processes use fast curing resin systems, and are subject to a fixed die length, a variable temperature,
and a variable speed at which the fiber and resin are pulled through the die. The resin must be at a sufficient state of cure leaving the die
for stability, which may vary depending on the materials and geometry. Ideally, the die length, temperature, and pull speed can be optimized
prior to investing in the tooling. However, the thermal mechanics of the resin system cannot be measured directly due to the fast cure times.

## Solution
Various methods have been developed by the chemistry community. The isoconversional techniques employed here do not require knowledge of the
precise chemical mechanisms. As the resin cures, it releases heat energy. The amount of heat released as a fraction of the total amount of 
heat that will get released during cure provides the current degree of cure. The rate at which the degree of cure advances is the information
we want, understanding that it depends only on the current degree of cure and the current temperature.

## More Details
Please see the presentation in the "images" folder, or the python code TimeAtTemperature.py for more information!

## Running it Yourself
All that is required is the TimeAtTemperature.py script and the contents in the data folder. To run this program using different data, the 
pandas.cvs_read functions will need to be changed to refer to your data. The data file format must have an index column, a time column in s, 
a temperature column in degrees C, and a measurement column in mW. See the existing data files for example formatting. The .sln file and 
.pyproj are overhead for Visual Studio, if you choose to run it that way.

## Reference
 * Sbirrazzuoli, Nicolas, et al. "Integral, differential and advanced isoconversional methods: complex
   mechanisms and isothermal predicted conversion–time curves." _Chemometrics and Intelligent Laboratory
   Systems 96.2 (2009)_: 219-226.