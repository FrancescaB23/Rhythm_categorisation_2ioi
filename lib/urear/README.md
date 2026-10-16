This is a wrapper code to comfortably run the auditory nerve and inferior colliculus models from the UR-EAR library. 

Requires [rnb_tools](https://github.com/TomasLenc/rnb_tools)

## Installation

If you've just downloaded the code onto a new machine, you need to compile the mex files. 

Open matlab, go into the `./UR_EAR2020b/source` folder, and run `mexANmodel` from matlab command line. 

This will genearte some mex files. You just need to move everything to the `UR_EAR2020b` folder. You can use this command:  

``` bash
mv ./UR_EAR2020b/source/*.mex* ./UR_EAR2020b
```

