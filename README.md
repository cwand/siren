# siren
SImple RENography analysis of dynamic SPECT or PET images

## Getting siren
Clone the repository to your computer using git:
```
> git clone https://github.com/cwand/siren
```

Enter the directory.
Make sure you are on the main branch:
```
> git checkout main
```

Create a new virtual python environment:
```
> python -m venv my_venv
```

Activate the virtual environment. Commands vary according to OS and shell 
(see [the venv documentation](https://docs.python.org/3/library/venv.html)), 
but in a Windows PowerShell:
```
> my_venv\Scripts\Activate.ps1
```

Install siren and required dependencies
```
> pip install .
```

If everything has gone right, you should be able to run tictac
```
> python -m siren
Starting SIREN 1.0.1
...
```

## Using siren
First and foremost: a help message is displayed when running siren with the 
```-h``` flag:
```
> python -m siren -h
```

To use siren you need
* A dynamic SPECT or PET dicom image series in a directory, say ```img_dir```
* A ROI file for each ROI you want to analyse, e.g. ```left_kidney.nii```, 
  ```left_cortex.nii```, ```right_kidney.nii``` and ```right_cortex.nii```.
  The ROI-files and the dynamic images must be in the same physical image space.
  (Any file-format that can be read by SimpleITK will do).

To get a renogram (time activity curve) for each region and a table printed
for all available kidney metrics, we run siren
```
> python -m siren -i img_dir -r left_kidney.nii left_cortex.nii right_kidney.nii right_cortex.nii
```

All regions of interest must be in a separate ROI-file, and in each ROI-file,
the structure is read as the value ```1``` in the ROI-image. Other values are
ignored.

After this call, a result file ```siren.txt``` will be generated and a renogram
shown.
