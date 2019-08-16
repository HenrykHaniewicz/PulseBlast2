# PulseBlast V2
A basic python package that deals with pulsar data manipulation in a variety of ways. Options include the ability to create Times-Of-Arrival (TOAs), excise RFI from pulse scans, calibrate data and easily create pulsar template profiles.

## **Mac Setup:**

After cloning both PulseBlast2 and PyPulse (and installing packages as needed), add PulseBlast2 and PyPulse to your PYTHONPATH in .bashrc as follows:

```bash
export PYTHONPATH=[directory_for_python3]:[directory_for_pypulse]:[directory_for_pulseblast2]:$PYTHONPATH
```

## **Windows Setup:**

Coming soon... maybe

## **Usage:**

Run `pulseblast2.py` in the terminal to start PulseBlast2. You will be met by the following command prompt:

```
Welcome to PulseBlast V2

Please type your command in here...
```

Type `help` to get a full list of features and input flags.

The first two commands must be the directories to search through (can be `None` for a default option) followed by the pulsar name as given in the PSRFITS files you wish to analyze.

Main block commands (-r, -c, -t) are interpreted in order, so `dirs psr_name -r -c -t` will do RFI excision followed by calibration and then will create TOAs, whereas `dirs psr_name -c -r -t` will calibrate, then do RFI excision before making TOAs. All main block flags are optional.

### **Timing**

Coming soon.

### **RFI Excision**

`-r/-rs/-rb/-rn [-i iterations] [-q "temp_dir"] [-g "saveddata_dir"] [-e]`

`-r`                          RFI mitigation (defaults to sigma-clipping)

`-rs`                         RFI mitigation by sigma-clipping

`-rb`                         RFI mitigation by Bayesian interference modulation (not-imp)

`-rn`                         RFI mitigation by deep-learning neural network 1D image recognition (not-imp)

`-i iterations`               Specifies number excision iterations {int} (optional)

`-q "temp_dir"`               Specifies the template directory {str} (optional)

`-g "saveddata_dir"`          Specifies the directory to save data to {str} (optional)

`-e`                          Specifies whether to average in the time domain before RFI excision (defaults to False if not given)

### **Calibration**

`-c -n "cont_name" [-w "cont_dir"] [-d "saveddata_dir"]`

`-c`                          Flux calibration

`-n "cont_name"`              Specifies continuum source name {str} (required)

`-w "cont_dir"`               Specifies continuum source cal file directory {str} (optional)

`-d "saveddata_dir"`          Specifies the directory to save data to {str} (optional)


### **Templates**

Builder:

`-m -f "frontend" -s subbands [-p "template directory"]`

`-m`                          Template builder

`-f "frontend"`               Specifies the observation frontend {str} (required)

`-s subbands`                 Specifies the number of sub-bands to make templates for {int} (required)

`-p "template directory"`     Specifies the template output directory {str} (optional)

Plotter:

Run `python template_plotter.py` as a separate program once you have a template.

## **Requires:**  

Python 3.X  

PyPulse  
numpy  
scipy  
astropy  
matplotlib



## **Notes:**
It is suggested to have the most up-to-date versions of the packages listed. To upgrade all python packages at once, try the following command:

```shell
pip list --outdated --format=freeze | grep -v '^\-e' | cut -d = -f 1  | xargs -n1 pip install -U
```
