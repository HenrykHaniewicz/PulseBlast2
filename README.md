# PulseBlast V2
A basic python package that deals with pulsar data manipulation in a variety of ways. Options include the ability to create Times-Of-Arrival (TOAs), excise RFI from pulse scans, calibrate data and easily create pulsar template profiles.

## **Mac Setup:**

After cloning both PulseBlast and PyPulse (as needed), add PulseBlast and PyPulse to your PYTHONPATH in .bashrc as follows:

```bash
export PYTHONPATH=[directory_for_python3]:[directory_for_pypulse]:[directory_for_pulseblast]:$PYTHONPATH
```

and add an update check in ~/.bash_profile (if not already set up):

```bash
test -f ~/.bashrc && source ~/.bashrc
```

## **Usage:**

### **Timing**

### **RFI Excision**

### **Calibration**

### **Templates**


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
