# PulseBlast
A basic python package that deals with pulsar data manipulation in a variety of ways. Options include the ability to create Times-Of-Arrival (TOAs), excise RFI from pulse scans and easily create pulsar template profiles.

## **Mac Setup:**

After cloning both PulseBlast and PyPulse (as needed), add PulseBlast and PyPulse to your PYTHONPATH in .bashrc as follows:

```bash
export PYTHONPATH=[directory_for_python3]:[directory_for_pypulse]:[directory_for_pulseblast]:$PYTHONPATH
```

and add an update check in .bash_profile (if not already set up):

```bash
test -f ~/.bashrc && source ~/.bashrc
```

## **Usage:**

### **Timing**

To get Times-of-Arrival (TOAs), run the following command in the terminal:

```shell
python main.py -x [text_files_containing_directories_and_/_or_files] -t [frequency_band] --temp [full_path_to_template] -s [sub-integrations_to_scrunch_to] -n [sub-bands_to_scrunch_to] -j [jump_after_fluxerr] -od [toa_output_file_directory] -o [toa_output_filename] -r [number_of_generations] -v
```

Text files should only contain directories (ending in a "/") or files (ending in a ".???"). E.g. `input.txt`:

```
/Users/User1/Pulsars/directory1/
/Users/User2/PSRs/directory/
/Users/User1/Pulsars/directory2/file1.fits
/Users/User1/Pulsars/directory2/file2.fits
/Users/User2/otherPSRs/profile.fits
```

Directories / files inside the text file need not be in any particular order.

The jump flag, `-j`, is used to add a string to the end of each line of the output TOA file. These are most likely to be jumps (hence the name) but the string itself can be anything supplied by the user.  

TOAs will be saved in TEMPO2 format to the output path and output filename given however, currently, the telescope name *will* need to be changed by the user to the newly used TEMPO2 codes until I can figure out how to implement this. If no output path is given, the current working directory will be used. If no filename is supplied, a default filename, `PSR_TOAs.toa`, will be used.  

The final two arguments denote the rejection and verbose flag respectively. The rejection flag plus its argument runs an RFI excision algorithm to the number of generations supplied by the user. If `-r` is not set, the timing will happen without any RFI excision. The verbose flag, `-v`, if set, will display more detailed information to the user about what is being loaded, how long tasks take, as well as many other features that might be useful for developers.  

### **Templates**

**Creating templates**

*Currently, GUI options provided by gooey ONLY work with* **[my fork of gooey](https://github.com/HenrykHaniewicz/Gooey "HenrykHaniewicz/Gooey")** *until my pull request goes through*

If you wish you use Graphical User Interface (GUI) features, `pip install gooey`.
The GUI is a simple field based argument handler to select files and directories from.  

If you are running Python in a [virtual environment](https://docs.python.org/3/tutorial/venv.html "Virtual environment documentation") within the shell, you will need to make you're running a **framework version** of Python. To do this, edit the [fwpy](https://github.com/HenrykHaniewicz/PulseBlast/blob/master/fwpy "Framework bash file") file as necessary and paste it into the directory where your Python executable is. Once this is done, make sure `fwpy` has executable permissions. You can set this as follows:

```shell
chmod a+x fwpy
```

To run the program in GUI mode, provided you have followed the setup above, type the command:

```shell
fwpy PSRTemplate.py
```

If you do not wish to use the GUI, the template program can also be used via Command Line Interface (CLI). To create templates using CLI, use the following command in the terminal:

```shell
python PSRTemplate.py -b [frequency_band] -d [directories_to_search_for_psrfits_files_in] -o [output_directory_and_filename]
```

Directories parsed to this command can either be local to the current working directory or absolute paths.

**Deleting templates**

If you need to delete a template in your code, you can run the `deleteTemplate()` method after initializing an instance of the Template class in your code (here called `templateObject`). This method takes in a required filename and **full path** directory where the template can be found.  
The syntax is as follows:

```python
templateObject.deleteTemplate( filename_of_template, full_path_to_directory )
```

**Warning**: The `deleteTemplate()` method *can* be used to delete anything so please use with care. The program will also warn the user before deleting.

## **Requires:**  

Python 3.X  

PyPulse  
numpy  
scipy  
astropy  
matplotlib  
filemagic  

**Required for GUI:**

gooey  



It is suggested to have the most up-to-date versions of the packages listed. To upgrade all python packages at once, try the following command:

```shell
pip list --outdated --format=freeze | grep -v '^\-e' | cut -d = -f 1  | xargs -n1 pip install -U
```
