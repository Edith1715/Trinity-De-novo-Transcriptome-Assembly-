# Import and run the RNA-Seq Workshop VM

Once you've installed [VirtualBox](https://www.virtualbox.org/wiki/Downloads) and downloaded the RNA-Seq workshop VM [Trinity2015.ova](ftp://ftp.broadinstitute.org/pub/Trinity/RNASEQ_WORKSHOP/Trinity2015.ova), complete the following steps to import and run the VM.

## Import the VM Appliance

Run the VirtualBox application.  From the VirtualBox menu, select the menu option 'File' -> 'Import Appliance...':

<img src="https://raw.githubusercontent.com/wiki/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/images/VM_install/import_appliance_menu.png" width=300 />

After which, you'll see the 'Appliance to Import' dialog':

<img src="https://raw.githubusercontent.com/wiki/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/images/VM_install/import_appliance_dialog.png" width=450 />

Click the folder icon next to the input box and browse your file system to locate where you previously downloaded the Trinity2015.ova file.  Select the file, and click 'Continue'.

<img src="https://raw.githubusercontent.com/wiki/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/images/VM_install/import_appliance_file_select.png" width=450 />

Now, in the main VirtualBox application window, you should see the 'Trinity2015' VM listed on the left. If you've never used VirtualBox before, then this may be your only option.  Otherwise, you'll see it among a list of the other VMs you have installed (such as is the case for me).

<img src="https://raw.githubusercontent.com/wiki/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/images/VM_install/appliance_loaded.png" width=450 />


## Start the VM

To start the VM, simply double-click the VM button.  It should load up a recent version of Ubuntu Linux, and you should see a nice clean desktop:

<img src="https://raw.githubusercontent.com/wiki/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/images/VM_install/appliance_started.png" width=450 />

To be ready to start the workshop, click the black '>_' icon on the left to open up a terminal window.  List the files via typing the 'ls' command in the terminal, and you should see the directory 'RNASEQ_Trinity_Tuxedo_Workshop/'.

<img src="https://raw.githubusercontent.com/wiki/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/images/VM_install/opened_terminal.png" width=450 />


Now you can begin one of the workshops.  Have fun!! :)
